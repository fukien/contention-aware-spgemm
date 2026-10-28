#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
DATASET_DIR="$ROOT/dataset"
mkdir -p "$DATASET_DIR"
WORK_DIR=$(mktemp -d "${TMPDIR:-/tmp}/join-aggregate-work.XXXXXX")
OUT_DIR="$ROOT/dataset"

mkdir -p "$OUT_DIR"
trap 'rm -rf -- "$WORK_DIR"' EXIT

need_command() {
	command -v "$1" >/dev/null 2>&1 || {
		echo "Required command not found: $1" >&2
		exit 1
	}
}

need_command curl
need_command python3

download() {
	local url=$1 destination=$2 kind=$3 partial="${2}.part"
	if [[ -s "$destination" ]]; then
		if [[ "$kind" == "gzip" ]] && gzip -t "$destination" >/dev/null 2>&1; then
			echo "Using existing $destination"
			return
		elif [[ "$kind" == "zip" ]] && python3 - "$destination" <<'PY'
import sys, zipfile
try:
	with zipfile.ZipFile(sys.argv[1]) as archive:
		ok = archive.testzip() is None
except Exception:
	ok = False
raise SystemExit(0 if ok else 1)
PY
		then
			echo "Using existing $destination"
			return
		fi
		echo "Existing download is incomplete or invalid; resuming it: $destination"
	fi
	echo "Downloading $url"
	curl --fail --location --retry 3 --continue-at - --output "$partial" "$url"
	if [[ "$kind" == "gzip" ]]; then
		gzip -t "$partial"
	else
		python3 - "$partial" <<'PY'
import sys, zipfile
with zipfile.ZipFile(sys.argv[1]) as archive:
	bad = archive.testzip()
if bad:
	raise SystemExit(f"Corrupt member in ZIP: {bad}")
PY
	fi
	mv -- "$partial" "$destination"
}

# ClickBench's official TSV is the source for the user–URL interaction sample.
download \
	https://datasets.clickhouse.com/hits_compatible/hits.tsv.gz \
	"$DATASET_DIR/clickbench-hits.tsv.gz" gzip

# GroupLens Book Genome ratings and IMDb's movie/principal tables.
download \
	https://files.grouplens.org/datasets/book-genome/book-genome.zip \
	"$DATASET_DIR/book-genome.zip" zip
download \
	https://datasets.imdbws.com/title.basics.tsv.gz \
	"$DATASET_DIR/title.basics.tsv.gz" gzip
download \
	https://datasets.imdbws.com/title.principals.tsv.gz \
	"$DATASET_DIR/title.principals.tsv.gz" gzip

echo "Building chbc_s_uni (CLH0), bkgn (BKGN), and imdb_movie (IMDB)."
echo "Output matrices will be written under dataset/ for this repo's SpGEMM runners."

python3 - \
	"$DATASET_DIR/clickbench-hits.tsv.gz" \
	"$DATASET_DIR/book-genome.zip" \
	"$DATASET_DIR/title.basics.tsv.gz" \
	"$DATASET_DIR/title.principals.tsv.gz" \
	"$OUT_DIR" \
	"$WORK_DIR" <<'PY'
import csv
import gzip
import json
import os
import random
import sqlite3
import sys
import zipfile

clickbench_path, bookgenome_path, basics_path, principals_path, out_dir, work_dir = sys.argv[1:]
SAMPLE_USERS = 850_000
SEED = 42
BATCH_SIZE = 20_000


def connect(path):
	conn = sqlite3.connect(path)
	conn.execute("PRAGMA journal_mode=OFF")
	conn.execute("PRAGMA synchronous=OFF")
	conn.execute("PRAGMA temp_store=FILE")
	conn.execute("PRAGMA cache_size=-262144")
	conn.execute("CREATE TABLE edges (r INTEGER NOT NULL, c INTEGER NOT NULL, PRIMARY KEY (r,c)) WITHOUT ROWID")
	return conn


def intern(mapping, value):
	idx = mapping.get(value)
	if idx is None:
		idx = len(mapping)
		mapping[value] = idx
	return idx


def flush_edges(conn, pending):
	if pending:
		conn.executemany("INSERT OR IGNORE INTO edges VALUES (?,?)", pending)
		conn.commit()
		pending.clear()


def export_matrix(conn, name, nrows, ncols):
	nnz = conn.execute("SELECT count(*) FROM edges").fetchone()[0]
	if not nrows or not ncols or not nnz:
		raise RuntimeError(f"{name}: empty matrix ({nrows} rows, {ncols} cols, {nnz} entries)")
	paths = {
		"org": os.path.join(out_dir, f"{name}.org.mtx"),
		"trs": os.path.join(out_dir, f"{name}.trs.mtx"),
	}
	for kind, path in paths.items():
		tmp = path + ".tmp"
		if kind == "org":
			query = "SELECT r,c FROM edges ORDER BY r,c"
			dimensions = (nrows, ncols)
		else:
			query = "SELECT c,r FROM edges ORDER BY c,r"
			dimensions = (ncols, nrows)
		with open(tmp, "w", encoding="ascii", newline="\n") as stream:
			stream.write("%%MatrixMarket matrix coordinate real general\n")
			stream.write("% binary incidence matrix prepared for contention-aware SpGEMM\n")
			stream.write(f"{dimensions[0]} {dimensions[1]} {nnz}\n")
			for row, col in conn.execute(query):
				stream.write(f"{row + 1} {col + 1} 1.0\n")
		os.replace(tmp, path)
	print(f"{name}: {nrows:,} x {ncols:,}, {nnz:,} entries")


def header_indexes(fieldnames, required, path):
	if not fieldnames:
		raise RuntimeError(f"Missing TSV header in {path}")
	lookup = {name.strip().lower(): name for name in fieldnames}
	missing = [name for name in required if name.lower() not in lookup]
	if missing:
		raise RuntimeError(f"{path} is missing columns {missing}; found {fieldnames}")
	return {name: lookup[name.lower()] for name in required}


def prepare_clickbench():
	print("[1/3] Processing ClickBench hits and drawing the CLH0 user sample...")
	dbpath = os.path.join(work_dir, "chbc_s_uni.sqlite")
	conn = sqlite3.connect(dbpath)
	conn.execute("PRAGMA journal_mode=OFF")
	conn.execute("PRAGMA synchronous=OFF")
	conn.execute("PRAGMA temp_store=FILE")
	conn.execute("PRAGMA cache_size=-262144")
	conn.execute("CREATE TABLE users (id TEXT PRIMARY KEY) WITHOUT ROWID")

	batch = []
	with gzip.open(clickbench_path, "rt", encoding="utf-8", errors="replace", newline="") as source:
		reader = csv.DictReader(source, delimiter="\t")
		cols = header_indexes(reader.fieldnames, ["UserID"], clickbench_path)
		for row in reader:
			user = row.get(cols["UserID"])
			if user and user != "\\N":
				batch.append((user,))
			if len(batch) >= BATCH_SIZE:
				conn.executemany("INSERT OR IGNORE INTO users VALUES (?)", batch)
				conn.commit()
				batch.clear()
	if batch:
		conn.executemany("INSERT OR IGNORE INTO users VALUES (?)", batch)
		conn.commit()

	# Assign dense indices in stable ID order, then sample users uniformly.
	conn.execute("CREATE TABLE dense_users (idx INTEGER PRIMARY KEY, id TEXT UNIQUE)")
	conn.execute("INSERT INTO dense_users SELECT row_number() OVER (ORDER BY id), id FROM users")
	conn.commit()
	user_count = conn.execute("SELECT count(*) FROM dense_users").fetchone()[0]
	if user_count < SAMPLE_USERS:
		raise RuntimeError(f"ClickBench has only {user_count} users; need {SAMPLE_USERS}")
	chosen = random.Random(SEED).sample(range(1, user_count + 1), SAMPLE_USERS)
	row_map = {}
	for start in range(0, len(chosen), 500):
		part = chosen[start:start + 500]
		placeholders = ",".join("?" for _ in part)
		for idx, user in conn.execute(
			f"SELECT idx,id FROM dense_users WHERE idx IN ({placeholders})", part
		):
			row_map[user] = idx - 1
	conn.close()

	conn = connect(dbpath)
	url_map = {}
	batch = []
	with gzip.open(clickbench_path, "rt", encoding="utf-8", errors="replace", newline="") as source:
		reader = csv.DictReader(source, delimiter="\t")
		cols = header_indexes(reader.fieldnames, ["UserID", "URL"], clickbench_path)
		for row in reader:
			user = row.get(cols["UserID"])
			url = row.get(cols["URL"])
			if user not in row_map or not url or url == "\\N":
				continue
			batch.append((row_map[user], intern(url_map, url)))
			if len(batch) >= BATCH_SIZE:
				flush_edges(conn, batch)
	flush_edges(conn, batch)
	export_matrix(conn, "chbc_s_uni", SAMPLE_USERS, len(url_map))
	conn.close()


def prepare_bookgenome():
	print("[2/3] Processing GroupLens Book Genome ratings...")
	dbpath = os.path.join(work_dir, "bkgn.sqlite")
	conn = connect(dbpath)
	users, books = {}, {}
	pending = []
	with zipfile.ZipFile(bookgenome_path) as archive:
		members = [p for p in archive.namelist() if p.endswith("book_dataset/raw/ratings.json")]
		if not members:
			raise RuntimeError("book-genome.zip has no book_dataset/raw/ratings.json")
		with archive.open(members[0]) as raw:
			for line_num, encoded in enumerate(raw, 1):
				record = json.loads(encoded)
				user = str(record["user_id"])
				book = str(record.get("book_id", record.get("item_id")))
				pending.append((intern(users, user), intern(books, book)))
				if len(pending) >= BATCH_SIZE:
					flush_edges(conn, pending)
				if line_num % 1_000_000 == 0:
					print(f"  read {line_num:,} ratings")
	flush_edges(conn, pending)
	export_matrix(conn, "bkgn", len(users), len(books))
	conn.close()


def prepare_imdb():
	print("[3/3] Processing IMDb movie principals...")
	dbpath = os.path.join(work_dir, "imdb_movie.sqlite")
	conn = connect(dbpath)
	movies = set()
	with gzip.open(basics_path, "rt", encoding="utf-8", errors="replace", newline="") as source:
		reader = csv.DictReader(source, delimiter="\t")
		for record in reader:
			if record.get("titleType") == "movie":
				movies.add(record["tconst"])
	print(f"  found {len(movies):,} movie titles")

	titles, people = {}, {}
	pending = []
	with gzip.open(principals_path, "rt", encoding="utf-8", errors="replace", newline="") as source:
		reader = csv.DictReader(source, delimiter="\t")
		for line_num, record in enumerate(reader, 1):
			title = record.get("tconst")
			person = record.get("nconst")
			if title in movies and person and person != "\\N":
				pending.append((intern(titles, title), intern(people, person)))
			if len(pending) >= BATCH_SIZE:
				flush_edges(conn, pending)
			if line_num % 10_000_000 == 0:
				print(f"  scanned {line_num:,} principal records")
	flush_edges(conn, pending)
	export_matrix(conn, "imdb_movie", len(titles), len(people))
	conn.close()


os.makedirs(out_dir, exist_ok=True)
prepare_clickbench()
prepare_bookgenome()
prepare_imdb()
print(f"Prepared all three original/transpose pairs in {out_dir}")
PY

echo "Done. Run from the repository root, e.g.:"
echo "  ./bin/ams_hybaccspgemm chbc_s_uni"
echo "  ./bin/ams_hybaccspgemm bkgn"
echo "  ./bin/ams_hybaccspgemm imdb_movie"
