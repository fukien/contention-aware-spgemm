#include "mm_utils.h"

#ifdef IN_VERIFY
int64_t over_predict = 0;
int64_t under_predict = 0;
#endif /* IN_VERIFY */

#ifdef IN_DEBUG
pthread_barrier_t debug_barrier;
int debug_rv = 0;
#endif /* IN_DEBUG */


// Initialize status array for 'num_statuses' entries

// Set the status of a specific index

// Get the status of a specific index

// Clear (reset to STATUS_0) the status of a specific index

// Free the allocated memory

// Debug function to print status array

void free_csr(csr_t *csr) {
	dealloc_memory( csr->rats, sizeof(my_rat_t) * csr->num_nnz );
	dealloc_memory( csr->cols, sizeof(my_cnt_t) * csr->num_nnz );
	dealloc_memory( csr->row_ptrs, sizeof(my_cnt_t) * (csr->num_row+1) );
	csr->num_row = 0;
	csr->num_col = 0;
	csr->num_nnz = 0;
}

void swap(triple_t *a, triple_t *b) {
	triple_t t = *a;
	*a = *b;
	*b = t;
}

int compare_mycnts(const void *a, const void *b) {
	my_cnt_t *cnt_a = (my_cnt_t *) a;
	my_cnt_t *cnt_b = (my_cnt_t *) b;
	return *cnt_a - *cnt_b;
}

int compare_triples(const void *a, const void *b) {
	triple_t *trpl_a = (triple_t *) a;
	triple_t *trpl_b = (triple_t *) b;
	if (trpl_a->user_id != trpl_b->user_id) {
		return (trpl_a->user_id - trpl_b->user_id);
	} else {
		return (trpl_a->item_id - trpl_b->item_id);
	}
}

int compare_keyvals_3(const void *a, const void *b) {	// val first (ascending), then key
	keyval_t *kv_a = (keyval_t *) a;
	keyval_t *kv_b = (keyval_t *) b;

	int val_cmp = (kv_a->val > kv_b->val) - (kv_a->val < kv_b->val);
	int key_cmp = (kv_a->key > kv_b->key) - (kv_a->key < kv_b->key);

	return val_cmp ? val_cmp : key_cmp;
}

void check_sorted_triples(triple_t *trpl, int64_t num_nnz) {
	for (int64_t idx = 0; idx < num_nnz -1; idx ++) {
		if (compare_triples(trpl+idx, trpl+idx+1) > 0) {
			printf(
				"check false: trpl[%ld]{%ld, %ld}, trpl[%ld]{%ld, %ld}\n",
				idx, trpl[idx].user_id, trpl[idx].item_id,
				idx+1, trpl[ idx+1].user_id, trpl[ idx+1].item_id
			);
			continue;
		}
	}
}


void check_sorted_keyvals(keyval_t *keyval, int64_t num_keyval, 
	int (*func_cmp)(const void*, const void*)) {
	for (int64_t idx = 0; idx < num_keyval -1; idx ++) {
		if (func_cmp(keyval+idx, keyval+idx+1) > 0) {
			printf(
				"check false: keyval[%ld]{%ld, %f}, keyval[%ld]{%ld, %f}\n",
				idx, keyval[idx].key, keyval[idx].val,
				idx+1, keyval[idx+1].key, keyval[idx+1].val
			);
			continue;
		}
	}
}

#if THREAD_NUM > 1
static void* pthread_qsort(void *param) {
	pthread_qsort_arg_t *arg = (pthread_qsort_arg_t *) param;
	arg->res = (triple_t *) alloc_memory(
		sizeof(triple_t) * arg->count, arg->numa_node_idx 
	);

	memcpy(
		arg->res, arg->trpl+arg->start, 
		sizeof(triple_t) * arg->count
	);

	qsort(
		arg->res,
		arg->count,
		sizeof(triple_t), 
		compare_triples
	);

	// check_sorted_triples(arg->res, arg->count);

	return NULL;
}


static void pthread_merge(void *param) {
	pthread_merge_arg_t *arg = (pthread_merge_arg_t *) param;
	int64_t index[2] = {0, 0};

	for (int64_t idx = 0; idx < arg->count0 + arg->count1; idx ++) {
		short flag = 0;

		if (index[1] >= arg->count1) {
			flag = 0;
		} else if (index[0] >= arg->count0) {
			flag = 1;
		} else if (compare_triples(arg->source0 + index[0], arg->source1 + index[1])<=0) {
			flag = 0;
		} else {
			flag = 1;
		}

		if (!flag) {
			arg->res[idx] = arg->source0[index[0]];
			index[0] ++;
		} else {
			arg->res[idx] = arg->source1[index[1]];
			index[1] ++;
		}
	}

	dealloc_memory(arg->source0, sizeof(triple_t)*arg->count0);
	dealloc_memory(arg->source1, sizeof(triple_t)*arg->count1);

	// check_sorted_triples(arg->res, arg->count0+arg->count1);
	
	// return NULL;
}

void p_qsortmerge(triple_t *trpl, my_cnt_t num_nnz, int numa_node_idx) {
	if (num_nnz < 1 << 10) {
		qsort(trpl, num_nnz, sizeof(triple_t), compare_triples);
	} else {
		pthread_t threadpool[THREAD_NUM];
		pthread_attr_t attr;
		pthread_attr_init(&attr);
		cpu_set_t set;
		int core_idx;
		int rv;

		pthread_qsort_arg_t qsort_args[THREAD_NUM];
		memset(qsort_args, 0, sizeof(pthread_qsort_arg_t)*THREAD_NUM);

		int64_t num_per_thr = num_nnz / THREAD_NUM;

		for (int idx = 0; idx < THREAD_NUM; idx++) {
			core_idx = get_cpu_id(idx);
			CPU_ZERO(&set);
			CPU_SET(core_idx, &set);
			pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);
			qsort_args[idx]._tid = idx;
			qsort_args[idx].numa_node_idx = numa_node_idx;
			qsort_args[idx].start = idx * num_per_thr;
			qsort_args[idx].count = idx == THREAD_NUM -1 ? 
				(num_nnz - qsort_args[idx].start) : num_per_thr;
			qsort_args[idx].trpl = trpl;
			rv = pthread_create(&threadpool[idx], &attr, pthread_qsort, qsort_args+idx);
			if (rv){
				printf("ERROR; return code from pthread_create() is %d\n", rv);
				exit(EXIT_FAILURE);
			}
		}
		for (int idx = 0; idx < THREAD_NUM; idx++) {
			rv = pthread_join(threadpool[idx], NULL);
			if (rv){
				printf("ERROR; return code from pthread_join() is %d\n", rv);
				exit(EXIT_FAILURE);
			}
		}

		merge_part_t merge_parts[THREAD_NUM];
		memset(merge_parts, 0, sizeof(merge_part_t)*THREAD_NUM);
		for (int idx = 0; idx < THREAD_NUM; idx++) {
			merge_parts[idx].count = qsort_args[idx].count;
			merge_parts[idx].res = qsort_args[idx].res;
		}

		int merge_num = THREAD_NUM;
		// assert(merge_num >=2);
		while (merge_num > 2){
			pthread_merge_arg_t merge_args[merge_num/2];
			memset(merge_args, 0, sizeof(pthread_merge_arg_t)*(merge_num/2));
			
			for (int idx = 0; idx < merge_num/2; idx ++) {
				core_idx = get_cpu_id(idx);
				CPU_ZERO(&set);
				CPU_SET(core_idx, &set);
				pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);
				merge_args[idx]._tid = idx;
				merge_args[idx].count0 = merge_parts[2 * idx].count;
				merge_args[idx].count1 = merge_parts[2 * idx + 1].count;
				merge_args[idx].source0 = merge_parts[2 * idx].res;
				merge_args[idx].source1 = merge_parts[2 * idx + 1].res;
				merge_args[idx].res = alloc_memory(
					sizeof(triple_t)*(merge_args[idx].count0+merge_args[idx].count1),
					numa_node_idx
				);
				rv = pthread_create(&threadpool[idx], &attr, (void *(*)(void *)) (uintptr_t) pthread_merge, merge_args+idx);
				if (rv){
					printf("ERROR; return code from pthread_create() is %d\n", rv);
					exit(EXIT_FAILURE);
				}
			}

			for (int idx = 0; idx <  merge_num/2; idx++) {
				rv = pthread_join(threadpool[idx], NULL);
				if (rv){
					printf("ERROR; return code from pthread_join() is %d\n", rv);
					exit(EXIT_FAILURE);
				}
				merge_parts[idx].count = merge_args[idx].count0 + merge_args[idx].count1;
				merge_parts[idx].res = merge_args[idx].res;
			}

			if (merge_num % 2 == 1){
				merge_parts[ (merge_num/2) ] = merge_parts[merge_num - 1];
				merge_num = merge_num / 2 + 1;
			} else{
				merge_num /= 2;
			}
		}

		pthread_merge_arg_t last_merge_arg = {
			._tid = 0,
			.count0 = merge_parts[0].count,
			.count1 = merge_parts[1].count,
			.source0 = merge_parts[0].res,
			.source1 = merge_parts[1].res,
			.res = trpl
		};

		// assert(merge_parts[0].count + merge_parts[1].count == num_nnz);
		pthread_merge(&last_merge_arg);

		// dealloc_memory(last_merge_arg.source0, sizeof(triple_t)*last_merge_arg.count0);
		// dealloc_memory(last_merge_arg.source1, sizeof(triple_t)*last_merge_arg.count1);

		pthread_attr_destroy(&attr);

		// check_sorted_triples(trpl, num_nnz);
	}
}


static void* pthread_qsort_keyval(void *param) {
	pthread_qsort_keyval_arg_t *arg = (pthread_qsort_keyval_arg_t *) param;
	arg->res = (keyval_t *) alloc2_memory(
		sizeof(keyval_t) * arg->count, arg->numa_cfg 
	);

	memcpy(
		arg->res, arg->keyval+arg->start, 
		sizeof(keyval_t) * arg->count
	);

	qsort(
		arg->res,
		arg->count,
		sizeof(keyval_t), 
		arg->func_cmp
	);

	// check_sorted_keyvals(arg->res, arg->count, arg->func_cmp);

	return NULL;
}


static void pthread_merge_keyval(void *param) {
	pthread_merge_keyval_arg_t *arg = (pthread_merge_keyval_arg_t *) param;
	int64_t index[2] = {0, 0};

	for (int64_t idx = 0; idx < arg->count0 + arg->count1; idx ++) {
		short flag = 0;

		if (index[1] >= arg->count1) {
			flag = 0;
		} else if (index[0] >= arg->count0) {
			flag = 1;
		} else if (arg->func_cmp(arg->source0 + index[0], arg->source1 + index[1])<=0) {
			flag = 0;
		} else {
			flag = 1;
		}

		if (!flag) {
			arg->res[idx] = arg->source0[index[0]];
			index[0] ++;
		} else {
			arg->res[idx] = arg->source1[index[1]];
			index[1] ++;
		}
	}

	dealloc_memory(arg->source0, sizeof(keyval_t)*arg->count0);
	dealloc_memory(arg->source1, sizeof(keyval_t)*arg->count1);

	// check_sorted_keyvals(arg->res, arg->count0+arg->count1, arg->func_cmp);

	// return NULL;
}


void p_qsortmerge_keyval(keyval_t *keyval, my_cnt_t num_keyval,
	int (*func_cmp)(const void*, const void*), numa_cfg_t *numa_cfg) {
	if (num_keyval < 1 << 10) {
		qsort(keyval, num_keyval, sizeof(keyval_t), func_cmp);
	} else {
		pthread_t threadpool[THREAD_NUM];
		pthread_attr_t attr;
		pthread_attr_init(&attr);
		cpu_set_t set;
		int core_idx;
		int rv;

		pthread_qsort_keyval_arg_t qsort_args[THREAD_NUM];
		memset(qsort_args, 0, sizeof(pthread_qsort_keyval_arg_t)*THREAD_NUM);

		int64_t num_per_thr = num_keyval / THREAD_NUM;

		for (int idx = 0; idx < THREAD_NUM; idx ++) {
			core_idx = get_cpu_id(idx);
			CPU_ZERO(&set);
			CPU_SET(core_idx, &set);
			pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);
			qsort_args[idx]._tid = idx;
			qsort_args[idx].numa_cfg = numa_cfg;
			qsort_args[idx].start = idx * num_per_thr;
			qsort_args[idx].count = idx == THREAD_NUM -1 ? 
				(num_keyval - qsort_args[idx].start) : num_per_thr;
			qsort_args[idx].keyval = keyval;
			qsort_args[idx].func_cmp = func_cmp;
			rv = pthread_create(&threadpool[idx], &attr, pthread_qsort_keyval, qsort_args+idx);
			if (rv) {
				printf("ERROR; return code from pthread_create() is %d\n", rv);
				exit(EXIT_FAILURE);
			}
		}

		for (int idx = 0; idx < THREAD_NUM; idx ++) {
			rv = pthread_join(threadpool[idx], NULL);
			if (rv) {
				printf("ERROR; return code from pthread_join() is %d\n", rv);
				exit(EXIT_FAILURE);
			}
		}


		merge_part_keyval_t merge_parts[THREAD_NUM];
		memset(merge_parts, 0, sizeof(merge_part_keyval_t)*(THREAD_NUM));

		for (int idx = 0; idx < THREAD_NUM; idx ++) {
			merge_parts[idx].count = qsort_args[idx].count;
			merge_parts[idx].res = qsort_args[idx].res;
		}

		int merge_num = THREAD_NUM;
		while (merge_num > 2) {
			pthread_merge_keyval_arg_t merge_args[merge_num/2];
			memset(merge_args, 0, sizeof(pthread_merge_keyval_arg_t)*(merge_num/2));

			for (int idx = 0; idx < merge_num/2; idx ++) {
				core_idx = get_cpu_id(idx);
				CPU_ZERO(&set);
				CPU_SET(core_idx, &set);
				pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);
				merge_args[idx]._tid = idx;
				merge_args[idx].count0 = merge_parts[2 * idx].count;
				merge_args[idx].count1 = merge_parts[2 * idx + 1].count;
				merge_args[idx].source0 = merge_parts[2 * idx].res;
				merge_args[idx].source1 = merge_parts[2 * idx + 1].res;
				merge_args[idx].res = alloc2_memory(
					sizeof(keyval_t)*(merge_args[idx].count0+merge_args[idx].count1),
					numa_cfg
				);
				merge_args[idx].func_cmp = func_cmp;
				rv = pthread_create(&threadpool[idx], &attr, (void *(*)(void *)) (uintptr_t) pthread_merge_keyval, merge_args+idx);
				if (rv) {
					printf("ERROR; return code from pthread_create() is %d\n", rv);
					exit(EXIT_FAILURE);
				}
			}

			for (int idx = 0; idx <  merge_num /2; idx++) {
				rv = pthread_join(threadpool[idx], NULL);
				if (rv){
					printf("ERROR; return code from pthread_join() is %d\n", rv);
					exit(EXIT_FAILURE);
				}
				merge_parts[idx].count = merge_args[idx].count0 + merge_args[idx].count1;
				merge_parts[idx].res = merge_args[idx].res;
			}

			if (merge_num % 2 == 1){
				merge_parts[ (merge_num/2) ] = merge_parts[merge_num - 1];
				merge_num = merge_num / 2 + 1;
			} else{
				merge_num /= 2;
			}
		}

		pthread_merge_keyval_arg_t last_merge_arg = {
			._tid = 0,
			.count0 = merge_parts[0].count,
			.count1 = merge_parts[1].count,
			.source0 = merge_parts[0].res,
			.source1 = merge_parts[1].res,
			.res = keyval,
			.func_cmp = func_cmp
		};

		// assert(merge_parts[0].count + merge_parts[1].count == num_keyval);

		// dealloc_memory(last_merge_arg.source0, sizeof(keyval_t)*last_merge_arg.count0);
		// dealloc_memory(last_merge_arg.source1, sizeof(keyval_t)*last_merge_arg.count1);

		pthread_merge_keyval(&last_merge_arg);
		pthread_attr_destroy(&attr);

		// check_sorted_keyvals(keyval, num_keyval, func_cmp);

	}
}
#endif /* THREAD_NUM > 1 */


bool mtx2trpl(triple_t **trpl, my_cnt_t* num_row, my_cnt_t* num_col,
	my_cnt_t* num_nnz, const char *filename, const int numa_node_idx) {
	bool isSymmetric = false;
	my_cnt_t tmp_num_nnz = 0, tmp_user_id = 0, tmp_item_id = 0;
	my_rat_t tmp_rating = 0.0, tmp_rating_double = 0.0;

	triple_t* tmp_trpl = NULL;

	// Open file using open system call with O_DIRECT
	int fd = open(filename, O_RDONLY | O_DIRECT);
	if (fd == -1) {
		fprintf(stderr, "Error opening file with O_DIRECT: %s\n", filename);
		exit(EXIT_FAILURE);
	}

	// Get file size and alignment
	struct stat sb;
	if (fstat(fd, &sb) == -1) {
		perror("Error getting file size");
		close(fd);
		exit(EXIT_FAILURE);
	}
	size_t alignment = sb.st_blksize; // Typically 4096 bytes
	size_t buffer_size = (sb.st_size + alignment - 1) & ~(alignment - 1);
	void *aligned_buffer;
	if (posix_memalign(&aligned_buffer, alignment, buffer_size)) {
		perror("Error allocating aligned memory");
		close(fd);
		exit(EXIT_FAILURE);
	}

	ssize_t bytes_read = 0;
	char *current_ptr = (char *)aligned_buffer;
	bool header_read = false;
	my_cnt_t cnt = 0;

	char leftover[4096];  // Temporary buffer for an incomplete line
	size_t leftover_len = 0;

	while ((bytes_read = read(fd, aligned_buffer, buffer_size)) > 0) {

		current_ptr = (char *)aligned_buffer;

		while (current_ptr < (char *)aligned_buffer + bytes_read) {
			char *line = current_ptr;
			// Find the end of the current line
			char *newline = memchr(line, '\n', (char *)aligned_buffer + bytes_read - current_ptr);
			if (!newline) {

				leftover_len = (char *)aligned_buffer + bytes_read - line;
				memset(leftover, 0, 4096);
				memcpy(leftover, line, leftover_len);
				break;  // Partial line, handle on the next read
			}

			*newline = '\0';  // Null-terminate the line
			current_ptr = newline + 1;

			// Skip comment lines starting with '%'
			if (line[0] == '%') {
				if (strstr(line, "symmetric")) {
					isSymmetric = true;
				}
				continue;
			}

			// Parse the header line
			if (!header_read) {
				if (sscanf(line, "%ld %ld %ld", num_row, num_col, &tmp_num_nnz) != 3) {
					fprintf(stderr, "Error reading header line: %s\n", line);
					free(aligned_buffer);
					close(fd);
					exit(EXIT_FAILURE);
				}
				header_read = true;

				if (isSymmetric) {
					tmp_num_nnz = tmp_num_nnz * 2;
				}

				tmp_trpl = (triple_t *)alloc_dram(sizeof(triple_t) * tmp_num_nnz, numa_node_idx);
				memset(tmp_trpl, 0, sizeof(triple_t) * tmp_num_nnz);
			} else {
				// Parse non-header lines
				int read_columns = sscanf(line, "%ld %ld %lf",
						&tmp_user_id, &tmp_item_id, &tmp_rating_double);

				// try to deal with leftover
				if (leftover_len > 0) { 
					memcpy(leftover+leftover_len, line, current_ptr-line);
					read_columns = sscanf(leftover, "%ld %ld %lf",
						&tmp_user_id, &tmp_item_id, &tmp_rating_double);
					leftover_len = 0;
					memset(leftover, 0, 4096);
				}

				if (read_columns == 2) {
					tmp_rating_double = 1.0;
				} else if (read_columns != 3) {
					fprintf(stderr, "Error reading line: %s\n", line);
					free(aligned_buffer);
					close(fd);
					exit(EXIT_FAILURE);
				}


				if (tmp_rating_double == 0.0) {
					continue;
				}

				tmp_rating = 1.0;
				tmp_trpl[cnt].user_id = tmp_user_id - 1;
				tmp_trpl[cnt].item_id = tmp_item_id - 1;
				tmp_trpl[cnt].rating = tmp_rating;
				cnt++;

				if (isSymmetric) {
					if (tmp_trpl[cnt - 1].user_id != tmp_trpl[cnt - 1].item_id) {
						tmp_trpl[cnt].user_id = tmp_trpl[cnt - 1].item_id;
						tmp_trpl[cnt].item_id = tmp_trpl[cnt - 1].user_id;
						tmp_trpl[cnt].rating = tmp_trpl[cnt - 1].rating;
						cnt++;
					}
				}
			}
		}
	}

	if (bytes_read == -1) {
		perror("Error reading file with O_DIRECT");
		free(aligned_buffer);
		close(fd);
		exit(EXIT_FAILURE);
	}

	free(aligned_buffer);
	close(fd);

	*trpl = (triple_t *) alloc_memory(
		sizeof(triple_t)*cnt, numa_node_idx
	);

	memset(*trpl, 0, sizeof(triple_t)*cnt);

	*num_nnz = cnt;

	memcpy(*trpl, tmp_trpl, sizeof(triple_t)*cnt);
	dealloc_memory(tmp_trpl, sizeof(triple_t)*tmp_num_nnz);

#if THREAD_NUM > 1
	p_qsortmerge(*trpl, *num_nnz, numa_node_idx);
#else	/* THREAD_NUM > 1 */
	qsort(*trpl, *num_nnz, sizeof(triple_t), compare_triples);
#endif	/* THREAD_NUM > 1 */

	return isSymmetric;
}


void trpl2mtx(const char *filename, triple_t *trpl, my_cnt_t num_row, 
	my_cnt_t num_col, my_cnt_t num_nnz, bool isSymmetric) {
	FILE *file = fopen(filename, "w");
	if (!file) {
		perror("Error opening file");
		exit(EXIT_FAILURE);
	}

	// Write the MTX header
	if (isSymmetric) {
		fprintf(file, "%%MatrixMarket matrix coordinate real symmetric\n");
		my_cnt_t actual_cnt = 0;
		for (my_cnt_t idx = 0; idx < num_nnz; idx++) {
			if (trpl[idx].user_id <= trpl[idx].item_id) {
				actual_cnt ++;
			}
		}
		fprintf(file, "%ld %ld %ld\n", num_row, num_col, actual_cnt);
		for (my_cnt_t idx = 0; idx < num_nnz; idx++) {
			if (trpl[idx].user_id <= trpl[idx].item_id) {
				fprintf(file, "%ld %ld %f\n",
					trpl[idx].user_id + 1, trpl[idx].item_id + 1, 1.0
				);
			}
		}
	} else {
		fprintf(file, "%%MatrixMarket matrix coordinate real general\n");
		fprintf(file, "%ld %ld %ld\n", num_row, num_col, num_nnz);
		for (my_cnt_t idx = 0; idx < num_nnz; idx++) {
			fprintf(file, "%ld %ld %f\n",
				trpl[idx].user_id + 1, trpl[idx].item_id + 1, 1.0
			);
		}

	}

	fclose(file);
}


void trpl2csr(csr_t *csr, triple_t *trpl, const my_cnt_t num_row, 
	const my_cnt_t num_col, const my_cnt_t num_nnz, numa_cfg_t* numa_cfg) {
	memset(csr, 0, sizeof(csr_t));
	csr->num_row = num_row;
	csr->num_col = num_col;
	csr->num_nnz = num_nnz;

	if (IS_POWER_OF_TWO(*(numa_cfg->numa_mask->maskp))) {
		csr->row_ptrs = (my_cnt_t *) alloc_dram(
			sizeof(my_cnt_t) * (num_row+1), numa_cfg->node_idx
		);
		csr->cols = (my_cnt_t *) alloc_dram(
			sizeof(my_cnt_t) * num_nnz, numa_cfg->node_idx
		);
		csr->rats = (my_rat_t *) alloc_dram( 
			sizeof(my_rat_t) * num_nnz, numa_cfg->node_idx
		);
	} else {
		csr->row_ptrs = (my_cnt_t *) alloc2_memory(
			sizeof(my_cnt_t) * (num_row+1), numa_cfg
		);
		csr->cols = (my_cnt_t *) alloc2_memory(
			sizeof(my_cnt_t) * num_nnz, numa_cfg
		);
		csr->rats = (my_rat_t *) alloc2_memory( 
			sizeof(my_rat_t) * num_nnz, numa_cfg
		);
	}

	memset(csr->row_ptrs, 0, sizeof(my_cnt_t)*(num_row+1));
	memset(csr->cols, 0, sizeof(my_cnt_t)*num_nnz);
	memset(csr->rats, 0, sizeof(my_rat_t)*num_nnz);

	for (my_cnt_t idx = 0; idx < num_nnz; idx ++) {
		csr->row_ptrs[ trpl[idx].user_id+1 ] ++;
		csr->cols[idx] = trpl[idx].item_id;
		csr->rats[idx] = trpl[idx].rating;
	}

	for (my_cnt_t idx = 1; idx <= num_row; idx ++ ) {
		csr->row_ptrs[idx] += csr->row_ptrs[idx-1];
	}
}


void csr2trpl(triple_t **trpl, my_cnt_t* num_row, my_cnt_t* num_col, 
	my_cnt_t* num_nnz, csr_t *csr, const int numa_node_idx) {
	*num_row = csr->num_row;
	*num_col = csr->num_col;
	*num_nnz = csr->num_nnz;
	*trpl = (triple_t *) alloc_memory(
		sizeof(triple_t) * (*num_nnz), numa_node_idx
	);
	memset(*trpl, 0, sizeof(triple_t) * (*num_nnz));

	my_cnt_t kdx = 0;
	for (my_cnt_t idx = 0; idx < *num_row; idx ++) {
		for (my_cnt_t jdx = csr->row_ptrs[idx]; jdx < csr->row_ptrs[idx+1]; jdx ++) {
			(*trpl)[kdx].user_id = idx;
			(*trpl)[kdx].item_id = csr->cols[jdx];
			(*trpl)[kdx].rating = csr->rats[jdx];
		}
	}

	// qsort(*trpl, *num_nnz, sizeof(triple_t), compare_triples);
}


void mtx2csr(csr_t *csr, const char *filename, numa_cfg_t* numa_cfg) {
	triple_t* trpl = NULL;
	my_cnt_t num_row, num_col, num_nnz;
	mtx2trpl(&trpl, &num_row, &num_col, &num_nnz, filename, numa_cfg->node_idx);
	trpl2csr(csr, trpl, num_row, num_col, num_nnz, numa_cfg);
	dealloc_memory(trpl, sizeof(triple_t) * num_nnz);
}





void raw2csr(csr_t *csr, const char *filename, numa_cfg_t* numa_cfg) {
	mtx2csr(csr, filename, numa_cfg);
}


void trpl_ones(triple_t* trpl, my_cnt_t num_nnz) {
	for (my_cnt_t idx = 0; idx < num_nnz; idx ++) {
		trpl[idx].rating = 1.0;
	}
}


void trpl_transpose(triple_t* trpl, my_cnt_t *num_row, 
	my_cnt_t *num_col, my_cnt_t num_nnz) {
	for (my_cnt_t idx = 0; idx < num_nnz; idx ++) {
		my_cnt_t tmp = trpl[idx].user_id;
		trpl[idx].user_id = trpl[idx].item_id;
		trpl[idx].item_id = tmp;
	}
	qsort(trpl, num_nnz, sizeof(triple_t), compare_triples);
	my_cnt_t tmp = *num_row;
	*num_row = *num_col;
	*num_col = tmp;
}




// Tailored specifically for PaRMAT-synthesized datasets
void txt2trpl(triple_t **trpl, my_cnt_t* num_row, my_cnt_t* num_col,
	my_cnt_t* num_nnz, const char *filename, const int numa_node_idx) {

	int scale = 0, edge_factor = 0;
	// Scan for the pattern 'er_' followed by two numbers separated by '_'
	if (sscanf(filename, RMAT_DIR"/%*[^_]_%d_%d.txt", &scale, &edge_factor) == 2) {
		printf("Extracted numbers: %d and %d\n", scale, edge_factor);
	} else {
		printf("Failed to extract numbers.\n");
	}
	*num_row = *num_col = (1 << scale);

	my_cnt_t tmp_user_id = 0, tmp_item_id = 0,\
		tmp_num_nnz = count_lines_num(filename) * 2;
	assert ( tmp_num_nnz == *num_row * edge_factor * 2 );

	my_rat_t tmp_rating = 0.0, tmp_rating_double = 0.0;

	triple_t* tmp_trpl = (triple_t*) alloc_dram(
		sizeof(triple_t) * tmp_num_nnz, numa_node_idx
	);
	memset(tmp_trpl, 0, sizeof(triple_t)*tmp_num_nnz);


	FILE* file = fopen(filename, "r");
	if (file == NULL) {
		fprintf(stderr, "Error opening file, %s\n", filename);
		exit(EXIT_FAILURE);
	}

	char line[CHAR_BUFFER_LEN];

	my_cnt_t cnt = 0;

	while (fgets(line, sizeof(line), file) != NULL) {
		if (cnt < tmp_num_nnz) {
			int read_columns = sscanf(line, "%ld %ld %lf", 
				&tmp_user_id, &tmp_item_id, &tmp_rating_double
			);

			if (read_columns == 2) {
				tmp_rating_double = 1.0;
			} else if (read_columns != 3) {
				fprintf(stderr, "Error reading line: %s\n", line);
				exit(EXIT_FAILURE);
			}

			if (tmp_rating_double == 0.0) {
				continue;
			}

			// tmp_rating = (float) tmp_rating_double;
			tmp_rating = 1.0;

			// Indexing starts from 0, so no need to subtract or adjust
			tmp_trpl[cnt].user_id = tmp_user_id;
			tmp_trpl[cnt].item_id = tmp_item_id;
			tmp_trpl[cnt].rating = tmp_rating;

			cnt++;


			assert (tmp_trpl[cnt-1].user_id != tmp_trpl[cnt-1].item_id);
			tmp_trpl[cnt].user_id = tmp_trpl[cnt-1].item_id;
			tmp_trpl[cnt].item_id = tmp_trpl[cnt-1].user_id;
			tmp_trpl[cnt].rating = tmp_trpl[cnt-1].rating;
			cnt ++;

		}
	}
	fclose(file);

	*trpl = (triple_t *) alloc_memory(
		sizeof(triple_t)*cnt, numa_node_idx
	);

	memset(*trpl, 0, sizeof(triple_t)*cnt);

	*num_nnz = cnt;

	memcpy(*trpl, tmp_trpl, sizeof(triple_t)*cnt);
	dealloc_memory(tmp_trpl, sizeof(triple_t)*tmp_num_nnz);

#if THREAD_NUM > 1
	p_qsortmerge(*trpl, *num_nnz, numa_node_idx);
#else	/* THREAD_NUM > 1 */
	qsort(*trpl, *num_nnz, sizeof(triple_t), compare_triples);
#endif	/* THREAD_NUM > 1 */
}


// Tailored specifically for PaRMAT-synthesized datasets
void txt_transpose(
	const char *filename_dest_trs, const char *filename_dest_org,
	const char *filename_src, const int numa_node_idx) {
	triple_t* trpl = NULL;
	my_cnt_t num_row, num_col, num_nnz;
	txt2trpl(&trpl, &num_row, &num_col, &num_nnz, filename_src, numa_node_idx);
	trpl_ones(trpl, num_nnz);
	trpl2mtx(filename_dest_org, trpl, num_row, num_col, num_nnz, true);
	trpl_transpose(trpl, &num_row, &num_col, num_nnz);
	trpl2mtx(filename_dest_trs, trpl, num_row, num_col, num_nnz, true);
	dealloc_triv_dram(trpl, sizeof(triple_t) * num_nnz);
}

// TODO: Parallelize this function

static void *pthread_scan(void *param) {
	pthread_scan_arg_t *arg = (pthread_scan_arg_t *) param;
	int _tid = arg->_tid;
	my_cnt_t num = arg->num;
	my_cnt_t *in = arg->in;
	my_cnt_t *out = arg->out;
	my_cnt_t *partial_sum = arg->partial_sum;

	my_cnt_t each_n = num / THREAD_NUM;
	my_cnt_t start_idx = _tid * each_n;
	my_cnt_t end_idx = (_tid == THREAD_NUM-1) ? num : start_idx + each_n;

	out[start_idx] = 0;
	for (my_cnt_t idx = start_idx; idx < end_idx - 1; idx++) {
		out[idx+1] = out[idx]+in[idx];
	}

	if (_tid < THREAD_NUM-1) {
		partial_sum[_tid] = out[end_idx-1] + in[end_idx-1];
	}

	int rv;
	BARRIER_ARRIVE(arg->barrier, rv);

	my_cnt_t offset = 0;
	for (int idx = 0; idx < _tid; ++idx) {
		offset += partial_sum[idx];
	}

	for (my_cnt_t idx = start_idx; idx < end_idx; ++ idx) {
		out[idx] += offset;
	}

	return NULL;
}


void scan(my_cnt_t *out, my_cnt_t *in, my_cnt_t num) {
	if ( num < (1 << 17) ) {							// single-thread execution
		out[0] = 0;
		for (my_cnt_t idx = 0; idx < num-1; idx ++) {
			out[idx+1] = out[idx] + in[idx];
		}
	} else {											// multi-thread execution
		pthread_t threadpool[THREAD_NUM];
		pthread_attr_t attr;
		pthread_attr_init(&attr);
		cpu_set_t set;
		int core_idx;
		int rv;

		pthread_scan_arg_t args[THREAD_NUM];
		memset(args, 0, sizeof(pthread_scan_arg_t)*THREAD_NUM);
		
		pthread_barrier_t barrier;
		rv = pthread_barrier_init(&barrier, NULL, THREAD_NUM);

		my_cnt_t partial_sum[THREAD_NUM];

		for (int idx = 0; idx < THREAD_NUM; idx ++) {
			core_idx = get_cpu_id(idx);
			CPU_ZERO(&set);
			CPU_SET(core_idx, &set);
			pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);

			args[idx]._tid = idx;
			args[idx].num = num;
			args[idx].in = in;
			args[idx].out = out;
			args[idx].partial_sum = partial_sum;
			args[idx].barrier = &barrier;

			rv = pthread_create(&threadpool[idx], &attr, pthread_scan, args+idx);
			if (rv){
				printf("ERROR; return code from pthread_create() is %d\n", rv);
				exit(EXIT_FAILURE);
			}

		}

		for (int idx = 0; idx < THREAD_NUM; idx ++) {
			rv = pthread_join(threadpool[idx], NULL);
			if (rv) {
				printf("ERROR; return code from pthread_join() is %d\n", rv);
				exit(EXIT_FAILURE);
			}
		}

		pthread_barrier_destroy(&barrier);
		pthread_attr_destroy(&attr);
	}
}





my_cnt_t* lower_bound(my_cnt_t *start, my_cnt_t *end, my_cnt_t val) {
	my_cnt_t* mid;
	while (start < end) {
		mid = start + (end - start) / 2;
		if (*mid < val) {
			start = mid + 1;
		} else {
			end = mid;
		}
	}
	return start;
}


my_rat_t func_mul(my_rat_t a, my_rat_t b) {
	return a * b;
}

static void* pthread_get_intprod(void *param) {
	pthread_intprod_arg_t *arg = (pthread_intprod_arg_t *)param;
	int _tid = arg->_tid;
	my_cnt_t *row_nz = arg->row_nz;
	csr_t *csr_a = arg->csr_a;
	csr_t *csr_b = arg->csr_b;

	my_cnt_t intprod = 0;
	my_cnt_t nz_per_row = 0;

	my_cnt_t start_idx = _tid * (csr_a->num_row/THREAD_NUM);
	my_cnt_t end_idx = (_tid == THREAD_NUM-1) ? csr_a->num_row : 
		(_tid+1) * (csr_a->num_row/THREAD_NUM);

	for (my_cnt_t idx = start_idx; idx < end_idx; idx ++) {
		nz_per_row = 0;
		for (my_cnt_t jdx = csr_a->row_ptrs[idx]; jdx < csr_a->row_ptrs[idx+1]; jdx ++) {
			nz_per_row += csr_b->row_ptrs[ csr_a->cols[jdx]+1 ] - csr_b->row_ptrs[ csr_a->cols[jdx] ];
		}
		row_nz[idx] = nz_per_row;
		intprod += nz_per_row;
	}

	arg->intprod = intprod;
	return NULL;
}

my_cnt_t get_intprod(csr_t* csr_a, csr_t* csr_b, my_cnt_t *row_nz) {
	my_cnt_t total_intprod = 0;

	pthread_t threadpool[THREAD_NUM];
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	cpu_set_t set;
	int core_idx;
	int rv;
	pthread_intprod_arg_t args[THREAD_NUM];
	memset(args, 0, sizeof(pthread_intprod_arg_t)*THREAD_NUM);

	for (int idx = 0; idx < THREAD_NUM; idx ++) {
		core_idx = get_cpu_id(idx);
		CPU_ZERO(&set);
		CPU_SET(core_idx, &set);
		pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);
		args[idx]._tid = idx;
		args[idx].intprod = 0;
		args[idx].row_nz = row_nz;
		args[idx].csr_a = csr_a;
		args[idx].csr_b = csr_b;
		rv = pthread_create(&threadpool[idx], &attr, pthread_get_intprod, args+idx);
		if (rv){
			printf("ERROR; return code from pthread_create() is %d\n", rv);
			exit(EXIT_FAILURE);
		}
	}

	for (int idx = 0; idx < THREAD_NUM; idx ++) {
		rv = pthread_join(threadpool[idx], NULL);
		if (rv){
			printf("ERROR; return code from pthread_join() is %d\n", rv);
			exit(EXIT_FAILURE);
		}
#ifdef IN_DEBUG
		// printf("\t%d_intprod: %ld\n", idx, args[idx].intprod);
		// if (idx == THREAD_NUM-1) {
		// 	printf("\n");
		// }
#endif /* IN_DEBUG */
		total_intprod += args[idx].intprod;
	}

	pthread_attr_destroy(&attr);

	return total_intprod;
}







static void* pthread_flp_populate(void *param) {
	pthread_intprod_arg_t *arg = (pthread_intprod_arg_t *)param;
	int _tid = arg->_tid;
	csr_t *csr_a = arg->csr_a;
	csr_t *csr_b = arg->csr_b;

	my_cnt_t *row_flp = arg->row_flp;

	my_cnt_t flp_per_row = 0;

	my_cnt_t start_idx = _tid * (csr_a->num_row/THREAD_NUM);
	my_cnt_t end_idx = (_tid == THREAD_NUM-1) ? csr_a->num_row : 
		(_tid+1) * (csr_a->num_row/THREAD_NUM);

	for (my_cnt_t idx = start_idx; idx < end_idx; idx ++) {
		flp_per_row = 0;
		for (my_cnt_t jdx = csr_a->row_ptrs[idx]; jdx < csr_a->row_ptrs[idx+1]; jdx ++) {
			flp_per_row += csr_b->row_ptrs[ csr_a->cols[jdx]+1 ] - csr_b->row_ptrs[ csr_a->cols[jdx] ];
		}
		row_flp[idx] = flp_per_row;
	}

	return NULL;
}


void flp_populate(csr_t* csr_a, csr_t* csr_b, my_cnt_t *row_flp) {
	pthread_t threadpool[THREAD_NUM];
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	cpu_set_t set;
	int core_idx;
	int rv;
	pthread_intprod_arg_t args[THREAD_NUM];
	memset(args, 0, sizeof(pthread_intprod_arg_t)*THREAD_NUM);

	for (int idx = 0; idx < THREAD_NUM; idx ++) {
		core_idx = get_cpu_id(idx);
		CPU_ZERO(&set);
		CPU_SET(core_idx, &set);
		pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);
		args[idx]._tid = idx;
		args[idx].csr_a = csr_a;
		args[idx].csr_b = csr_b;

		args[idx].row_flp = row_flp;
		rv = pthread_create(&threadpool[idx], &attr, pthread_flp_populate, args+idx);
		if (rv){
			printf("ERROR; return code from pthread_create() is %d\n", rv);
			exit(EXIT_FAILURE);
		}
	}

	for (int idx = 0; idx < THREAD_NUM; idx ++) {
		rv = pthread_join(threadpool[idx], NULL);
		if (rv){
			printf("ERROR; return code from pthread_join() is %d\n", rv);
			exit(EXIT_FAILURE);
		}
	}

	pthread_attr_destroy(&attr);
}



static void* pthread_set_rows_offset(void *param) {
	pthread_set_rows_offset_arg_t* arg = (pthread_set_rows_offset_arg_t*) param;
	int _tid = arg->_tid;
	my_cnt_t num_row = arg->num_row;
	my_cnt_t average_intprod = arg->average_intprod;
	my_cnt_t* rows_offset = arg->rows_offset;
	my_cnt_t* row_start_idcs = arg->row_start_idcs;

	rows_offset[_tid+1] = lower_bound(
		row_start_idcs, 
		row_start_idcs + num_row + 1, 
		average_intprod * (_tid + 1)
	) - row_start_idcs;

	return NULL;
}


void set_rows_offset(my_cnt_t num_row, my_cnt_t total_intprod, 
	my_cnt_t* rows_offset, my_cnt_t* row_start_idcs) {
	pthread_t threadpool[THREAD_NUM];
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	cpu_set_t set;
	int core_idx;
	int rv;
	pthread_set_rows_offset_arg_t args[THREAD_NUM];
	memset(args, 0, sizeof(pthread_set_rows_offset_arg_t)*THREAD_NUM);

	my_cnt_t average_intprod = (total_intprod + THREAD_NUM - 1) / THREAD_NUM;

	for (int idx = 0; idx < THREAD_NUM; idx ++) {
		core_idx = get_cpu_id(idx);
		CPU_ZERO(&set);
		CPU_SET(core_idx, &set);
		pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);
		args[idx]._tid = idx;
		args[idx].num_row = num_row;
		args[idx].average_intprod = average_intprod;
		args[idx].rows_offset = rows_offset;
		args[idx].row_start_idcs = row_start_idcs;
		rv = pthread_create(&threadpool[idx], &attr, pthread_set_rows_offset, args+idx);
		if (rv){
			printf("ERROR; return code from pthread_create() is %d\n", rv);
			exit(EXIT_FAILURE);
		}
	}

	for (int idx = 0; idx < THREAD_NUM; idx ++) {
		rv = pthread_join(threadpool[idx], NULL);
		if (rv){
			printf("ERROR; return code from pthread_join() is %d\n", rv);
			exit(EXIT_FAILURE);
		}
#ifdef IN_DEBUG
		// my_cnt_t flop_estimated = row_start_idcs[rows_offset[idx+1]]
			// - row_start_idcs[rows_offset[idx]];
		if (idx == THREAD_NUM -1 ) {
			rows_offset[THREAD_NUM] = num_row;
			// flop_estimated = row_start_idcs[rows_offset[idx+1]]
				// - row_start_idcs[rows_offset[idx]];
		}
		// printf("\t%d_flop_estimated: %ld\n", idx, flop_estimated);
		// if (idx == THREAD_NUM-1) {
		// 	printf("\n");
		// }
#endif /* IN_DEBUG */
	}

	rows_offset[THREAD_NUM] = num_row;

	pthread_attr_destroy(&attr);
}
















void write_non_temporal_batched(
	my_cnt_t *data_cols,
	my_rat_t *data_vals,
	my_cnt_t offset,
	my_cnt_t length,
	my_cnt_t *ht_check,
	my_rat_t *ht_value
) {
	my_cnt_t write_pos = offset;
	my_cnt_t index = 0;

	// Cacheline-sized temporary buffers for `cols` and `vals`
	my_cnt_t cols_buffer[NUM64_PER_CACHELINE] __attribute__((aligned(CACHELINE_SIZE)));
	my_rat_t vals_buffer[NUM64_PER_CACHELINE] __attribute__((aligned(CACHELINE_SIZE)));;

	uintptr_t start_address_cols = (uintptr_t)(data_cols + offset);  // Starting address for cols
	uintptr_t aligned_start_cols = (start_address_cols + CACHELINE_SIZE - 1) & ~(CACHELINE_SIZE - 1);  // Align to next cacheline

	my_cnt_t unaligned_prefix_entries = (aligned_start_cols - start_address_cols) / NUM64_T_SIZE;  // Entries in unaligned prefix

	// Handle unaligned prefix
	for (my_cnt_t i = 0; i < unaligned_prefix_entries && index < length; index++) {
		if (ht_check[index] != -1) {
			data_cols[write_pos] = ht_check[index];
			data_vals[write_pos] = ht_value[index];
			write_pos++;
		}
	}

	// Handle cacheline-sized aligned entries
	my_cnt_t buffer_pos = 0;
	for (; index < length; index++) {
		if (ht_check[index] != -1) {
			cols_buffer[buffer_pos] = ht_check[index];  // Write to temp buffer
			vals_buffer[buffer_pos] = ht_value[index];
			buffer_pos++;

			// If buffer is full (one cacheline), perform non-temporal write
			if (buffer_pos == NUM64_PER_CACHELINE) {
				WRITE_NT_64(&data_cols[write_pos], cols_buffer);
				WRITE_NT_64(&data_vals[write_pos], vals_buffer);
				write_pos += NUM64_PER_CACHELINE;
				buffer_pos = 0;  // Reset buffer
			}
		}
	}

	// Handle unaligned suffix (remaining entries in the buffer)
	for (my_cnt_t i = 0; i < buffer_pos; i++) {
		data_cols[write_pos] = cols_buffer[i];
		data_vals[write_pos] = vals_buffer[i];
		write_pos++;
	}
}

keyval_t* radix_sort_keyval_toggling_unrolling(keyval_t *arr, keyval_t *aux, size_t n) {
	if (n < 2) {
		return arr;
	}

	/* ------------------------------------------------------------------
	 * 1) Quick check if all keys are the same. If so, no need to sort.
	 * ----------------------------------------------------------------*/
	bool all_same = true;
	for (size_t i = 1; i < n; i++) {
		if (arr[i].key != arr[0].key) {
			all_same = false;
			break;
		}
	}
	if (all_same) {
		return arr; // already "sorted" since all keys identical
	}

	/* ------------------------------------------------------------------
	 * 2) Determine how many bytes we actually need to sort
	 *    by finding the maximum key in 'arr'.
	 * ----------------------------------------------------------------*/
	my_cnt_t max_key = 0;
	for (size_t i = 0; i < n; i++) {
		if (arr[i].key > max_key) {
			max_key = arr[i].key;
		}
	}

	// Count how many bytes are needed:
	size_t needed_passes = 0;
	while (max_key > 0) {
		max_key >>= 8; // shift by one byte
		needed_passes++;
	}
	if (needed_passes == 0) {
		needed_passes = 1; // at least 1 pass if all keys happen to be 0
	}
	// Safety check, in case something is off (e.g., 8-byte keys => max 8 passes)
	if (needed_passes > sizeof(my_cnt_t)) {
		needed_passes = sizeof(my_cnt_t);
	}

	/* ------------------------------------------------------------------
	 * 3) LSD Radix Sort with toggling buffers
	 * ----------------------------------------------------------------*/
	keyval_t *in  = arr;
	keyval_t *out = aux;

	for (size_t byte_index = 0; byte_index < needed_passes; byte_index++) {
		// 3a) Build histogram (256 buckets)
		size_t count[256] = {0};

		// Unrolled counting in chunks of 4
		size_t i = 0;
		for (; i + 3 < n; i += 4) {
			unsigned char d0 = ((unsigned char *)&in[i+0].key)[byte_index];
			unsigned char d1 = ((unsigned char *)&in[i+1].key)[byte_index];
			unsigned char d2 = ((unsigned char *)&in[i+2].key)[byte_index];
			unsigned char d3 = ((unsigned char *)&in[i+3].key)[byte_index];
			count[d0]++;
			count[d1]++;
			count[d2]++;
			count[d3]++;
		}
		for (; i < n; i++) {
			unsigned char d = ((unsigned char *)&in[i].key)[byte_index];
			count[d]++;
		}

		// 3b) Prefix sums for stable positions
		for (int d = 1; d < 256; d++) {
			count[d] += count[d - 1];
		}

		// 3c) Distribute elements (stable) into 'out'
		long j = (long)n - 1;
		for (; j - 3 >= 0; j -= 4) {
			unsigned char d0 = ((unsigned char *)&in[j  ].key)[byte_index];
			unsigned char d1 = ((unsigned char *)&in[j-1].key)[byte_index];
			unsigned char d2 = ((unsigned char *)&in[j-2].key)[byte_index];
			unsigned char d3 = ((unsigned char *)&in[j-3].key)[byte_index];
			out[ --count[d0] ] = in[j  ];
			out[ --count[d1] ] = in[j-1];
			out[ --count[d2] ] = in[j-2];
			out[ --count[d3] ] = in[j-3];
		}
		// handle leftover elements
		for (; j >= 0; j--) {
			unsigned char d = ((unsigned char *)&in[j].key)[byte_index];
			out[ --count[d] ] = in[j];
		}

		// 3d) Swap buffers (toggle)
		keyval_t *tmp = in;
		in  = out;
		out = tmp;
	}

	// Return pointer to the final sorted data (could be arr or aux)
	return in;
}




// Deepseek
// Partition function with branch prediction hints
static inline size_t omp_qsort_partition(my_cnt_t *arr, size_t low, size_t high) {
    // Median-of-three pivot selection to avoid worst-case scenarios
    size_t mid = low + (high - low) / 2;
    
    // Sort low, mid, high to get median
    if (__builtin_expect(arr[low] > arr[mid], 0)) {
        my_cnt_t tmp = arr[low];
        arr[low] = arr[mid];
        arr[mid] = tmp;
    }
    if (__builtin_expect(arr[low] > arr[high], 0)) {
        my_cnt_t tmp = arr[low];
        arr[low] = arr[high];
        arr[high] = tmp;
    }
    if (__builtin_expect(arr[mid] > arr[high], 0)) {
        my_cnt_t tmp = arr[mid];
        arr[mid] = arr[high];
        arr[high] = tmp;
    }
    
    // Place pivot at high-1
    my_cnt_t pivot = arr[mid];
    arr[mid] = arr[high - 1];
    arr[high - 1] = pivot;
    
    size_t i = low;
    size_t j = high - 1;
    
    for (;;) {
        // Move i forward while elements are less than pivot
        while (__builtin_expect(arr[++i] < pivot, 1)) {}
        
        // Move j backward while elements are greater than pivot
        while (__builtin_expect(arr[--j] > pivot, 1)) {}
        
        if (i >= j) break;
        
        // Swap elements
        my_cnt_t tmp = arr[i];
        arr[i] = arr[j];
        arr[j] = tmp;
    }
    
    // Restore pivot
    arr[high - 1] = arr[i];
    arr[i] = pivot;
    
    return i;
}

// Insertion sort for small subarrays
static inline void insertion_sort(my_cnt_t *arr, size_t low, size_t high) {
    for (size_t i = low + 1; i <= high; i++) {
        my_cnt_t key = arr[i];
        size_t j = i;
        
        // Move elements that are greater than key
        while (j > low && __builtin_expect(arr[j - 1] > key, 0)) {
            arr[j] = arr[j - 1];
            j--;
        }
        arr[j] = key;
    }
}

// Stack frame for iterative quicksort
typedef struct {
    size_t low;
    size_t high;
} stack_frame;

// Parallel quicksort function
void qsort_omp_parallel(my_cnt_t *arr, size_t n, numa_cfg_t *numa_cfg) {

	size_t insertion_sort_threshold = 32;

    if (n <= 1) return;
    
    // Use a stack to avoid recursion
    stack_frame *stack = (stack_frame *) alloc2_memory(
    	sizeof(stack_frame) * (n / 2 + 1), numa_cfg
    );
    if (!stack) return; // Handle allocation failure
    
    size_t stack_top = 0;
    stack[stack_top].low = 0;
    stack[stack_top].high = n - 1;
    stack_top++;
    
    // Determine the number of threads to use
    int thread_num = THREAD_NUM;
    
    #pragma omp parallel num_threads(thread_num)
    {
        while (1) {
            size_t current_low, current_high;
            int has_work = 0;
            
            // Critical section to get work from stack
            #pragma omp critical
            {
                if (stack_top > 0) {
                    stack_top--;
                    current_low = stack[stack_top].low;
                    current_high = stack[stack_top].high;
                    has_work = 1;
                }
            }
            
            if (!has_work) break;
            
            while (current_high > current_low) {
                if (current_high - current_low < insertion_sort_threshold) {
                    insertion_sort(arr, current_low, current_high);
                    break;
                }
                
                size_t p = omp_qsort_partition(arr, current_low, current_high);
                
                // Push the larger partition onto the stack
                #pragma omp critical
                {
                    if (p - current_low > current_high - p) {
                        if (p > current_low) {
                            stack[stack_top].low = current_low;
                            stack[stack_top].high = p - 1;
                            stack_top++;
                        }
                        current_low = p + 1;
                    } else {
                        if (current_high > p) {
                            stack[stack_top].low = p + 1;
                            stack[stack_top].high = current_high;
                            stack_top++;
                        }
                        current_high = p - 1;
                    }
                }
            }
        }
    }
    
    dealloc_memory(stack, sizeof(stack_frame) * (n / 2 + 1));
}





#ifdef IN_DEBUG
void debug_init_barrier() {
	debug_rv = pthread_barrier_init(&debug_barrier, NULL, THREAD_NUM);
	// BARRIER_ARRIVE(&debug_barrier, debug_rv);
}
#endif /* IN_DEBUG */


void cfg_print() {
	printf("CONFIGURATION:");
#ifdef USE_HUGE
	printf("\tUSE_HUGE");
#endif /* USE_HUGE */
#ifdef CSRC_PREPOPULATED
	printf("\tCSRC_PREPOPULATED");
#endif /* CSRC_PREPOPULATED */
#ifdef CSRC_WRITE_OPTIMIZED
	printf("\tCSRC_WRITE_OPTIMIZED");
#endif /* CSRC_WRITE_OPTIMIZED */
#ifdef USE_HYPERTHREADING
	printf("\tUSE_HYPERTHREADING");
#endif /* USE_HYPERTHREADING */
#ifdef USE_WEIGHTED_INTERLEAVING
	printf("\tUSE_WEIGHTED_INTERLEAVING");
#endif /* USE_WEIGHTED_INTERLEAVING */
#ifdef USE_INTERLEAVING
	printf("\tUSE_INTERLEAVING");
#endif /* USE_INTERLEAVING */
#ifdef MEMSV
	printf("\tMEMSV");
#endif /* MEMSV */
#ifdef MEMSVB
	printf("\tMEMSVB");
#endif /* MEMSVB */
#ifdef MEM_MON
	printf("\tMEM_MON");
#endif /* MEM_MON */
	printf("\tHASH_CONSTANT: %ld", HASH_CONSTANT);
	printf("\tTHREAD_NUM: %d", THREAD_NUM);
	printf("\tAMS_HIGH_THREAD_NUM: %d", AMS_HIGH_THREAD_NUM);
	printf("\tAMS_SIZE: %ld", AMS_SIZE);
	printf("\tAMS_MIN_NUM: %ld", AMS_MIN_NUM);
#ifdef AMS_LPT
	printf("\tAMS_LPT");
#elif AMS_SIMP
	printf("\tAMS_SIMP");
#else /* AMS_LPT */
	printf("\tAMS_PRT: %.6f", AMS_PRT);
#endif /* AMS_LPT */
#ifdef AMS_MAXFLP
	printf("\tAMS_MAXFLP");
#endif /* AMS_MAXFLP */
#ifdef AMS_DNR
	printf("\tAMS_DNR");
#endif /* AMS_DNR */
#ifdef AMS_HYOPT
	printf("\tAMS_HYOPT");
	printf("\tAMS_SYMB_THRES: %.6f", AMS_SYMB_THRES);
	printf("\tAMS_NUMC_THRES: %.6f", AMS_NUMC_THRES);
#endif /* AMS_HYOPT */
	printf("\tAMS_ACC_THRES: %.6f", AMS_ACC_THRES);
#ifdef SMTACT
	printf("\tSMTACT");
	printf("\tSMTACT_REAR");
	printf("\tSMTACT_QUANTILE: %.6f", SMTACT_QUANTILE);
	printf("\tSMTACT_SYMB_THRES: %ld", SMTACT_SYMB_THRES);
	printf("\tSMTACT_NUMC_THRES: %ld", SMTACT_NUMC_THRES);
#endif /* SMTACT */
	printf("\tDN_MIN: %ld", DN_MIN);
	printf("\tGUD_ALPHA: %.6f", GUD_ALPHA);
	printf("\tGUD_MIN: %ld", GUD_MIN);
	printf("\tHYB_THRES: %.6f", HYB_THRES);
	printf("\n");
}













