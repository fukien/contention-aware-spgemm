#include "common.h"

#ifdef IN_VERIFY
extern int64_t over_predict;
extern int64_t under_predict;
#endif /* IN_VERIFY */


static void* comm_pthread_init_local(void *param) {
	comm_ptrl_init_arg_t *arg = (comm_ptrl_init_arg_t *) param;

	int _tid = arg->_tid;
	csr_t *csr_a = arg->csr_a;
	csr_t *csr_b = arg->csr_b;
	my_cnt_t *rows_offset = arg->rows_offset;

	my_cnt_t local_nnz = 0, local_nnz_max = 0;

	for (my_cnt_t idx = rows_offset[_tid]; idx < rows_offset[_tid+1]; idx ++) {
		local_nnz = 0;		
		for (my_cnt_t jdx = csr_a->row_ptrs[idx]; jdx < csr_a->row_ptrs[idx+1]; jdx ++) {
			my_cnt_t t_acol = csr_a->cols[jdx];

			local_nnz += csr_b->row_ptrs[t_acol+1] - csr_b->row_ptrs[t_acol];
		}
		if (local_nnz_max < local_nnz) {
			local_nnz_max = local_nnz;
		}
	}

	*(arg->local_nnz_max) = local_nnz_max;

	return NULL;
}


void comm_init_local(csr_t *csr_a, csr_t *csr_b, my_cnt_t *rows_offset, 
	my_cnt_t*local_nnz_max, numa_cfg_t* numa_cfg) {

	pthread_t threadpool[THREAD_NUM];
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	cpu_set_t set;
	int core_idx;
	int rv;

	comm_ptrl_init_arg_t args[THREAD_NUM];
	memset(args, 0, sizeof(comm_ptrl_init_arg_t)*THREAD_NUM);

	for (int  idx = 0 ; idx < THREAD_NUM; idx ++) {
		core_idx = get_cpu_id(idx);
		CPU_ZERO(&set);
		CPU_SET(core_idx, &set);
		pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);
		args[idx]._tid = idx;
		args[idx].csr_a = csr_a;
		args[idx].csr_b = csr_b;
		args[idx].rows_offset = rows_offset;
		args[idx].local_nnz_max = local_nnz_max+idx;
		args[idx].numa_cfg = numa_cfg;
		rv = pthread_create(&threadpool[idx], &attr, comm_pthread_init_local, args+idx);
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


static void* pthread_esc_symbolic(void *param) {
#ifdef IN_DEBUG
	struct timespec time_start, time_end;
	clock_gettime(CLOCK_REALTIME, &time_start);
#endif /* IN_DEBUG */

	esc_symbolic_arg_t *arg = (esc_symbolic_arg_t *) param;
	int _tid = arg->_tid;
	csr_t *csr_a = arg->csr_a;
	csr_t *csr_b = arg->csr_b;
	my_cnt_t *rows_offset = arg->rows_offset;
	my_cnt_t *row_nz = arg->row_nz;
	numa_cfg_t * numa_cfg = arg->numa_cfg;

	my_cnt_t local_nnz_max = arg->local_nnz_max;

	if (local_nnz_max > csr_b->num_col) {
		local_nnz_max = csr_b->num_col;
	}

	local_nnz_max = NEXT_POW_2(csr_b->num_col);

	my_cnt_t *check = (my_cnt_t*) alloc2_memory(
		sizeof(my_cnt_t) * local_nnz_max, numa_cfg
	);


	for (my_cnt_t idx = rows_offset[_tid]; idx < rows_offset[_tid + 1]; idx ++) {

		register my_cnt_t nz = 0;

		if (row_nz[idx] > 0) {
			my_cnt_t ht_size = NEXT_POW_2(row_nz[idx]);
			if (ht_size > local_nnz_max) {
				ht_size = local_nnz_max;
			}
			memset(check, -1, sizeof(my_cnt_t)*ht_size);

			for (my_cnt_t jdx = csr_a->row_ptrs[idx]; jdx < csr_a->row_ptrs[idx+1]; jdx ++) {
				my_cnt_t t_acol = csr_a->cols[jdx];
				for (my_cnt_t kdx = csr_b->row_ptrs[t_acol]; kdx < csr_b->row_ptrs[t_acol+1]; kdx ++) {
					my_cnt_t key = csr_b->cols[kdx];
					my_cnt_t hash = (key * HASH_CONSTANT) & (ht_size - 1);

					while (true) { // Loop for hash probing
						if (check[hash] == key) { // if the key is already inserted, it's ok
							break;
						} else if (check[hash] == -1) { // if the key has not been inserted yet, then it's added.
							check[hash] = key;
							nz++;
							break;
						} else { // linear probing: check next entry
							hash = (hash + 1) & (ht_size - 1); //hash = (hash + 1) % ht_size
						}
					}

				}
			}

		}



#ifdef IN_VERIFY
		if (row_nz[idx] > nz) {
			(arg->over_predict) ++;
		}

		if (row_nz[idx] < nz) {
			(arg->under_predict) ++;
		}
#endif /* IN_VERIFY */

		row_nz[idx] = nz;

	}

	dealloc_memory(check, sizeof(my_cnt_t) * local_nnz_max);


#ifdef IN_DEBUG
	clock_gettime(CLOCK_REALTIME, &time_end);
	arg->_symb_time = diff_sec(time_start, time_end);
#endif /* IN_DEBUG */

	return NULL;
}


void esc_symbolic(csr_t* csr_a, csr_t* csr_b,
	my_cnt_t* rows_offset, my_cnt_t* row_nz, 
	my_cnt_t* local_nnz_max, numa_cfg_t* numa_cfg) {
	pthread_t threadpool[THREAD_NUM];
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	cpu_set_t set;
	int core_idx;
	int rv;

	esc_symbolic_arg_t args[THREAD_NUM];
	memset(args, 0, sizeof(esc_symbolic_arg_t) * THREAD_NUM);

	for (int idx = 0; idx < THREAD_NUM; idx ++) {
		core_idx = get_cpu_id(idx);
		CPU_ZERO(&set);
		CPU_SET(core_idx, &set);
		pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);
		args[idx]._tid = idx;
		args[idx].csr_a = csr_a;
		args[idx].csr_b = csr_b;
		args[idx].rows_offset = rows_offset;
		args[idx].row_nz = row_nz;
		args[idx].local_nnz_max = local_nnz_max[idx];
		args[idx].numa_cfg = numa_cfg;
		rv = pthread_create(&threadpool[idx], &attr, pthread_esc_symbolic, args+idx);
		if (rv){
			printf("ERROR; return code from pthread_create() is %d\n", rv);
			exit(EXIT_FAILURE);
		}
	}

#ifdef IN_VERIFY
	over_predict = 0;
	under_predict = 0;
#endif /* IN_VERIFY */

	for (int idx = 0; idx < THREAD_NUM; idx ++) {
		rv = pthread_join(threadpool[idx], NULL);

#ifdef IN_DEBUG
		printf("\t%d_symb_time: %.9f\n", idx, args[idx]._symb_time);
		if (idx == THREAD_NUM-1) {
			printf("\n");
		}
#endif /* IN_DEBUG */

#ifdef IN_VERIFY
		over_predict += args[idx].over_predict;
		under_predict += args[idx].under_predict;
#endif /* IN_VERIFY */
		if (rv) {
			printf("ERROR; return code from pthread_create() is %d\n", rv);
			exit(EXIT_FAILURE);
		}
	}

	pthread_attr_destroy(&attr);

#ifdef IN_VERIFY
	printf("VERIFYPREDICT:\t%ld\t%ld\n", over_predict, under_predict);
#endif /* IN_VERIFY */

}


void init_bin(bin_t* bin, my_cnt_t num_row, my_cnt_t ht_size, numa_cfg_t *numa_cfg) {
	memset(bin, 0, sizeof(bin_t));

	bin->row_nz = (my_cnt_t*) alloc2_memory(
		sizeof(my_cnt_t) * num_row, numa_cfg
	);

	bin->rows_offset = (my_cnt_t*) alloc2_memory(
		sizeof(my_cnt_t) * (THREAD_NUM+1), numa_cfg
	);

	bin->bin_id = (char*) alloc2_memory(
		sizeof(char) * num_row, numa_cfg
	);

	bin->local_hash_table_id = (my_cnt_t **) alloc2_memory(
		sizeof(my_cnt_t*) * THREAD_NUM, numa_cfg
	);

	bin->local_hash_table_val = (my_rat_t **) alloc2_memory(
		sizeof(my_rat_t*) * THREAD_NUM, numa_cfg
	);

	// memset(bin->row_nz, 0, sizeof(my_cnt_t) * num_row);
	// memset(bin->rows_offset, 0, sizeof(my_cnt_t) * (THREAD_NUM+1));
	// memset(bin->bin_id, 0, sizeof(char) * num_row);
	// memset(bin->local_hash_table_id, 0, sizeof(my_cnt_t*) * THREAD_NUM);
	// memset(bin->local_hash_table_val, 0, sizeof(my_rat_t*) * THREAD_NUM);

	bin->num_row = num_row;
	bin->min_ht_size = ht_size;
}


void free_bin(bin_t *bin) {
	dealloc_memory(bin->row_nz, sizeof(my_cnt_t)*bin->num_row);
	dealloc_memory(bin->rows_offset, sizeof(my_cnt_t)*(THREAD_NUM+1));
	dealloc_memory(bin->bin_id, sizeof(char)*bin->num_row);

	for (int idx = 0; idx < THREAD_NUM; idx ++) {
		dealloc_memory(bin->local_hash_table_id[idx], sizeof(my_cnt_t)*bin->ht_size_array[idx]);
		dealloc_memory(bin->local_hash_table_val[idx], sizeof(my_rat_t)*bin->ht_size_array[idx]);
	}

	dealloc_memory(bin->local_hash_table_id, sizeof(my_cnt_t*)*THREAD_NUM);
	dealloc_memory(bin->local_hash_table_val, sizeof(my_rat_t*)*THREAD_NUM);
}


static void *pthread_set_intprod_num(void *param) {
	pthread_set_intprod_num_arg_t *arg = (pthread_set_intprod_num_arg_t *) param;
	int _tid = arg->_tid;
	bin_t *bin = arg->bin;
	csr_t *csr_a = arg->csr_a;
	csr_t *csr_b = arg->csr_b;

	my_cnt_t each_intprod = 0;
	my_cnt_t nz_per_row = 0;

	my_cnt_t start_idx = _tid * (csr_a->num_row/THREAD_NUM);
	my_cnt_t end_idx = (_tid == THREAD_NUM-1) ? csr_a->num_row : 
		(_tid+1) * (csr_a->num_row/THREAD_NUM);

	for (my_cnt_t idx = start_idx; idx < end_idx; idx++) {
		nz_per_row = 0;
		for (my_cnt_t jdx = csr_a->row_ptrs[idx]; jdx < csr_a->row_ptrs[idx+1]; jdx ++) {
			nz_per_row += csr_b->row_ptrs[ csr_a->cols[jdx]+1 ] - csr_b->row_ptrs[ csr_a->cols[jdx] ];
		}
		bin->row_nz[idx] = nz_per_row;
		each_intprod += nz_per_row;
	}

	arg->intprod = each_intprod;
	// __sync_fetch_and_add(&bin->total_intprod, each_intprod);
	return NULL;
}


void set_intprod_num(bin_t *bin, csr_t *csr_a, csr_t *csr_b) {
	pthread_t threadpool[THREAD_NUM];
	pthread_attr_t attr;
	pthread_attr_init(&attr);

	pthread_set_intprod_num_arg_t args[THREAD_NUM];
	memset(args, 0, sizeof(pthread_set_intprod_num_arg_t)*THREAD_NUM);
	cpu_set_t set;
	int core_idx;
	int rv;
	for (int idx = 0; idx < THREAD_NUM; idx ++) {
		core_idx = get_cpu_id(idx);
		CPU_ZERO(&set);
		CPU_SET(core_idx, &set);
		pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);
		args[idx]._tid = idx;
		args[idx].intprod = 0;
		args[idx].bin = bin;
		args[idx].csr_a = csr_a;
		args[idx].csr_b = csr_b;
		rv = pthread_create(&threadpool[idx], &attr, pthread_set_intprod_num, args+idx);
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
		printf("\t%d_intprod: %ld\n", idx, args[idx].intprod);
		if (idx == THREAD_NUM-1) {
			printf("\n");
		}
#endif /* IN_DEBUG */
		bin->total_intprod += args[idx].intprod;
	}

	pthread_attr_destroy(&attr);
}


void *pthread_hash_set_rows_offset(void *param) {
	pthread_hash_set_rows_offset_arg_t *arg = (pthread_hash_set_rows_offset_arg_t *) param;
	int _tid = arg->_tid;
	my_cnt_t *ps_row_nz = arg->ps_row_nz;
	bin_t *bin = arg->bin;
	my_cnt_t num_row = arg->num_row;
	my_cnt_t average_intprod = arg->average_intprod;

	bin->rows_offset[_tid+1] = lower_bound(
		ps_row_nz, 
		ps_row_nz + num_row + 1, 
		average_intprod * (_tid + 1)
	) - ps_row_nz;

	return NULL;
}


void hash_set_rows_offset(bin_t *bin, my_cnt_t num_row, numa_cfg_t* numa_cfg) {
	my_cnt_t *ps_row_nz = (my_cnt_t *) alloc2_memory(
		sizeof(my_cnt_t) * (num_row+1), numa_cfg
	);
	// memset(ps_row_nz, 0, sizeof(my_cnt_t) * num_row);
	scan(ps_row_nz, bin->row_nz, num_row+1);

	my_cnt_t average_intprod = (bin->total_intprod + THREAD_NUM - 1) / THREAD_NUM;

	/* Search end point of each range */
	bin->rows_offset[0] = 0;

	pthread_t threadpool[THREAD_NUM];
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	cpu_set_t set;
	int core_idx;
	int rv;

	pthread_hash_set_rows_offset_arg_t args[THREAD_NUM];
	memset(args, 0, sizeof(pthread_hash_set_rows_offset_arg_t)*THREAD_NUM);

	for (int idx = 0; idx < THREAD_NUM; idx ++) {
		core_idx = get_cpu_id(idx);
		CPU_ZERO(&set);
		CPU_SET(core_idx, &set);
		pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);
		args[idx]._tid = idx;
		args[idx].num_row = num_row;
		args[idx].average_intprod = average_intprod;
		args[idx].ps_row_nz = ps_row_nz;
		args[idx].bin = bin;
		rv = pthread_create(&threadpool[idx], &attr, pthread_hash_set_rows_offset, args+idx);
		if (rv){
			printf("ERROR; return code from pthread_create() is %d\n", rv);
			exit(EXIT_FAILURE);
		}
	}

	for (int idx = 0; idx < THREAD_NUM; idx ++) {
		rv = pthread_join(threadpool[idx], NULL);
#ifdef IN_DEBUG
		// my_cnt_t flop_estimated = ps_row_nz[bin->rows_offset[idx+1]]
		// 	- ps_row_nz[bin->rows_offset[idx]];
		if (idx == THREAD_NUM -1 ) {
			bin->rows_offset[THREAD_NUM] = num_row;
			// flop_estimated = ps_row_nz[bin->rows_offset[idx+1]]
				// - ps_row_nz[bin->rows_offset[idx]];
		}
		// printf("\t%d_flop_estimated: %ld\n", idx, flop_estimated);
		// if (idx == THREAD_NUM-1) {
		// 	printf("\n");
		// }
#endif /* IN_DEBUG */

		if (rv) {
			printf("ERROR; return code from pthread_join() is %d\n", rv);
			exit(EXIT_FAILURE);
		}
	}

	bin->rows_offset[THREAD_NUM] = num_row;

	dealloc_memory(ps_row_nz, sizeof(my_cnt_t) * (num_row+1));

	pthread_attr_destroy(&attr);
}


static void *pthread_set_bin_id(void *param) {
	pthread_set_bin_id_arg_t *arg = (pthread_set_bin_id_arg_t *) param;
	int _tid = arg->_tid;
	my_cnt_t num_row = arg->num_row;
	my_cnt_t num_col = arg->num_col;
	bin_t *bin = arg->bin;

	my_cnt_t each_num_row = num_row/THREAD_NUM;
	my_cnt_t start_idx = _tid * each_num_row;
	my_cnt_t end_idx = (_tid == THREAD_NUM-1) ? num_row : (_tid+1) * each_num_row;

	for (my_cnt_t idx = start_idx; idx < end_idx; idx ++) {
		my_cnt_t nz_per_row = bin->row_nz[idx];

		if (nz_per_row > num_col) {
			nz_per_row = num_col;
		}

		if (nz_per_row == 0) {
			bin->bin_id[idx] = 0;
		} else {
			my_cnt_t jdx = 0;
			while (nz_per_row > (bin->min_ht_size << jdx)) {
				jdx++;
			}
			bin->bin_id[idx] = jdx + 1;
		}
	}

	return NULL;
}


void set_bin_id(bin_t *bin, my_cnt_t num_row, my_cnt_t num_col) {
	pthread_t threadpool[THREAD_NUM];
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	cpu_set_t set;
	int core_idx;
	int rv;

	pthread_set_bin_id_arg_t args[THREAD_NUM];
	memset(args, 0, sizeof(pthread_set_bin_id_arg_t)*THREAD_NUM);

	for (int idx = 0; idx < THREAD_NUM; idx ++) {
		core_idx = get_cpu_id(idx);
		CPU_ZERO(&set);
		CPU_SET(core_idx, &set);
		pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);
		args[idx]._tid = idx;
		args[idx].num_row = num_row;
		args[idx].num_col = num_col;
		args[idx].bin = bin;
		rv = pthread_create(&threadpool[idx], &attr, pthread_set_bin_id, args+idx);
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

	pthread_attr_destroy(&attr);	
}


// TODO: set ht_size to num_col to skip this computation cost (may not need to, as num_col is several orders larger)
static void *pthread_bin_create_local_hash_table(void *param) {
	pthread_bin_create_local_hash_table_arg_t *arg = \
		(pthread_bin_create_local_hash_table_arg_t *) param;\
	int _tid = arg->_tid;
	my_cnt_t num_col = arg->num_col;
	bin_t *bin = arg->bin;

	my_cnt_t ht_size = 0;

	for (my_cnt_t idx = bin->rows_offset[_tid]; idx < bin->rows_offset[_tid+1]; idx++) {
		if (ht_size < bin->row_nz[idx]) {
			ht_size = bin->row_nz[idx];
		}
	}

	/* the size of hash table is aligned as 2^n */
	if (ht_size > 0) {
		if (ht_size > num_col) {
			ht_size = num_col;
		}

		my_cnt_t k = bin->min_ht_size;
		while (k < ht_size) {
			k <<= 1;
		}
		ht_size = k;
	}

	bin->local_hash_table_id[_tid] = (my_cnt_t*) alloc2_memory(
		sizeof(my_cnt_t)*ht_size, arg->numa_cfg
	);
	// memset(bin->local_hash_table_id[_tid], 0, sizeof(my_cnt_t)*ht_size);

	bin->local_hash_table_val[_tid] = (my_rat_t*) alloc2_memory(
		sizeof(my_rat_t)*ht_size, arg->numa_cfg
	);
	// memset(bin->local_hash_table_val[_tid], 0, sizeof(my_cnt_t)*ht_size);

	bin->ht_size_array[_tid] = ht_size;

	return NULL;
}


void bin_create_local_hash_table(bin_t *bin, my_cnt_t num_col, numa_cfg_t* numa_cfg) {
	pthread_t threadpool[THREAD_NUM];
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	cpu_set_t set;
	int core_idx;
	int rv;

	pthread_bin_create_local_hash_table_arg_t args[THREAD_NUM];
	memset(args, 0, sizeof(pthread_bin_create_local_hash_table_arg_t)*THREAD_NUM);

	for (int idx = 0; idx < THREAD_NUM; idx ++) {
		core_idx = get_cpu_id(idx);
		CPU_ZERO(&set);
		CPU_SET(core_idx, &set);
		pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);
		args[idx]._tid = idx;
		args[idx].numa_cfg = numa_cfg;
		args[idx].num_col = num_col;
		args[idx].bin = bin;
		rv = pthread_create(&threadpool[idx], &attr, 
			pthread_bin_create_local_hash_table, args+idx);
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

	pthread_attr_destroy(&attr);
}


static void *pthread_hash_symbolic(void *param) {
#ifdef IN_DEBUG
	struct timespec time_start, time_end;
	clock_gettime(CLOCK_REALTIME, &time_start);
#endif /* IN_DEBUG */

	hash_symbolic_arg_t *arg = (hash_symbolic_arg_t *) param;
	int _tid = arg->_tid;
	csr_t *csr_a = arg->csr_a;
	csr_t *csr_b = arg->csr_b;
	bin_t *bin = arg->bin;

	my_cnt_t start_row = bin->rows_offset[_tid];
	my_cnt_t end_row = bin->rows_offset[_tid+1];

	my_cnt_t *check = bin->local_hash_table_id[_tid];

	for (my_cnt_t idx = start_row; idx < end_row; idx++) {
		my_cnt_t nz = 0;
		my_cnt_t bid = bin->bin_id[idx];

		if (bid > 0) {
			my_cnt_t ht_size = MIN_HT_S << (bid - 1); 	// determine hash table size for i-th row
			memset(check, -1, sizeof(my_cnt_t)*ht_size);

			for (my_cnt_t jdx = csr_a->row_ptrs[idx]; jdx < csr_a->row_ptrs[idx+1]; jdx++) {
				my_cnt_t t_acol = csr_a->cols[jdx];
				for (my_cnt_t kdx = csr_b->row_ptrs[t_acol]; kdx < csr_b->row_ptrs[t_acol+1]; kdx++) {

					my_cnt_t key = csr_b->cols[kdx];
					my_cnt_t hash = (key * HASH_CONSTANT) & (ht_size - 1);

					while (true) { // Loop for hash probing
						if (check[hash] == key) { // if the key is already inserted, it's ok
							break;
						} else if (check[hash] == -1) { // if the key has not been inserted yet, then it's added.
							check[hash] = key;
							nz++;
							break;
						} else { // linear probing: check next entry
							hash = (hash + 1) & (ht_size - 1); //hash = (hash + 1) % ht_size
						}
					}
				}
			}
		}

#ifdef IN_VERIFY
		if (bin->row_nz[idx] > nz) {
			(arg->over_predict) ++;
		}

		if (bin->row_nz[idx] < nz) {
			(arg->under_predict) ++;
		}
#endif /* IN_VERIFY */
		bin->row_nz[idx] = nz;

	}

#ifdef IN_DEBUG
	clock_gettime(CLOCK_REALTIME, &time_end);
	arg->_symb_time = diff_sec(time_start, time_end);
#endif /* IN_DEBUG */

	return NULL;
}

void hash_symbolic(csr_t *csr_a, csr_t *csr_b, bin_t *bin) {
	pthread_t threadpool[THREAD_NUM];
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	cpu_set_t set;
	int core_idx;
	int rv;

	hash_symbolic_arg_t args[THREAD_NUM];
	memset(args, 0, sizeof(hash_symbolic_arg_t)*THREAD_NUM);

	for (int idx = 0; idx < THREAD_NUM; idx ++) {
		core_idx = get_cpu_id(idx);
		CPU_ZERO(&set);
		CPU_SET(core_idx, &set);
		pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);
		args[idx]._tid = idx;
		args[idx].csr_a = csr_a;
		args[idx].csr_b = csr_b;
		args[idx].bin = bin;
		rv = pthread_create(&threadpool[idx], &attr, pthread_hash_symbolic, args+idx);
		if (rv){
			printf("ERROR; return code from pthread_create() is %d\n", rv);
			exit(EXIT_FAILURE);
		}
	}

#ifdef IN_VERIFY
	over_predict = 0;
	under_predict = 0;
#endif /* IN_VERIFY */

	for (int idx = 0; idx < THREAD_NUM; idx ++) {
		rv = pthread_join(threadpool[idx], NULL);

#ifdef IN_DEBUG
		printf("\t%d_symb_time: %.9f\n", idx, args[idx]._symb_time);
		if (idx == THREAD_NUM-1) {
			printf("\n");
		}
#endif /* IN_DEBUG */

#ifdef IN_VERIFY
		over_predict += args[idx].over_predict;
		under_predict += args[idx].under_predict;
#endif /* IN_VERIFY */
		if (rv) {
			printf("ERROR; return code from pthread_create() is %d\n", rv);
			exit(EXIT_FAILURE);
		}
	}

	pthread_attr_destroy(&attr);


#ifdef IN_VERIFY
	printf("VERIFYPREDICT:\t%ld\t%ld\n", over_predict, under_predict);
#endif /* IN_VERIFY */
}



