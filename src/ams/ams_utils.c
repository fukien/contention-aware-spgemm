#include "ams_utils.h"

#ifdef AMS_HYOPT
my_rat_t global_cur_phase_threshold;
my_rat_t global_b_arr_cur_val[THREAD_NUM];
#endif /* AMS_HYOPT */

#ifdef AMS_DNR
my_cnt_t global_row_idx_a;
#endif /* AMS_DNR */

/* target hash-table occupancy for the ams_* kernels -- see the note on
 * mm_ht_load in src/mm/common.c. 1.0 is the identity, so every existing
 * binary keeps its previous sizing; only ams_lfhashspgemm sets it.
 *
 * NOT under #ifdef SMTACT: ams_utils.h declares both unconditionally, the four
 * call sites below are unconditional, and ams_hashspgemm.c assigns ams_ht_load
 * unconditionally. Guarding the definitions left every non-SMTACT build with an
 * undefined reference to ams_ht_slots_for in libAb.so. */
my_rat_t ams_ht_load = 1.0;

my_cnt_t ams_ht_slots_for(my_cnt_t nz) {
	if (ams_ht_load >= 1.0 || nz <= 0) {
		return nz;
	}
	return (my_cnt_t) (nz / ams_ht_load);
}

#ifdef SMTACT
// current_phase_max_flp_threshold
my_cnt_t cph_mf_th;
bool ht_dis_flag;
#endif /* SMTACT */



#ifdef IN_EXAMINE
double _ams_set_time[THREAD_NUM];
double _ams_comp_time[THREAD_NUM];
double _ams_write_time[THREAD_NUM];
my_cnt_t _ams_row_cnt[THREAD_NUM];
my_cnt_t _ams_flop_cnt[THREAD_NUM];
my_cnt_t _ams_nnzc_cnt[THREAD_NUM];
my_cnt_t _ams_nnza_cnt[THREAD_NUM];
#endif /* IN_EXAMINE */

#ifdef IN_INSP
double *_ams_numc_time_array;
#endif /* IN_INSP */


static inline int lt_my_cnt_t(my_cnt_t a, my_cnt_t b) {
	return a < b;
}

static inline int le_my_cnt_t(my_cnt_t a, my_cnt_t b) {
	return a <= b;
}

static inline int gt_my_cnt_t(my_cnt_t a, my_cnt_t b) {
	return a > b;
}

static inline int ge_my_cnt_t(my_cnt_t a, my_cnt_t b) {
	return a >= b;
}

static inline my_cnt_t atomic_fetch_add(volatile my_cnt_t *ptr, int value) {
	// return __atomic_fetch_add(ptr, value, __ATOMIC_RELAXED);	// __ATOMIC_RELAXED/__ATOMIC_SEQ_CST
	return __sync_fetch_and_add(ptr, value);
}

static inline my_cnt_t atomic_fetch_sub(volatile my_cnt_t *ptr, int value) {
	// return __atomic_fetch_sub(ptr, value, __ATOMIC_RELAXED);	// __ATOMIC_RELAXED/__ATOMIC_SEQ_CST
	return __sync_fetch_and_sub(ptr, value);
}


static inline my_cnt_t find_max(const my_cnt_t *arr, my_cnt_t num) {
	// if (num == 0) {
	// 	fprintf(stderr, "Error: Array num is zero.\n");
	// 	return LONG_MIN;
	// }

	my_cnt_t max_value = arr[0];

	// Use branchless conditional move to avoid misprediction
	for (my_cnt_t i = 1; i < num; i++) {
		max_value = (arr[i] > max_value) ? arr[i] : max_value;
	}

	return max_value;
}


static inline int get_hyp_tid(int tid) {
	return tid < SOCKET_CORE_NUM ? tid + SOCKET_CORE_NUM : tid - SOCKET_CORE_NUM;
}


static inline void ams_hash_symb_krnl(csr_t *csr_a, csr_t *csr_b, ams_bin_t *ams_bin, 
	my_cnt_t *thr_ht_check, my_cnt_t *head_idx, my_cnt_t *tail_idx, int self_update, 
	int (*cmp_my_cnt)(my_cnt_t, my_cnt_t)) {

	ams_pack_t *ams_pack = ams_bin->ams_pack;

	while(1) {
		my_cnt_t current_head, current_tail, new_head;
		do {
			current_head = *head_idx;
			current_tail = *tail_idx;
			if (cmp_my_cnt(current_head, current_tail)) {
				// CAREFUL
				return;
			}
			new_head = current_head + self_update;

		} while (__sync_val_compare_and_swap(head_idx, current_head, new_head) != current_head);

		my_cnt_t idx2roff = ams_pack->b_arr[current_head].key;
		my_cnt_t start_row_offset = ams_pack->b_roff_maxc[idx2roff].key;
		my_cnt_t end_row_offset = ams_pack->b_roff_maxc[idx2roff+1].key;

		for (my_cnt_t row_idx_a = start_row_offset; row_idx_a < end_row_offset; row_idx_a ++) {
			my_cnt_t ht_size = ams_bin->row_flp[row_idx_a];

			my_cnt_t nnz_c = 0;

			if (ht_size > 0) {
				my_cnt_t ht_cap = ams_ht_slots_for(csr_b->num_col);
				ht_size = ams_ht_slots_for(ht_size);
				ht_size = ht_size > ht_cap ? ht_cap : ht_size;
				ht_size = NEXT_POW_2(ht_size);
				memset(thr_ht_check, -1, sizeof(my_cnt_t)*ht_size);

				for (my_cnt_t col_offset_a = csr_a->row_ptrs[row_idx_a];
					col_offset_a < csr_a->row_ptrs[row_idx_a+1]; col_offset_a++) {
					my_cnt_t col_idx_a = csr_a->cols[col_offset_a];

					for (my_cnt_t col_offset_b = csr_b->row_ptrs[col_idx_a];
						col_offset_b < csr_b->row_ptrs[col_idx_a+1]; col_offset_b++) {
						my_cnt_t key = csr_b->cols[col_offset_b];
						my_cnt_t hash = (key * HASH_CONSTANT) & (ht_size - 1);

						while (true) { // loop for hash probing
							if (thr_ht_check[hash] == key) { // if the key is already inserted, it's ok
								break;
							} else if (thr_ht_check[hash] == -1) { // if the key has not been inserted yet, then it's added.
								thr_ht_check[hash] = key;
								nnz_c ++;
								break;
							} else { // linear probing: check next entry
								hash = (hash + 1) & (ht_size - 1); //hash = (hash + 1) % ht_size
							}
						}

					}
				}
			}
			ams_bin->row_nz[row_idx_a] = nnz_c;
		}
	}

}



static inline void ams_hash_numc_krnl(csr_t *csr_a, csr_t *csr_b, csr_t *csr_c,
	ams_bin_t *ams_bin, my_cnt_t *thr_ht_check, my_rat_t *thr_ht_value,
	my_cnt_t *head_idx, my_cnt_t *tail_idx, int self_update, 
	int (*cmp_my_cnt)(my_cnt_t, my_cnt_t)) {

	int __attribute__((unused)) _tid = omp_get_thread_num();

	ams_pack_t *ams_pack = ams_bin->ams_pack;


#ifdef AMS_HYOPT
	bool hyopt_flag = (_tid < SOCKET_CORE_NUM) && (omp_get_num_threads() > SOCKET_CORE_NUM);
	int _hyp_tid = get_hyp_tid(_tid);
#endif /* AMS_HYOPT */

#ifdef IN_EXAMINE
	struct timespec examine_time_start, examine_time_end;
	double _set_time = 0.0, _comp_time = 0.0, _write_time = 0.0;
	my_cnt_t _row_cnt = 0, _flop_cnt = 0, _nnzc_cnt = 0, _nnza_cnt = 0;
#endif /* IN_EXAMINE */


	while(1) {
		my_cnt_t current_head, current_tail, new_head;
		do {
			current_head = *head_idx;
			current_tail = *tail_idx;

#ifdef AMS_HYOPT
			hyopt_flag &= global_b_arr_cur_val[_hyp_tid] > global_cur_phase_threshold;
			if (cmp_my_cnt(current_head, current_tail) || hyopt_flag) {
#else /* AMS_HYOPT */
			if (cmp_my_cnt(current_head, current_tail)) {
#endif /* AMS_HYOPT */

#ifdef IN_EXAMINE
				_ams_set_time[_tid] += _set_time;
				_ams_comp_time[_tid] += _comp_time;
				_ams_write_time[_tid] += _write_time;
				_ams_row_cnt[_tid] += _row_cnt;
				_ams_flop_cnt[_tid] += _flop_cnt;
				_ams_nnzc_cnt[_tid] += _nnzc_cnt;
				_ams_nnza_cnt[_tid] += _nnza_cnt;
#endif /* IN_EXAMINE */
				return;
			}
			new_head = current_head + self_update;

		} while (__sync_val_compare_and_swap(head_idx, current_head, new_head) != current_head);

#ifdef AMS_HYOPT
		global_b_arr_cur_val[_tid] = ams_pack->b_arr[current_head].val;
#endif /* AMS_HYOPT */

		my_cnt_t idx2roff = ams_pack->b_arr[current_head].key;
		my_cnt_t start_row_offset = ams_pack->b_roff_maxc[idx2roff].key;
		my_cnt_t end_row_offset = ams_pack->b_roff_maxc[idx2roff+1].key;

		for (my_cnt_t row_idx_a = start_row_offset; row_idx_a < end_row_offset; row_idx_a ++) {

#ifdef IN_INSP
			struct timespec insp_time_start, insp_time_end;
			clock_gettime(CLOCK_REALTIME, &insp_time_start);
#endif /* IN_INSP */

			my_cnt_t ht_size = ams_bin->row_flp[row_idx_a];

			if (ht_size > 0) {

#ifdef IN_EXAMINE
				_row_cnt ++;
				_nnza_cnt += csr_a->row_ptrs[row_idx_a+1] - csr_a->row_ptrs[row_idx_a];
				clock_gettime(CLOCK_REALTIME, &examine_time_start);
#endif /* IN_EXAMINE */

				my_cnt_t ht_cap = ams_ht_slots_for(csr_b->num_col);
				ht_size = ams_ht_slots_for(ht_size);
				ht_size = ht_size > ht_cap ? ht_cap : ht_size;
				ht_size = NEXT_POW_2(ht_size);
				memset(thr_ht_check, -1, sizeof(my_cnt_t)*ht_size);

#ifdef IN_EXAMINE
				clock_gettime(CLOCK_REALTIME, &examine_time_end);
				_set_time += diff_sec(examine_time_start, examine_time_end);
				clock_gettime(CLOCK_REALTIME, &examine_time_start);
#endif /* IN_EXAMINE */

				for (my_cnt_t col_offset_a = csr_a->row_ptrs[row_idx_a];
					col_offset_a < csr_a->row_ptrs[row_idx_a+1]; col_offset_a++) {
					my_cnt_t col_idx_a = csr_a->cols[col_offset_a];
					my_rat_t t_aval = csr_a->rats[col_offset_a];

#ifdef IN_EXAMINE
					_flop_cnt += csr_b->row_ptrs[col_idx_a+1] - csr_b->row_ptrs[col_idx_a];
#endif /* IN_EXAMINE */

					for (my_cnt_t col_offset_b = csr_b->row_ptrs[col_idx_a];
						col_offset_b < csr_b->row_ptrs[col_idx_a+1]; col_offset_b++) {
						my_cnt_t key = csr_b->cols[col_offset_b];
						my_cnt_t hash = (key * HASH_CONSTANT) & (ht_size - 1);

						my_rat_t t_val = func_mul(t_aval, csr_b->rats[col_offset_b]);

						while (true) { // loop for hash probing
							if (thr_ht_check[hash] == key) { // if the key is already inserted, it's ok
								thr_ht_value[hash] += t_val;
								break;
							} else if (thr_ht_check[hash] == -1) { // if the key has not been inserted yet, then it's added.
								thr_ht_check[hash] = key;
								thr_ht_value[hash] = t_val;
								break;
							} else { // linear probing: check next entry
								hash = (hash + 1) & (ht_size - 1); //hash = (hash + 1) % ht_size
							}
						}

					}
				}

#ifdef IN_EXAMINE
				clock_gettime(CLOCK_REALTIME, &examine_time_end);
				_comp_time += diff_sec(examine_time_start, examine_time_end);
				clock_gettime(CLOCK_REALTIME, &examine_time_start);
#endif /* IN_EXAMINE */

				my_cnt_t offset = csr_c->row_ptrs[row_idx_a];
#ifdef CSRC_WRITE_OPTIMIZED
				write_non_temporal_batched(
					csr_c->cols, csr_c->rats, 
					offset, ht_size,
					thr_ht_check, thr_ht_value
				);
#else /* CSRC_WRITE_OPTIMIZED */
				my_cnt_t index = 0;

				for (my_cnt_t jdx = 0; jdx < ht_size; jdx++) {
					if (thr_ht_check[jdx] != -1) {
						csr_c->cols[offset+index] = thr_ht_check[jdx];
						csr_c->rats[offset+index] = thr_ht_value[jdx];
						index++;
					}
				}
#endif /* CSRC_WRITE_OPTIMIZED */

#ifdef IN_EXAMINE
				clock_gettime(CLOCK_REALTIME, &examine_time_end);
				_write_time += diff_sec(examine_time_start, examine_time_end);
				_nnzc_cnt += index;
#endif /* IN_EXAMINE */

			}


#ifdef IN_INSP
			clock_gettime(CLOCK_REALTIME, &insp_time_end);
			_ams_numc_time_array[row_idx_a] += diff_sec(insp_time_start, insp_time_end);
#endif /* IN_INSP */

		}
	}
}






static inline void ams_hsersc_numc_krnl(
	csr_t *csr_a, csr_t *csr_b, csr_t *csr_c,
	ams_bin_t *ams_bin, keyval_t *keyval_org, keyval_t *aux, 
	my_cnt_t *head_idx, my_cnt_t *tail_idx, int self_update, 
	int (*cmp_my_cnt)(my_cnt_t, my_cnt_t)) {

	int __attribute__((unused)) _tid = omp_get_thread_num();

	ams_pack_t *ams_pack = ams_bin->ams_pack;


#ifdef AMS_HYOPT
	bool hyopt_flag = (_tid < SOCKET_CORE_NUM) && (omp_get_num_threads() > SOCKET_CORE_NUM);
	int _hyp_tid = get_hyp_tid(_tid);
#endif /* AMS_HYOPT */

#ifdef IN_EXAMINE
	struct timespec examine_time_start, examine_time_end;
	double _set_time = 0.0, _comp_time = 0.0, _write_time = 0.0;
	my_cnt_t _row_cnt = 0, _flop_cnt = 0, _nnzc_cnt = 0, _nnza_cnt = 0;
#endif /* IN_EXAMINE */

	while(1) {
		my_cnt_t current_head, current_tail, new_head;
		do {
			current_head = *head_idx;
			current_tail = *tail_idx;

#ifdef AMS_HYOPT
			hyopt_flag &= global_b_arr_cur_val[_hyp_tid] > global_cur_phase_threshold;
			if (cmp_my_cnt(current_head, current_tail) || hyopt_flag) {
#else /* AMS_HYOPT */
			if (cmp_my_cnt(current_head, current_tail)) {
#endif /* AMS_HYOPT */

#ifdef IN_EXAMINE
				_ams_set_time[_tid] += _set_time;
				_ams_comp_time[_tid] += _comp_time;
				_ams_write_time[_tid] += _write_time;
				_ams_row_cnt[_tid] += _row_cnt;
				_ams_flop_cnt[_tid] += _flop_cnt;
				_ams_nnzc_cnt[_tid] += _nnzc_cnt;
				_ams_nnza_cnt[_tid] += _nnza_cnt;
#endif /* IN_EXAMINE */
				return;
			}
			new_head = current_head + self_update;

		} while (__sync_val_compare_and_swap(head_idx, current_head, new_head) != current_head);

#ifdef AMS_HYOPT
		global_b_arr_cur_val[_tid] = ams_pack->b_arr[current_head].val;
#endif /* AMS_HYOPT */

		my_cnt_t idx2roff = ams_pack->b_arr[current_head].key;
		my_cnt_t start_row_offset = ams_pack->b_roff_maxc[idx2roff].key;
		my_cnt_t end_row_offset = ams_pack->b_roff_maxc[idx2roff+1].key;

		for (my_cnt_t row_idx_a = start_row_offset; row_idx_a < end_row_offset; row_idx_a ++) {

#ifdef IN_INSP
			struct timespec insp_time_start, insp_time_end;
			clock_gettime(CLOCK_REALTIME, &insp_time_start);
#endif /* IN_INSP */


			if (ams_bin->row_flp[row_idx_a] > 0) {
				my_cnt_t index_0 = 0;
				keyval_t *keyval = keyval_org;

#ifdef IN_EXAMINE
				_row_cnt ++;
				_nnza_cnt += csr_a->row_ptrs[row_idx_a+1] - csr_a->row_ptrs[row_idx_a];
				clock_gettime(CLOCK_REALTIME, &examine_time_start);
#endif /* IN_EXAMINE */


				for (my_cnt_t col_offset_a = csr_a->row_ptrs[row_idx_a];
					col_offset_a < csr_a->row_ptrs[row_idx_a+1]; col_offset_a++) {
					my_cnt_t col_idx_a = csr_a->cols[col_offset_a];
					my_rat_t t_aval = csr_a->rats[col_offset_a];

#ifdef IN_EXAMINE
					_flop_cnt += csr_b->row_ptrs[col_idx_a+1] - csr_b->row_ptrs[col_idx_a];
#endif /* IN_EXAMINE */

					for (my_cnt_t col_offset_b = csr_b->row_ptrs[col_idx_a];
						col_offset_b < csr_b->row_ptrs[col_idx_a+1]; col_offset_b++) {

						keyval[index_0].key = csr_b->cols[col_offset_b];
						keyval[index_0].val = func_mul(t_aval, csr_b->rats[col_offset_b]);

						index_0 ++;
					}
				}

				keyval = radix_sort_keyval_toggling_unrolling(keyval_org, aux, index_0);

				my_cnt_t offset = csr_c->row_ptrs[row_idx_a];
				my_cnt_t index_1 = 0;

				keyval_t* aux_1 = keyval == keyval_org ? aux : keyval_org;
				aux_1[index_1] = keyval[0];

				for (my_cnt_t ldx = 1; ldx < index_0; ldx ++) {
					bool flag = keyval[ldx].key == keyval[ldx-1].key;
					index_1 += !flag;
					aux_1[index_1].key = keyval[ldx].key;
					aux_1[index_1].val = aux_1[index_1].val * flag + keyval[ldx].val;
				}

#ifdef IN_EXAMINE
				clock_gettime(CLOCK_REALTIME, &examine_time_end);
				_comp_time += diff_sec(examine_time_start, examine_time_end);
				clock_gettime(CLOCK_REALTIME, &examine_time_start);
#endif /* IN_EXAMINE */


				for (my_cnt_t ldx = 0; ldx < index_1+1; ldx ++) {
					csr_c->cols[ldx+offset] = aux_1[ldx].key;
					csr_c->rats[ldx+offset] = aux_1[ldx].val;
				}


#ifdef IN_EXAMINE
				clock_gettime(CLOCK_REALTIME, &examine_time_end);
				_write_time += diff_sec(examine_time_start, examine_time_end);
				_nnzc_cnt += index;
#endif /* IN_EXAMINE */

			}


#ifdef IN_INSP
			clock_gettime(CLOCK_REALTIME, &insp_time_end);
			_ams_numc_time_array[row_idx_a] += diff_sec(insp_time_start, insp_time_end);
#endif /* IN_INSP */

		}
	}
}







void ams_hash_symb(csr_t *csr_a, csr_t* csr_b, my_cnt_t *ht_check, ams_bin_t *ams_bin) {
	typedef int (*cmp_my_cnt_t)(my_cnt_t, my_cnt_t);
	ams_pack_t *ams_pack = ams_bin->ams_pack;

	int thread_num = THREAD_NUM;
#if defined(SMTACT) && (THREAD_NUM > SOCKET_CORE_NUM)
	if (ht_dis_flag) {
		thread_num = SOCKET_CORE_NUM;
	}
#endif /* defined(SMTACT) && (THREAD_NUM > SOCKET_CORE_NUM) */

	int __attribute__((unused)) ams_high_thread_num = AMS_HIGH_THREAD_NUM;

#ifdef IN_DEBUG
	double symb_time[thread_num];
#endif /* IN_DEBUG */

	#pragma omp parallel num_threads(thread_num)
	{
#ifdef IN_DEBUG
		struct timespec time_start, time_end;
		clock_gettime(CLOCK_REALTIME, &time_start);
#endif /* IN_DEBUG */
		// --- Bind each thread to a physical core (Linux example) ---
#ifdef __linux__
		int tid = omp_get_thread_num();
		cpu_set_t cpuset;
		CPU_ZERO(&cpuset);
		// Here we assume that the OpenMP thread number maps directly to a physical core.
		CPU_SET(get_cpu_id(tid), &cpuset);
		if (pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset) != 0) {
			perror("pthread_setaffinity_np");
			// Optionally handle the error or exit.
		}
#endif	/* __linux__ */

		my_cnt_t *thr_ht_check = ht_check + tid * NEXT_POW_2(csr_b->num_col);

		my_cnt_t *head_idx, *tail_idx;
		int self_update;
		cmp_my_cnt_t cmp_my_cnt;


#ifdef AMS_LPT
		head_idx = &(ams_pack->head_idx);
		tail_idx = &(ams_pack->tail_idx);
		self_update = -1;
		cmp_my_cnt = le_my_cnt_t;

		ams_hash_symb_krnl(csr_a, csr_b, ams_bin, thr_ht_check,
			head_idx, tail_idx, self_update, cmp_my_cnt);
#elif AMS_SIMP	/* AMS_LPT */
		if (tid >= thread_num - ams_high_thread_num) {
			head_idx = &(ams_pack->head_idx);
			tail_idx = &(ams_pack->tail_idx);
			self_update = -1;
			cmp_my_cnt = le_my_cnt_t;
		} else {
			head_idx = &(ams_pack->tail_idx);
			tail_idx = &(ams_pack->head_idx);
			self_update = 1;
			cmp_my_cnt = gt_my_cnt_t;
		}

		ams_hash_symb_krnl(csr_a, csr_b, ams_bin, thr_ht_check,
			head_idx, tail_idx, self_update, cmp_my_cnt);
#else	/* AMS_LPT */
		if (tid >= thread_num - ams_high_thread_num/2 ) {				// region a: start from high_head_idx
			head_idx = &(ams_pack->high_head_idx);
			tail_idx = &(ams_pack->high_tail_idx);
			self_update = -1;
			cmp_my_cnt = le_my_cnt_t;
		} else if (tid >= thread_num - ams_high_thread_num) {			// region b: start from high_tail_idx
			head_idx = &(ams_pack->high_tail_idx);
			tail_idx = &(ams_pack->high_head_idx);
			self_update = 1;
			cmp_my_cnt = gt_my_cnt_t;
		} else {														// region d: start from low_head_idx
			head_idx = &(ams_pack->low_head_idx);
			tail_idx = &(ams_pack->low_tail_idx);
			self_update = 1;
			cmp_my_cnt = gt_my_cnt_t;
		}


		ams_hash_symb_krnl(csr_a, csr_b, ams_bin, thr_ht_check,
			head_idx, tail_idx, self_update, cmp_my_cnt);

		// reassign head and tail indices to check for unfinished jobs in the opposite region
		if (tid >= thread_num - ams_high_thread_num) {			// region a/b switch to d
			head_idx = &(ams_pack->low_head_idx);
			tail_idx = &(ams_pack->low_tail_idx);
			self_update = 1;
			cmp_my_cnt = gt_my_cnt_t;
		} else {												// region c/d switch to b
			head_idx = &(ams_pack->high_tail_idx);
			tail_idx = &(ams_pack->high_head_idx);
			self_update = 1;
			cmp_my_cnt = gt_my_cnt_t;
		}

		ams_hash_symb_krnl(csr_a, csr_b, ams_bin, thr_ht_check,
			head_idx, tail_idx, self_update, cmp_my_cnt);
#endif	/* AMS_LPT */

#ifdef IN_DEBUG
		clock_gettime(CLOCK_REALTIME, &time_end);
		symb_time[tid] = diff_sec(time_start, time_end);
#endif /* IN_DEBUG */
	}

#ifdef IN_DEBUG
	for (int tid = 0; tid < thread_num; tid ++) {
		printf("\t%d_symb_time: %.9f\n", tid, symb_time[tid]);
		if (tid == thread_num-1) {
			printf("\n");
		}
	}
#endif /* IN_DEBUG */
}





void ams_hash_numc(csr_t *csr_a, csr_t* csr_b, csr_t *csr_c, 
	my_cnt_t *ht_check, my_rat_t *ht_value, ams_bin_t *ams_bin) {
	typedef int (*cmp_my_cnt_t)(my_cnt_t, my_cnt_t);
	ams_pack_t *ams_pack = ams_bin->ams_pack;

	int thread_num = THREAD_NUM;
#if defined(SMTACT) && (THREAD_NUM > SOCKET_CORE_NUM)
	if (ht_dis_flag) {
		thread_num = SOCKET_CORE_NUM;
	}
#endif /* defined(SMTACT) && (THREAD_NUM > SOCKET_CORE_NUM) */

	int __attribute__((unused)) ams_high_thread_num = AMS_HIGH_THREAD_NUM;

#ifdef IN_DEBUG
	double numc_time[thread_num];
#endif /* IN_DEBUG */

#ifdef IN_EXAMINE
	memset(_ams_set_time, 0, sizeof(double)*THREAD_NUM);
	memset(_ams_comp_time, 0, sizeof(double)*THREAD_NUM);
	memset(_ams_write_time, 0, sizeof(double)*THREAD_NUM);
	memset(_ams_row_cnt, 0, sizeof(my_cnt_t)*THREAD_NUM);
	memset(_ams_flop_cnt, 0, sizeof(my_cnt_t)*THREAD_NUM);
	memset(_ams_nnzc_cnt, 0, sizeof(my_cnt_t)*THREAD_NUM);
	memset(_ams_nnza_cnt, 0, sizeof(my_cnt_t)*THREAD_NUM);
#endif /* IN_EXAMINE */


#ifdef AMS_HYOPT
	memset(global_b_arr_cur_val, 0, sizeof(my_rat_t)*THREAD_NUM);
	global_cur_phase_threshold = AMS_NUMC_THRES;
#endif /* AMS_HYOPT */

	#pragma omp parallel num_threads(thread_num)
	{
#ifdef IN_DEBUG
		struct timespec time_start, time_end;
		clock_gettime(CLOCK_REALTIME, &time_start);
#endif /* IN_DEBUG */
		// --- Bind each thread to a physical core (Linux example) ---
#ifdef __linux__
		int tid = omp_get_thread_num();
		cpu_set_t cpuset;
		CPU_ZERO(&cpuset);
		// Here we assume that the OpenMP thread number maps directly to a physical core.
		CPU_SET(get_cpu_id(tid), &cpuset);
		if (pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset) != 0) {
			perror("pthread_setaffinity_np");
			// Optionally handle the error or exit.
		}
#endif	/* __linux__ */


		my_cnt_t *thr_ht_check = ht_check + tid * NEXT_POW_2(csr_b->num_col);
		my_rat_t *thr_ht_value = ht_value + tid * NEXT_POW_2(csr_c->num_col);

		my_cnt_t *head_idx, *tail_idx;
		int self_update;
		cmp_my_cnt_t cmp_my_cnt;


#ifdef AMS_LPT
		head_idx = &(ams_pack->head_idx);
		tail_idx = &(ams_pack->tail_idx);
		self_update = -1;
		cmp_my_cnt = le_my_cnt_t;

		ams_hash_numc_krnl(csr_a, csr_b, csr_c, 
			ams_bin, thr_ht_check, thr_ht_value,
			head_idx, tail_idx, self_update, cmp_my_cnt
		);
#elif AMS_SIMP	/* AMS_LPT */
		if (tid >= thread_num - ams_high_thread_num) {
			head_idx = &(ams_pack->head_idx);
			tail_idx = &(ams_pack->tail_idx);
			self_update = -1;
			cmp_my_cnt = le_my_cnt_t;
		} else {
			head_idx = &(ams_pack->tail_idx);
			tail_idx = &(ams_pack->head_idx);
			self_update = 1;
			cmp_my_cnt = gt_my_cnt_t;
		}
		ams_hash_numc_krnl(csr_a, csr_b, csr_c, 
			ams_bin, thr_ht_check, thr_ht_value,
			head_idx, tail_idx, self_update, cmp_my_cnt
		);
#else	/* AMS_LPT */
		if (tid >= thread_num - ams_high_thread_num/2 ) {				// region a: start from high_head_idx
			head_idx = &(ams_pack->high_head_idx);
			tail_idx = &(ams_pack->high_tail_idx);
			self_update = -1;
			cmp_my_cnt = le_my_cnt_t;
		} else if (tid >= thread_num - ams_high_thread_num) {			// region b: start from high_tail_idx
			head_idx = &(ams_pack->high_tail_idx);
			tail_idx = &(ams_pack->high_head_idx);
			self_update = 1;
			cmp_my_cnt = gt_my_cnt_t;
		} else {														// region d: start from low_head_idx
			head_idx = &(ams_pack->low_head_idx);
			tail_idx = &(ams_pack->low_tail_idx);
			self_update = 1;
			cmp_my_cnt = gt_my_cnt_t;
		}

		ams_hash_numc_krnl(csr_a, csr_b, csr_c, 
			ams_bin, thr_ht_check, thr_ht_value,
			head_idx, tail_idx, self_update, cmp_my_cnt
		);


		// reassign head and tail indices to check for unfinished jobs in the opposite region
		if (tid >= thread_num - ams_high_thread_num) {			// region a/b switch to d
			head_idx = &(ams_pack->low_head_idx);
			tail_idx = &(ams_pack->low_tail_idx);
			self_update = 1;
			cmp_my_cnt = gt_my_cnt_t;
		} else {												// region c/d switch to b
			head_idx = &(ams_pack->high_tail_idx);
			tail_idx = &(ams_pack->high_head_idx);
			self_update = 1;
			cmp_my_cnt = gt_my_cnt_t;
		}

		ams_hash_numc_krnl(csr_a, csr_b, csr_c, 
			ams_bin, thr_ht_check, thr_ht_value,
			head_idx, tail_idx, self_update, cmp_my_cnt
		);
#endif	/* AMS_LPT */

#ifdef IN_DEBUG
		clock_gettime(CLOCK_REALTIME, &time_end);
		numc_time[tid] = diff_sec(time_start, time_end);
#endif /* IN_DEBUG */

	}



#ifdef IN_DEBUG
	for (int tid = 0; tid < thread_num; tid ++) {
		printf("\t%d_numc_time: %.9f", tid, numc_time[tid]);

#ifdef IN_EXAMINE
		printf("\t_set_time: %.9f\t_comp_time: %.9f\t_write_time: %.9f", 
			_ams_set_time[tid], _ams_comp_time[tid], _ams_write_time[tid]
		);
		printf("\t_row_cnt: %ld\t_flop_cnt: %ld\t_nnzc_cnt: %ld",
			_ams_row_cnt[tid], _ams_flop_cnt[tid], _ams_nnzc_cnt[tid]
		);
		printf("\t_nna_cnt: %ld\t_cfi: %.6f\t_cfo: %.6f\n", 
			_ams_nnza_cnt[tid], 
			(double)_ams_flop_cnt[tid]/_ams_nnza_cnt[tid], 
			(double)_ams_flop_cnt[tid]/_ams_nnzc_cnt[tid]
		);
#endif /* IN_EXAMINE */
		printf("\n");
		if (tid == thread_num-1) {
			printf("\n");
		}
	}


#endif /* IN_DEBUG */

}






void ams_hsersc_numc(
	csr_t *csr_a, csr_t* csr_b, csr_t *csr_c, 
	keyval_t *keyval_org, keyval_t *aux, 
	my_cnt_t max_row_flp, ams_bin_t *ams_bin) {
	typedef int (*cmp_my_cnt_t)(my_cnt_t, my_cnt_t);
	ams_pack_t *ams_pack = ams_bin->ams_pack;

	int thread_num = THREAD_NUM;
#if defined(SMTACT) && (THREAD_NUM > SOCKET_CORE_NUM)
	if (ht_dis_flag) {
		thread_num = SOCKET_CORE_NUM;
	}
#endif /* defined(SMTACT) && (THREAD_NUM > SOCKET_CORE_NUM) */

	int __attribute__((unused)) ams_high_thread_num = AMS_HIGH_THREAD_NUM;

#ifdef IN_DEBUG
	double numc_time[thread_num];
#endif /* IN_DEBUG */

#ifdef IN_EXAMINE
	memset(_ams_set_time, 0, sizeof(double)*THREAD_NUM);
	memset(_ams_comp_time, 0, sizeof(double)*THREAD_NUM);
	memset(_ams_write_time, 0, sizeof(double)*THREAD_NUM);
	memset(_ams_row_cnt, 0, sizeof(my_cnt_t)*THREAD_NUM);
	memset(_ams_flop_cnt, 0, sizeof(my_cnt_t)*THREAD_NUM);
	memset(_ams_nnzc_cnt, 0, sizeof(my_cnt_t)*THREAD_NUM);
	memset(_ams_nnza_cnt, 0, sizeof(my_cnt_t)*THREAD_NUM);
#endif /* IN_EXAMINE */


#ifdef AMS_HYOPT
	memset(global_b_arr_cur_val, 0, sizeof(my_rat_t)*THREAD_NUM);
	global_cur_phase_threshold = AMS_NUMC_THRES;
#endif /* AMS_HYOPT */

	#pragma omp parallel num_threads(thread_num)
	{
#ifdef IN_DEBUG
		struct timespec time_start, time_end;
		clock_gettime(CLOCK_REALTIME, &time_start);
#endif /* IN_DEBUG */
		// --- Bind each thread to a physical core (Linux example) ---
#ifdef __linux__
		int tid = omp_get_thread_num();
		cpu_set_t cpuset;
		CPU_ZERO(&cpuset);
		// Here we assume that the OpenMP thread number maps directly to a physical core.
		CPU_SET(get_cpu_id(tid), &cpuset);
		if (pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset) != 0) {
			perror("pthread_setaffinity_np");
			// Optionally handle the error or exit.
		}
#endif	/* __linux__ */


		keyval_t *thr_keyval_org = keyval_org + tid * max_row_flp;
		keyval_t *thr_aux = aux + tid * max_row_flp;

		my_cnt_t *head_idx, *tail_idx;
		int self_update;
		cmp_my_cnt_t cmp_my_cnt;


#ifdef AMS_LPT
		head_idx = &(ams_pack->head_idx);
		tail_idx = &(ams_pack->tail_idx);
		self_update = -1;
		cmp_my_cnt = le_my_cnt_t;

		ams_hsersc_numc_krnl(csr_a, csr_b, csr_c, 
			ams_bin, thr_keyval_org, thr_aux,
			head_idx, tail_idx, self_update, cmp_my_cnt
		);
#elif AMS_SIMP	/* AMS_LPT */
		if (tid >= thread_num - ams_high_thread_num) {
			head_idx = &(ams_pack->head_idx);
			tail_idx = &(ams_pack->tail_idx);
			self_update = -1;
			cmp_my_cnt = le_my_cnt_t;	
		} else {
			head_idx = &(ams_pack->tail_idx);
			tail_idx = &(ams_pack->head_idx);
			self_update = 1;
			cmp_my_cnt = gt_my_cnt_t;	
		}

		ams_hsersc_numc_krnl(csr_a, csr_b, csr_c, 
			ams_bin, thr_keyval_org, thr_aux,
			head_idx, tail_idx, self_update, cmp_my_cnt
		);
#else	/* AMS_LPT */
		if (tid >= thread_num - ams_high_thread_num/2 ) {				// region a: start from high_head_idx
			head_idx = &(ams_pack->high_head_idx);
			tail_idx = &(ams_pack->high_tail_idx);
			self_update = -1;
			cmp_my_cnt = le_my_cnt_t;
		} else if (tid >= thread_num - ams_high_thread_num) {			// region b: start from high_tail_idx
			head_idx = &(ams_pack->high_tail_idx);
			tail_idx = &(ams_pack->high_head_idx);
			self_update = 1;
			cmp_my_cnt = gt_my_cnt_t;
		} else {														// region d: start from low_head_idx
			head_idx = &(ams_pack->low_head_idx);
			tail_idx = &(ams_pack->low_tail_idx);
			self_update = 1;
			cmp_my_cnt = gt_my_cnt_t;
		}

		ams_hsersc_numc_krnl(csr_a, csr_b, csr_c, 
			ams_bin, thr_keyval_org, thr_aux,
			head_idx, tail_idx, self_update, cmp_my_cnt
		);

		// reassign head and tail indices to check for unfinished jobs in the opposite region
		if (tid >= thread_num - ams_high_thread_num) {			// region a/b switch to d
			head_idx = &(ams_pack->low_head_idx);
			tail_idx = &(ams_pack->low_tail_idx);
			self_update = 1;
			cmp_my_cnt = gt_my_cnt_t;
		} else {												// region c/d switch to b
			head_idx = &(ams_pack->high_tail_idx);
			tail_idx = &(ams_pack->high_head_idx);
			self_update = 1;
			cmp_my_cnt = gt_my_cnt_t;
		}

		ams_hsersc_numc_krnl(csr_a, csr_b, csr_c, 
			ams_bin, thr_keyval_org, thr_aux,
			head_idx, tail_idx, self_update, cmp_my_cnt
		);
#endif	/* AMS_LPT */

#ifdef IN_DEBUG
		clock_gettime(CLOCK_REALTIME, &time_end);
		numc_time[tid] = diff_sec(time_start, time_end);
#endif /* IN_DEBUG */

	}



#ifdef IN_DEBUG
	for (int tid = 0; tid < thread_num; tid ++) {
		printf("\t%d_numc_time: %.9f", tid, numc_time[tid]);

#ifdef IN_EXAMINE
		printf("\t_set_time: %.9f\t_comp_time: %.9f\t_write_time: %.9f", 
			_ams_set_time[tid], _ams_comp_time[tid], _ams_write_time[tid]
		);
		printf("\t_row_cnt: %ld\t_flop_cnt: %ld\t_nnzc_cnt: %ld",
			_ams_row_cnt[tid], _ams_flop_cnt[tid], _ams_nnzc_cnt[tid]
		);
		printf("\t_nna_cnt: %ld\t_cfi: %.6f\t_cfo: %.6f\n", 
			_ams_nnza_cnt[tid], 
			(double)_ams_flop_cnt[tid]/_ams_nnza_cnt[tid], 
			(double)_ams_flop_cnt[tid]/_ams_nnzc_cnt[tid]
		);
#endif /* IN_EXAMINE */
		printf("\n");
		if (tid == thread_num-1) {
			printf("\n");
		}
	}


#endif /* IN_DEBUG */

}






void hybacc_hash_krnl(
	my_cnt_t row_idx_a, 
	hybacc_krnl_arg_t *hybacc_krnl_arg, 
	my_cnt_t __attribute__((unused)) ht_size) {
	csr_t *csr_a = hybacc_krnl_arg->csr_a;
	csr_t *csr_b = hybacc_krnl_arg->csr_b;
	csr_t *csr_c = hybacc_krnl_arg->csr_c;

	my_cnt_t *thr_ht_check = hybacc_krnl_arg->thr_ht_check;
	my_rat_t *thr_ht_value = hybacc_krnl_arg->thr_ht_value;

	ht_size = ht_size > csr_b->num_col ? csr_b->num_col : ht_size;
	ht_size = NEXT_POW_2(ht_size);
	memset(thr_ht_check, -1, sizeof(my_cnt_t)*ht_size);

	for (my_cnt_t col_offset_a = csr_a->row_ptrs[row_idx_a];
		col_offset_a < csr_a->row_ptrs[row_idx_a+1]; col_offset_a++) {
		my_cnt_t col_idx_a = csr_a->cols[col_offset_a];
		my_rat_t t_aval = csr_a->rats[col_offset_a];

		for (my_cnt_t col_offset_b = csr_b->row_ptrs[col_idx_a];
			col_offset_b < csr_b->row_ptrs[col_idx_a+1]; col_offset_b++) {
			my_cnt_t key = csr_b->cols[col_offset_b];
			my_cnt_t hash = (key * HASH_CONSTANT) & (ht_size - 1);

			my_rat_t t_val = func_mul(t_aval, csr_b->rats[col_offset_b]);

			while (true) { // loop for hash probing
				if (thr_ht_check[hash] == key) { // if the key is already inserted, it's ok
					thr_ht_value[hash] += t_val;
					break;
				} else if (thr_ht_check[hash] == -1) { // if the key has not been inserted yet, then it's added.
					thr_ht_check[hash] = key;
					thr_ht_value[hash] = t_val;
					break;
				} else { // linear probing: check next entry
					hash = (hash + 1) & (ht_size - 1); //hash = (hash + 1) % ht_size
				}
			}

		}
	}

	my_cnt_t offset = csr_c->row_ptrs[row_idx_a];

#ifdef CSRC_WRITE_OPTIMIZED
	write_non_temporal_batched(
		csr_c->cols, csr_c->rats, 
		offset, ht_size,
		thr_ht_check, thr_ht_value
	);
#else /* CSRC_WRITE_OPTIMIZED */
	my_cnt_t index = 0;

	for (my_cnt_t jdx = 0; jdx < ht_size; jdx++) {
		if (thr_ht_check[jdx] != -1) {
			csr_c->cols[offset+index] = thr_ht_check[jdx];
			csr_c->rats[offset+index] = thr_ht_value[jdx];
			index++;
		}
	}
#endif /* CSRC_WRITE_OPTIMIZED */

}



void hybacc_hsersc_krnl(
	my_cnt_t row_idx_a, 
	hybacc_krnl_arg_t *hybacc_krnl_arg, 
	my_cnt_t __attribute__((unused)) ht_size) {

	csr_t *csr_a = hybacc_krnl_arg->csr_a;
	csr_t *csr_b = hybacc_krnl_arg->csr_b;
	csr_t *csr_c = hybacc_krnl_arg->csr_c;

	keyval_t *keyval_org = hybacc_krnl_arg->keyval_org;
	keyval_t *aux = hybacc_krnl_arg->aux;

	my_cnt_t index_0 = 0;
	keyval_t *keyval = keyval_org;


	for (my_cnt_t col_offset_a = csr_a->row_ptrs[row_idx_a];
		col_offset_a < csr_a->row_ptrs[row_idx_a+1]; col_offset_a++) {
		my_cnt_t col_idx_a = csr_a->cols[col_offset_a];
		my_rat_t t_aval = csr_a->rats[col_offset_a];

		for (my_cnt_t col_offset_b = csr_b->row_ptrs[col_idx_a];
			col_offset_b < csr_b->row_ptrs[col_idx_a+1]; col_offset_b++) {

			keyval[index_0].key = csr_b->cols[col_offset_b];
			keyval[index_0].val = func_mul(t_aval, csr_b->rats[col_offset_b]);

			index_0 ++;
		}
	}

	keyval = radix_sort_keyval_toggling_unrolling(keyval_org, aux, index_0);

	my_cnt_t offset = csr_c->row_ptrs[row_idx_a];
	my_cnt_t index_1 = 0;

	keyval_t* aux_1 = keyval == keyval_org ? aux : keyval_org;
	aux_1[index_1] = keyval[0];

	for (my_cnt_t ldx = 1; ldx < index_0; ldx ++) {
		bool flag = keyval[ldx].key == keyval[ldx-1].key;
		index_1 += !flag;
		aux_1[index_1].key = keyval[ldx].key;
		aux_1[index_1].val = aux_1[index_1].val * flag + keyval[ldx].val;
	}


	for (my_cnt_t ldx = 0; ldx < index_1+1; ldx ++) {
		csr_c->cols[ldx+offset] = aux_1[ldx].key;
		csr_c->rats[ldx+offset] = aux_1[ldx].val;
	}

}




static inline void ams_hybacc_numc_krnl(
	csr_t *csr_a, csr_t *csr_b, csr_t *csr_c, ams_bin_t *ams_bin, 
	keyval_t *keyval_org, keyval_t *aux, 
	my_cnt_t *thr_ht_check, my_rat_t *thr_ht_value,
	my_cnt_t *head_idx, my_cnt_t *tail_idx, int self_update, 
	int (*cmp_my_cnt)(my_cnt_t, my_cnt_t)) {

	int __attribute__((unused)) _tid = omp_get_thread_num();

	ams_pack_t *ams_pack = ams_bin->ams_pack;


#ifdef AMS_HYOPT
	bool hyopt_flag = (_tid < SOCKET_CORE_NUM) && (omp_get_num_threads() > SOCKET_CORE_NUM);
	int _hyp_tid = get_hyp_tid(_tid);
#endif /* AMS_HYOPT */


	while(1) {
		my_cnt_t current_head, current_tail, new_head;
		do {
			current_head = *head_idx;
			current_tail = *tail_idx;

#ifdef AMS_HYOPT
			hyopt_flag &= global_b_arr_cur_val[_hyp_tid] > global_cur_phase_threshold;
			if (cmp_my_cnt(current_head, current_tail) || hyopt_flag) {
#else /* AMS_HYOPT */
			if (cmp_my_cnt(current_head, current_tail)) {
#endif /* AMS_HYOPT */

				return;
			}
			new_head = current_head + self_update;

		} while (__sync_val_compare_and_swap(head_idx, current_head, new_head) != current_head);

#ifdef AMS_HYOPT
		global_b_arr_cur_val[_tid] = ams_pack->b_arr[current_head].val;
#endif /* AMS_HYOPT */

		my_cnt_t idx2roff = ams_pack->b_arr[current_head].key;
		my_cnt_t start_row_offset = ams_pack->b_roff_maxc[idx2roff].key;
		my_cnt_t end_row_offset = ams_pack->b_roff_maxc[idx2roff+1].key;

		hybacc_krnl_arg_t hybacc_krnl_arg;
		hybacc_krnl_arg.csr_a = csr_a;
		hybacc_krnl_arg.csr_b = csr_b;
		hybacc_krnl_arg.csr_c = csr_c;
		hybacc_krnl_arg.keyval_org = keyval_org;
		hybacc_krnl_arg.aux = aux;
		hybacc_krnl_arg.thr_ht_check = thr_ht_check;
		hybacc_krnl_arg.thr_ht_value = thr_ht_value;

		typedef void (*hybacc_kernel_func_t)(my_cnt_t, hybacc_krnl_arg_t *, my_cnt_t);
		hybacc_kernel_func_t hybacc_kernel_func;

		for (my_cnt_t row_idx_a = start_row_offset; row_idx_a < end_row_offset; row_idx_a ++) {

			if ( ams_bin->row_flp[row_idx_a] > 0) {

				my_rat_t row_idx_cf = (my_rat_t) ams_bin->row_flp[row_idx_a] 
						/ 
					( csr_c->row_ptrs[row_idx_a+1] - csr_c->row_ptrs[row_idx_a] ) ;

				hybacc_kernel_func = (row_idx_cf < AMS_ACC_THRES) 
					? hybacc_hsersc_krnl 
					: hybacc_hash_krnl;

				hybacc_kernel_func(row_idx_a, &hybacc_krnl_arg, ams_bin->row_flp[row_idx_a]);

			}

		}
	}
}



void ams_hybacc_numc(
	csr_t *csr_a, csr_t* csr_b, csr_t *csr_c, 
	keyval_t *keyval_org, keyval_t *aux, my_cnt_t max_row_flp,
	my_cnt_t *ht_check, my_rat_t *ht_value, ams_bin_t *ams_bin) {
	typedef int (*cmp_my_cnt_t)(my_cnt_t, my_cnt_t);
	ams_pack_t *ams_pack = ams_bin->ams_pack;

	int thread_num = THREAD_NUM;
#if defined(SMTACT) && (THREAD_NUM > SOCKET_CORE_NUM)
	if (ht_dis_flag) {
		thread_num = SOCKET_CORE_NUM;
	}
#endif /* defined(SMTACT) && (THREAD_NUM > SOCKET_CORE_NUM) */

	int __attribute__((unused)) ams_high_thread_num = AMS_HIGH_THREAD_NUM;

#ifdef IN_DEBUG
	double numc_time[thread_num];
#endif /* IN_DEBUG */

#ifdef IN_EXAMINE
	memset(_ams_set_time, 0, sizeof(double)*THREAD_NUM);
	memset(_ams_comp_time, 0, sizeof(double)*THREAD_NUM);
	memset(_ams_write_time, 0, sizeof(double)*THREAD_NUM);
	memset(_ams_row_cnt, 0, sizeof(my_cnt_t)*THREAD_NUM);
	memset(_ams_flop_cnt, 0, sizeof(my_cnt_t)*THREAD_NUM);
	memset(_ams_nnzc_cnt, 0, sizeof(my_cnt_t)*THREAD_NUM);
	memset(_ams_nnza_cnt, 0, sizeof(my_cnt_t)*THREAD_NUM);
#endif /* IN_EXAMINE */


#ifdef AMS_HYOPT
	memset(global_b_arr_cur_val, 0, sizeof(my_rat_t)*THREAD_NUM);
	global_cur_phase_threshold = AMS_NUMC_THRES;
#endif /* AMS_HYOPT */

	#pragma omp parallel num_threads(thread_num)
	{
#ifdef IN_DEBUG
		struct timespec time_start, time_end;
		clock_gettime(CLOCK_REALTIME, &time_start);
#endif /* IN_DEBUG */
		// --- Bind each thread to a physical core (Linux example) ---
#ifdef __linux__
		int tid = omp_get_thread_num();
		cpu_set_t cpuset;
		CPU_ZERO(&cpuset);
		// Here we assume that the OpenMP thread number maps directly to a physical core.
		CPU_SET(get_cpu_id(tid), &cpuset);
		if (pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset) != 0) {
			perror("pthread_setaffinity_np");
			// Optionally handle the error or exit.
		}
#endif	/* __linux__ */


		keyval_t *thr_keyval_org = keyval_org + tid * max_row_flp;
		keyval_t *thr_aux = aux + tid * max_row_flp;

		my_cnt_t *thr_ht_check = ht_check + tid * NEXT_POW_2(csr_b->num_col);
		my_rat_t *thr_ht_value = ht_value + tid * NEXT_POW_2(csr_c->num_col);

		my_cnt_t *head_idx, *tail_idx;
		int self_update;
		cmp_my_cnt_t cmp_my_cnt;


#ifdef AMS_LPT
		head_idx = &(ams_pack->head_idx);
		tail_idx = &(ams_pack->tail_idx);
		self_update = -1;
		cmp_my_cnt = le_my_cnt_t;

		ams_hybacc_numc_krnl(
			csr_a, csr_b, csr_c, ams_bin, 
			thr_keyval_org, thr_aux, 
			thr_ht_check, thr_ht_value, 
			head_idx, tail_idx, self_update, cmp_my_cnt
		);
#elif AMS_SIMP	/* AMS_LPT */
		if (tid >= thread_num - ams_high_thread_num) {
			head_idx = &(ams_pack->head_idx);
			tail_idx = &(ams_pack->tail_idx);
			self_update = -1;
			cmp_my_cnt = le_my_cnt_t;	
		} else {
			head_idx = &(ams_pack->tail_idx);
			tail_idx = &(ams_pack->head_idx);
			self_update = 1;
			cmp_my_cnt = gt_my_cnt_t;	
		}

		ams_hybacc_numc_krnl(
			csr_a, csr_b, csr_c, ams_bin, 
			thr_keyval_org, thr_aux, 
			thr_ht_check, thr_ht_value, 
			head_idx, tail_idx, self_update, cmp_my_cnt
		);
#else	/* AMS_LPT */
		if (tid >= thread_num - ams_high_thread_num/2 ) {				// region a: start from high_head_idx
			head_idx = &(ams_pack->high_head_idx);
			tail_idx = &(ams_pack->high_tail_idx);
			self_update = -1;
			cmp_my_cnt = le_my_cnt_t;
		} else if (tid >= thread_num - ams_high_thread_num) {			// region b: start from high_tail_idx
			head_idx = &(ams_pack->high_tail_idx);
			tail_idx = &(ams_pack->high_head_idx);
			self_update = 1;
			cmp_my_cnt = gt_my_cnt_t;
		} else {														// region d: start from low_head_idx
			head_idx = &(ams_pack->low_head_idx);
			tail_idx = &(ams_pack->low_tail_idx);
			self_update = 1;
			cmp_my_cnt = gt_my_cnt_t;
		}

		ams_hybacc_numc_krnl(
			csr_a, csr_b, csr_c, ams_bin, 
			thr_keyval_org, thr_aux, 
			thr_ht_check, thr_ht_value, 
			head_idx, tail_idx, self_update, cmp_my_cnt
		);

		// reassign head and tail indices to check for unfinished jobs in the opposite region
		if (tid >= thread_num - ams_high_thread_num) {			// region a/b switch to d
			head_idx = &(ams_pack->low_head_idx);
			tail_idx = &(ams_pack->low_tail_idx);
			self_update = 1;
			cmp_my_cnt = gt_my_cnt_t;
		} else {												// region c/d switch to b
			head_idx = &(ams_pack->high_tail_idx);
			tail_idx = &(ams_pack->high_head_idx);
			self_update = 1;
			cmp_my_cnt = gt_my_cnt_t;
		}

		ams_hybacc_numc_krnl(
			csr_a, csr_b, csr_c, ams_bin, 
			thr_keyval_org, thr_aux, 
			thr_ht_check, thr_ht_value, 
			head_idx, tail_idx, self_update, cmp_my_cnt
		);
#endif	/* AMS_LPT */

#ifdef IN_DEBUG
		clock_gettime(CLOCK_REALTIME, &time_end);
		numc_time[tid] = diff_sec(time_start, time_end);
#endif /* IN_DEBUG */

	}



#ifdef IN_DEBUG
	for (int tid = 0; tid < thread_num; tid ++) {
		printf("\t%d_numc_time: %.9f", tid, numc_time[tid]);

#ifdef IN_EXAMINE
		printf("\t_set_time: %.9f\t_comp_time: %.9f\t_write_time: %.9f", 
			_ams_set_time[tid], _ams_comp_time[tid], _ams_write_time[tid]
		);
		printf("\t_row_cnt: %ld\t_flop_cnt: %ld\t_nnzc_cnt: %ld",
			_ams_row_cnt[tid], _ams_flop_cnt[tid], _ams_nnzc_cnt[tid]
		);
		printf("\t_nna_cnt: %ld\t_cfi: %.6f\t_cfo: %.6f\n", 
			_ams_nnza_cnt[tid], 
			(double)_ams_flop_cnt[tid]/_ams_nnza_cnt[tid], 
			(double)_ams_flop_cnt[tid]/_ams_nnzc_cnt[tid]
		);
#endif /* IN_EXAMINE */
		printf("\n");
		if (tid == thread_num-1) {
			printf("\n");
		}
	}


#endif /* IN_DEBUG */

}






void free_ams_pack(ams_pack_t *ams_pack) {
	if (ams_pack->b_roff_maxc) {
		dealloc_memory(ams_pack->b_roff_maxc, sizeof(keyval_t) * (ams_pack->total_cnt+1));
		ams_pack->b_roff_maxc = NULL;
	}

	if (ams_pack->b_arr) {
		dealloc_memory(ams_pack->b_arr, sizeof(keyval_t) * ams_pack->total_cnt);
		ams_pack->b_arr = NULL;
	}

	ams_pack->total_cnt = 0;
	// ams_pack->high_cnt = 0;
	// ams_pack->low_cnt = 0;

#if defined(AMS_LPT) || defined(AMS_SIMP)
	ams_pack->head_idx = 0;
	ams_pack->tail_idx = 0;
#else	/* defined(AMS_LPT) || defined(AMS_SIMP) */
	ams_pack->high_head_idx = 0;
	ams_pack->high_tail_idx = 0;
	ams_pack->low_head_idx = 0;
	ams_pack->low_tail_idx = 0;
#endif	/* defined(AMS_LPT) || defined(AMS_SIMP) */
}


void init_ams_bin(ams_bin_t *ams_bin, ams_pack_t *ams_pack,
	csr_t *csr, numa_cfg_t *numa_cfg) {
	ams_bin->ams_pack = ams_pack;
	ams_bin->csr = csr;
	ams_bin->row_flp = (my_cnt_t *) alloc2_memory(
		sizeof(my_cnt_t) * csr->num_row, numa_cfg
	);
	ams_bin->ps_row_flp = (my_cnt_t *) alloc2_memory(
		sizeof(my_cnt_t) * (csr->num_row+1), numa_cfg
	);
	ams_bin->row_nz = (my_cnt_t *) alloc2_memory(
		sizeof(my_cnt_t) * csr->num_row, numa_cfg
	);
}


void free_ams_bin(ams_bin_t *ams_bin) {
	if (ams_bin->row_flp) {
		dealloc_memory(ams_bin->row_flp, sizeof(my_cnt_t) * ams_bin->csr->num_row);
		ams_bin->row_flp = NULL;
	}
	if (ams_bin->ps_row_flp) {
		dealloc_memory(ams_bin->ps_row_flp, sizeof(my_cnt_t) * (ams_bin->csr->num_row + 1));
		ams_bin->ps_row_flp = NULL;
	}
	if (ams_bin->row_nz) {
		dealloc_memory(ams_bin->row_nz, sizeof(my_cnt_t) * ams_bin->csr->num_row);
		ams_bin->row_nz = NULL;
	}
	ams_bin->csr = NULL;
	ams_bin->ams_pack = NULL;
	// if (ams_bin->ams_pack->total_cnt) {
	// 	free_ams_pack(ams_bin->ams_pack);
	// 	ams_bin->ams_pack = NULL;
	// }
}


my_cnt_t pack_ams_pack(my_cnt_t *ps_row_nz, ams_bin_t *ams_bin, numa_cfg_t *numa_cfg) {

	my_rat_t max_row_flp = 0.0;	// will cast to my_cnt_t in the end of the function

	ams_pack_t *ams_pack = ams_bin->ams_pack;

	csr_t * csr = ams_bin->csr;

	my_cnt_t ams_nnz_cnt = MIN(AMS_MIN_NUM, AMS_SIZE/sizeof(my_cnt_t));

	// CAREFUL HERE, SHOULD THINK ABOUT WHAT TO BE DIVIDED
	ams_pack->total_cnt = ps_row_nz[csr->num_row]/ams_nnz_cnt;

	ams_pack->b_arr = (keyval_t *) alloc2_memory(
		sizeof(keyval_t) * ams_pack->total_cnt, numa_cfg
	);

	ams_pack->b_roff_maxc = (keyval_t *) alloc2_memory(
		sizeof(keyval_t) * (ams_pack->total_cnt+1), numa_cfg
	);

	ams_pack->b_roff_maxc[0].key = 0;

#ifdef SMTACT
#ifndef SMTACT_REAR 
	my_cnt_t *batch_max_flp = (my_cnt_t *) alloc2_memory(sizeof(my_cnt_t)* ams_pack->total_cnt, numa_cfg);
#endif /* SMTACT_REAR */
#endif /* SMTACT */

	int thread_num = THREAD_NUM;
	my_cnt_t static_scheduled_batch_size = ams_pack->total_cnt/thread_num;
	static_scheduled_batch_size = static_scheduled_batch_size > 0 
		? static_scheduled_batch_size : 1;

	#pragma omp parallel num_threads(thread_num)
	{
		// --- Bind each thread to a physical core (Linux example) ---
	#ifdef __linux__
		int tid = omp_get_thread_num();
		cpu_set_t cpuset;
		CPU_ZERO(&cpuset);
		// Here we assume that the OpenMP thread number maps directly to a physical core.
		CPU_SET(get_cpu_id(tid), &cpuset);
		if (pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset) != 0) {
			perror("pthread_setaffinity_np");
			// Optionally handle the error or exit.
		}
	#endif	/* __linux__ */

		#pragma omp for schedule(static, static_scheduled_batch_size)
		for (my_cnt_t ams_idx = 0; ams_idx < ams_pack->total_cnt; ams_idx ++) {
			ams_pack->b_roff_maxc[ams_idx+1].key = lower_bound(
					ps_row_nz,
					ps_row_nz + csr->num_row + 1,
					ams_nnz_cnt * (ams_idx + 1)
				) - ps_row_nz;

			ams_pack->b_arr[ams_idx].key = ams_idx;

		} // <-- implicit barrier here

		#pragma omp single
		{
			ams_pack->b_roff_maxc[ams_pack->total_cnt].key = csr->num_row;	// set the end address
		}

		my_cnt_t thr_invalid_cnt = 0;
		my_rat_t thr_max_row_flp = 0;	// will cast to my_cnt_t in the end of the function
		#pragma omp for schedule(static, static_scheduled_batch_size)
		for (my_cnt_t ams_idx = 0; ams_idx < ams_pack->total_cnt; ams_idx ++) {

			my_cnt_t start_row_offset = ams_pack->b_roff_maxc[ams_idx].key;
			my_cnt_t end_row_offset = ams_pack->b_roff_maxc[ams_idx+1].key;

			my_cnt_t offset_diff = end_row_offset - start_row_offset;

			if (offset_diff > 0) {

#ifdef AMS_MAXFLP
				ams_pack->b_arr[ams_idx].val = 
					ams_pack->b_roff_maxc[ams_idx].val =
						find_max(
							ams_bin->row_flp+start_row_offset, 
							end_row_offset-start_row_offset
						);
#else /* AMS_MAXFLP */
				ams_pack->b_arr[ams_idx].val = (my_rat_t) 
					(
						ams_bin->ps_row_flp[end_row_offset] -
						ams_bin->ps_row_flp[start_row_offset]
					) / (
						ps_row_nz[end_row_offset] -
						ps_row_nz[start_row_offset]
					);
				ams_pack->b_roff_maxc[ams_idx].val = find_max(
					ams_bin->row_flp+start_row_offset, 
					end_row_offset-start_row_offset
				);
#endif /* AMS_MAXFLP */

				thr_max_row_flp = 
					ams_pack->b_roff_maxc[ams_idx].val > thr_max_row_flp ?
					ams_pack->b_roff_maxc[ams_idx].val : thr_max_row_flp;
			} else {
				ams_pack->b_roff_maxc[ams_idx].val = 0.0;
				ams_pack->b_arr[ams_idx].val = 0.0;
				thr_invalid_cnt ++;
			}

#ifdef SMTACT
#ifndef SMTACT_REAR
			batch_max_flp[ams_idx] = (my_cnt_t) ams_pack->b_roff_maxc[ams_idx].val;
#endif /* SMTACT_REAR */
#endif /* SMTACT */

		} // <-- implicit barrier here

		#pragma omp atomic
#if defined(AMS_LPT) || defined(AMS_SIMP)
		ams_pack->tail_idx += thr_invalid_cnt;
#else	/* defined(AMS_LPT) || defined(AMS_SIMP) */
		ams_pack->low_head_idx += thr_invalid_cnt;  // atomically update global variable
#endif	/* defined(AMS_LPT) || defined(AMS_SIMP) */

		#pragma omp critical
		max_row_flp = max_row_flp > thr_max_row_flp ? max_row_flp : thr_max_row_flp;
	}

#ifdef SMTACT
#ifndef SMTACT_REAR
	qsort_omp_parallel(batch_max_flp, ams_pack->total_cnt, numa_cfg);
	if (batch_max_flp[ (my_cnt_t) ( (ams_pack->total_cnt-1) * SMTACT_QUANTILE) ] > cph_mf_th) {
		ht_dis_flag = true;
	}
	dealloc_memory(batch_max_flp, sizeof(my_cnt_t) * ams_pack->total_cnt);
#else /* SMTACT_REAR */
	if (max_row_flp > cph_mf_th) {
		ht_dis_flag = true;
	}
#endif /* SMTACT_REAR */
#endif /* SMTACT */


	// compare_keyvals_3, val first (ascending), then key
	p_qsortmerge_keyval(ams_pack->b_arr, ams_pack->total_cnt, compare_keyvals_3, numa_cfg);
	// ams_pack->tail_idx = 0;
	// ams_pack->low_head_idx = 0;

#if defined(AMS_LPT) || defined(AMS_SIMP)
	ams_pack->head_idx = ams_pack->total_cnt-1;
#else	/* defined(AMS_LPT) || defined(AMS_SIMP) */
	ams_pack->high_head_idx = ams_pack->total_cnt-1;
	ams_pack->high_tail_idx = (ams_pack->total_cnt - ams_pack->low_head_idx) * AMS_PRT + ams_pack->low_head_idx;
	ams_pack->low_tail_idx = ams_pack->high_tail_idx - 1;
	// assert( ams_pack->low_tail_idx > ams_pack->low_head_idx );
#endif	/* defined(AMS_LPT) || defined(AMS_SIMP) */


#ifdef IN_DEBUG
	// for (my_cnt_t ams_idx = ams_pack->low_head_idx; ams_idx < ams_pack->total_cnt; ams_idx ++) {
	// 	printf("\t%ld_ams_pack: \tbatch_index: %ld\tbatch_cf: %.6f\tbatch_max_flp: %.2f\n",
	// 		ams_idx, ams_pack->b_arr[ams_idx].key, ams_pack->b_arr[ams_idx].val,
	// 		ams_pack->b_roff_maxc[ams_pack->b_arr[ams_idx].key].val
	// 	);
	// }
	// // fflush(stdout);

	// printf("\tlow_head_idx: %ld\t"
	// 	"low_tail_idx: %ld\t"
	// 	"high_tail_idx: %ld\t"
	// 	"high_head_idx: %ld\n",
	// 	ams_pack->low_head_idx,
	// 	ams_pack->low_tail_idx, 
	// 	ams_pack->high_tail_idx,
	// 	ams_pack->high_head_idx
	// );

	// printf("ams_nnz_cnt: %ld\ttotal_cnt: %ld\n", ams_nnz_cnt, ams_pack->total_cnt);
#endif /* IN_DEBUG */

	// printf("%s:%d:\tmax_cf: %f\n", 
	// 	__FILE__, __LINE__,
	// 	ams_pack->b_arr[ams_pack->total_cnt-1].val
	// );

	return (my_cnt_t) max_row_flp;
}



#ifdef AMS_DNR

typedef struct {
	int _tid;
	csr_t *csr_a;
	csr_t *csr_b;
	my_cnt_t *row_flp;
	my_cnt_t *row_nz;
	my_cnt_t *ht_check;
#ifdef IN_DEBUG
	double _symb_time;
#endif /* IN_DEBUG */
} ams_hash_symb_dnr_arg_t;


static inline gud_rng_t get_next_rng(my_cnt_t global_end_idx) {
	gud_rng_t gud_rng;
	my_cnt_t expected, new_idx;
	do {
		// Read the current starting index atomically.
		expected = global_row_idx_a;
		// Compute remaining work from the current position.
		new_idx = expected + DN_MIN;
		new_idx = new_idx > global_end_idx ? global_end_idx : new_idx;
	} while (__sync_val_compare_and_swap(&global_row_idx_a, expected, new_idx) != expected);
	gud_rng.start_row_idx = expected;
	gud_rng.end_row_idx = new_idx;
	return gud_rng;
}

static void* pthread_ams_hash_symb_dnr(void *param) {
#ifdef IN_DEBUG
	struct timespec time_start, time_end;
	clock_gettime(CLOCK_REALTIME, &time_start);
#endif /* IN_DEBUG */

	ams_hash_symb_dnr_arg_t *arg = (ams_hash_symb_dnr_arg_t *) param;
	int _tid = arg->_tid;
	csr_t *csr_a = arg->csr_a;
	csr_t *csr_b = arg->csr_b;
	my_cnt_t *row_flp = arg->row_flp;
	my_cnt_t *row_nz = arg->row_nz;
	my_cnt_t *thr_ht_check = arg->ht_check + _tid * NEXT_POW_2(csr_b->num_col);

	gud_rng_t gud_rng = get_next_rng(csr_a->num_row);
	while (gud_rng.start_row_idx < csr_a->num_row) {

		for (my_cnt_t row_idx_a = gud_rng.start_row_idx; row_idx_a < gud_rng.end_row_idx; row_idx_a++) {
			my_cnt_t ht_size = row_flp[row_idx_a];

			my_cnt_t nnz_c = 0;

			if (ht_size > 0) {
				my_cnt_t ht_cap = ams_ht_slots_for(csr_b->num_col);
				ht_size = ams_ht_slots_for(ht_size);
				ht_size = ht_size > ht_cap ? ht_cap : ht_size;
				ht_size = NEXT_POW_2(ht_size);
				memset(thr_ht_check, -1, sizeof(my_cnt_t)*ht_size);

				for (my_cnt_t col_offset_a = csr_a->row_ptrs[row_idx_a];
					col_offset_a < csr_a->row_ptrs[row_idx_a+1]; col_offset_a++) {
					my_cnt_t col_idx_a = csr_a->cols[col_offset_a];

					for (my_cnt_t col_offset_b = csr_b->row_ptrs[col_idx_a];
						col_offset_b < csr_b->row_ptrs[col_idx_a+1]; col_offset_b++) {
						my_cnt_t key = csr_b->cols[col_offset_b];
						my_cnt_t hash = (key * HASH_CONSTANT) & (ht_size - 1);

						while (true) { // loop for hash probing
							if (thr_ht_check[hash] == key) { // if the key is already inserted, it's ok
								break;
							} else if (thr_ht_check[hash] == -1) { // if the key has not been inserted yet, then it's added.
								thr_ht_check[hash] = key;
								nnz_c ++;
								break;
							} else { // linear probing: check next entry
								hash = (hash + 1) & (ht_size - 1); //hash = (hash + 1) % ht_size
							}
						}

					}
				}
			}

			row_nz[row_idx_a] = nnz_c;

		}

		gud_rng = get_next_rng(csr_a->num_row);
	}



#ifdef IN_DEBUG
	clock_gettime(CLOCK_REALTIME, &time_end);
	arg->_symb_time = diff_sec(time_start, time_end);
#endif /* IN_DEBUG */

	return NULL;
}

void ams_hash_symb_dnr(csr_t *csr_a, csr_t* csr_b, my_cnt_t *ht_check, ams_bin_t *ams_bin) {
	int thread_num = THREAD_NUM;
#if defined(SMTACT) && (THREAD_NUM > SOCKET_CORE_NUM)
	if (ht_dis_flag) {
		thread_num = SOCKET_CORE_NUM;
	}
#endif /* defined(SMTACT) && (THREAD_NUM > SOCKET_CORE_NUM) */

	pthread_t threadpool[thread_num];
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	cpu_set_t set;
	int core_idx;
	int rv;
	ams_hash_symb_dnr_arg_t symb_args[thread_num];
	memset(symb_args, 0, sizeof(ams_hash_symb_dnr_arg_t)*thread_num);
	global_row_idx_a = 0;

	for (int idx = 0; idx < thread_num; idx ++) {
		core_idx = get_cpu_id(idx);
		CPU_ZERO(&set);
		CPU_SET(core_idx, &set);
		pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);
		symb_args[idx]._tid = idx;
		symb_args[idx].csr_a = csr_a;
		symb_args[idx].csr_b = csr_b;
		symb_args[idx].row_flp = ams_bin->row_flp;
		symb_args[idx].row_nz = ams_bin->row_nz;
		symb_args[idx].ht_check = ht_check;
		rv = pthread_create(&threadpool[idx], &attr, pthread_ams_hash_symb_dnr, symb_args+idx);
		if (rv){
			printf("ERROR; return code from pthread_create() is %d\n", rv);
			exit(EXIT_FAILURE);
		}
	}

	for (int idx = 0; idx < thread_num; idx ++) {
		rv = pthread_join(threadpool[idx], NULL);

#ifdef IN_DEBUG
		printf("\t%d_symb_time: %.9f\n", idx, symb_args[idx]._symb_time);
		if (idx == thread_num-1) {
			printf("\n");
		}
#endif /* IN_DEBUG */

		if (rv) {
			printf("ERROR; return code from pthread_create() is %d\n", rv);
			exit(EXIT_FAILURE);
		}
	}
}

void ams_hash_symb_omp_dnr(csr_t *csr_a, csr_t* csr_b, my_cnt_t *ht_check, ams_bin_t *ams_bin) {
	int thread_num = THREAD_NUM;
#if defined(SMTACT) && (THREAD_NUM > SOCKET_CORE_NUM)
	if (ht_dis_flag) {
		thread_num = SOCKET_CORE_NUM;
	}
#endif /* defined(SMTACT) && (THREAD_NUM > SOCKET_CORE_NUM) */

#ifdef IN_DEBUG
	double symb_time[thread_num];
#endif /* IN_DEBUG */

	#pragma omp parallel num_threads(thread_num)
	{
#ifdef IN_DEBUG
		struct timespec time_start, time_end;
		clock_gettime(CLOCK_REALTIME, &time_start);
#endif /* IN_DEBUG */
		// --- Bind each thread to a physical core (Linux example) ---
#ifdef __linux__
		int tid = omp_get_thread_num();
		cpu_set_t cpuset;
		CPU_ZERO(&cpuset);
		// Here we assume that the OpenMP thread number maps directly to a physical core.
		CPU_SET(get_cpu_id(tid), &cpuset);
		if (pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset) != 0) {
			perror("pthread_setaffinity_np");
			// Optionally handle the error or exit.
		}
#endif	/* __linux__ */


		my_cnt_t *thr_ht_check = ht_check + tid * NEXT_POW_2(csr_b->num_col);

		#pragma omp for schedule(dynamic,DN_MIN) nowait
		for (my_cnt_t row_idx_a = 0; row_idx_a < csr_a->num_row; row_idx_a ++) {
			my_cnt_t ht_size = ams_bin->row_flp[row_idx_a];

			my_cnt_t nnz_c = 0;

			if (ht_size > 0) {
				my_cnt_t ht_cap = ams_ht_slots_for(csr_b->num_col);
				ht_size = ams_ht_slots_for(ht_size);
				ht_size = ht_size > ht_cap ? ht_cap : ht_size;
				ht_size = NEXT_POW_2(ht_size);
				memset(thr_ht_check, -1, sizeof(my_cnt_t)*ht_size);

				for (my_cnt_t col_offset_a = csr_a->row_ptrs[row_idx_a];
					col_offset_a < csr_a->row_ptrs[row_idx_a+1]; col_offset_a++) {
					my_cnt_t col_idx_a = csr_a->cols[col_offset_a];

					for (my_cnt_t col_offset_b = csr_b->row_ptrs[col_idx_a];
						col_offset_b < csr_b->row_ptrs[col_idx_a+1]; col_offset_b++) {
						my_cnt_t key = csr_b->cols[col_offset_b];
						my_cnt_t hash = (key * HASH_CONSTANT) & (ht_size - 1);

						while (true) { // loop for hash probing
							if (thr_ht_check[hash] == key) { // if the key is already inserted, it's ok
								break;
							} else if (thr_ht_check[hash] == -1) { // if the key has not been inserted yet, then it's added.
								thr_ht_check[hash] = key;
								nnz_c ++;
								break;
							} else { // linear probing: check next entry
								hash = (hash + 1) & (ht_size - 1); //hash = (hash + 1) % ht_size
							}
						}

					}
				}
			}

			ams_bin->row_nz[row_idx_a] = nnz_c;
		}


#ifdef IN_DEBUG
		clock_gettime(CLOCK_REALTIME, &time_end);
		symb_time[tid] = diff_sec(time_start, time_end);
#endif /* IN_DEBUG */
	}


#ifdef IN_DEBUG
	for (int tid = 0; tid < thread_num; tid ++) {
		printf("\t%d_symb_time: %.9f\n", tid, symb_time[tid]);
		if (tid == thread_num-1) {
			printf("\n");
		}
	}
#endif /* IN_DEBUG */
}

#endif /* AMS_DNR */



my_cnt_t pack_ams_pack_mmmv(csr_t *csr, ams_pack_t* ams_pack, 
	my_cnt_t *ps_row_coef, my_cnt_t *nnz_per_row, 
	my_cnt_t coef, numa_cfg_t *numa_cfg) {

	my_rat_t max_row_flp = 0.0;	// will cast to my_cnt_t in the end of the function

	my_cnt_t ams_nnz_cnt = MIN(AMS_MIN_NUM, AMS_SIZE/sizeof(my_cnt_t));
	ams_pack->total_cnt = csr->num_nnz*coef/ams_nnz_cnt;


	ams_pack->b_arr = (keyval_t *) alloc2_memory(
		sizeof(keyval_t) * ams_pack->total_cnt, numa_cfg
	);
	ams_pack->b_roff_maxc = (keyval_t *) alloc2_memory(
		sizeof(keyval_t) * (ams_pack->total_cnt+1), numa_cfg
	);
	ams_pack->b_roff_maxc[0].key = 0;

#ifdef SMTACT
#ifndef SMTACT_REAR 
	my_cnt_t *batch_max_flp = (my_cnt_t *) alloc2_memory(sizeof(my_cnt_t)* ams_pack->total_cnt, numa_cfg);
#endif /* SMTACT_REAR */
#endif /* SMTACT */

	int thread_num = THREAD_NUM;
	my_cnt_t static_scheduled_batch_size = ams_pack->total_cnt/thread_num;
	static_scheduled_batch_size = static_scheduled_batch_size > 0 
		? static_scheduled_batch_size : 1;

	#pragma omp parallel num_threads(THREAD_NUM)
	{
		// --- Bind each thread to a physical core (Linux example) ---
	#ifdef __linux__
		int tid = omp_get_thread_num();
		cpu_set_t cpuset;
		CPU_ZERO(&cpuset);
		// Here we assume that the OpenMP thread number maps directly to a physical core.
		CPU_SET(get_cpu_id(tid), &cpuset);
		if (pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset) != 0) {
			perror("pthread_setaffinity_np");
			// Optionally handle the error or exit.
		}
	#endif	/* __linux__ */

		#pragma omp for schedule(static, static_scheduled_batch_size)
		for (my_cnt_t ams_idx = 0; ams_idx < ams_pack->total_cnt; ams_idx ++) {
			ams_pack->b_roff_maxc[ams_idx+1].key = lower_bound(
					ps_row_coef,
					ps_row_coef + csr->num_row + 1,
					ams_nnz_cnt * (ams_idx + 1)
				) - ps_row_coef;

			ams_pack->b_arr[ams_idx].key = ams_idx;

		} // <-- implicit barrier here

		#pragma omp single
		{
			ams_pack->b_roff_maxc[ams_pack->total_cnt].key = csr->num_row;	// set the end address
		}

		my_cnt_t thr_invalid_cnt = 0;
		my_rat_t thr_max_row_flp = 0;	// will cast to my_cnt_t in the end of the function
		#pragma omp for schedule(static, static_scheduled_batch_size)
		for (my_cnt_t ams_idx = 0; ams_idx < ams_pack->total_cnt; ams_idx ++) {

			my_cnt_t start_row_offset = ams_pack->b_roff_maxc[ams_idx].key;
			my_cnt_t end_row_offset = ams_pack->b_roff_maxc[ams_idx+1].key;

			my_cnt_t offset_diff = end_row_offset - start_row_offset;

			if (offset_diff > 0) {

#ifdef AMS_MAXFLP
				ams_pack->b_arr[ams_idx].val = 
					ams_pack->b_roff_maxc[ams_idx].val =
						find_max(
							nnz_per_row+start_row_offset, 
							end_row_offset-start_row_offset
						);
#else /* AMS_MAXFLP */
				// no need to multiplied by coefficient here
				ams_pack->b_arr[ams_idx].val = (my_rat_t) 
					(
						csr->row_ptrs[end_row_offset] -
						csr->row_ptrs[start_row_offset]
					) / (end_row_offset - start_row_offset);

				ams_pack->b_roff_maxc[ams_idx].val = find_max(
					nnz_per_row+start_row_offset, 
					end_row_offset-start_row_offset
				);
#endif /* AMS_MAXFLP */

				thr_max_row_flp = 
					ams_pack->b_roff_maxc[ams_idx].val > thr_max_row_flp ?
					ams_pack->b_roff_maxc[ams_idx].val : thr_max_row_flp;
			} else {
				ams_pack->b_roff_maxc[ams_idx].val = 0.0;
				ams_pack->b_arr[ams_idx].val = 0.0;
				thr_invalid_cnt ++;
			}

#ifdef SMTACT
#ifndef SMTACT_REAR
			batch_max_flp[ams_idx] = (my_cnt_t) ams_pack->b_roff_maxc[ams_idx].val;
#endif /* SMTACT_REAR */
#endif /* SMTACT */

		} // <-- implicit barrier here

		#pragma omp atomic
#if defined(AMS_LPT) || defined(AMS_SIMP)
		ams_pack->tail_idx += thr_invalid_cnt;
#else	/* defined(AMS_LPT) || defined(AMS_SIMP) */
		ams_pack->low_head_idx += thr_invalid_cnt;  // atomically update global variable
#endif	/* defined(AMS_LPT) || defined(AMS_SIMP) */

		#pragma omp critical
		max_row_flp = max_row_flp > thr_max_row_flp ? max_row_flp : thr_max_row_flp;
	}

#ifdef SMTACT
#ifndef SMTACT_REAR
	qsort_omp_parallel(batch_max_flp, ams_pack->total_cnt, numa_cfg);
	if (batch_max_flp[ (my_cnt_t) ( (ams_pack->total_cnt-1) * SMTACT_QUANTILE) ] > cph_mf_th) {
		ht_dis_flag = true;
	}
	dealloc_memory(batch_max_flp, sizeof(my_cnt_t) * ams_pack->total_cnt);
#else /* SMTACT_REAR */
	if (max_row_flp > cph_mf_th) {
		ht_dis_flag = true;
	}
#endif /* SMTACT_REAR */
#endif /* SMTACT */


	// compare_keyvals_3, val first (ascending), then key
	p_qsortmerge_keyval(ams_pack->b_arr, ams_pack->total_cnt, compare_keyvals_3, numa_cfg);
	// ams_pack->tail_idx = 0;
	// ams_pack->low_head_idx = 0;

#if defined(AMS_LPT) || defined(AMS_SIMP)
	ams_pack->head_idx = ams_pack->total_cnt-1;
#else	/* defined(AMS_LPT) || defined(AMS_SIMP) */
	ams_pack->high_head_idx = ams_pack->total_cnt-1;
	ams_pack->high_tail_idx = (ams_pack->total_cnt - ams_pack->low_head_idx) * AMS_PRT + ams_pack->low_head_idx;
	ams_pack->low_tail_idx = ams_pack->high_tail_idx - 1;
	// assert( ams_pack->low_tail_idx > ams_pack->low_head_idx );
#endif	/* defined(AMS_LPT) || defined(AMS_SIMP) */


#ifdef IN_DEBUG
	// for (my_cnt_t ams_idx = ams_pack->low_head_idx; ams_idx < ams_pack->total_cnt; ams_idx ++) {
	// 	printf("\t%ld_ams_pack: \tbatch_index: %ld\tbatch_cf: %.6f\tbatch_max_flp: %.2f\n",
	// 		ams_idx, ams_pack->b_arr[ams_idx].key, ams_pack->b_arr[ams_idx].val,
	// 		ams_pack->b_roff_maxc[ams_pack->b_arr[ams_idx].key].val
	// 	);
	// }
	// // fflush(stdout);

	// printf("\tlow_head_idx: %ld\t"
	// 	"low_tail_idx: %ld\t"
	// 	"high_tail_idx: %ld\t"
	// 	"high_head_idx: %ld\n",
	// 	ams_pack->low_head_idx,
	// 	ams_pack->low_tail_idx, 
	// 	ams_pack->high_tail_idx,
	// 	ams_pack->high_head_idx
	// );

	// printf("ams_nnz_cnt: %ld\ttotal_cnt: %ld\n", ams_nnz_cnt, ams_pack->total_cnt);
#endif /* IN_DEBUG */

	return (my_cnt_t) max_row_flp;
}


