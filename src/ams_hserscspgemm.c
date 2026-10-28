#include "ams/ams_utils.h"

extern numa_cfg_t numa_cfg;
extern int64_t file_idx;

#ifdef SMTACT
extern my_cnt_t cph_mf_th;
extern bool ht_dis_flag;
#endif /* SMTACT */


#ifdef IN_DEBUG
extern pthread_barrier_t debug_barrier;
extern int debug_rv;
#endif /* IN_DEBUG */

int local_numa_idx;
int cxl_numa_idx;
int cur_numa_idx;

#ifdef MEMSV
numa_cfg_t cxl_numa_cfg;
#endif /* MEMSV */


#ifdef IN_VERIFY
extern int64_t over_predict;
extern int64_t under_predict;
#endif /* IN_VERIFY */

#ifdef IN_INSP
extern double *_ams_numc_time_array;
#endif /* IN_INSP */


static runtime_t hserscspgemm(csr_t* csr_a, csr_t* csr_b, csr_t* csr_c, numa_cfg_t* numa_cfg) {
	struct timespec time_start, time_end;
	runtime_t runtime;

	// INIT
	clock_gettime(CLOCK_REALTIME, &time_start);

	csr_c->num_row = csr_a->num_row;
	csr_c->num_col = csr_b->num_col;
	csr_c->row_ptrs = (my_cnt_t*) alloc2_memory(
		sizeof(my_cnt_t)*(csr_c->num_row+1), numa_cfg
	);
	memset(csr_c->row_ptrs, 0, sizeof(my_cnt_t)*(csr_c->num_row+1));


	ams_pack_t ams_pack;
	memset(&ams_pack, 0, sizeof(ams_pack_t));
	ams_bin_t ams_bin;
	init_ams_bin(&ams_bin, &ams_pack, csr_a, numa_cfg);
	my_cnt_t __attribute__((unused)) max_row_flp = 0;

	flp_populate(csr_a, csr_b, ams_bin.row_flp);

	scan(ams_bin.ps_row_flp, ams_bin.row_flp, csr_a->num_row+1);


	typedef void (*symb_func_t)(csr_t *, csr_t *, my_cnt_t*, ams_bin_t*);
	symb_func_t symb_func;

#ifdef SMTACT
	ht_dis_flag = false;
	cph_mf_th = SMTACT_SYMB_THRES;
#endif /* SMTACT */

#ifdef AMS_DNR
	symb_func = ams_hash_symb_dnr;
#ifdef SMTACT
	max_row_flp = pack_ams_pack(csr_a->row_ptrs, &ams_bin, numa_cfg);
#endif /* SMTACT */
#else /* AMS_DNR */
	free_ams_pack(&ams_pack);
	max_row_flp = pack_ams_pack(csr_a->row_ptrs, &ams_bin, numa_cfg);
	symb_func = ams_hash_symb;
#endif /* AMS_DNR */

	clock_gettime(CLOCK_REALTIME, &time_end);
	runtime.init_time = diff_sec(time_start, time_end);

	// SYMBOLIC
	clock_gettime(CLOCK_REALTIME, &time_start);

	my_cnt_t *ht_check = (my_cnt_t *) alloc2_memory(
		sizeof(my_cnt_t)*THREAD_NUM*NEXT_POW_2(csr_c->num_col), numa_cfg
	);


	// 2nd pass, csr_a->row_ptrs, csr_a->cols, csr_b->row_ptrs, csr_b->cols
#if defined(MEMSV) && defined(MEMSVB)
	numa_cfg_t intr_numa_cfg;
	numa_cfg_init(&intr_numa_cfg);
	numa_bitmask_setbit(intr_numa_cfg.numa_mask, LOCAL_NUMA_NODE_IDX-1);
	numa_bitmask_setbit(intr_numa_cfg.numa_mask, CXL_NUMA_NODE_IDX-1);
	csr_t intr_csr_b;
	intr_csr_b.num_row = csr_b->num_row;
	intr_csr_b.num_col = csr_b->num_col;
	intr_csr_b.num_nnz = csr_b->num_nnz;
	intr_csr_b.row_ptrs = (my_cnt_t *) alloc_weighted(
		sizeof(my_cnt_t) * (intr_csr_b.num_row+1), &intr_numa_cfg
	);
	intr_csr_b.cols = (my_cnt_t *) alloc_weighted(
		sizeof(my_cnt_t) * intr_csr_b.num_nnz, &intr_numa_cfg
	);
	parallel_memcpy(
		intr_csr_b.row_ptrs, csr_b->row_ptrs, 
		sizeof(my_cnt_t)*(csr_b->num_row+1), SOCKET_CORE_NUM
	);
	parallel_memcpy(
		intr_csr_b.cols, csr_b->cols, 
		sizeof(my_cnt_t)*csr_b->num_nnz, SOCKET_CORE_NUM
	);
	symb_func(csr_a, &intr_csr_b, ht_check, &ams_bin);
	dealloc_memory(intr_csr_b.row_ptrs, sizeof(my_cnt_t)*(intr_csr_b.num_row+1));
	dealloc_memory(intr_csr_b.cols, sizeof(my_cnt_t)*intr_csr_b.num_nnz);
	intr_csr_b.num_row = intr_csr_b.num_col = intr_csr_b.num_nnz = 0;
	numa_cfg_free(&intr_numa_cfg);
#else	/* defined(MEMSV) && defined(MEMSVB) */
	symb_func(csr_a, csr_b, ht_check, &ams_bin);
#endif /* defined(MEMSV) && defined(MEMSVB) */ 

	scan(csr_c->row_ptrs, ams_bin.row_nz, csr_a->num_row+1);

	csr_c->num_nnz = csr_c->row_ptrs[csr_c->num_row];

	clock_gettime(CLOCK_REALTIME, &time_end);
	runtime.symb_time = diff_sec(time_start, time_end);

#ifdef IN_STATS
	printf("\t%ld\t%ld\t%ld\t%ld\t%ld\t%f\n", 
		csr_a->num_row, csr_a->num_col, csr_a->num_nnz,
		csr_c->num_nnz, ams_bin.ps_row_flp[csr_c->num_row], 
		(double) ams_bin.ps_row_flp[csr_c->num_row]/csr_c->num_nnz
	);
	exit(1);
#endif /* IN_STATS */

	dealloc_memory(ht_check, sizeof(my_cnt_t)*THREAD_NUM*NEXT_POW_2(csr_c->num_col));

	// MANAGEMENT
	clock_gettime(CLOCK_REALTIME, &time_start);

#ifdef MEMSV
	csr_c->cols = (my_cnt_t*) alloc_dram(
		sizeof(my_cnt_t)*(csr_c->num_nnz), cxl_numa_cfg.node_idx
	);
	csr_c->rats = (my_rat_t*) alloc_dram(
		sizeof(my_rat_t)*csr_c->num_nnz, cxl_numa_cfg.node_idx
	);
#else /* MEMSV */
	csr_c->cols = (my_cnt_t*) alloc2_memory(
		sizeof(my_cnt_t)*(csr_c->num_nnz), numa_cfg
	);
	csr_c->rats = (my_rat_t*) alloc2_memory(
		sizeof(my_rat_t)*csr_c->num_nnz, numa_cfg
	);
#endif /* MEMSV */


#ifdef CSRC_PREPOPULATED
	// memset(csr_c->cols, 0, sizeof(my_cnt_t)*csr_c->num_nnz);
	// memset(csr_c->rats, 0, sizeof(my_rat_t)*csr_c->num_nnz);
	parallel_memset(csr_c->cols, 0, sizeof(my_cnt_t)*csr_c->num_nnz, SOCKET_CORE_NUM);
	parallel_memset(csr_c->rats, 0, sizeof(my_rat_t)*csr_c->num_nnz, SOCKET_CORE_NUM);
#endif /* CSRC_PREPOPULATED */

	// CAREFUL HERE
	free_ams_pack(&ams_pack);

#ifdef SMTACT
	ht_dis_flag = false;
	cph_mf_th = SMTACT_NUMC_THRES;
#endif /* SMTACT */
	max_row_flp = pack_ams_pack(csr_c->row_ptrs, &ams_bin, numa_cfg);
	// printf("%s:%d:\tcph_mf_th: %ld\tht_dis_flag: %d\n", 
	// 	__FILE__, __LINE__,
	// 	cph_mf_th, ht_dis_flag
	// );

	clock_gettime(CLOCK_REALTIME, &time_end);
	runtime.mgmt_time = diff_sec(time_start, time_end);

	// NUMERIC
	clock_gettime(CLOCK_REALTIME, &time_start);


	keyval_t *keyval_org = (keyval_t *) alloc2_memory(
		sizeof(keyval_t) * THREAD_NUM * max_row_flp, numa_cfg
	);
	keyval_t *aux = (keyval_t *) alloc2_memory(
		sizeof(keyval_t) * THREAD_NUM * max_row_flp, numa_cfg
	);

	// 3rd pass, csr_a->row_ptrs, csr_a->cols, csr_a->rats, csr_b->row_ptrs, csr_b->cols, csr_b->rats
	ams_hsersc_numc(csr_a, csr_b, csr_c, keyval_org, aux, max_row_flp, &ams_bin);

	dealloc_memory(keyval_org, sizeof(keyval_t) * THREAD_NUM * max_row_flp);
	dealloc_memory(aux, sizeof(keyval_t) * THREAD_NUM * max_row_flp);


	free_ams_pack(&ams_pack);
	free_ams_bin(&ams_bin);

	clock_gettime(CLOCK_REALTIME, &time_end);
	runtime.numc_time = diff_sec(time_start, time_end);

	runtime.total_time = runtime.init_time+runtime.symb_time+runtime.mgmt_time+runtime.numc_time;

	return runtime;
}


int main(int argc, char** argv) {
	srand(32768);

#ifdef IN_DEBUG
	debug_init_barrier();
#endif /* IN_DEBUG */

	char matrix_org_path[CHAR_BUFFER_LEN];
	char matrix_trs_path[CHAR_BUFFER_LEN];
	assert(argc == 2);
	char *dataset = argv[1];
	get_matrix_paths(matrix_org_path, matrix_trs_path, dataset);

	cfg_print();
#ifdef IN_STATS
	printf("%s", dataset);
#endif /* IN_STATS */

	/*** NUMA initialization ***/
	numa_cfg_init(&numa_cfg);

	local_numa_idx = numa_cfg.local_cores[0]/SOCKET_CORE_NUM;
	cxl_numa_idx = 2;
	cur_numa_idx = local_numa_idx;
	cur_numa_idx = cxl_numa_idx;
#if IS_POWER_OF_TWO(NUMA_MASK_VAL)
	numa_cfg.node_idx = cur_numa_idx = POWER_OF_TWO_TO_POWER(NUMA_MASK_VAL);
#endif /* IS_POWER_OF_TWO(NUMA_MASK_VAL) */
	numa_cfg.node_idx = cur_numa_idx;

#if defined(USE_WEIGHTED_INTERLEAVING) && LINUX_VERSION_CODE >= KERNEL_VERSION(6, 9, 0)
	numa_bitmask_clearall(numa_cfg.numa_mask);
	for (int node_idx = 0; node_idx < numa_cfg.nr_nodes; node_idx ++) {
		if ( GETBIT( NUMA_MASK_VAL, node_idx ) ) {
			numa_bitmask_setbit(numa_cfg.numa_mask, node_idx);
		}
	}
	// numa_bitmask_setbit(numa_cfg.numa_mask, 1);
	// numa_bitmask_setbit(numa_cfg.numa_mask, 2);
	numa_cfg.weights[0] = 6;
	numa_cfg.weights[1] = 3;
	numa_cfg.weights[2] = 2;
	// set_interleaving_weight(&numa_cfg);
#elif !defined(USE_INTERLEAVING)	/* defined(USE_WEIGHTED_INTERLEAVING) && LINUX_VERSION_CODE >= KERNEL_VERSION(6, 9, 0) */
	numa_bitmask_clearall(numa_cfg.numa_mask);
	numa_bitmask_setbit(numa_cfg.numa_mask, cur_numa_idx);
	numa_cfg.node_idx = cur_numa_idx;
#endif /* defined(USE_WEIGHTED_INTERLEAVING) && LINUX_VERSION_CODE >= KERNEL_VERSION(6, 9, 0) */


#ifdef MEMSV
	numa_cfg_init(&cxl_numa_cfg);
	cxl_numa_cfg.node_idx = cxl_numa_idx;
	numa_bitmask_clearall(cxl_numa_cfg.numa_mask);
	numa_bitmask_setbit(cxl_numa_cfg.numa_mask, cxl_numa_cfg.node_idx);
	// printf("cxl_numa_cfg.numa_mask: %ld\n", 
	// 	*(cxl_numa_cfg.numa_mask->maskp)
	// );
#endif /* MEMSV */

#ifdef USE_NVM
	newdir(NVM_LOCAL_PREFIX);
	newdir(NVM_REMOTE_PREFIX);
#endif /* USE_NVM */

	csr_t csr_a, csr_b, csr_c;
	memset(&csr_a, 0, sizeof(csr_t));
	memset(&csr_b, 0, sizeof(csr_t));
	memset(&csr_c, 0, sizeof(csr_t));

#ifdef MEMSV
	raw2csr(&csr_a, matrix_org_path, &cxl_numa_cfg);
#else /* MEMSV */
	raw2csr(&csr_a, matrix_org_path, &numa_cfg);
#endif /* MEMSV */
	raw2csr(&csr_b, matrix_trs_path, &numa_cfg);


#ifdef IN_INSP
	_ams_numc_time_array = (double *) alloc2_memory(sizeof(double) * csr_a.num_row, &numa_cfg);
	memset(_ams_numc_time_array, 0, sizeof(double) * csr_a.num_row);
#endif /* IN_INSP */

#ifdef MEM_MON
	char command[CHAR_BUFFER_LEN] = {'\0'};
	snprintf(
		command, sizeof(command),
		"bash %s/mem_mon.sh %d %s < /dev/null &",
		PROJECT_PATH, getpid(), dataset
	);
	system(command);
#endif /* MEM_MON */

	runtime_t global_runtime = {0.0, 0.0, 0.0, 0.0, 0.0};
	for (int idx = 0; idx < RUN_NUM; idx ++) {
		runtime_t runtime = hserscspgemm(&csr_a, &csr_b, &csr_c, &numa_cfg);

#ifdef IN_VERIFY
		if (idx > 0) {
			int64_t sum_col = 0; double sum_rat = 0.0;
			for (my_cnt_t jdx = 0; jdx < csr_c.num_nnz; jdx ++) {
				sum_col += csr_c.cols[jdx];
				sum_rat += csr_c.rats[jdx];
			}
			printf("VERIFYRESULT:\t%ld\t%ld\t%f\n", csr_c.num_nnz, sum_col, sum_rat);
		}
#endif /* IN_VERIFY */


#ifdef IN_INSP
		if (idx > 0 && idx == RUN_NUM - 1) {
			for (my_cnt_t row_idx = 0; row_idx < csr_a.num_row; row_idx ++) {
				my_cnt_t flop = 0;
				for (my_cnt_t col_offset = csr_a.row_ptrs[row_idx]; col_offset < csr_a.row_ptrs[row_idx+1]; col_offset ++) {
					my_cnt_t col_idx = csr_a.cols[col_offset];
					flop += csr_b.row_ptrs[col_idx+1] - csr_b.row_ptrs[col_idx];
				}

				printf("\t\trow_idx: %ld\tnumc_time: %f\ta_nnz: %ld\tc_nnz: %ld\tflop: %ld\n", 
					// note that should be divided by RUN_NUM for RUN_0 also accumulated
					row_idx, _ams_numc_time_array[row_idx]/RUN_NUM,
					csr_a.row_ptrs[row_idx+1] - csr_a.row_ptrs[row_idx],
					csr_c.row_ptrs[row_idx+1] - csr_c.row_ptrs[row_idx],
					flop
				);
			}
		}
#endif /* IN_INSP */


		free_csr(&csr_c);
		if (idx > 0) {
			global_runtime.init_time += runtime.init_time;
			global_runtime.symb_time += runtime.symb_time;
			global_runtime.mgmt_time += runtime.mgmt_time;
			global_runtime.numc_time += runtime.numc_time;
			global_runtime.total_time += runtime.total_time;
		}

		cyan();
		printf("%dth run\tinit_time: %.9f\tsymb_time: %.9f\tmgmt_time: %.9f\tnumc_time: %.9f\ttotal_time: %.9f\n", 
			idx, runtime.init_time, runtime.symb_time, runtime.mgmt_time, runtime.numc_time, runtime.total_time
		);
		reset();
		printf("\n");
	}

#if RUN_NUM > 1
	my_cnt_t *row_flop = (my_cnt_t *) alloc2_memory(sizeof(my_cnt_t)*csr_a.num_row, &numa_cfg);
	my_cnt_t __attribute__((unused)) total_flop_cnt = get_intprod(&csr_a, &csr_b, row_flop);
	dealloc_memory(row_flop, sizeof(my_cnt_t)*csr_a.num_row);

	green();
	printf("dataset: %s\tnum_trials: %d\tthread_num: %d\tnuma_mask: %ld\t"
		"init_time: %.9f\tsymb_time: %.9f\tmgmt_time: %.9f\tnumc_time: %.9f\ttotal_time: %.9f\t"
		"flop_cnt: %ld\tnumc_perf: %.9f\tperf: %.9f\n",
		dataset, RUN_NUM-1, THREAD_NUM, *(numa_cfg.numa_mask->maskp),
		global_runtime.init_time/(RUN_NUM-1),
		global_runtime.symb_time/(RUN_NUM-1),
		global_runtime.mgmt_time/(RUN_NUM-1),
		global_runtime.numc_time/(RUN_NUM-1),
		global_runtime.total_time/(RUN_NUM-1),
		total_flop_cnt,
		total_flop_cnt/(global_runtime.numc_time/(RUN_NUM-1)),
		total_flop_cnt/(global_runtime.total_time/(RUN_NUM-1))
	);
	reset();
#endif /* RUN_NUM > 1 */


#ifdef IN_INSP
	dealloc_memory(_ams_numc_time_array, sizeof(double) * csr_a.num_row);
#endif /* IN_INSP */


	free_csr(&csr_a);
	free_csr(&csr_b);

#ifdef MEMSV
	numa_cfg_free(&cxl_numa_cfg);
#endif /* MEMSV */
	numa_cfg_free(&numa_cfg);

	return 0;
}











