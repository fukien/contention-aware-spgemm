#ifndef COMMON_H
#define COMMON_H

#ifdef __cplusplus
extern "C" {
#endif

#include "mm_utils.h"



typedef struct {
	int _tid;
	my_cnt_t *rows_offset;
	csr_t *csr_a;
	csr_t *csr_b;
	my_cnt_t *local_nnz_max;
	numa_cfg_t *numa_cfg;
} comm_ptrl_init_arg_t;


typedef struct {
	int _tid;
	csr_t *csr_a;
	csr_t *csr_b;
	my_cnt_t *rows_offset;
	my_cnt_t *row_nz;
	my_cnt_t local_nnz_max;
	numa_cfg_t *numa_cfg;
#ifdef IN_VERIFY
	int64_t over_predict;
	int64_t under_predict;
#endif /* IN_VERIFY */
#ifdef IN_DEBUG
	double _symb_time;
#endif /* IN_DEBUG */
} esc_symbolic_arg_t;


typedef struct {
	// int thread_num;
	my_cnt_t total_intprod;

	my_cnt_t num_row;

	my_cnt_t max_intprod;
	my_cnt_t min_ht_size;

	// row-level
	char *bin_id;
	my_cnt_t *row_nz;

	// thread-level
	my_cnt_t *rows_offset;
	my_cnt_t ht_size_array[THREAD_NUM];
	my_cnt_t** local_hash_table_id;
	my_rat_t** local_hash_table_val;
} bin_t;


typedef struct {
	int _tid;
	my_cnt_t intprod;
	bin_t *bin;
	csr_t *csr_a;
	csr_t *csr_b;

	numa_cfg_t * __attribute__((unused)) numa_cfg;
} pthread_set_intprod_num_arg_t;


typedef struct {
	int _tid;
	my_cnt_t num_row;
	my_cnt_t average_intprod;
	my_cnt_t *ps_row_nz;
	bin_t *bin;
} pthread_hash_set_rows_offset_arg_t;


typedef struct {
	int _tid;
	my_cnt_t num_row;
	my_cnt_t num_col;
	bin_t *bin;
} pthread_set_bin_id_arg_t;


typedef struct {
	int _tid;
	numa_cfg_t* numa_cfg;
	my_cnt_t num_col;
	bin_t *bin;
} pthread_bin_create_local_hash_table_arg_t;


typedef struct {
	int _tid;
	csr_t *csr_a;
	csr_t *csr_b;
	bin_t *bin;

	// gud
	void * __attribute__((unused)) asst;

	// gud & hyb
	numa_cfg_t * __attribute__((unused)) numa_cfg;


#ifdef IN_VERIFY
	int64_t over_predict;
	int64_t under_predict;
#endif /* IN_VERIFY */
#ifdef IN_DEBUG
	double _symb_time;
#endif /* IN_DEBUG */
} hash_symbolic_arg_t;

typedef struct {
	int _tid;
	func_op_t func_op;
	csr_t *csr_a;
	csr_t *csr_b;
	csr_t *csr_c;
	bin_t *bin;

	// gud
	void * __attribute__((unused)) asst;

	// gud & hyb
	numa_cfg_t *numa_cfg;


#ifdef IN_DEBUG
	double _numc_time;
#endif /* IN_DEBUG */
#ifdef IN_EXAMINE
	double _set_time;
	double _comp_time;
	double _write_time;
	my_cnt_t _flop_cnt;
#endif /* IN_EXAMINE */
} hash_numeric_arg_t;



void comm_init_local(csr_t *csr_a, csr_t *csr_b, my_cnt_t *rows_offset, 
	my_cnt_t*local_nnz_max, numa_cfg_t* numa_cfg);

void esc_symbolic(csr_t* csr_a, csr_t* csr_b,
	my_cnt_t* rows_offset, my_cnt_t* row_nz, 
	my_cnt_t* local_nnz_max, numa_cfg_t* numa_cfg);


void init_bin(bin_t* bin, my_cnt_t num_row, my_cnt_t ht_size, numa_cfg_t *numa_cfg);
void free_bin(bin_t *bin);
void set_intprod_num(bin_t *bin, csr_t *csr_a, csr_t *csr_b);
void *pthread_hash_set_rows_offset(void *param);
void hash_set_rows_offset(bin_t *bin, my_cnt_t num_row, numa_cfg_t* numa_cfg);
void set_bin_id(bin_t *bin, my_cnt_t num_row, my_cnt_t num_col);
void bin_create_local_hash_table(bin_t *bin, my_cnt_t num_col, numa_cfg_t* numa_cfg);
void hash_symbolic(csr_t *csr_a, csr_t *csr_b, bin_t *bin);


#ifdef __cplusplus
}
#endif

#endif /* COMMON_H */















