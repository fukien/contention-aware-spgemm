#ifndef AMS_UTILS_H
#define AMS_UTILS_H

#ifdef __cplusplus
extern "C" {
#endif

#include "../mm/common.h"


typedef struct {
	keyval_t *b_roff_maxc;	// key: start_row_offset, val: max_flp 
	keyval_t *b_arr;		// key: ptr_2_b_roff_maxc, val: ams_cf

	my_cnt_t total_cnt;
	// my_cnt_t high_cnt;
	// my_cnt_t low_cnt;

#if defined(AMS_LPT) || defined(AMS_SIMP)
	/**								
		|________________|-----------------------------|
		invalid_range	tail_idx				head_idx
	*/
	my_cnt_t head_idx;
	my_cnt_t tail_idx;
#else	/* defined(AMS_LPT) || defined(AMS_SIMP) */
	/**
	  										region d | region c 			 	region b / region a	
		|________________|--------------------------------------------------------|--------------|
		invalid_range	low_head_idx							low_tail_idx/high_tail_idx		high_head_idx
	*/
	my_cnt_t high_head_idx;
	my_cnt_t high_tail_idx;
	my_cnt_t low_head_idx;
	my_cnt_t low_tail_idx;
#endif	/* defined(AMS_LPT) || defined(AMS_SIMP) */
} ams_pack_t;


typedef struct {
	csr_t *csr;
	ams_pack_t *ams_pack;
	my_cnt_t *row_flp;
	my_cnt_t *ps_row_flp;
	my_cnt_t *row_nz;
} ams_bin_t;


typedef struct {
	csr_t *csr_a;
	csr_t *csr_b;
	csr_t *csr_c;
	keyval_t *keyval_org;
	keyval_t *aux; 
	my_cnt_t *thr_ht_check;
	my_rat_t *thr_ht_value;
} hybacc_krnl_arg_t;


void free_ams_pack(ams_pack_t *ams_pack);
void init_ams_bin(ams_bin_t *ams_bin, ams_pack_t *ams_pack, 
	csr_t *csr, numa_cfg_t *numa_cfg);
void free_ams_bin(ams_bin_t *ams_bin);
my_cnt_t pack_ams_pack(my_cnt_t *ps_row_nz, ams_bin_t *ams_bin, numa_cfg_t *numa_cfg);

/* target hash-table occupancy for the ams_* kernels; 1.0 == original sizing */
extern my_rat_t ams_ht_load;
my_cnt_t ams_ht_slots_for(my_cnt_t nz);

void ams_hash_symb(csr_t *csr_a, csr_t* csr_b, 
	my_cnt_t *ht_check, ams_bin_t *ams_bin);
void ams_hash_numc(csr_t *csr_a, csr_t* csr_b, csr_t *csr_c, 
	my_cnt_t *ht_check, my_rat_t *ht_value, ams_bin_t *ams_bin);
void ams_hsersc_numc(
	csr_t *csr_a, csr_t* csr_b, csr_t *csr_c, 
	keyval_t *keyval_org, keyval_t *aux, 
	my_cnt_t max_row_flp, ams_bin_t *ams_bin
);
void ams_hybacc_numc(
	csr_t *csr_a, csr_t* csr_b, csr_t *csr_c, 
	keyval_t *keyval_org, keyval_t *aux, my_cnt_t max_row_flp,
	my_cnt_t *ht_check, my_rat_t *ht_value, ams_bin_t *ams_bin
);


#ifdef AMS_DNR
void ams_hash_symb_dnr(csr_t *csr_a, csr_t* csr_b, my_cnt_t *ht_check, ams_bin_t *ams_bin);
void ams_hash_symb_omp_dnr(csr_t *csr_a, csr_t* csr_b, my_cnt_t *ht_check, ams_bin_t *ams_bin);
#endif /* AMS_DNR */


my_cnt_t pack_ams_pack_mmmv(csr_t *csr, ams_pack_t* ams_pack, 
	my_cnt_t *ps_row_coef, my_cnt_t *nnz_per_row, 
	my_cnt_t coef, numa_cfg_t *numa_cfg);


#ifdef __cplusplus
}
#endif

#endif /* AMS_UTILS_H */
