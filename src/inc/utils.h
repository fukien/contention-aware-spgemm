#pragma once


#ifdef __cplusplus
extern "C" {
#endif

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif /* _GNU_SOURCE */
#include <assert.h>
#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <dirent.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <linux/version.h>

#include <numa.h>
#include <numaif.h>

#include <libconfig.h>

#include <emmintrin.h>
#include <immintrin.h>

#include <libpmem.h>

#ifndef __USE_GNU
#define __USE_GNU
#endif /* __USE_GNU */
#include <sched.h>
#include <pthread.h>

#include "store.h"

// Function-like macro to convert the value to a string
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#ifndef BARRIER_ARRIVE
/** barrier wait macro */
#define BARRIER_ARRIVE(B, RV)							\
	RV = pthread_barrier_wait(B);					   	\
	if(RV !=0 && RV != PTHREAD_BARRIER_SERIAL_THREAD){  \
		printf("Couldn't wait on barrier\n");		   	\
		exit(EXIT_FAILURE);							 	\
	}
#endif /* BARRIER_ARRIVE */

#ifndef ERROR_RETURN
#define ERROR_RETURN(retval) { fprintf(stderr, "Error %d %s:line %d: \n", retval,__FILE__,__LINE__);  exit(retval); }
#endif /* ERROR_RETURN */


#ifndef MAX
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#endif /* MAX */

#ifndef MIN
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#endif /* MIN */

#define IS_POWER_OF_TWO(x) (((x) != 0) && (((x) & ((x) - 1)) == 0))
	
#define POWER_OF_TWO_TO_POWER(x) ({ \
	int _power = 0; \
	int _val = (x); /* Create a temporary variable to hold x */ \
	if (IS_POWER_OF_TWO(_val)) { \
		while (_val > 1) { \
			_val >>= 1; \
			_power++; \
		} \
	} else { \
		_power = -1; /* Return -1 if x is not a power of two */ \
	} \
	_power; \
})

#define NEXT_POW_2(x) ( (my_cnt_t) pow(2, ceil(log((x)) / log(2))))


#define SETBIT(x, y) x |= (1 << y)		// set the yth bit of x is 1
#define CLRBIT(x, y) x &= ~(1 << y)		// set the yth bit of x is 0
#define GETBIT(x, y) ((x) >> (y)&1)		// get the yth bit of x


typedef struct {
	int nr_nodes;
	int node_idx;
	struct bitmask *numa_mask;
	int *weights;
	int local_cores[SOCKET_CORE_NUM];
	int remote_cores[SOCKET_CORE_NUM];
} numa_cfg_t;


typedef struct {
	int _tid;
	void *addr;
	int ch;
	size_t size;
} pthread_par_memset_arg_t;


typedef struct {
	int _tid;
	void *dest;
	void *src;
	size_t size;
} pthread_par_memcpy_arg_t;





void green();
void cyan();



void reset();
void newdir(const char * dirpath);
size_t count_lines_num(const char * filepath);	// count based on "\n"
char * get_command_output_short(const char * command);
double diff_sec(const struct timespec start, const struct timespec end);
void mc_cfg(int* local_cores, int* remote_cores);
int get_cpu_id(int thread_id);
void * create_pmem_buffer(const char *filepath, const size_t buffer_size);
void * alloc_dram(const size_t buffer_size, const int numa_node_idx);
void * alloc_memory(const size_t buffer_size, const int numa_node_idx);
void dealloc_memory(void * addr, const size_t size);
void parallel_memset(void *addr, int ch, size_t size, int thread_num);
void parallel_memcpy(void *dest, void *src, size_t size, int thread_num);


static inline void sfence(void) {asm volatile("sfence");}
static inline void lfence(void) {asm volatile("lfence");}
static inline void mfence(void) {asm volatile("mfence");}

inline void clflush(void* addr, size_t size) {
	char *p = (char*) addr;
	uintptr_t endline = ((uintptr_t)addr + size - 1) | (CACHELINE_SIZE-1);
	// set the offset-within-line address bits to get the last byte 
	// of the cacheline containing the end of the struct.

	do {
		// flush while p is in a cache line that contains any of the struct
		_mm_clflush(p);
		p += CACHELINE_SIZE;
	} while(p <= (const char*)endline);
}

inline void clflushopt(void* addr, size_t size) {
	char *p = (char*) addr;
	uintptr_t endline = ((uintptr_t)addr + size - 1) | (CACHELINE_SIZE-1);
	// set the offset-within-line address bits to get the last byte 
	// of the cacheline containing the end of the struct.

	do {
		// flush while p is in a cache line that contains any of the struct
		_mm_clflushopt(p);
		p += CACHELINE_SIZE;
	} while(p <= (const char*)endline);
}


#define pmem_clwb(addr) asm volatile(".byte 0x66; xsaveopt %0" : "+m"(*(volatile char *)(addr)));
inline void clwb(void* addr, size_t size) {
	char *p = (char*) addr;
	uintptr_t endline = ((uintptr_t)addr + size - 1) | (CACHELINE_SIZE-1);
	// set the offset-within-line address bits to get the last byte 
	// of the cacheline containing the end of the struct.

	do {
		// flush while p is in a cache line that contains any of the struct
		// _mm_clwb(p);
		pmem_clwb(p);
		p += CACHELINE_SIZE;
	} while(p <= (const char*)endline);
}



void numa_cfg_init(numa_cfg_t *numa_cfg);
void numa_cfg_free(numa_cfg_t *numa_cfg);

#if LINUX_VERSION_CODE >= KERNEL_VERSION(6, 9, 0)
void set_node_interleaving_weight(int weight, int node_idx);
void set_interleaving_weight(numa_cfg_t* numa_cfg);
void * alloc_weighted(const size_t buffer_size, numa_cfg_t *numa_cfg);
#endif	/* LINUX_VERSION_CODE >= KERNEL_VERSION(6, 9, 0) */


void * alloc_triv_dram(const size_t buffer_size);
void dealloc_triv_dram(void *addr, const size_t buffer_size);

void * alloc2_memory(const size_t buffer_size, numa_cfg_t *numa_cfg);


#ifdef __cplusplus
}
#endif


