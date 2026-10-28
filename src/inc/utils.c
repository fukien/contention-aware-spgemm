#include "utils.h"

/* Linux 6.9 UAPI value; older numactl headers may not expose this mode.
 * An older kernel rejects it via mbind when weighted interleaving is used. */
#ifndef MPOL_WEIGHTED_INTERLEAVE
#define MPOL_WEIGHTED_INTERLEAVE 6
#endif

numa_cfg_t numa_cfg;
size_t file_idx = 0;

void green() {
	printf("\033[0;32m");
}

void cyan() {
	printf("\033[0;36m");
}

void reset() {
	printf("\033[0m");
}

void newdir(const char * dirpath) {
	char command[CHAR_BUFFER_LEN] = {'\0'};
	strcat(strcpy(command, "mkdir -p "), dirpath);
	int val = system(command);
	if (val != 0) {
		printf("\"%s\" failed, errno: %d\n", command, errno);
		perror(command);
	}
}

size_t count_lines_num(const char * filepath) {		// count based on "\n"
	size_t num = 0;
	int fd;
	char *buffer;
	ssize_t bytes_read, i;

	// Open the file with O_DIRECT flag
	fd = open(filepath, O_RDONLY | O_DIRECT);
	if (fd == -1) {
		perror("Error opening file");
		return 0;
	}

	size_t ALIGNMENT = 4096, BUFFER_SIZE = 4096;
	// Allocate aligned memory for O_DIRECT
	if (posix_memalign((void **)&buffer, ALIGNMENT, BUFFER_SIZE)) {
		perror("posix_memalign failed");
		close(fd);
		return 0;
	}

	// Read file and count newline characters
	while ((bytes_read = read(fd, buffer, BUFFER_SIZE)) > 0) {
		for (i = 0; i < bytes_read; i++) {
			if (buffer[i] == '\n') {
				num++;
			}
		}
	}

	if (bytes_read == -1) {
		perror("Error reading file");
	}

	// Cleanup
	free(buffer);
	close(fd);
	return num;
}

char *get_command_output_short(const char * command) {
	FILE *fp;
	fp = popen(command, "r");
	if (fp == NULL) {
		printf("Failed to run command: %s, errno: %d\n", command, errno);
		exit(1);
	}

	char path[CHAR_BUFFER_LEN] = {'\0'};	
	/* Read the output a line at a time - output it. */
	while (fgets(path, sizeof(path), fp) == NULL) {
		printf("get command \"%s\" output failed, errno: %d", command, errno);
	}

	/* close */
	pclose(fp);

	size_t malloc_len = strlen(path) + 1;	// +1 to accomodate for the null terminator
	char *res = (char*) malloc( (malloc_len) * sizeof(char) );
	memset(res, '\0', (malloc_len) * sizeof(char) );
	strcpy(res, path);

	return res;
}


void mc_cfg(int *local_cores, int *remote_cores) {
	config_t cfg;
	config_init(&cfg);

	if(! config_read_file(&cfg, CFG_PATH)) {
		fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg),
			config_error_line(&cfg), config_error_text(&cfg));
		config_destroy(&cfg);
		exit(1);
	}

	int length;
	config_setting_t* setting;

	setting = config_lookup(&cfg, "cores.local_core_id_arrays");
	length = config_setting_length(setting);
	if(setting != NULL) {
		for (int idx = 0; idx < length; idx ++) {
			local_cores[idx] = config_setting_get_int_elem(setting, idx);
		}
	}

	setting = config_lookup(&cfg, "cores.remote_core_id_arrays");
	length = config_setting_length(setting);
	if(setting != NULL) {
		for (int idx = 0; idx < length; idx ++) {
			remote_cores[idx] = config_setting_get_int_elem(setting, idx);
		}
	}

	config_destroy(&cfg);
}


int get_cpu_id(int thread_id) {
#ifdef USE_NUMA
	if (thread_id < SOCKET_CORE_NUM) {
		return numa_cfg.remote_cores[thread_id%SOCKET_CORE_NUM];
	} else {
#ifdef USE_HYPERTHREADING
		return numa_cfg.remote_cores[thread_id%SOCKET_CORE_NUM]+CORE_NUM; 
#else /* USE_HYPERTHREADING */
		return numa_cfg.remote_cores[thread_id%SOCKET_CORE_NUM];
#endif  /* USE_HYPERTHREADING */
#else /* USE_NUMA */
	if (thread_id < SOCKET_CORE_NUM) {
		return numa_cfg.local_cores[thread_id%SOCKET_CORE_NUM];
	} else {
#ifdef USE_HYPERTHREADING
		return numa_cfg.local_cores[thread_id%SOCKET_CORE_NUM]+CORE_NUM; 
#else /* USE_HYPERTHREADING */
		return numa_cfg.local_cores[thread_id%SOCKET_CORE_NUM];
#endif  /* USE_HYPERTHREADING */
#endif /* USE_NUMA */
	}
}


double diff_sec(const struct timespec start, const struct timespec end) {
	return (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) * 1e-9;
}

void * create_pmem_buffer(const char *filepath, const size_t buffer_size) {
	void *pmemaddr;
	int is_pmem = 0;
	size_t mapped_len = 0;

	int res = access(filepath, F_OK);
	if (res == 0) {
		pmemaddr = (char*) pmem_map_file(filepath, buffer_size, 
			PMEM_FILE_CREATE, 0666, &mapped_len, &is_pmem);
	} else {
		// printf("create mapped filepath: %s\n", filepath);
		pmemaddr = (char*) pmem_map_file(filepath, buffer_size,
			PMEM_FILE_CREATE|PMEM_FILE_EXCL, 0666, &mapped_len, &is_pmem);
	}

	if (pmemaddr == NULL) {
		printf("%s pmem_map_file failed, errno: %d\n", filepath, errno);
		printf("arguments, buffer_size: %zu, mapped_len: %zu, is_pmem: %d\n", buffer_size, mapped_len, is_pmem);
		perror("pmem_map_file failed");
		exit(1);
	}

	// if (is_pmem) {
	// 	// printf("successfully mapped\n");
	// } else {
	// 	printf("%s error, not pmem, is_pmem: %d", filepath, is_pmem);
	// 	exit(1);
	// }

	if (!is_pmem) {
		// printf("%s error, not pmem, is_pmem: %d", filepath, is_pmem);
		exit(1);
	}

	return pmemaddr;
}


void * alloc_dram(const size_t buffer_size, const int numa_node_idx) { 
	struct bitmask *local_numa_mask = numa_bitmask_alloc(numa_cfg.nr_nodes);
	numa_bitmask_clearall(local_numa_mask);
	numa_bitmask_setbit(local_numa_mask, numa_node_idx);
	void *addr = mmap(NULL, buffer_size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if (addr == MAP_FAILED) {
		perror("mmap");
		printf("errno: %d\n", errno);
		printf("buffer_size: %zu\n", buffer_size);
		exit (1);
	}
	int rv = mbind(addr, buffer_size, MPOL_BIND, local_numa_mask->maskp, local_numa_mask->size+1, 0);

	if (rv < 0) {
		perror("mbind");
		exit (1);
	}

#ifdef USE_HUGE
	rv = madvise(addr, buffer_size, MADV_HUGEPAGE);
#else /* USE_HUGE */
	rv = madvise(addr, buffer_size, MADV_NOHUGEPAGE);
#endif /* USE_HUGE */

	if (rv < 0) {
		perror("madvise");
		exit (1);
	}

	numa_bitmask_free(local_numa_mask);

	return addr;
}




void * alloc_memory(const size_t buffer_size, const int numa_node_idx) {
#ifdef USE_NVM
	char file_path[ CHAR_BUFFER_LEN * 2 ];
	if (numa_node_idx == 1) {
		snprintf(file_path, CHAR_BUFFER_LEN * 2, "%s/%lu_%zu.bin", NVM_LOCAL_PREFIX, pthread_self(), file_idx ++);
	} else if (numa_node_idx == 0) {
		snprintf(file_path, CHAR_BUFFER_LEN * 2, "%s/%lu_%zu.bin", NVM_REMOTE_PREFIX, pthread_self(), file_idx ++);
	} else {
		exit(1);
	}
	return create_pmem_buffer(file_path, buffer_size);
#else /* USE_NVM */
	return alloc_dram(buffer_size, numa_node_idx);
#endif /* USE_NVM */
}


void dealloc_memory(void * addr, const size_t size) {
#ifdef USE_NVM
	pmem_unmap(addr, size);
#else /* USE_NVM */
	munmap(addr, size);
#endif /* USE_NVM */
}

void *pthread_par_memset(void *param) {
	pthread_par_memset_arg_t *arg = (pthread_par_memset_arg_t *) param;
	memset(arg->addr, arg->ch, arg->size);
	return NULL;
}


void parallel_memset(void *addr, int ch, size_t size, int thread_num) {
	pthread_t threadpool[thread_num];
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	cpu_set_t set;
	int core_idx;
	int rv;

	pthread_par_memset_arg_t args[thread_num];
	memset(args, 0, sizeof(pthread_par_memset_arg_t)*thread_num);

	for (int idx = 0; idx < thread_num; idx ++) {
		core_idx = get_cpu_id(idx);
		CPU_ZERO(&set);
		CPU_SET(core_idx, &set);
		pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);
		args[idx]._tid = idx;
		args[idx].addr = addr + idx * (size/thread_num);
		args[idx].ch = ch;
		args[idx].size = size/thread_num;
		if (idx == thread_num -1) {
			args[idx].size = size - size/thread_num * idx;
		}
		rv = pthread_create(&threadpool[idx], &attr, pthread_par_memset, args+idx);
		if (rv){
			printf("ERROR; return code from pthread_create() is %d\n", rv);
			exit(EXIT_FAILURE);
		}
	}
	for (int idx = 0; idx < thread_num; idx ++) {
		rv = pthread_join(threadpool[idx], NULL);
		if (rv) {
			printf("ERROR; return code from pthread_join() is %d\n", rv);
			exit(EXIT_FAILURE);
		}
	}
	pthread_attr_destroy(&attr);
}


void *pthread_par_memcpy(void *param) {
	pthread_par_memcpy_arg_t *arg = (pthread_par_memcpy_arg_t *) param;
	memcpy(arg->dest, arg->src, arg->size);
	return NULL;
}


void parallel_memcpy(void *dest, void *src, size_t size, int thread_num) {
	pthread_t threadpool[thread_num];
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	cpu_set_t set;
	int core_idx;
	int rv;

	pthread_par_memcpy_arg_t args[thread_num];
	memset(args, 0, sizeof(pthread_par_memcpy_arg_t)*thread_num);

	for (int idx = 0; idx < thread_num; idx ++) {
		core_idx = get_cpu_id(idx);
		CPU_ZERO(&set);
		CPU_SET(core_idx, &set);
		pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);
		args[idx]._tid = idx;
		args[idx].dest = dest + idx * (size/thread_num);
		args[idx].src = src + idx * (size/thread_num);;
		args[idx].size = size/thread_num;
		if (idx == thread_num -1) {
			args[idx].size = size - size/thread_num * idx;
		}
		rv = pthread_create(&threadpool[idx], &attr, pthread_par_memcpy, args+idx);
		if (rv){
			printf("ERROR; return code from pthread_create() is %d\n", rv);
			exit(EXIT_FAILURE);
		}
	}
	for (int idx = 0; idx < thread_num; idx ++) {
		rv = pthread_join(threadpool[idx], NULL);
		if (rv) {
			printf("ERROR; return code from pthread_join() is %d\n", rv);
			exit(EXIT_FAILURE);
		}
	}
	pthread_attr_destroy(&attr);
}


void numa_cfg_init(numa_cfg_t *numa_cfg) {
	if ( numa_available() < 0 ) {
		fprintf(stderr, "NUMA not available on this system.\n");
		exit(1);
	}

	memset(numa_cfg, 0, sizeof(numa_cfg_t));

	numa_cfg->nr_nodes = numa_num_configured_nodes();
	numa_cfg->numa_mask = numa_bitmask_alloc(numa_cfg->nr_nodes);
	numa_bitmask_clearall(numa_cfg->numa_mask);

	numa_cfg->node_idx = LOCAL_NUMA_NODE_IDX - 1;	// by default, we set it to local memory node

	numa_cfg->weights = malloc(numa_cfg->nr_nodes * sizeof(float));
	memset(numa_cfg->weights, 0, sizeof(int)*numa_cfg->nr_nodes);

// 	for (int node_idx = 0; node_idx < numa_cfg->nr_nodes; node_idx++) {
// 		numa_cfg->weights[node_idx] = 1; 	// equal distribution by default (set as 1, by default)
// #if defined(USE_WEIGHTED_INTERLEAVING) && LINUX_VERSION_CODE >= KERNEL_VERSION(6, 9, 0)
// 		set_node_interleaving_weight(1, node_idx);
// #endif	/* defined(USE_WEIGHTED_INTERLEAVING) && LINUX_VERSION_CODE >= KERNEL_VERSION(6, 9, 0) */
// 	}

	mc_cfg(numa_cfg->local_cores, numa_cfg->remote_cores);
}


void numa_cfg_free(numa_cfg_t *numa_cfg) {
// #if defined(USE_WEIGHTED_INTERLEAVING) && LINUX_VERSION_CODE >= KERNEL_VERSION(6, 9, 0)
// 	for (int node_idx = 0; node_idx < numa_cfg->nr_nodes; node_idx++) {
// 		set_node_interleaving_weight(1, node_idx);
// 	}
// #endif	/* defined(USE_WEIGHTED_INTERLEAVING) && LINUX_VERSION_CODE >= KERNEL_VERSION(6, 9, 0) */

	numa_bitmask_free(numa_cfg->numa_mask);
	free(numa_cfg->weights);
} 


#if LINUX_VERSION_CODE >= KERNEL_VERSION(6, 9, 0)

void set_node_interleaving_weight(int weight, int node_idx) {
	/**
	char command[CHAR_BUFFER_LEN];
	snprintf(command, CHAR_BUFFER_LEN, "%s%d%s%d%s", 
		"sh -c 'echo ", weight, 
		" > /sys/kernel/mm/mempolicy/weighted_interleave/node", node_idx, "'"
	 );
	printf("%s\n", command);

	int return_value = 0;;
	return_value = system(command);

	// Check the result of the system call
	if (return_value == -1) {
		printf("Error executing command, set node %d with interleaving weight %d\n", node_idx, weight);
	}
	*/

	char node_path[CHAR_BUFFER_LEN];
	snprintf(node_path, CHAR_BUFFER_LEN, "%s%d", "/sys/kernel/mm/mempolicy/weighted_interleave/node", node_idx);

	int fd;
	char buffer[16];  // Buffer to hold the integer as a string
	ssize_t bytes_written;
	int length;

	// Open the file with write-only permissions
	fd = open(node_path, O_WRONLY);
	if (fd == -1) {
		perror("Error opening the file");
		exit(1);
	}

	// Convert the integer weight to a string
	length = snprintf(buffer, sizeof(buffer), "%d\n", weight);
	if (length < 0) {
		perror("Error converting integer to string");
		close(fd);
		exit(1);
	}

	// Write the string to the file
	bytes_written = write(fd, buffer, length);
	if (bytes_written == -1) {
		perror("Error writing to the file");
		close(fd);
		exit(1);
	}

	// Close the file descriptor
	if (close(fd) == -1) {
		perror("Error closing the file");
		exit(1);
	}
}


void set_interleaving_weight(numa_cfg_t* numa_cfg) {
	for (int node_idx = 0; node_idx < numa_cfg->nr_nodes; node_idx ++) {
		if ( numa_bitmask_isbitset(numa_cfg->numa_mask, node_idx) ) {
			set_node_interleaving_weight(numa_cfg->weights[node_idx], node_idx);
		}
	}
}


void * alloc_weighted(const size_t buffer_size, numa_cfg_t *numa_cfg) {

#if LINUX_VERSION_CODE >= KERNEL_VERSION(6, 9, 0)
	numa_cfg_t local_numa_cfg = *numa_cfg;

	void *addr = mmap(NULL, buffer_size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if (addr == MAP_FAILED) {
		perror("mmap");
		printf("errno: %d\n", errno);
		printf("buffer_size: %zu\n", buffer_size);
		exit (1);
	}

	int rv = mbind(
		addr, 
		buffer_size, 
		MPOL_WEIGHTED_INTERLEAVE, 
		local_numa_cfg.numa_mask->maskp, 
		local_numa_cfg.numa_mask->size+1, 
		0
	);
	if (rv < 0) {
		perror("mbind");
		exit (1);
	}

#ifdef USE_HUGE
	rv = madvise(addr, buffer_size, MADV_HUGEPAGE);
#else /* USE_HUGE */
	rv = madvise(addr, buffer_size, MADV_NOHUGEPAGE);
#endif /* USE_HUGE */

	if (rv < 0) {
		perror("madvise");
		exit (1);
	}

	return addr;
#else	/* LINUX_VERSION_CODE >= KERNEL_VERSION(6, 9, 0) */
	return alloc_triv_dram(buffer_size);
#endif	/* LINUX_VERSION_CODE >= KERNEL_VERSION(6, 9, 0) */
}

#endif	/* LINUX_VERSION_CODE >= KERNEL_VERSION(6, 9, 0) */

void * alloc_triv_dram(const size_t buffer_size) {
	void * addr = mmap(NULL, buffer_size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
#ifdef USE_HUGE
	madvise(addr, buffer_size, MADV_HUGEPAGE);
#else /* USE_HUGE */
	madvise(addr, buffer_size, MADV_NOHUGEPAGE);
#endif /* USE_HUGE */
	return addr;
}

void dealloc_triv_dram(void *addr, const size_t buffer_size) {
	munmap(addr, buffer_size);
}

void * alloc2_memory(const size_t buffer_size, numa_cfg_t *numa_cfg __attribute__((unused)) ) {
#ifdef USE_INTERLEAVING
	return alloc_triv_dram(buffer_size);
#elif USE_WEIGHTED_INTERLEAVING && LINUX_VERSION_CODE >= KERNEL_VERSION(6, 9, 0)
	return alloc_weighted(buffer_size, numa_cfg);
#else	/* USE_INTERLEAVING */
	return alloc_memory(buffer_size, numa_cfg->node_idx);
#endif	/* USE_INTERLEAVING */
}








