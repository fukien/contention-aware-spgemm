# Algorithm and memory-placement options.
if(USE_HUGE)
	add_definitions(-DUSE_HUGE)
endif()

# Declare options with default values
option(USE_INTERLEAVING "Enable interleaving" OFF)
option(USE_WEIGHTED_INTERLEAVING "Enable weighted interleaving" OFF)
if(USE_INTERLEAVING)
	add_definitions(-DUSE_INTERLEAVING)
elseif(USE_WEIGHTED_INTERLEAVING)
	add_definitions(-DUSE_WEIGHTED_INTERLEAVING)
endif()

if(USE_NUMA)
	add_definitions(-DUSE_NUMA)
endif()

# NUMA_MASK_VAL starts from 1
if(NUMA_MASK_VAL)
	add_definitions(-DNUMA_MASK_VAL=${NUMA_MASK_VAL})
else()
	# use local numa node by default
	add_definitions(-DNUMA_MASK_VAL=1)
endif()

if (MEMSV)
	add_definitions(-DMEMSV)
endif()

if (MEMSVB)
	add_definitions(-DMEMSVB)
endif()

if (MEM_MON)
	add_definitions(-DMEM_MON)
endif()

if(USE_NVM)
	add_definitions(-DUSE_NVM)
	add_definitions(-DNVM_LOCAL_PREFIX="/pmemfs1/huang/hmem-io")
	add_definitions(-DNVM_REMOTE_PREFIX="/pmemfs0/huang/hmem-io")
endif()


add_definitions(-DRMAT_DIR="dataset")

if(LOCAL_NUMA_NODE_IDX)
	add_definitions(-DLOCAL_NUMA_NODE_IDX=${LOCAL_NUMA_NODE_IDX})
else()
	add_definitions(-DLOCAL_NUMA_NODE_IDX=2)
	# add_definitions(-DLOCAL_NUMA_NODE_IDX=8)
endif()




if(THREAD_NUM)
	add_definitions(-DTHREAD_NUM=${THREAD_NUM})
else()
	add_definitions(-DTHREAD_NUM=32)
endif()

if(USE_HYPERTHREADING)
	add_definitions(-DUSE_HYPERTHREADING)
endif()

if(CSRC_WRITE_OPTIMIZED)
	add_definitions(-DCSRC_WRITE_OPTIMIZED)
endif()

if(CSRC_PREPOPULATED)
	add_definitions(-DCSRC_PREPOPULATED)
endif()

if(MIN_HT_S)
	add_definitions(-DMIN_HT_S=${MIN_HT_S})
else()
	add_definitions(-DMIN_HT_S=8)
endif()

if(HASH_CONSTANT)
	add_definitions(-DHASH_CONSTANT=${HASH_CONSTANT})
else()
	add_definitions(-DHASH_CONSTANT=2654435761)	# for 32-bit
	# add_definitions(-DHASH_CONSTANT=11400714819323198485)	# for 64-bit
	# add_definitions(-DHASH_CONSTANT=107L)	# mtspgemmlib
	# add_definitions(-DHASH_CONSTANT=3L)
	# add_definitions(-DHASH_CONSTANT=1L)
endif()



if (AMS_HIGH_THREAD_NUM)
	add_definitions(-DAMS_HIGH_THREAD_NUM=${AMS_HIGH_THREAD_NUM})
else()
	add_definitions(-DAMS_HIGH_THREAD_NUM=4)
endif()
if(AMS_HIGH_THREAD_NUM GREATER THREAD_NUM)
	message(FATAL_ERROR "AMS_HIGH_THREAD_NUM (${AMS_HIGH_THREAD_NUM}) cannot be greater than THREAD_NUM (${THREAD_NUM})")
endif()


if(AMS_SIZE)
	add_definitions(-DAMS_SIZE=${AMS_SIZE}L)
else()
	add_definitions(-DAMS_SIZE=524288L) # 512K
endif()

if(AMS_MIN_NUM)
	add_definitions(-DAMS_MIN_NUM=${AMS_MIN_NUM}L)
else()
	add_definitions(-DAMS_MIN_NUM=2048L) # 512K
endif()

if (AMS_LPT)
	add_definitions(-DAMS_LPT)
endif()

if (AMS_SIMP)
	add_definitions(-DAMS_SIMP)
endif()

if(AMS_PRT)
	add_definitions(-DAMS_PRT=${AMS_PRT})
else()
	add_definitions(-DAMS_PRT=0.80)
endif()

if(AMS_MAXFLP)
	add_definitions(-DAMS_MAXFLP)
endif()

if(AMS_DNR)
	add_definitions(-DAMS_DNR=true)
endif()

if(AMS_HYOPT)
	add_definitions(-DAMS_HYOPT)
endif()

if (AMS_ACC_THRES)
	add_definitions(-DAMS_ACC_THRES=${AMS_ACC_THRES})
else()
	add_definitions(-DAMS_ACC_THRES=5.00)
endif()

if(AMS_SYMB_THRES)
	add_definitions(-DAMS_SYMB_THRES=${AMS_SYMB_THRES})
else()
	add_definitions(-DAMS_SYMB_THRES=66.0)
endif()

if(AMS_NUMC_THRES)
	add_definitions(-DAMS_NUMC_THRES=${AMS_NUMC_THRES})
else()
	add_definitions(-DAMS_NUMC_THRES=1.5)
endif()


if(SMTACT)
	add_definitions(-DSMTACT)
endif()

if(SMTACT_REAR)
	add_definitions(-DSMTACT_REAR)
endif()

if(SMTACT_QUANTILE)
	add_definitions(-DSMTACT_QUANTILE=${SMTACT_QUANTILE})
else()
	add_definitions(-DSMTACT_QUANTILE=0.90)
endif()

if(SMTACT_SYMB_THRES)
	add_definitions(-DSMTACT_SYMB_THRES=${SMTACT_SYMB_THRES}L)
else()
	add_definitions(-DSMTACT_SYMB_THRES=8000000L) # default: 8M
endif()

if(SMTACT_NUMC_THRES)
	add_definitions(-DSMTACT_NUMC_THRES=${SMTACT_NUMC_THRES}L)
else()
	add_definitions(-DSMTACT_NUMC_THRES=800000L) # default: 800K
endif()

if(DN_MIN)
	add_definitions(-DDN_MIN=${DN_MIN}L)
else()
	add_definitions(-DDN_MIN=1L)	
endif()

if(GUD_ALPHA)
	add_definitions(-DGUD_ALPHA=${GUD_ALPHA})
else()
	add_definitions(-DGUD_ALPHA=0.5)
endif()

if(GUD_MIN)
	add_definitions(-DGUD_MIN=${GUD_MIN}L)
else()
	add_definitions(-DGUD_MIN=1L)
endif()

if(HYB_THRES)
	add_definitions(-DHYB_THRES=${HYB_THRES})
else()
	add_definitions(-DHYB_THRES=0.8)
endif()

if(CF_MARK_0)
	add_definitions(-DCF_MARK_0=${CF_MARK_0})
else()
	add_definitions(-DCF_MARK_0=1.1)
endif()

if(CF_MARK_1)
	add_definitions(-DCF_MARK_1=${CF_MARK_1})
else()
	add_definitions(-DCF_MARK_1=2.1)
endif()

if(RUN_NUM)
	add_definitions(-DRUN_NUM=${RUN_NUM})
else()
	add_definitions(-DRUN_NUM=50)
endif()


if(IN_STATS)
	add_definitions(-DIN_STATS)
endif()

if(IN_VERIFY)
	add_definitions(-DIN_VERIFY)
endif()

if(IN_INSP)
	add_definitions(-DIN_INSP)
endif()

if(IN_DEBUG)
	add_definitions(-DIN_DEBUG)
endif()

if(IN_EXAMINE)
	add_definitions(-DIN_EXAMINE)
endif()

