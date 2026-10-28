cmake_minimum_required(VERSION 3.18)
project(contention_aware_spgemm LANGUAGES C CXX)

if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")
    message(FATAL_ERROR "These NUMA/pthread kernels require Linux (x86-64).")
endif()

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)

foreach(dep numa config pmem)
    find_library(${dep}_LIBRARY NAMES ${dep} REQUIRED)
endforeach()

find_path(NUMA_INCLUDE_DIR numa.h REQUIRED)
find_path(CONFIG_INCLUDE_DIR libconfig.h REQUIRED)
find_path(PMEM_INCLUDE_DIR libpmem.h REQUIRED)

set(THREAD_NUM 32 CACHE STRING "Worker threads")
set(AMS_HIGH_THREAD_NUM 4 CACHE STRING "Workers for high-contention bins")
set(RUN_NUM 2 CACHE STRING "Iterations including warmup")

if(THREAD_NUM LESS 2 OR AMS_HIGH_THREAD_NUM LESS 1 OR AMS_HIGH_THREAD_NUM GREATER_EQUAL THREAD_NUM)
    message(FATAL_ERROR "Require THREAD_NUM >= 2 and 1 <= AMS_HIGH_THREAD_NUM < THREAD_NUM")
endif()

if(RUN_NUM LESS 2)
    message(FATAL_ERROR "RUN_NUM must include a warmup and at least one measured iteration")
endif()

set(CORE_NUM 64 CACHE STRING "Physical cores across sockets")
set(SOCKET_CORE_NUM 32 CACHE STRING "Physical cores per socket")
set(CXL_NUMA_NODE_IDX 3 CACHE STRING "CXL node ID plus one")
set(LAST_NUMA_NODE_IDX 3 CACHE STRING "Last node ID plus one")
set(CFG_PATH "${PROJECT_SOURCE_DIR}/config/mc/mxc.cfg" CACHE FILEPATH "CPU topology config")

get_filename_component(CFG_PATH "${CFG_PATH}" ABSOLUTE BASE_DIR "${PROJECT_SOURCE_DIR}")
if(NOT EXISTS "${CFG_PATH}")
    message(FATAL_ERROR "CPU topology config does not exist: ${CFG_PATH}")
endif()

add_compile_options(-Wall -Wextra -march=native -UNDEBUG)
add_definitions(
    -DPROJECT_PATH="${PROJECT_SOURCE_DIR}"
    -DCFG_PATH="${CFG_PATH}"
    -DCHAR_BUFFER_LEN=256
    -DCORE_NUM=${CORE_NUM}
    -DSOCKET_CORE_NUM=${SOCKET_CORE_NUM}
    -DCXL_NUMA_NODE_IDX=${CXL_NUMA_NODE_IDX}
    -DLAST_NUMA_NODE_IDX=${LAST_NUMA_NODE_IDX}
)

include(cmake/Options.cmake)

add_library(Inc STATIC src/inc/utils.c)
target_include_directories(Inc PUBLIC
    ${NUMA_INCLUDE_DIR}
    ${CONFIG_INCLUDE_DIR}
    ${PMEM_INCLUDE_DIR}
)
target_link_libraries(Inc PUBLIC
    Threads::Threads
    ${numa_LIBRARY}
    ${config_LIBRARY}
    ${pmem_LIBRARY}
    m
)

add_library(Mm STATIC src/mm/common.c src/mm/mm_utils.c)
target_link_libraries(Mm PUBLIC Inc OpenMP::OpenMP_C)

add_library(Ams STATIC src/ams/ams_utils.c)
target_link_libraries(Ams PUBLIC Mm)

foreach(algo ams_hashspgemm ams_hserscspgemm ams_hybaccspgemm)
    add_executable(${algo} src/${algo}.c)

    target_link_libraries(${algo} PRIVATE Ams)
endforeach()

foreach(tool rmat_load rmat_ones)
    add_executable(${tool} src/${tool}.c)
    target_link_libraries(${tool} PRIVATE Mm)
endforeach()

# PaRMAT remains independently buildable through its original Release/makefile.
option(BUILD_PARMAT "Also build the bundled R-MAT generator" OFF)

if(BUILD_PARMAT)
    add_executable(PaRMAT
        PaRMAT/src/PaRMAT.cpp
        PaRMAT/src/Edge.cpp
        PaRMAT/src/Square.cpp
        PaRMAT/src/GraphGen_sorted.cpp
        PaRMAT/src/GraphGen_notSorted.cpp
        PaRMAT/src/utils.cpp
    )
    target_link_libraries(PaRMAT PRIVATE Threads::Threads)
endif()
