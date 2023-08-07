cmake_minimum_required(VERSION 3.12.4)
cmake_policy(VERSION 3.12.4)
cmake_policy(SET CMP0003 NEW)
cmake_policy(SET CMP0074 NEW)
project(cosmictiger CXX C Fortran ASM)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
find_package(HPX REQUIRED NO_CMAKE_PACKAGE_REGISTRY)
find_package(Boost REQUIRED)
option(COSMICTIGER_WITH_CUDA "Enable CUDA" ON)
if( COSMICTIGER_WITH_CUDA ) 
	enable_language(CUDA)
	set(CUDA_PROPAGATE_HOST_FLAGS FALSE)
	execute_process(COMMAND nvcc -lcuda ${PROJECT_SOURCE_DIR}/src/cuda_detect.cu -o cuda_detect)
	execute_process(COMMAND ./cuda_detect OUTPUT_VARIABLE CUDA_ARCH)
	message(STATUS "CUDA Architecture: ${CUDA_ARCH}")
endif()
set(CMAKE_CXX_FLAGS "-Wno-attributes --fast-math")
set(COSMICTIGER_WITH_ORDER "7" )
message(STATUS "cosmictiger expansion order P+1 = ${COSMICTIGER_WITH_ORDER}")

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -fno-pie -no-pie")
endif()

set(fft_source_files
    src/fft/fft.cpp
    src/fft/twiddles.cpp
)

set(cosmictiger_source_files
    src/cosmology.cpp
    src/drift.cpp
    src/kick.cpp
    src/main.cpp
    src/memory.cpp
    src/options.cpp
    src/particle.cpp
    src/util.cpp
    ${fft_source_files}
)      
                       
set(cosmictiger_header_files
    include/constants.hpp
    include/cosmictiger.hpp
    include/cosmology.hpp
    include/drift.hpp
    include/fft.hpp
    include/kick.hpp
    include/options.hpp
    include/particle.hpp
    include/util.hpp
)

add_hpx_executable( cosmictiger SOURCES ${cosmictiger_source_files} HEADERS ${cosmictiger_header_files} )
  
add_hpx_executable( fft_test SOURCES ${fft_source_files} src/util.cpp src/fft/fft_test.cpp HEADERS cosmictiger/fft.hpp )
  
if( COSMICTIGER_WITH_CUDA )
	target_link_libraries(cosmictiger PUBLIC cudart)
	target_compile_options(cosmictiger PRIVATE $<$<COMPILE_LANGUAGE:CUDA>: -DUSE_CUDA --use_fast_math -lineinfo ;>)
    target_compile_options(cosmictiger PRIVATE -DUSE_CUDA -DBOOST_STACKTRACE_USE_ADDR2LINE )
	set_property(TARGET cosmictiger PROPERTY CUDA_ARCHITECTURES ${CUDA_ARCH})
	set_property(TARGET cosmictiger PROPERTY CUDA_SEPARABLE_COMPILATION ON)
endif()
 
target_link_libraries(cosmictiger PUBLIC sfmm simd)
target_include_directories(cosmictiger PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(cosmictiger PUBLIC ${PROJECT_SOURCE_DIR}/spherical_fmm)
target_compile_definitions(cosmictiger PUBLIC ORDER=${COSMICTIGER_WITH_ORDER})

target_include_directories(fft_test PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(fft_test PUBLIC ${PROJECT_SOURCE_DIR}/spherical_fmm)
target_link_libraries(fft_test PUBLIC fftw3)



