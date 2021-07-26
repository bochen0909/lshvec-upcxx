cmake_minimum_required(VERSION 3.10)

project(LSHVEC-UPCXX VERSION 0.1.0 LANGUAGES CXX)

#specify the C++ standard
set(CMAKE_STANDARD 14)
set(CMAKE_STANDARD_REQUIRED True)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g -ggdb -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -s")

option(BUILD_LSHVEC_UPCXX "build UPCXX version of LSHVEC" ON)



set(LOCAL_EXT_PREFIX_DIR ${PROJECT_SOURCE_DIR}/extlib/local)       

add_custom_target(hsds_upcxx ALL)

include(cmake/CMakefile.external)


############# upcxx programs ###########################
if(BUILD_LSHVEC_UPCXX)
	############# check UPCXX ##############################
	set (UPCXX_VERBOSE ON)
	include(cmake/UPCXXConfig.cmake)
	if (NOT EXISTS "${UPCXX_META_EXECUTABLE}")
	   message( FATAL_ERROR "cannot find upcxx" )
	endif()
	############# end check UPCXX ##############################
	
	find_package(MPI 2.0 REQUIRED COMPONENTS CXX)
	find_package(ZLIB REQUIRED)

	link_directories( "${LOCAL_EXT_PREFIX_DIR}/lib" )
	link_directories( "${LOCAL_EXT_PREFIX_DIR}/lib64" )

 
	add_executable(kgraph_upcxx
		src/trainer.cpp 
		src/log.cpp src/utils.cpp
		)
 

						
	foreach(EXEC  kgraph_upcxx )
		add_dependencies(${EXEC} spdloglib gzstreamlib)
		target_include_directories(${EXEC} PUBLIC ${UPCXX_INCLUDE_DIRS} "${LOCAL_EXT_PREFIX_DIR}/include" )
		target_compile_definitions(${EXEC} PUBLIC ${UPCXX_DEFINITIONS})
		target_compile_definitions(${EXEC} PUBLIC "-DUSE_UPCXXCLIENT")
        target_link_libraries(${EXEC} PUBLIC ${MPI_CXX_LIBRARIES})
		target_link_libraries(${EXEC} PUBLIC spdlog gzstream ZLIB::ZLIB ${UPCXX_LIBRARIES})
	    set_target_properties(${EXEC} PROPERTIES EXCLUDE_FROM_ALL TRUE)
	    add_dependencies(hsds_upcxx ${EXEC})
	endforeach(EXEC)
	
	
endif()


find_package(OpenMP REQUIRED)
find_package(ZLIB REQUIRED)
include_directories("${OpenMP_CXX_FLAGS}")
add_executable(hash src/hasher_main.cpp src/CRandProj.cpp src/utils.cpp src/kmer.cpp src/io.cpp)
target_link_libraries(hash PUBLIC spdlog gzstream ZLIB::ZLIB  OpenMP::OpenMP_CXX)


############# unit test ###########################
include_directories( "${LOCAL_EXT_PREFIX_DIR}/include")
add_executable(test_utils src/test_utils.cpp src/utils.cpp)
add_executable(test_kmer src/test_kmer.cpp src/utils.cpp src/kmer.cpp)

add_executable(test_read_fasta src/test_read_fasta.cpp src/utils.cpp src/io.cpp)
add_dependencies(test_read_fasta spdloglib gzstreamlib)
target_include_directories(test_read_fasta PUBLIC   "${LOCAL_EXT_PREFIX_DIR}/include" )
target_link_libraries(test_read_fasta PUBLIC spdlog gzstream ZLIB::ZLIB  )
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data/knucleotide.fasta DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data/knucleotide.fasta.gz  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


find_package(OpenMP REQUIRED)
include_directories("${OpenMP_CXX_FLAGS}")
add_executable(test_hash src/test_hash.cpp src/CRandProj.cpp src/utils.cpp )
target_link_libraries( test_hash  OpenMP::OpenMP_CXX )