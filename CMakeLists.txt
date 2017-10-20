CMAKE_MINIMUM_REQUIRED(VERSION 3.3)
PROJECT (PLAPACK VERSION 3.0 LANGUAGES C Fortran)

option(BUILD_TESTS "${PROJECT_NAME} - Build tests" ON)
set(PACKAGE_CONFIG_DIR "lib/cmake" CACHE STRING "Directory for package config files (relative to CMAKE_INSTALL_PREFIX)")
option(WITH_OPENMP "${PROJECT_NAME} - Use OpenMP" ON)

# Find BLAS separately so we can use 'blas' target.
find_package(BLAS ${BLAS_VERSION} REQUIRED)
find_package(LAPACK ${LAPACK_VERSION} REQUIRED)
find_package(MPI REQUIRED)

if (WIN32)
    set(MACHINE_TYPE 600)
else()
    set(MACHINE_TYPE 500)
endif()
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DMACHINE_TYPE=${MACHINE_TYPE} -DMANUFACTURE=50")
SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 -DDEBUG=0")
if (CMAKE_COMPILER_IS_GNUFortran)
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-f2c -Q")
endif()

if (WITH_OPENMP)
  find_package(OpenMP QUIET)
  if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  else()
    message(WARNING "OpenMP requested but not found. Skipping.")
  endif()
endif()

# create mtocpp
file(GLOB template TEMPLATE/*.c)
file(GLOB api API/*.c)
file(GLOB blas1 BLAS1/*.c)
file(GLOB blas2 BLAS2/*.c)
file(GLOB blas3 BLAS3/*.c)
file(GLOB copy COPY/*.c)
SET(lblas 
    Local_BLAS/PLA_Local_BLAS1.c 
	Local_BLAS/PLA_Local_BLAS2.c 
	Local_BLAS/PLA_Local_BLAS1_enter_exit.c 
	Local_BLAS/PLA_Local_BLAS2_enter_exit.c 
	Local_BLAS/PLA_Local_BLAS3_enter_exit.c 
	Local_BLAS/PLA_Local_iamax.c 
	Local_BLAS/PLA_Local_syr.c 
	Local_BLAS/PLA_Local_syr2.c 
    Local_BLAS/PLA_Local_gemm.c 
    Local_BLAS/PLA_Local_symm.c 
	Local_BLAS/PLA_Local_syrk.c 
	Local_BLAS/PLA_Local_herk.c 
	Local_BLAS/PLA_Local_syr2k.c 
	Local_BLAS/PLA_Local_trsm.c 
	Local_BLAS/PLA_Local_trmm.c 
	Local_BLAS/PLA_Local_inquiry.c 
	Local_BLAS/PLA_Local_misc.c 
	Local_BLAS/PLA_Local_absolute_max.c 
	Local_BLAS/PLA_Conjugate.c 
	Local_BLAS/PLA_OMP_gemm.c
)
file(GLOB obj OBJ/*.c)
file(GLOB reduce REDUCE/*.c)
file(GLOB util UTIL/*.c)
file(GLOB io IO/*.c)
file(GLOB slv_c SOLVERS/Cholesky/*.c)
file(GLOB slv_e SOLVERS/Eigensolvers/*.c SOLVERS/Eigensolvers/*.f)
file(GLOB slv_i SOLVERS/Invert/*.c)
file(GLOB slv_l SOLVERS/LINPACK/*.c)
file(GLOB slv_lu SOLVERS/LU/*.c)
file(GLOB slv_qr SOLVERS/QR/*.c)
SET(slv_r
    SOLVERS/Reduction/PLA_Apply_sym_House.c 
    SOLVERS/Reduction/PLA_Tri_red.c 
    SOLVERS/Reduction/PLA_Tri_red_enter_exit.c 
    SOLVERS/Reduction/PLA_Create_sym_eigenproblem.c 
    SOLVERS/Reduction/PLA_Bi_red_simple.c 
    SOLVERS/Reduction/PLA_Bi_red_enter_exit.c
    SOLVERS/Reduction/PLA_Create_svd_problem.c
)
file(GLOB slv_s SOLVERS/SVD/*.c SOLVERS/SVD/*.f)
SET(slv ${slv_c} ${slv_e} ${slv_i} ${slv_l} ${slv_lu} ${slv_qr} ${slv_r} ${slv_s})
file(GLOB pla POOCLAPACK/*.c)
file(GLOB fortran FORTRAN_interface/*.c FORTRAN_interface/*.f)
if (WIN32)
    set(rand-fix rand-fix.c)
endif()
add_library(plapack ${template} ${api} ${blas1} ${blas2} ${blas3} ${lblas} ${copy} ${obj} ${reduce} ${util} ${io} ${pla} ${fortran} ${slv} ${rand-fix})
if (DEFINED FORTRAN_MANGLING)
    target_compile_definitions(plapack PUBLIC $<BUILD_INTERFACE:${FORTRAN_MANGLING}>)
endif()
target_include_directories(plapack PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/INCLUDE> 
)
target_link_libraries(plapack PUBLIC blas lapack mpi-c)

install(TARGETS plapack EXPORT plapack-config
    DESTINATION lib
    INCLUDES DESTINATION include)
install(EXPORT plapack-config DESTINATION ${PACKAGE_CONFIG_DIR})
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/INCLUDE/ DESTINATION include
        FILES_MATCHING PATTERN "*.h")
include(CMakePackageConfigHelpers)
WRITE_BASIC_PACKAGE_VERSION_FILE(${CMAKE_CURRENT_BINARY_DIR}/plapack-config-version.cmake COMPATIBILITY AnyNewerVersion)
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/plapack-config-dependencies.cmake 
    "include(CMakeFindDependencyMacro)\r\n"
    	"set(${PROJECT_NAME}_IMPORT_PREFIX \${_IMPORT_PREFIX})\r\n"
)
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/plapack-config-dependencies.cmake "find_dependency(LAPACK ${LAPACK_VERSION})\r\n")
if (WITH_OPENMP)
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/plapack-config-dependencies.cmake "find_dependency(OpenMP)\r\n")
endif()
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/plapack-config-dependencies.cmake "set(_IMPORT_PREFIX \${${PROJECT_NAME}_IMPORT_PREFIX})")
install(FILES 
    ${CMAKE_CURRENT_BINARY_DIR}/plapack-config-version.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/plapack-config-dependencies.cmake
    DESTINATION ${PACKAGE_CONFIG_DIR})

if (BUILD_TESTS)
    enable_testing()
    add_subdirectory(TEST)
endif()
