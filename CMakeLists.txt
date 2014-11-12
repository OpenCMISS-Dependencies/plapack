CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT (PLAPACK VERSION 3.0 LANGUAGES C Fortran)

option(BUILD_TESTS "${PROJECT_NAME} - Build tests" ON)

FIND_PACKAGE(BLAS CONFIG REQUIRED)
FIND_PACKAGE(LAPACK CONFIG REQUIRED)
FIND_PACKAGE(MPI REQUIRED)

SET(CMAKE_C_COMPILER ${MPI_C_COMPILER})
SET(CMAKE_C_FLAGS "-DMACHINE_TYPE=500 -DMANUFACTURE=50")
SET(CMAKE_C_FLAGS_RELEASE "-O3 -DDEBUG=0")
SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
SET(CMAKE_Fortran_FLAGS "-fno-f2c -O3 -Q")

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
file(GLOB slv_e SOLVERS/Eigensolvers/*.c)
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
file(GLOB slv_s SOLVERS/SVD/*.c)
SET(slv ${slv_c} ${slv_e} ${slv_i} ${slv_l} ${slv_lu} ${slv_qr} ${slv_r} ${slv_s})
file(GLOB pla POOCLAPACK/*.c)
file(GLOB fortran FORTRAN_interface/*.c FORTRAN_interface/*.f)
add_library(plapack STATIC ${template} ${api} ${blas1} ${blas2} ${blas3} ${lblas} ${copy} ${obj} ${reduce} ${util} ${io} ${pla} ${fortran} ${slv})
target_link_libraries(plapack PUBLIC blas lapack ${MPI_C_LIBRARIES} ${MPI_Fortran_LIBRARIES})

# need also mpi for plapack (c/fortran)
TARGET_INCLUDE_DIRECTORIES(plapack
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/INCLUDE> 
    PRIVATE ${MPI_C_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH}
)
install(TARGETS plapack EXPORT plapack-config
    DESTINATION lib
    INCLUDES DESTINATION include)
install(EXPORT plapack-config DESTINATION lib/cmake)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/INCLUDE/ DESTINATION include
        FILES_MATCHING PATTERN "*.h")
include(CMakePackageConfigHelpers)
WRITE_BASIC_PACKAGE_VERSION_FILE(${CMAKE_CURRENT_BINARY_DIR}/plapack-config-version.cmake COMPATIBILITY AnyNewerVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/plapack-config-version.cmake DESTINATION lib/cmake)

if (BUILD_TESTS)
    enable_testing()
    add_subdirectory(TEST)
endif()