set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PLAPACK_BINARY_DIR}/TESTING)

INCLUDE( ProcessorCount )
ProcessorCount(SYSTEM_PROCESSOR_COUNT)
	
MACRO(CREATE_TEST NAME)
    add_executable(${NAME} ${NAME}.c ${util})
    target_include_directories(${NAME} PRIVATE ../INCLUDE)
    target_link_libraries(${NAME} plapack pthread) #blas lapack
    if(UNIX)
      target_link_libraries(${NAME} g2c m)
    endif()
    add_test(${NAME} ${MPIEXEC} 
        ${MPIEXEC_NUMPROC_FLAG} ${SYSTEM_PROCESSOR_COUNT}
        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${NAME})
ENDMACRO()

SET(util create_problem.c create_problem_mv.c)
CREATE_TEST(test_lu)
CREATE_TEST(test_elem_mult)
CREATE_TEST(test_inv)
CREATE_TEST(test_tri_inv)
CREATE_TEST(test_sym_inv)
#CREATE_TEST(test_spectral_decomp)
#CREATE_TEST(test_svd)
#CREATE_TEST(test_tri_red) 
CREATE_TEST(test_shift)
CREATE_TEST(test_matrix_one_norm)
CREATE_TEST(test_matrix_inf_norm)
CREATE_TEST(test_iamax)
CREATE_TEST(test_syrk)
CREATE_TEST(test_gemm)
CREATE_TEST(test_trmm)
CREATE_TEST(test_trsm)
CREATE_TEST(test_chol)
CREATE_TEST(test_gemv)
#CREATE_TEST(test_trmv) # not implemented
CREATE_TEST(test_diag)

SET(util ${util} create_problem_ooc.c)
CREATE_TEST(test_ooc)
# Note: The sources for PLA_OOC_Syrk are not included in this PLAPACK version and are only referenced to in this test.
# So i guess this is a bug/artifact?
#CREATE_TEST(test_ooc_syrk)