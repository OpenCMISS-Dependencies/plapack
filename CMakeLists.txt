PROJECT (PLAPACK C Fortran)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)

# TODO machine types & manufacturer
# enable make test
#ENABLE_TESTING()

FIND_PACKAGE(BLAS REQUIRED)
FIND_PACKAGE(LAPACK REQUIRED)
FIND_PACKAGE(MPI REQUIRED)

SET(CMAKE_C_COMPILER ${MPI_C_COMPILER})
SET(CMAKE_C_FLAGS "-O3 -DMACHINE_TYPE=500 -DMANUFACTURE=50 -DDEBUG=0")
SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
SET(CMAKE_Fortran_FLAGS "-fno-f2c -O3 -Q")


# pathes where to look for header files
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/INCLUDE ${MPI_C_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH})

# create mtocpp
file(GLOB template TEMPLATE/*.c)
file(GLOB api API/*.c)
file(GLOB blas1 BLAS1/*.c)
file(GLOB blas2 BLAS2/*.c)
file(GLOB blas3 BLAS3/*.c)
file(GLOB copy COPY/*.c)
file(GLOB lblas Local_BLAS/*.c)
file(GLOB obj OBJ/*.c)
file(GLOB reduce REDUCE/*.c)
file(GLOB util UTIL/*.c)
file(GLOB io IO/*.c)
file(GLOB_RECURSE slv SOLVERS/*.c)
file(GLOB pla POOCLAPACK/*.c)
file(GLOB fortran FORTRAN_interface/*.c FORTRAN_interface/*.f)
ADD_LIBRARY(PLAPACK STATIC ${template} ${api} ${blas1} ${blas2} ${blas3} ${copy} ${lblas} ${obj} ${reduce} ${util} ${io} ${pla} ${fortran})