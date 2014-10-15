include /scratch/home/pb24/PLAPACKR32/Make.include

OBJS = create_problem.o driver.o process_answer.o \
       PLA_LU_linpack.o

#OBJS = create_problem1.o driver.o process_answer.o \
#       PLA_LU_linpack.o

LINPACK : $(OBJS) $(PLAPACKLIB) 
	rm -f LU_driver*
	$(LINKER) -o linpack.x $(LFLAGS) $(OBJS) $(PLAPACKLIB) $(LIB)
	touch LU_driver

clean :
	rm -f LU_driver* *.o core *~

