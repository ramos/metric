F90 = f90 -CB -check all
F90LD = f90 -CB -check all

all:
	@make geometry.o bases.o metric.out 

bases.o: bases.f90
	$(F90) -c -I/home/alberto/fortran/lib bases.f90

geometry.o: geometry.f90
	$(F90) -c -I/home/alberto/fortran/lib geometry.f90

metric.o: metric.f90
	$(F90) -c -I/home/alberto/fortran/lib metric.f90

metric.out: bases.o geometry.o metric.o
	$(F90LD) -o metric.out bases.o geometry.o metric.o \
-L/home/alberto/fortran/lib -lf90

cleanall:
	rm *.o *.out *.mod
