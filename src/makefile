FC = gfortran #ifort
OPT = #-no-wrap-margin

EXE = adj.xx

SOURCE = bioptimod_memory.o fm34.o SLAE.o adj_new.o main_adj.o

all: $(SOURCE)
	$(FC) $(OPT) $(SOURCE) -o $(EXE)

%.o: %.f90
	$(FC) $(OPT) -c $< -o $*.o

%.o: %.F90
	$(FC) $(OPT) -c $< -o $*.o

clear:
	rm *.o

clean: clear

SAVE = *.f90 *.txt *.nml makefile

tar: 
	tar -zcf backup_`date +%F`.tar.gz $(SAVE)
	echo 'backup_'`date +%F`'.tar.gz created'

