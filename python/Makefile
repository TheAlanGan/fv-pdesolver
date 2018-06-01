##--------------------- Makefile for parallel -------------------------#
##usage:  make compile ; make run  or  make pbs
##---------------------------------------------------------------------#
 
##============================ set MPI, compiler ======================#
##-----> set appropriate compiler_wrapper: mpif77 mpif90 mpicc mpic++
###  COMP   = $(MPI)/bin/mpif77
##-----> set appropriate extension: f  c  cpp
###  EXT    = f
###  LFLAGs = 
#for C:  LFLAGs = -lm
##-------------------------- for all:
###  FLAGs  = -g  $(MPIlnk) 
# FLAGs  = -O3 $(MPIlnk) 
###  MPIlnk = -I$(MPI)/include  -L$(MPI)/lib
##---------------------> set path to openMPI on local:
###  MPI  = /opt/open-mpi
##---------------------> parallel on cluster:  for open-mpi with intel
# MPI  = /opt/open-mpi/tcp-intel11
##-------------------- parallel on cluster:  for open-mpi with gnu44
# MPI  = /opt/open-mpi/tcp-gnu44
#
##========================== set source code  =========================#
##--------------->set names for your PROGram and std I/O files: 
# PROG = code1D.x 
  PROG = code2D.x
  INPUT  = ./dat 
  OUTPUT = out 
##--------------------> set code components: 
  CODE_o = *.pyc

##-------------------- not needed for most 'make':
# CODE  = z.main.$(EXT) z.mainMR.$(EXT) z.mainWR.$(EXT) z.io.$(EXT)  \
#         z.setup.$(EXT) z.update.$(EXT) z.messaging.$(EXT)
 
##====================== create executable: make compile ============# 
#$(CODE_o):%.o: %.$(EXT)
#	$(COMP) $(FLAGs) -c $< -o $@

compile:
	python -m compileall .
	zip -r $(PROG) $(CODE_o)
	rm -f *.pyc
##====================== execute: make run | make pbs ==========# 
run: 
	mpiexec -n  3 python $(PROG)  $(INPUT)  $(OUTPUT) 
#or:	$(MPI)/bin/mpirun  -np 3  $(PROG).x < $(INPUT) > $(OUTPUT) 
pbs: 
	@ vi PBSscript 
	make clean 
	qsub  PBSscript 

clean: 
	rm -f *.pyc
	rm -f  o.*  DONE  watch 
##---------------------------------------------------------------------#