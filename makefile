#==================================
# Makefile for micro.exe
#==================================
.SUFFIXES: .f .c .cpp .f90

.f90.o:         $<
	$(MPIF) -c   $<
.f.o:           $<
	$(MPIF) -c  $<

FILES = ZQMPI.f90 \
	Modules.f90 \
	MPIMODULE.f90 \
	FunctionMD.f90 \
	TreeModule.f90 \
	Communication.f90 \
	CommonUsed.f90 \
	ComputeStrain.f90 \
	EliminateSmall.f90 \
	Output.f90 \
	ReadWritePast.f90 \
	LoopRearrange.f90 \
	NewAnnihilation.f90 \
	SplitBuild.f90 \
	Crossslip.f90\
	Dynamics.f90 \
	Update.f90 \
	SimInitialization.f90 \
	Simulation.f90 \
	Microplasticity.f90

OBJSA  = $(FILES:.f90=.o)
OBJS   = $(OBJSA:.f=.o)

TARGET = micro
OPT     =-g 
MPIF     = mpif90 $(OPT) $(MPI_COMPILE_FLAGS)


all:    $(OBJS)
	$(MPIF) -o $(TARGET) $(OBJS) $(MPI_LD_FLAGS) -lmpi

clean:
	touch nofile.o
	rm -f $(OBJS) $(TARGET) MOD*  core *.mod *.o 
