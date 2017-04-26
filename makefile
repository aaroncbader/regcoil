# makefile for NERSC Edison and Cori
# You must first load the cray-netcdf module and python module:
#   module load cray-netcdf python
# It is convenient to run
#   module unload cray-libsci
# to avoid warning messages about libsci during compiling.


ifdef NERSC_HOST
        HOSTNAME = $(NERSC_HOST)
else ifeq ($(HOST),abader.cptc.wisc.edu)
  HOSTNAME="abader"
else
  HOSTNAME="laptop"
endif


ifeq ($(HOSTNAME),edison)
	FC = ftn
	## NERSC documentation recommends against specifying -O3
	## -mkl MUST APPEAR AT THE END!!
	EXTRA_COMPILE_FLAGS = -openmp -mkl
	EXTRA_LINK_FLAGS =  -openmp -mkl -Wl,-ydgemm_
	# Above, the link flag "-Wl,-ydgemm_" causes the linker to report which version of DGEMM (the BLAS3 matrix-matrix-multiplication subroutine) is used.
	# For batch systems, set the following variable to the command used to run jobs. This variable is used by 'make test'.
	REGCOIL_COMMAND_TO_SUBMIT_JOB = srun -n 1 -c 24

else ifeq ($(HOSTNAME),cori)
	FC = ftn
	## NERSC documentation recommends against specifying -O3
	## -mkl MUST APPEAR AT THE END!!
	EXTRA_COMPILE_FLAGS = -qopenmp -mkl
	EXTRA_LINK_FLAGS =  -qopenmp -mkl -Wl,-ydgemm_
	# Above, the link flag "-Wl,-ydgemm_" causes the linker to report which version of DGEMM (the BLAS3 matrix-matrix-multiplication subroutine) is used.
	# For batch systems, set the following variable to the command used to run jobs. This variable is used by 'make test'.
	REGCOIL_COMMAND_TO_SUBMIT_JOB = srun -n 1 -c 32
else ifeq ($(HOSTNAME),"abader") 
  FC = mpif90
  EXTRA_COMPILE_FLAGS = -O2 -fopenmp -I/home/abader/netcdf/include  -I/usr/include -ffree-line-length-none
  EXTRA_LINK_FLAGS = -O2 -fopenmp -L/home/abader/netcdf/lib  -lnetcdff  -lnetcdf -lblas -llapack

else
  FC = mpif90
  #EXTRA_COMPILE_FLAGS = -fopenmp -I/opt/local/include -ffree-line-length-none -cpp
  EXTRA_COMPILE_FLAGS = -fopenmp -I/opt/local/include -ffree-line-length-none
  EXTRA_LINK_FLAGS =  -fopenmp -L/opt/local/lib -lnetcdff  -lnetcdf -framework Accelerate

  # For batch systems, set the following variable to the command used to run jobs. This variable is used by 'make test'.
	REGCOIL_COMMAND_TO_SUBMIT_JOB =
  
endif


# End of system-dependent variable assignments

LIBSTELL_DIR = mini_libstell
TARGET = regcoil

export

.PHONY: all clean

all: $(TARGET)

include makefile.depend

%.o: %.f90 $(LIBSTELL_DIR)/mini_libstell.a
	$(FC) $(EXTRA_COMPILE_FLAGS) -I $(LIBSTELL_DIR) -c $<

%.o: %.f $(LIBSTELL_DIR)/mini_libstell.a
	$(FC) $(EXTRA_COMPILE_FLAGS) -I $(LIBSTELL_DIR) -c $<

$(TARGET): $(OBJ_FILES) $(LIBSTELL_DIR)/mini_libstell.a
	$(FC) -o $(TARGET) $(OBJ_FILES) $(LIBSTELL_DIR)/mini_libstell.a $(EXTRA_LINK_FLAGS)
#	$(FC) -o $(TARGET) $(OBJ_FILES) $(LIBSTELL_DIR)/libstell.a $(EXTRA_LINK_FLAGS)

$(LIBSTELL_DIR)/mini_libstell.a:
	$(MAKE) -C mini_libstell

clean:
	rm -f *.o *.mod *.MOD *~ $(TARGET)
	cd $(LIBSTELL_DIR); rm -f *.o *.mod *.MOD *.a

test: $(TARGET)
	@echo "Beginning functional tests." && cd examples && export REGCOIL_RETEST=no && ./runExamples.py

retest: $(TARGET)
	@echo "Testing existing output files for examples without re-running then." && cd examples && export REGCOIL_RETEST=yes && ./runExamples.py

test_make:
	@echo HOSTNAME is $(HOSTNAME)
	@echo FC is $(FC)
	@echo EXTRA_COMPILE_FLAGS is $(EXTRA_COMPILE_FLAGS)
	@echo EXTRA_LINK_FLAGS is $(EXTRA_LINK_FLAGS)
