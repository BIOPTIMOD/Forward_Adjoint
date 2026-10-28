# Compiler settings for the build inside an OPA chain.
# Compilers (Intel oneAPI ifx/icx with Intel MPI) are taken from the spack
# environment of the chain, in $(OPA_HOME)/HOST/$(OPA_HOSTNAME)/spack_env;
# no module has to be loaded.

ifeq ($(strip $(OPA_HOME)),)
$(error OPA_HOME is not defined: this file can be used only inside an OPA chain)
endif
ifeq ($(strip $(OPA_HOSTNAME)),)
$(error OPA_HOSTNAME is not defined: this file can be used only inside an OPA chain)
endif

SPACK_VIEW = $(OPA_HOME)/HOST/$(OPA_HOSTNAME)/spack_env/.spack-env/view

# The Intel MPI wrappers find their installation through I_MPI_ROOT and the
# compilers (ifx, icx) through PATH
export I_MPI_ROOT := $(SPACK_VIEW)/mpi/latest
export PATH := $(SPACK_VIEW)/oneapi-compilers/compiler/latest/bin:$(I_MPI_ROOT)/bin:$(PATH)

# -static-intel links the Intel runtime libraries statically: otherwise the
# executable would not start in a clean environment, because libimf.so needs
# libintlc.so.5 and the RUNPATH of the executable does not apply to it.
# The Makefile links with FC (without LDFLAGS), so the option is set here; it
# has no effect when compiling.
FC      = mpiifx $(OPENMP_FLAG) -static-intel
CC      = mpiicx
F90     = mpiifx $(OPENMP_FLAG)
F77     = mpiifx
CPP     = mpiifx -cpp -E
LD      = $(FC) $(OPENMP_FLAG)
