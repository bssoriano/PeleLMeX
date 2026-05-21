# AMReX
DIM = 3
COMP = gnu
PRECISION = DOUBLE
USE_EB = TRUE
USE_HYPRE = FALSE

# Profiling
PROFILE = FALSE
TINY_PROFILE = FALSE
COMM_PROFILE = FALSE
TRACE_PROFILE = FALSE
MEM_PROFILE = FALSE
USE_GPROF = FALSE

# Performance
USE_MPI = TRUE
USE_OMP = FALSE
USE_CUDA = FALSE
USE_HIP = FALSE
USE_SYCL = FALSE

# Debugging
DEBUG = FALSE
FSANITIZER = FALSE
THREAD_SANITIZER = FALSE

# PeleLMeX
USE_PLASMA = FALSE
# CEXE_headers += EBUserDefined.H

# PelePhysics
# Set USE_NONREACT=TRUE for a non-reacting test (e.g.: make -j USE_NONREACT=TRUE)
USE_NONREACT ?= TRUE
ifeq ($(USE_NONREACT), TRUE)
  USE_MANIFOLD    = FALSE
  Chemistry_Model = Null
  Eos_Model       = GammaLaw
  Transport_Model = Constant
else
  USE_MANIFOLD = TRUE
endif

ifeq ($(USE_MANIFOLD), TRUE)
  Eos_Model       = Manifold
  Chemistry_Model = Null
  Transport_Model = Manifold
  Manifold_Dim = 3
  Manifold_Type = Table
endif

# PELE_HOME ?= /home/bsouzas/PeleLMeX/
PELE_PHYSICS_HOME ?= /Users/soriano/Documents/codes/PelePhysics/
PELE_HOME ?= /Users/soriano/Documents/codes/PeleLMeX/
include $(PELE_HOME)/Exec/Make.PeleLMeX
