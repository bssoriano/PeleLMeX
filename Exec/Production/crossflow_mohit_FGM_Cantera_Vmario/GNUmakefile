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
USE_MANIFOLD = TRUE
ifeq ($(USE_MANIFOLD), TRUE)
  Eos_Model       = Manifold
  Chemistry_Model = Null
  Transport_Model = Manifold
  Manifold_Dim = 3
  Manifold_Type = Table
else
  Chemistry_Model = drm19
  Eos_Model = Fuego
  Transport_Model = Simple
endif

# PELE_HOME ?= /home/bsouzas/PeleLMeX/
PELE_HOME ?= /home/sabino-mariovane/Documentos/running-codes/PeleLMeX_for_tests
#PELE_PHYSICS_HOME ?= /home/bsouzas/my_LMeX/PelePhysics/
#PELE_HOME ?= /home/bsouzas/my_LMeX/PeleLMeX/
include $(PELE_HOME)/Exec/Make.PeleLMeX
