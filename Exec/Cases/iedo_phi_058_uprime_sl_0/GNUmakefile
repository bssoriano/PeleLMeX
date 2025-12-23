# AMReX
DIM = 3
COMP = llvm
PRECISION = DOUBLE
USE_EB = FALSE
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
USE_HIP = TRUE
USE_SYCL = FALSE

# Debugging
DEBUG = FALSE
FSANITIZER = FALSE
THREAD_SANITIZER = FALSE

# PeleLMeX

# PelePhysics
# Chemistry_Model = air
Chemistry_Model = H2_Dario_final
Eos_Model = Fuego
Transport_Model = Simple
# PELE_CVODE_FORCE_YCORDER = TRUE
PELE_USE_MAGMA=TRUE

# PELE_HOME ?= /home/bsouzas/ITV-RWTH/PeleLMeX/
PELE_HOME      ?= /ccs/home/bsoriano/my_LMeX/PeleLMeX
PELE_PHYSICS_HOME ?= /ccs/home/bsoriano/my_LMeX/PelePhysics
include $(PELE_HOME)/Exec/Make.PeleLMeX
