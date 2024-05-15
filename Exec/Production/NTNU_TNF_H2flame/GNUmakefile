PELE_HOME = /home/bsouzas/my_LMeX/PeleLMeX

# AMReX
DIM             = 3
DEBUG           = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE
TINY_PROFILE    = TRUE
USE_EB          = TRUE

# Compilation
COMP            = gnu
USE_MPI         = TRUE
USE_OMP         = FALSE
USE_CUDA        = FALSE
USE_HIP         = FALSE

# PeleLMeX
USE_EFIELD      = FALSE
USE_HYPRE       = FALSE

# PelePhysics
Chemistry_Model = LiDryer
# Chemistry_Model = H2_ammonia
# Chemistry_Model = air
Eos_Model       = Fuego
Transport_Model = Simple

USE_MF          = TRUE
NUM_MF          = 2
USE_AGE         = TRUE

CEXE_headers    += EBUserDefined.H

include $(PELE_HOME)/Exec/Make.PeleLMeX

ifeq ($(USE_CUDA),TRUE)
DEFINES += -DAMREX_GPUS_PER_NODE=2
endif
