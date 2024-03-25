PELE_HOME      ?= /home/bsouzas/my_LMeX/PeleLMeX/

# AMReX
DIM             = 3
DEBUG           = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE
TINY_PROFILE    = FALSE
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
Chemistry_Model = jeta_soot
# Chemistry_Model = air
Eos_Model       = Fuego
Transport_Model = Simple
PELE_USE_MAGMA = FALSE

USE_PARTICLES = TRUE
SPRAY_FUEL_NUM = 1

USE_SOOT = TRUE
NUM_SOOT_MOMENTS = 3
USE_RADIATION = TRUE

ifeq ($(USE_RADIATION), TRUE)
ifeq ($(USE_HIP), TRUE)
LIBRARIES += -lstdc++fs
endif
endif

include $(PELE_HOME)/Exec/Make.PeleLMeX

