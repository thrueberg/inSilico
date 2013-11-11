# determine mode of compilation
DEBUG  = 
# choose solver
SOLVER = 
# name the compilation targets
TARGET = bruteForce 

# include configuration file
include $(INSILICOROOT)/config/convenience.mk

# conditional to call CG solver of Eigen
ifeq ($(strip $(SOLVER)),EIGENCG)
	CPPFLAGS += -DCALL_EIGENCG
endif

