# determine mode of compilation
DEBUG  = YES
# choose solver
SOLVER = SUPERLU
# name the compilation targets
TARGET = linearElastic compressible incompressible
# test make files
TESTMK = compressibleTest.mk linearElasticTest.mk

# spatial dimension
SPACEDIM=2
# pass onwards
APPFLAGS=-DSPACEDIM=$(SPACEDIM)
# include configuration file
include $(INSILICOROOT)/config/convenience.mk


