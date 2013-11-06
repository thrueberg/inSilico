# determine mode of compilation
DEBUG  = NO
# name the compilation targets
TARGET = unitCubeSMF unitCubeSGF unitCubeSMFsimplex

# pass value to application flags
APPFLAGS=

BIN=$(INSILICOROOT)/tools/bin

# include configuration file
include $(INSILICOROOT)/config/convenience.mk


#######################################
# new rules

# unitCubeSMF with simplex elements
unitCubeSMFsimplex:
	$(CXX) $(CPPFLAGS) -DSIMPLEX $(CXXFLAGS) -o unitCubeSMFsimplex unitCubeSMF.cpp


all: $(TARGET) install
