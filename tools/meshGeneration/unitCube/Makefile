# determine mode of compilation
DEBUG  = NO
# name the compilation targets
TARGET = unitCubeSGF unitCubeSMF unitCubeSMFsimplex unitCubeSMF2 #unitCubeSMFsimplex2

# pass value to application flags
APPFLAGS=

BIN=$(INSILICOROOT)/tools/bin

# include configuration file
include $(INSILICOROOT)/config/convenience.mk


#######################################
# new rules

# unitCubeSMF with simplex elements
unitCubeSMF2:
	$(CXX) $(CPPFLAGS) -DQUADRATIC $(CXXFLAGS) -o unitCubeSMF2 unitCubeSMF.cpp

unitCubeSMFsimplex:
	$(CXX) $(CPPFLAGS) -DSIMPLEX $(CXXFLAGS) -o unitCubeSMFsimplex unitCubeSMF.cpp

# not implemented yet
unitCubeSMFsimplex2:
	$(CXX) $(CPPFLAGS) -DSIMPLEX -DQUADRATIC $(CXXFLAGS) -o unitCubeSMFsimplex2 unitCubeSMF.cpp


all: $(TARGET) install
