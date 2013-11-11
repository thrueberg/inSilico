APPNAME=nitsche.interface
LOG=$(APPNAME).log
MESHFILES=$(wildcard mesh?D.???.smf)
VTKFILES= $(wildcard mesh?D.???.vtk)
TARGET=    interface2D interface3D

all: pre clean build gridgen run test

pre:
	@echo -e "\n\n--------------------------------------------------------"
	@echo    "Testing application $(APPNAME) in $(PWD)"
	@echo        "--------------------------------------------------------"
clean:
	rm -f $(LOG)
	@echo "--- Clean invoked"         >> $(LOG)
	make -f Makefile clean            >> $(LOG) 
	rm -f $(VTKFILES) $(MESHFILES)    >> $(LOG)
	rm -f $(TARGET) interface?D?.dat  >> $(LOG)
	@echo "--- Clean successful"     >> $(LOG)

# the executables
interface2D:
	(make -f Makefile DEBUG=NO SPACEDIM=2 interface -B && mv interface interface2D) >> $(LOG)

interface3D:
	(make -f Makefile DEBUG=NO SPACEDIM=3 interface -B && mv interface interface3D) >> $(LOG)


# build process
build: $(TARGET)
	@echo "--- Build successful"             >> $(LOG)

# generate the grids
gridgen: 
	@echo "--- Grid generation invoked"       >> $(LOG)
	./gridGen.sh 2 U
	./gridGen.sh 3 U
	@echo "--- Grid generation successful"

# execute
run: build gridgen
	@echo "--- Run invoked"       >> $(LOG)
	./exec.sh interface 2 U
	./exec.sh interface 3 U
	@echo "--- Run successful"    >> $(LOG)

# perform tests
test: 
	@echo "--- Test invoked"       >> $(LOG)
	numdiff -X1:4- -X2:4- -r1.e-6 interface2D.dat interface2D.ref.dat >> $(LOG)
	numdiff -X1:4- -X2:4- -r1.e-6 interface3D.dat interface3D.ref.dat >> $(LOG)
	@echo "--- Test successful"   >> $(LOG)

plot: 
	gnuplot interface.gp

triumph: clean
	@echo "(II) Testing application $(APPNAME) successful"
	rm -f $(LOG)

defeat:
	@echo "(EE) Application $(APPNAME) failed" >&2
	cat $(LOG) >&2

