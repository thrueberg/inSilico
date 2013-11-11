APPNAME=nitsche.box
LOG=$(APPNAME).log
MESHFILES= $(wildcard grid?D.???.sgf) $(wildcard mesh?D.???.smf)
VTKFILES=  $(wildcard grid?D.???.vtk) $(wildcard mesh?D.???.vtk)
TARGET=    box2DS box3DS box2DU box3DU


all: pre clean build gridgen run test

pre:
	@echo -e "\n\n--------------------------------------------------------"
	@echo    "Testing application $(APPNAME) in $(PWD)" 
	@echo        "--------------------------------------------------------"

clean:
	rm -f $(LOG)
	@echo "--- Clean invoked"       >> $(LOG)
	make -f Makefile clean          >> $(LOG)  
	rm -f $(VTKFILES) $(MESHFILES)  >> $(LOG)
	rm -f $(TARGET)  box?D?.dat     >> $(LOG)
	@echo "--- Clean successful"    >> $(LOG)

# the executables
box2DS:
	(make -f Makefile DEBUG=NO SPACEDIM=2 STRUCTURED=YES box -B && mv box box2DS) >> $(LOG)
box3DS:
	(make -f Makefile DEBUG=NO SPACEDIM=3 STRUCTURED=YES box -B && mv box box3DS) >> $(LOG)
box2DU:
	(make -f Makefile DEBUG=NO SPACEDIM=2 STRUCTURED=NO  box -B && mv box box2DU) >> $(LOG)
box3DU:
	(make -f Makefile DEBUG=NO SPACEDIM=3 STRUCTURED=NO  box -B && mv box box3DU) >> $(LOG)

# build process
build: $(TARGET)
	@echo "--- Build successful"              >> $(LOG)

# generate the grids
gridgen: 
	@echo "--- Grid generation invoked"       >> $(LOG)
	./gridGen.sh 2 S
	./gridGen.sh 3 S
	./gridGen.sh 2 U
	./gridGen.sh 3 U
	@echo "--- Grid generation successful"

# execute
run: build gridgen
	@echo "--- Run invoked"       >> $(LOG)
	./exec.sh box 2 S
	./exec.sh box 3 S
	./exec.sh box 2 U
	./exec.sh box 3 U
	@echo "--- Run successful"

# perform tests
test: 
	@echo "--- Test invoked"       >> $(LOG)
	numdiff -X1:4- -X2:4- -r1.e-6 box2DS.dat box2DS.ref.dat  >> $(LOG)
	numdiff -X1:4- -X2:4- -r1.e-6 box3DS.dat box3DS.ref.dat  >> $(LOG)
	numdiff -X1:4- -X2:4- -r1.e-6 box2DU.dat box2DU.ref.dat  >> $(LOG)
	numdiff -X1:4- -X2:4- -r1.e-6 box3DU.dat box3DU.ref.dat  >> $(LOG)
	@echo "--- Test successful"

plot: run
	gnuplot box.gp

triumph: clean
	@echo "(II) Testing application $(APPNAME) successful"	
	rm -f $(LOG)

defeat:
	@echo "(EE) Application $(APPNAME) failed" >&2
	cat $(LOG) >&2

