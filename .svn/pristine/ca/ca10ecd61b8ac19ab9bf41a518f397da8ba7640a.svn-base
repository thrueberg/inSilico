APPNAME=mesh
LOG=$(APPNAME).log

BOUNDSMF= squareLT_boundary.smf cubeLT_boundary.smf squareGrid_boundary.smf cubeGrid_boundary.smf
VTKFILES= squareLT.vtk cubeLT.vtk squareGrid.vtk cubeGrid.vtk

all: pre clean build run test

pre:
	@echo -e "\n\n--------------------------------------------------------"
	@echo    "Testing application $(APPNAME) in $(PWD)"
	@echo        "--------------------------------------------------------"

clean:
	rm -f $(LOG)
	@echo "--- Clean invoked"     >> $(LOG)
	make -f Makefile clean        >> $(LOG)
	rm -f $(BOUNDSMF) $(VTKFILES) >> $(LOG)
	rm -f unstructured?D structured?D
	@echo "--- Clean successful"  >> $(LOG)

build:
	@echo "--- Build invoked"      >> $(LOG)
	(make -f Makefile SPACEDIM=2 unstructured -B && mv unstructured unstructured2D) >> $(LOG)
	(make -f Makefile SPACEDIM=3 unstructured -B && mv unstructured unstructured3D) >> $(LOG)
	(make -f Makefile SPACEDIM=2   structured -B && mv structured structured2D    ) >> $(LOG)
	(make -f Makefile SPACEDIM=3   structured -B && mv structured structured3D    ) >> $(LOG) 
	@echo "--- Build successful"   >> $(LOG)

run:
	@echo "--- Run invoked"         >> $(LOG)
	./unstructured2D squareLT.smf   >> $(LOG) 
	./unstructured3D cubeLT.smf     >> $(LOG) 
	./structured2D   squareGrid.sgf >> $(LOG) 
	./structured3D   cubeGrid.sgf   >> $(LOG) 
	@echo "--- Run successful"      >> $(LOG)

test:
	@echo "--- Test invoked"         >> $(LOG)
	numdiff squareLT_boundary.smf   squareLT_boundary.ref.smf     >> $(LOG)
	numdiff cubeLT_boundary.smf     cubeLT_boundary.ref.smf       >> $(LOG)
	numdiff squareGrid_boundary.smf squareGrid_boundary.ref.smf   >> $(LOG)
	numdiff cubeGrid_boundary.smf   cubeGrid_boundary.ref.smf     >> $(LOG)
	numdiff squareLT.vtk   squareLT.ref.vtk                       >> $(LOG)
	numdiff cubeLT.vtk     cubeLT.ref.vtk                         >> $(LOG)
	numdiff squareGrid.vtk squareGrid.ref.vtk                     >> $(LOG)
	numdiff cubeGrid.vtk   cubeGrid.ref.vtk                       >> $(LOG)
	@echo "--- Test successful"       >> $(LOG)

triumph: clean
	@echo "(II) Testing application $(APPNAME) successful"
	rm -f $(LOG)

defeat:
	@echo "(EE) Application $(APPNAME) failed" >&2
	cat $(LOG) >&2
