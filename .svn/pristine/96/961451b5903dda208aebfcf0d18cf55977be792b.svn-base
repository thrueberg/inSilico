APPNAME=compressible
LOG=$(APPNAME).log

VTKFILES= quad.020.00??.vtk quad.020.?sol.vtk
DATFILES= compSolOut?.dat


all: pre clean build run test

pre:
	@echo -e "\n\n--------------------------------------------------------"
	@echo    "Testing application $(APPNAME) in $(PWD)"
	@echo        "--------------------------------------------------------"

clean:
	rm -f $(LOG)
	@echo "--- Clean invoked"     >> $(LOG)
	make -f Makefile clean        >> $(LOG)
	rm -f $(VTKFILES) $(DATFILES) >> $(LOG)
	@echo "--- Clean successful"  >> $(LOG)

build:
	@echo "--- Build invoked"                               >> $(LOG)
	(make -f Makefile DEBUG=NO SPACEDIM=2 compressible -B ) >> $(LOG)
	@echo "--- Build successful"                            >> $(LOG)

run:
	@echo "--- Run invoked"         >> $(LOG)
	(./compressible quad.020.smf inputCompRefD.dat \
	&& cp quad.020.0005.vtk quad.020.Dsol.vtk)   >> compSolOutD.dat
	(./compressible quad.020.smf inputCompRefN.dat \
	&& cp quad.020.0005.vtk quad.020.Nsol.vtk)   >> compSolOutN.dat
	@echo "--- Run successful"      >> $(LOG)

test:
	@echo "--- Test invoked"         >> $(LOG)
	numdiff quad.020.Dsol.vtk quad.020.Dref.vtk >> $(LOG)
	numdiff quad.020.Nsol.vtk quad.020.Nref.vtk >> $(LOG)
	numdiff compSolOutD.dat   compRefOutD.dat   >> $(LOG)
	numdiff compSolOutN.dat   compRefOutN.dat   >> $(LOG)
	@echo "--- Test successful"       >> $(LOG)

triumph: clean
	@echo "(II) Testing application $(APPNAME) successful"
	rm -f $(LOG)

defeat:
	@echo "(EE) Application $(APPNAME) failed" >&2
	cat $(LOG) >&2
