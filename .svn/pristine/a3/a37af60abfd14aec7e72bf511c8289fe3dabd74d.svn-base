APPNAME=linearElastic
LOG=$(APPNAME).log

VTKFILES= quad.???.vtk cube.???.vtk
DATFILES= linearElastic?D.dat


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
	rm -f linearElastic2D linearElastic3D >> $(LOG)
	@echo "--- Clean successful"  >> $(LOG)

build:
	@echo "--- Build invoked"                                >> $(LOG)
	(make -f Makefile DEBUG=NO SPACEDIM=2 linearElastic -B \
	 && cp linearElastic linearElastic2D)                    >> $(LOG)
	(make -f Makefile DEBUG=NO SPACEDIM=3 linearElastic -B \
	 && cp linearElastic linearElastic3D)                    >> $(LOG)
	@echo "--- Build successful"                             >> $(LOG)

run:
	@echo "--- Run invoked"         >> $(LOG)
	./exec.sh 2 
	./exec.sh 3
	@echo "--- Run successful"      >> $(LOG)

test:
	@echo "--- Test invoked"         >> $(LOG)
	numdiff linearElastic2D.dat linearElastic2D.ref.dat >> $(LOG)
	numdiff linearElastic3D.dat linearElastic3D.ref.dat >> $(LOG)
	@echo "--- Test successful"       >> $(LOG)

triumph: clean
	@echo "(II) Testing application $(APPNAME) successful"
	rm -f $(LOG)

defeat:
	@echo "(EE) Application $(APPNAME) failed" >&2
	cat $(LOG) >&2
