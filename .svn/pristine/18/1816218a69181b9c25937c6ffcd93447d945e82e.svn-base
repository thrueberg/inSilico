APPNAME=drivenCavity
LOG=$(APPNAME).log

VTKFILES= cavity.vtk

all: pre clean build run test

pre:
	@echo -e "\n\n--------------------------------------------------------"
	@echo    "Testing application $(APPNAME) in $(PWD)"
	@echo        "--------------------------------------------------------"

clean:
	rm -f $(LOG)
	@echo "--- Clean invoked"     >> $(LOG)
	make -f Makefile clean        >> $(LOG)
	rm -f $(VTKFILES)             >> $(LOG)
	@echo "--- Clean successful"  >> $(LOG)

build:
	@echo "--- Build invoked"                 >> $(LOG)
	(make -f Makefile drivenCavity DEBUG=NO ) >> $(LOG)
	@echo "--- Build successful"              >> $(LOG)

run:
	@echo "--- Run invoked"         >> $(LOG)
	./drivenCavity inputRef.dat     >> $(LOG)
	@echo "--- Run successful"      >> $(LOG)

test:
	@echo "--- Test invoked"           >> $(LOG)
	numdiff cavity.ref.vtk cavity.vtk  >> $(LOG)
	@echo "--- Test successful"        >> $(LOG)

triumph: clean
	@echo "(II) Testing application $(APPNAME) successful"
	rm -f $(LOG)

defeat:
	@echo "(EE) Application $(APPNAME) failed" >&2
	cat $(LOG) >&2
