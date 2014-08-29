APPNAME=convection
LOG=$(APPNAME).log

DATFILE=convectionOut
REFFILE=convectionOutRef

all: pre clean build run test

pre:
	@echo -e "\n\n--------------------------------------------------------"
	@echo    "Testing application $(APPNAME) in $(PWD)"
	@echo        "--------------------------------------------------------"

clean:
	rm -f $(LOG)
	@echo "--- Clean invoked"     >> $(LOG)
	rm -f $(DATFILE)              >> $(LOG)
	rm -f $(APPNAME)              >> $(LOG)
	@echo "--- Clean successful"  >> $(LOG)

build:
	@echo "--- Build invoked"      >> $(LOG)
	make -f Makefile $(APPNAME) -B DEBUG=NO  >> $(LOG)
	@echo "--- Build successful"   >> $(LOG)

run:
	@echo "--- Run invoked"       >> $(LOG)
	./$(APPNAME) convectionInput.dat > $(DATFILE)
	@echo "--- Run successful"    >> $(LOG)

test:
	@echo "--- Test invoked"                   >> $(LOG)
	numdiff $(DATFILE) $(REFFILE)              >> $(LOG)
	@echo "--- Test successful"                >> $(LOG)

triumph: clean
	@echo "(II) Testing application $(APPNAME) successful"
	rm -f $(LOG)

defeat:
	@echo "(EE) Application $(APPNAME) failed" >&2
	cat $(LOG) >&2
