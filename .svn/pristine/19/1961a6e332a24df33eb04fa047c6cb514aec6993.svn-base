APPNAME=areaVolume
LOG=$(APPNAME).log

all: pre clean build run test

pre:
	@echo -e "\n\n--------------------------------------------------------"
	@echo    "Testing application $(APPNAME) in $(PWD)"
	@echo        "--------------------------------------------------------"

clean:
	rm -f $(LOG)
	@echo "--- Clean invoked"     >> $(LOG)
	make -f Makefile clean        >> $(LOG)
	rm -f measure.dat             >> $(LOG)
	@echo "--- Clean successful"  >> $(LOG)

build:
	@echo "--- Build invoked"      >> $(LOG)
	(make -f Makefile ) >> $(LOG)
	@echo "--- Build successful"   >> $(LOG)

run:
	@echo "--- Run invoked"         >> $(LOG)
	./areaVolume                    >> measure.dat
	@echo "--- Run successful"      >> $(LOG)

test:
	@echo "--- Test invoked"            >> $(LOG)
	numdiff measure.dat measure.ref.dat >> $(LOG)
	@echo "--- Test successful"       >> $(LOG)

triumph: clean
	@echo "(II) Testing application $(APPNAME) successful"
	rm -f $(LOG)

defeat:
	@echo "(EE) Application $(APPNAME) failed" >&2
	cat $(LOG) >&2
