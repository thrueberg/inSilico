APPNAME=doFHandler
LOG=$(APPNAME).log

DATFILES= sparsity.1.dat sparsity.2.dat sparsity.3.dat

all: pre clean build run test

pre:
	@echo -e "\n\n--------------------------------------------------------"
	@echo    "Testing application $(APPNAME) in $(PWD)"
	@echo        "--------------------------------------------------------"

clean:
	rm -f $(LOG)
	@echo "--- Clean invoked"     >> $(LOG)
	make -f Makefile clean        >> $(LOG)
	rm -f $(DATFILES)             >> $(LOG)
	rm -f doFHandler?             >> $(LOG)
	@echo "--- Clean successful"  >> $(LOG)

build:
	@echo "--- Build invoked"      >> $(LOG)
	(make -f Makefile DEGREE=1 doFHandler -B && mv doFHandler doFHandler1) >> $(LOG)
	(make -f Makefile DEGREE=2 doFHandler -B && mv doFHandler doFHandler2) >> $(LOG)
	(make -f Makefile DEGREE=3 doFHandler -B && mv doFHandler doFHandler3) >> $(LOG)
	@echo "--- Build successful"   >> $(LOG)

run:
	@echo "--- Run invoked"       >> $(LOG)
	./doFHandler1  square_20.smf  >> $(LOG) 
	./doFHandler2  square_20.smf  >> $(LOG) 
	./doFHandler3  square_20.smf  >> $(LOG) 
	@echo "--- Run successful"    >> $(LOG)

test:
	@echo "--- Test invoked"                   >> $(LOG)
	numdiff sparsity.1.dat sparsity.1.ref.dat  >> $(LOG)
	numdiff sparsity.2.dat sparsity.2.ref.dat  >> $(LOG)
	numdiff sparsity.3.dat sparsity.3.ref.dat  >> $(LOG)
	@echo "--- Test successful"                >> $(LOG)

triumph: clean
	@echo "(II) Testing application $(APPNAME) successful"
	rm -f $(LOG)

defeat:
	@echo "(EE) Application $(APPNAME) failed" >&2
	cat $(LOG) >&2
