# the list of directories 
SUBDIRS := ${sort ${dir ${wildcard ??-*/}}}	

# need tools to be in a good state
tools:
	cd $(INSILICOROOT)/tools && $(MAKE) all


# generate the tutorials
tutorials: 
	for i in $(SUBDIRS); do \
	$(MAKE) -C $$i;\
	done


# run test applications
appTest: tools
	for i in $(SUBDIRS); do \
	$(MAKE) -C $$i appTest;\
	done


