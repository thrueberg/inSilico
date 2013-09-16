# the list of directories 
SUBDIRS := ${sort ${dir ${wildcard */}}}	

# generate the tutorials
tutorials: 
	for i in $(SUBDIRS); do \
	$(MAKE) -C $$i;\
  done

