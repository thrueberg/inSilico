# the list of directories 
SUBDIRS := ${sort ${dir ${wildcard */}}}	

# generate the conversion tools
converter: 
	for i in $(SUBDIRS); do \
	$(MAKE) -C $$i;\
  done

