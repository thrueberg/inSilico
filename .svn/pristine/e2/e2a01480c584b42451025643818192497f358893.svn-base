# the list of directories 
SUBDIRS := ${sort ${dir ${wildcard */}}}	

# generate the tutorials
meshgeneration:
	for i in $(SUBDIRS); do \
	$(MAKE) -C $$i;\
	done



all: meshgeneration

