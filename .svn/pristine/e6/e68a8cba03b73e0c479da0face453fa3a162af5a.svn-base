# the list of directories 
SUBDIRS := meshGeneration converter

# generate the tutorials
tools:
	for i in $(SUBDIRS); do \
	$(MAKE) -C $$i;\
	done



all: tools

