# the command to execute doxygen
DOXY   = doxygen
# the pictures directory
PICDIR = pictures
# the list of directories inside picture
SUBDIRS := ${sort ${dir ${wildcard $(PICDIR)/*/ }}}	
# directory of references with Makefile for more pictures
REFDIR  := $(INSILICOROOT)/reference

all: doxy

# generate the pictures
pics: 
	for i in $(SUBDIRS); do \
	$(MAKE) -C $$i;\
	done;\
	$(MAKE) -C $(REFDIR) pics;\

# generate the doxygen documentation
doxy: pics
	$(DOXY)
