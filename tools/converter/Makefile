# the list of directories 
SUBDIRS := gmsh2smf sgf2vtk sgfRefine smf2gp smf2vtk smfMap smfReverse

# generate the conversion tools
converter: 
	for i in $(SUBDIRS); do \
	$(MAKE) -C $$i;\
  done

all: converter
