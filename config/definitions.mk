# determine the hostname
HOST=$(shell hostname -s)

# directory of local definitions
LOCALDEFSDIR  = $(INSILICOROOT)/config

# filename of local definitions file
LOCALDEFSFILE = $(LOCALDEFSDIR)/$(addsuffix .$(HOST), config)

# default file
LOCALDEFSDEFAULT = config.default

# check existence of definitions file
# "test -f $(LOCALDEFSFILE)" has exit status 1 (failed) and or 0 (successful)
# thus "&& echo 'true'" is only executed, when file is found,
# otherwise call branches to "|| echo 'false'"
ifeq ($(shell test -f $(LOCALDEFSFILE) && echo 'true' || echo 'false'),false)
# put default filename since 
	LOCALDEFSFILE = $(LOCALDEFSDIR)/$(LOCALDEFSDEFAULT)
endif


# Compilation environment (default=GNU)
COMPENV = GNU
# Hack to check existence of INSILICOCOMPENV as variable
ifeq ($(shell test -z $(INSILICOCOMPENV) && echo 'unset' || echo 'set'),set)
    # check the value of the variable (currently only intel is interesting)
    ifeq ($(INSILICOCOMPENV),INTEL)
	COMPENV=INTEL
    endif
endif

# inform user
$(info Using definitions file $(LOCALDEFSFILE) and compiler environment $(COMPENV))

# include local definitions into this file
include $(LOCALDEFSFILE)

