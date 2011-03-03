################################################################################
# 
# This Makefile can be used to rebuild one or more object directories at a time.
#
# Usage:
#   'make'                 # rebuilds all obj*/piernik executables calling make
#                            for all obj* directories
#   'make obj_abc obj_xyz' # rebuilds just the two specified directories
#   'make resetup'         # calls setup for all obj* directories, relying on 
#                            obj*/.setup.call and current .setuprc file.
#   'RS=1 make obj_abc'    # calls setup only for obj_abc directory
#   'make obj_abc RS=1'    # same as above
#
# Resetup will also call make for the object directories, unless you've 
# specified --nocompile either in your .setuprc* files or it was stored in 
# .setup.call file.
#
# ToDo: drop support for reading setup call from obj*/env.dat as soon as 
#       everyone adopts .setup.call-aware setup script.
#
################################################################################

MAKEFLAGS += -s -l8

ALLOBJ = $(wildcard obj*)

ECHO ?= /bin/echo

.PHONY: $(ALLOBJ)

all: $(ALLOBJ)

*:
ifeq ("$(RS)","1")
	@if [ -e $@/.setup.call ] ; then \
		eval `grep -v "^#" $@/.setup.call`; \
	elif [ -e $@/env.dat ] ; then \
		eval `head -n 1 $@/env.dat`; \
	else \
		$(ECHO) -e "\033[31;1mDon't know how to resetup '"$@"'\033[0m"; \
	fi
else
ifeq ("$(CL)","1")
	$(MAKE) -k -C $@ clean && $(ECHO) -e "\033[32;1m"Cleaned $@"\033[0m" || echo -e "\033[31;1m"Unable to clean $@"\033[0m"
else
	@$(MAKE) -k -C $@ && $(ECHO) -e "\033[32;1m"$@ ready"\033[0m" || echo -e "\033[31;1m"$@ failed"\033[0m"
endif
endif

resetup:
	@RS=1 $(MAKE) all

clean:
	@CL=1 $(MAKE) all
