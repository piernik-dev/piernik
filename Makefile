################################################################################
# 
# This Makefile can be used to rebuild one or more object directories at a time.
#
# Usage:
#
#   'make'                 # rebuilds all obj*/piernik executables calling make
#                            for all obj* directories
#   'make obj_abc obj_xyz' # rebuilds just the two specified directories
#   'make resetup'         # calls setup for all obj* directories, relying on 
#                            obj*/.setup.call and current .setuprc file.
#   'RS=1 make obj_abc'    # calls setup only for obj_abc directory
#   'make obj_abc RS=1'    # same as above
#   'make clean'           # calls make clean for all obj* directories
#   'make allsetup'        # creates object directories for all valid problems,
#                            but does not compile them
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
	@$(MAKE) -k -C $@ PNAME=$@ && $(ECHO) -e "\033[32;1m"$@ ready"\033[0m" || echo -e "\033[31;1m"$@ failed"\033[0m"
endif
endif

resetup:
	@RS=1 $(MAKE) -k all

clean:
	@CL=1 $(MAKE) -k all

allsetup:
	for i in problems/* ../problems/* ; do \
		if [ ! -e $$i/OBSOLETE ] ; then \
			if [ $$( dirname $$( dirname $$i ) ) == "." ] ; then \
				nm=$$( basename $$i ); \
			else \
				nm="../"$$i; \
			fi; \
			./setup $$nm -o "_"$$( basename $$i )"_" --nocompile && sed -i 's/ --nocompile//' "obj__"$$( basename $$i )"_/"{.setup.call,Makefile,env.dat,version.F90} & \
		fi; \
	done
