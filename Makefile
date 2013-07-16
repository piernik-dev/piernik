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
#   'make ctags'           # recreate ctags for {src,problems}
#   'make dep [P=problem]  # create and show dependency graph
#                            $P defaults to mcrwind
#
# Resetup will also call make for the object directories, unless you've 
# specified --nocompile either in your .setuprc* files or it was stored in 
# .setup.call file.
#
################################################################################

MAKEFLAGS += -s -l8

ALLOBJ = $(wildcard obj*)

ECHO ?= /bin/echo

.PHONY: $(ALLOBJ) check dep

all: $(ALLOBJ)

*:
ifeq ("$(RS)","1")
	@if [ -e $@/.setup.call ] ; then \
		eval `grep -v "^#" $@/.setup.call`; \
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

ctags:
	ctags -R {src,problems} --fortran-kinds=+iL

resetup:
	@RS=1 $(MAKE) -k all

clean:
	@CL=1 $(MAKE) -k all

allsetup:
	for i in {,../}problems/* ; do \
		if [ ! -e $$i/OBSOLETE ] ; then \
			if [ $$( dirname $$( dirname $$i ) ) == "." ] ; then \
				nm=$$( basename $$i ); \
			else \
				nm="../"$$i; \
			fi; \
			./setup $$nm -o "A_"$$( basename $$i ) --nocompile && sed -i 's/ --nocompile//' "obj_A_"$$( basename $$i )"/"{.setup.call,Makefile,env.dat,version.F90} & \
		fi; \
	done

check:
	TMPDIR=$$(mktemp -d /dev/shm/test_XXXXXX);\
	bitten-slave -d . --build-dir $$TMPDIR -k bitten/trunk.mcrtest.xml ;\
	rm -rf $$TMPDIR

ifndef P
P = "mcrwind"
endif
dep:
	TMPDIR=$$(mktemp XXXXXX) ;\
	rm $$TMPDIR ;\
	OTMPDIR="obj_"$$TMPDIR ;\
	PROBLEM=$P ;\
	GRAPH="dep.png" ;\
	rm $$GRAPH 2> /dev/null ;\
	./setup $$PROBLEM -n -o $$TMPDIR | grep -v "skipped" ;\
	if [ -e $$OTMPDIR ] ; then \
		$(MAKE) -k -C $$OTMPDIR $$GRAPH ;\
		mv $$OTMPDIR"/"$$GRAPH . ;\
		rm -r $$OTMPDIR "runs/"$${PROBLEM}"_"$$TMPDIR ;\
	fi ;\
	which display > /dev/null 2> /dev/null && [ -e $$GRAPH ] && display $$GRAPH ;\
	if [ -e $$GRAPH ] ; then \
		$(ECHO) "Dependency graph for the "$$PROBLEM" problem stored in "$$GRAPH ;\
	else \
		$(ECHO) "Cannot create dependency graph" ;\
	fi
