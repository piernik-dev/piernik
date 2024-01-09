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
#   'make obj_abc CL=1'    # calls make clean for specified directories
#   'make allsetup'        # creates object directories for all valid problems,
#                            but does not compile them
#   'make ctags'           # recreate ctags for {src,problems}
#   'make dep [P=problem]' # create and show dependency graph
#                            $P defaults to mcrwind
#   'make qa'              # run qa.py on all F90 files in src and problems
#                            directories
#   'make pep8'            # run pep8 on all Python scripts, ignore long lines
#                            (obsoleted by pycodestyle)
#   'make pycodestyle'     # run pycodestyle on all Python scripts, ignore long lines
#   'make chk_err_msg'     # check filenames in error messages
#   'make doxy'            # generate/updare Doxygen documentation
#   'make gold'            # run the gold tests from ./jenkins directory
#   'make gold-serial'     # run the gold tests from ./jenkins directory in serial mode
#   'make gold-clean'      # remove files after gold test
#
# Resetup will also call make for the object directories, unless you've
# specified --nocompile either in your .setuprc* files or it was stored in
# .setup.call file.
#
################################################################################

MAKEFLAGS += -s

ALLOBJ = $(wildcard obj*)

ECHO ?= /bin/echo

.PHONY: $(ALLOBJ) check dep qa pep8 pycodestyle doxy chk_err_msg gold gold-serial gold-clean

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
	( for i in $$( find problems/* -type d ) ; do \
		if [ ! -e $$i/.skipauto ] ; then \
			pnm=$$( echo $$i | sed 's-^problems/--' ); \
			onm=$${pnm//\//___}; \
			./setup $$pnm -o "A_"$$onm --nocompile && \
				for f in .setup.call Makefile env.dat version.F90 ; do \
					of="obj_A_"$$onm"/"$$f; \
					[ -f $$of ] && sed -i 's/ --nocompile//' $$of; \
				done & \
			sleep .1; \
		fi; \
	done; \
	wait )

qa:
	./bin/qa.py $$( git ls-files | grep -vE "^(compilers/tests|doc/general)" | grep "\.F90$$" )

QA:
	make -k  chk_err_msg chk_lic_hdr pycodestyle qa

pep8: pycodestyle

pycodestyle:
	echo 'Pycodestyle check (--ignore=E501,E722,W504,W605)'
	pycodestyle `git ls-files | grep '\.py$$'` bin/gdf_distance bin/ask_jenkins --ignore=E501,E722,W504,W605

chk_err_msg:
	echo Check filenames in error messages
	./bin/checkmessages.sh

chk_lic_hdr:
	./bin/check_license_headers.sh

gold:
	./jenkins/gold_test_list.sh

gold-serial:
	SERIAL=1 ./jenkins/gold_test_list.sh

gold-clean:
	\rm -rf jenkins/workspace/*

doxy:
	doxygen piernik.doxy

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
