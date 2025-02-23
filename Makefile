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
#                            $P defaults to testing_and_debuging/chimaera
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
#   'make CI'              # run all checks locally (replaces Jenkins to some extent)
#
# Resetup will also call make for the object directories, unless you've
# specified --nocompile either in your .setuprc* files or it was stored in
# .setup.call file.
#
################################################################################

MAKEFLAGS += -s

ALLOBJ = $(wildcard obj*)
GREEN = "\033[32;1m"
RED = "\033[31;1m"
RESET = "\033[0m"
PASSED = $(GREEN)passed$(RESET)
FAILED = $(RED)failed$(RESET)
ARTIFACTS = "./jenkins/workspace/artifacts"

ECHO ?= /bin/echo

.PHONY: $(ALLOBJ) check dep qa pep8 pycodestyle doxy chk_err_msg gold gold-serial gold-clean CI allgold artifact_tests gold_CRESP gold_mcrtest gold_mcrwind gold_MHDsedovAMR gold_resist gold_streaming_instability view_dep py3 noHDF5

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
	./bin/qa.py -q $$( git ls-files | grep -vE "^(compilers/tests|doc/general)" | grep "\.F90$$" ) && \
	echo -e "QA checks "$(PASSED) || \
	( echo -e "QA checks "$(FAILED) && exit 1 )

QA:
	make -k chk_err_msg chk_lic_hdr pycodestyle qa

pep8: pycodestyle

pycodestyle:
	TSTNAME="Pycodestyle check "; \
	REMARK=" (--ignore=E501,E722,W504,W605)"; \
	pycodestyle `git ls-files | grep '\.py$$'` bin/gdf_distance bin/ask_jenkins --ignore=E501,E722,W504,W605 && \
		echo -e "$$TSTNAME"$(PASSED)"$$REMARK" ||\
		( echo -e "$$TSTNAME"$(FAILED)"$$REMARK" && exit 1 )

chk_err_msg:
	./bin/checkmessages.sh && \
		echo -e "Message checks "$(PASSED) || \
		( echo -e "Message checks "$(FAILED)": incorrect file references found" && exit 1 )

chk_lic_hdr:
	./bin/check_license_headers.sh && \
		echo -e "License header checks "$(PASSED) || \
		( echo -e "License header checks "$(FAILED)": exceptional license headers found" && exit 1 )

gold:
	./jenkins/gold_test_list.sh

gold-serial:
	SERIAL=1 ./jenkins/gold_test_list.sh

gold-clean:
	\rm -rf jenkins/workspace/*

doxy:
	doxygen piernik.doxy

ifndef P
P = "testing_and_debuging/chimaera"
endif
dep:
	OTMPDIR=$$(mktemp -d obj_XXXXXX) ;\
	PROBLEM=$P ;\
	GRAPH="dep.png" ;\
	[ ! -e $(ARTIFACTS) ] && mkdir -p $(ARTIFACTS) ;\
	[ -e $(ARTIFACTS)"/$$GRAPH" ] && rm $(ARTIFACTS)"/$$GRAPH" || true ;\
	./setup $$PROBLEM -n -o $${OTMPDIR//obj_/} > setup.stdout ;\
	if [ -e $$OTMPDIR ] ; then \
		mv setup.stdout $$OTMPDIR ;\
		$(MAKE) -k -C $$OTMPDIR $$GRAPH ;\
		mv $$OTMPDIR"/"$$GRAPH $(ARTIFACTS)/ &&\
			rm -r $$OTMPDIR "runs/"$$( basename $${PROBLEM} )"_"$${OTMPDIR//obj_/} ;\
	fi ;\
	if [ -e $(ARTIFACTS)"/$$GRAPH" ] ; then \
		$(ECHO) -e "Dependency graph for the "$$PROBLEM" problem stored in $(ARTIFACTS)/"$$GRAPH" . Test "$(PASSED) ;\
	else \
		[ -e $$OTMPDIR/setup.stdout ] && ( cat $$OTMPDIR/setup.stdout ; rm -r $$OTMPDIR ) ;\
		$(ECHO) -e "Dependency graph creation "$(FAILED) && exit 1 ;\
	fi

view_dep: dep
	which display > /dev/null 2> /dev/null && [ -e $(ARTIFACTS)"/dep.png" ] && display $(ARTIFACTS)"/dep.png" || true

py3:
	OTMPDIR=$$(mktemp -d obj_XXXXXX) ;\
	python3 ./python/piernik_setup.py maclaurin -n -o $${OTMPDIR//obj_/} > $(ARTIFACTS)/py3.setup.stdout && \
		( rm -r $${OTMPDIR} runs/maclaurin_$${OTMPDIR//obj_/} ; echo -e "Python 3 test "$(PASSED) ) || \
		( rm -r $${OTMPDIR} ; echo -e "Python 3 test "$(FAILED) && exit 1 )

noHDF5:
	OTMPDIR=$$(mktemp -d obj_XXXXXX) ;\
	./setup crtest --param problem.par.build -d I_KNOW_WHAT_I_AM_DOING -o $${OTMPDIR//obj_/} > $(ARTIFACTS)/noHDF5.setup.stdout && \
		( rm -r $${OTMPDIR} runs/crtest_$${OTMPDIR//obj_/} ; echo -e "No HDF5 test "$(PASSED) ) || \
		( rm -r $${OTMPDIR} ; echo -e "No HDF5 test "$(FAILED) && exit 1 )

artifact_tests: dep py3 noHDF5

gold_CRESP:
	./jenkins/gold_test.sh ./jenkins/gold_configs/mcrtest_CRESP.config > ./jenkins/workspace/CRESP.gold_stdout 2> ./jenkins/workspace/CRESP.gold_stderr && \
		echo -e "CRESP test "$(PASSED) || \
		( echo -e "CRESP test "$(FAILED)" (more details in ./jenkins/workspace/CRESP.gold_std*)" && exit 1 )
	# Fails on Fedora 41

gold_mcrtest:
	./jenkins/gold_test.sh ./jenkins/gold_configs/mcrtest_new.config > ./jenkins/workspace/mcrtest.gold_stdout 2> ./jenkins/workspace/mcrtest.gold_stderr && \
		echo -e "mcrtest "$(PASSED) || \
		( echo -e "mcrtest "$(FAILED)" (more details in ./jenkins/workspace/mcrtest.gold_std*)" && exit 1 )

gold_mcrwind:
	./jenkins/gold_test.sh ./jenkins/gold_configs/mcrwind.config > ./jenkins/workspace/mcrwind.gold_stdout 2> ./jenkins/workspace/mcrwind.gold_stderr && \
		echo -e "mcrwind "$(PASSED) || \
		( echo -e "mcrwind "$(FAILED)" (more details in ./jenkins/workspace/mcrwind.gold_std*)" && exit 1 )

gold_MHDsedovAMR:
	./jenkins/gold_test.sh ./jenkins/gold_configs/MHDsedovAMR.config > ./jenkins/workspace/MHDsedovAMR.gold_stdout 2> ./jenkins/workspace/MHDsedovAMR.gold_stderr && \
		echo -e "MHDsedovAMR test "$(PASSED) || \
		( echo -e "MHDsedovAMR test "$(FAILED)" (more details in ./jenkins/workspace/MHDsedovAMR.gold_std*)" && exit 1 )

gold_resist:
	./jenkins/gold_test.sh ./jenkins/gold_configs/resist.config > ./jenkins/workspace/resist.gold_stdout 2> ./jenkins/workspace/resist.gold_stderr && \
		echo -e "Resistivity test "$(PASSED) || \
		( echo -e "Resistivity test "$(FAILED)" (more details in ./jenkins/workspace/resist.gold_std*)" && exit 1 )

gold_streaming_instability:
	./jenkins/gold_test.sh ./jenkins/gold_configs/streaming_instability.config > ./jenkins/workspace/streaming_instability.gold_stdout 2> ./jenkins/workspace/streaming_instability.gold_stderr && \
		echo -e "Streaming instability test "$(PASSED) || \
		( echo -e "Streaming instability test "$(FAILED)" (more details in ./jenkins/workspace/streaming_instability.gold_std*)" && exit 1 )

allgold: gold_CRESP gold_mcrtest gold_mcrwind gold_MHDsedovAMR gold_resist gold_streaming_instability

CI: QA
	[ -e $(ARTIFACTS) ] && rm -rf $(ARTIFACTS) || true
	[ ! -e $(ARTIFACTS) ] && mkdir -p $(ARTIFACTS)
	$(MAKE) -k artifact_tests
	$(MAKE) -k allgold
