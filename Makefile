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
BLUE = "\033[34;1m"
RESET = "\033[0m"
PASSED = $(GREEN)passed$(RESET)
FAILED = $(RED)failed$(RESET)
ARTIFACTS = "./jenkins/artifacts"
GOLDSPACE = "./jenkins/workspace"

ECHO ?= /bin/echo
RM ?= /bin/rm

.PHONY: $(ALLOBJ) check dep qa pep8 pycodestyle doxy chk_err_msg gold gold-serial gold-clean CI allgold artifact_tests gold_CRESP gold_mcrtest gold_mcrwind gold_MHDsedovAMR gold_resist gold_streaming_instability view_dep py3 noHDF5 I64 IOv2

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
	$(MAKE) -k -C $@ clean && $(ECHO) -e "\033[32;1m"Cleaned $@"\033[0m" || $(ECHO) -e "\033[31;1m"Unable to clean $@"\033[0m"
else
	@$(MAKE) -k -C $@ PNAME=$@ && $(ECHO) -e "\033[32;1m"$@ ready"\033[0m" || $(ECHO) -e "\033[31;1m"$@ failed"\033[0m"
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
	echo -e "  qa.py checks "$(PASSED) || \
	( $(ECHO) -e "  qa.py checks "$(FAILED) && exit 1 )

QA:
	$(ECHO) -e $(BLUE)"Starting QA checks ..."$(RESET)
	$(MAKE) -k chk_err_msg chk_lic_hdr pycodestyle qa
	$(ECHO) -e $(BLUE)"QA checks "$(PASSED)"."

pep8: pycodestyle

pycodestyle:
	TSTNAME="  Pycodestyle check "; \
	REMARK=" (with --ignore=E501,E722,W504,W605)"; \
	pycodestyle `git ls-files | grep '\.py$$'` bin/gdf_distance bin/ask_jenkins --ignore=E501,E722,W504,W605 && \
		$(ECHO) -e "$$TSTNAME"$(PASSED)"$$REMARK" ||\
		( $(ECHO) -e "$$TSTNAME"$(FAILED)"$$REMARK" && exit 1 )

chk_err_msg:
	./bin/checkmessages.sh && \
		$(ECHO) -e "  Message checks "$(PASSED) || \
		( $(ECHO) -e "  Message checks "$(FAILED)": incorrect file references found" && exit 1 )

chk_lic_hdr:
	./bin/check_license_headers.sh && \
		$(ECHO) -e "  License header checks "$(PASSED) || \
		( $(ECHO) -e "  License header checks "$(FAILED)": exceptional license headers found" && exit 1 )

gold:
	./jenkins/gold_test_list.sh

gold-serial:
	SERIAL=1 ./jenkins/gold_test_list.sh

gold-clean:
	$(RM) -rf $(GOLDSPACE)/* $(ARTIFACTS)

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
		$(ECHO) -e "  Dependency test "$(PASSED)". The graph for the $$PROBLEM problem was stored as $(ARTIFACTS)/$$GRAPH" ;\
	else \
		[ -e $$OTMPDIR/setup.stdout ] && ( cat $$OTMPDIR/setup.stdout ; rm -r $$OTMPDIR ) ;\
		$(ECHO) -e "  Dependency graph creation "$(FAILED) && exit 1 ;\
	fi

view_dep: dep
	which display > /dev/null 2> /dev/null && [ -e $(ARTIFACTS)"/dep.png" ] && display $(ARTIFACTS)"/dep.png" || true

py3:
	OTMPDIR=$$(mktemp -d obj_XXXXXX) ;\
	python3 ./python/piernik_setup.py maclaurin -n -o $${OTMPDIR//obj_/} > $(ARTIFACTS)/py3.setup.stdout && \
		( rm -r $${OTMPDIR} runs/maclaurin_$${OTMPDIR//obj_/} ; $(ECHO) -e "  Python 3 test "$(PASSED) ) || \
		( rm -r $${OTMPDIR} ; $(ECHO) -e "  Python 3 test "$(FAILED) && exit 1 )

noHDF5:
	OTMPDIR=$$(mktemp -d obj_XXXXXX) ;\
	./setup crtest --param problem.par.build -d I_KNOW_WHAT_I_AM_DOING -o $${OTMPDIR//obj_/} > $(ARTIFACTS)/noHDF5.setup.stdout && \
		( rm -r $${OTMPDIR} runs/crtest_$${OTMPDIR//obj_/} ; $(ECHO) -e "  No HDF5 test "$(PASSED) ) || \
		( rm -r $${OTMPDIR} ; $(ECHO) -e "  No HDF5 test "$(FAILED) && exit 1 )

I64:
	OTMPDIR=$$(mktemp -d obj_XXXXXX) ;\
	./setup $(P) -o $${OTMPDIR//obj_/} --f90flags="-fdefault-integer-8 -Werror=conversion" > $(ARTIFACTS)/I64.setup.stdout && \
		( rm -r $${OTMPDIR} runs/chimaera_$${OTMPDIR//obj_/} ; $(ECHO) -e "  I64 test "$(PASSED) ) || \
		( rm -r $${OTMPDIR} ; $(ECHO) -e "  I64 test "$(FAILED) && exit 1 )

IOv2:
	OTMPDIR=$$(mktemp -d obj_XXXXXX) ;\
	RUNDIR=runs/advection_test_$${OTMPDIR//obj_/} ;\
	./setup advection_test -p problem.par.restart_test_v2 -o $${OTMPDIR//obj_/} > $(ARTIFACTS)/IOv2.setup.stdout && \
		( cd $${RUNDIR} ;\
			$(ECHO) -e "run_id = ts1\n  Start:   t = 0.0, nproc = 1" ;\
			mpiexec -n 1 ./piernik > ts1.out ;\
			$(ECHO) -e "  Restart: t = 1.0, nproc = 1" ;\
			mpiexec -n 1 ./piernik -n '&END_CONTROL tend = 2.0 /' >> ts1.out ;\
			$(ECHO) -e "  Finish:  t = 2.0\n\nrun_id = ts2\n  Start:   t = 0.0, nproc = 5" ;\
			mpiexec -n 5 ./piernik -n '&OUTPUT_CONTROL run_id = "ts2" /' > ts2.out ;\
			$(ECHO) -e "  Restart: t = 1.0, nproc = 3" ;\
			mpiexec -n 3 ./piernik -n '&END_CONTROL tend = 2.0 /' -n '&OUTPUT_CONTROL run_id = "ts2" /' >> ts2.out ;\
			$(ECHO) -e "  Finish:  t = 2.0" ;\
		) >> $(ARTIFACTS)/IOv2.setup.stdout && \
		./bin/gdf_distance $${RUNDIR}/moving_pulse_ts{1,2}_0002.h5 | tee $${RUNDIR}/compare.log >> $(ARTIFACTS)/IOv2.setup.stdout && \
		( [ $$( grep "^Total difference between" $${RUNDIR}/compare.log | awk '{print $$NF}' ) == 0 ] || exit 1 ) && \
		( rm -r $${OTMPDIR} $${RUNDIR} ; $(ECHO) -e "  IO v2 test "$(PASSED) ) || \
		( rm -r $${OTMPDIR} ; $(ECHO) -e "  IO v2 test "$(FAILED) && exit 1 )

artifact_tests:
	$(ECHO) -e $(BLUE)"Starting artifact tests ..."$(RESET)
	[ -e $(ARTIFACTS) ] && rm -rf $(ARTIFACTS) || true
	[ ! -e $(ARTIFACTS) ] && mkdir -p $(ARTIFACTS)
	$(MAKE) -k dep py3 noHDF5 I64 IOv2
	$(ECHO) -e $(BLUE)"All artifact tests "$(PASSED)". Details can be found in "$(ARTIFACTS)" directory."

define run_gold_test
	./jenkins/gold_test.sh ./jenkins/gold_configs/$(1).config > $(GOLDSPACE)/$(1).gold_stdout 2> $(GOLDSPACE)/$(1).gold_stderr && \
		$(ECHO) -e "  $(1) test "$(PASSED) || \
		( $(ECHO) -e "  $(1) test "$(FAILED)" (more details in $(GOLDSPACE)/$(1).gold_std*)" && exit 1 )
endef

gold_CRESP:
	$(call run_gold_test,mcrtest_CRESP)

gold_mcrtest:
	$(call run_gold_test,mcrtest_new)

gold_mcrwind:
	$(call run_gold_test,mcrwind)

gold_MHDsedovAMR:
	$(call run_gold_test,MHDsedovAMR)

gold_resist:
	$(call run_gold_test,resist)

gold_streaming_instability:
	$(call run_gold_test,streaming_instability)

allgold:
	$(ECHO) -e $(BLUE)"Starting gold tests ..."$(RESET)
	$(MAKE) -k gold_CRESP gold_mcrtest gold_mcrwind gold_MHDsedovAMR gold_resist gold_streaming_instability
	$(ECHO) -e $(BLUE)"All gold tests "$(PASSED)". Details can be found in $(GOLDSPACE) directory."

# This set of tests is meant to be run locally, when Jenkins server is not available or one wants to test things before pushing to the repository.
# It is not meant to replace Jenkins, but to provide a way to run all tests locally.
# New Jenkins configuration should use them instead of custom scripts.
CI: QA
	$(MAKE) -k artifact_tests
	$(MAKE) -k allgold
