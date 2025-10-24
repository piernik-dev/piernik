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

# Define variables for colors and paths
ALLOBJ = $(wildcard obj*)

# Terminal colors and formatting
RED    = "\033[31;1m"
GREEN  = "\033[32;1m"
YELLOW = "\033[33;1m"
BLUE   = "\033[34;1m"
PURPLE = "\033[35;1m"
CYAN   = "\033[36;1m"
WHITE  = "\033[37;1m"
GRAY   = "\033[90;1m"
RESET  = "\033[0m"
PASSED = $(GREEN)passed$(RESET)
FAILED = $(RED)failed$(RESET)
OK = $(CYAN)OK$(RESET)
DIFF = $(YELLOW)Differs!$(RESET)
COLORTEST = $(RED) RED $(GREEN) GREEN $(YELLOW) YELLOW $(BLUE) BLUE $(PURPLE) PURPLE $(CYAN) CYAN $(WHITE) WHITE $(GRAY) GRAY $(RESET) back to normal

# Directory structure
PROBLEMS_DIR = problems
RUNS_DIR = runs
SRC_DIR = src
BIN_DIR = bin
PYTHON_DIR = python
JENKINS_DIR = ./jenkins
ARTIFACTS = "$(JENKINS_DIR)/artifacts"
GOLDSPACE = "$(JENKINS_DIR)/workspace"
CONFIG_DIR = $(JENKINS_DIR)/gold_configs

# System commands
ECHO ?= /bin/echo
RM ?= /bin/rm
MPIEXEC ?= mpiexec
DISPLAY ?= display
PYTHON ?= python3
SHELL := /bin/bash

export LC_ALL = C

# Test environment
SETUP = ./setup
GOLD_TEST_SCRIPT = $(JENKINS_DIR)/gold_test.sh


# Sets of tests
QA_TESTS = chk_err_msg chk_lic_hdr pycodestyle qa
ARTIFACT_TESTS = dep py3 noHDF5 I64 IOv2 Jeans Maclaurin 3body CPAW
GOLD_TESTS = gold_CRESP gold_mcrtest gold_mcrwind gold_MHDsedovAMR gold_resist gold_streaming_instability custom_gold

# Common command sequences
define cleanup_tmpdir
	rm -r $${OTMPDIR} $${RUNDIR}
endef

# Define phony targets
.PHONY: $(ALLOBJ) ctags resetup clean check allsetup doxy help \
	oldgold gold-serial gold-clean pep8 view_dep colortest \
	CI QA artifacts gold $(QA_TESTS) $(ARTIFACT_TESTS) $(GOLD_TESTS)

# Default target to build all object directories
all: $(ALLOBJ)

# Target to build a specific object directory
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

# Target to recreate ctags
ctags:
	ctags -R {$(SRC_DIR),$(PROBLEMS_DIR)} --fortran-kinds=+iL

# Target to resetup all object directories
resetup:
	@RS=1 $(MAKE) -k all

# Target to clean all object directories
clean:
	@CL=1 $(MAKE) -k all

# Function to setup a problem (used in allsetup target)
define setup_problem
	pnm=$$( echo $(1) | sed 's-^$(PROBLEMS_DIR)/--' ); \
	onm=$${pnm//\//___}; \
	$(SETUP) $$pnm -o "A_"$$onm --nocompile && \
	for f in .setup.call Makefile env.dat version.F90 ; do \
		of="obj_A_"$$onm"/"$$f; \
		[ -f $$of ] && sed -i 's/ --nocompile//' $$of; \
	done
endef

# Target to setup all problems (compile them e.g. with "make obj_A_* -j")
allsetup:
	( for i in $$( find $(PROBLEMS_DIR)/* -type d ) ; do \
		if [ ! -e $$i/.skipauto ] ; then \
			$(call setup_problem,$$i) & \
			sleep .1; \
		fi; \
	done; \
	wait )

# Target to generate Doxygen documentation
doxy:
	doxygen piernik.doxy

# Target to run qa.py checks
qa:
	$(BIN_DIR)/qa.py -q $$( git ls-files | grep -vE "^(compilers/tests|doc/general)" | grep "\.F90$$" ) && \
	echo -e "  qa.py checks "$(PASSED) || \
	( $(ECHO) -e "  qa.py checks "$(FAILED) && exit 1 )

# Targets to run pycodestyle checks
pep8: pycodestyle

pycodestyle:
	TSTNAME="  Pycodestyle check "; \
	REMARK=" (with --ignore=E501,E722,W504,W605)"; \
	pycodestyle `git ls-files | grep '\.py$$'` $(BIN_DIR)/gdf_distance --ignore=E501,E722,W504,W605 && \
		$(ECHO) -e "$$TSTNAME"$(PASSED)"$$REMARK" ||\
		( $(ECHO) -e "$$TSTNAME"$(FAILED)"$$REMARK" && exit 1 )

# Target to check error messages in the Fortran source files
chk_err_msg:
	$(BIN_DIR)/checkmessages.sh && \
		$(ECHO) -e "  Message checks "$(PASSED) || \
		( $(ECHO) -e "  Message checks "$(FAILED)": incorrect file references found" && exit 1 )

# Target to check consistency of license headers
chk_lic_hdr:
	$(BIN_DIR)/check_license_headers.sh && \
		$(ECHO) -e "  License header checks "$(PASSED) || \
		( $(ECHO) -e "  License header checks "$(FAILED)": exceptional license headers found" && exit 1 )

# Help target
help:
	@$(ECHO) "Testing targets:"
	@$(ECHO) "  QA        - Run all QA checks       ($(QA_TESTS))"
	@$(ECHO) "  artifacts - Run all artifact tests  ($(ARTIFACT_TESTS))"
	@$(ECHO) "  gold      - Run all gold tests      ($(GOLD_TESTS))"
	@$(ECHO) "  CI        - Run all CI checks       (QA artifacts gold)"
	@$(ECHO) "  help      - Show this help"

# Testing ANSI colors in current terminal
colortest:
	$(ECHO) -e "Color test:"$(COLORTEST)

# Target to run all QA checks
QA:
	$(ECHO) -e $(BLUE)"Starting QA checks ..."$(RESET)
	$(MAKE) -k $(QA_TESTS)
	$(ECHO) -e $(BLUE)"QA checks "$(PASSED)"."

# Target to create dependency graph
ifndef P
P = "testing_and_debuging/chimaera"
endif
dep:
	OTMPDIR=$$(mktemp -d obj_XXXXXX) ;\
	PROBLEM=$P ;\
	RUNDIR="runs/"$$( basename $${PROBLEM} )"_"$${OTMPDIR//obj_/} ;\
	GRAPH="dep.png" ;\
	[ ! -e $(ARTIFACTS) ] && mkdir -p $(ARTIFACTS) ;\
	[ -e $(ARTIFACTS)"/$$GRAPH" ] && rm $(ARTIFACTS)"/$$GRAPH" || true ;\
	$(SETUP) $$PROBLEM -n -o $${OTMPDIR//obj_/} > setup.stdout ;\
	if [ -e $$OTMPDIR ] ; then \
		mv setup.stdout $$OTMPDIR ;\
		$(MAKE) -k -C $$OTMPDIR $$GRAPH ;\
		mv $$OTMPDIR"/"$$GRAPH $(ARTIFACTS)/ && $(cleanup_tmpdir) ;\
	fi ;\
	if [ -e $(ARTIFACTS)"/$$GRAPH" ] ; then \
		$(ECHO) -e "  Dependency test "$(PASSED)", the graph for the $$PROBLEM problem was stored as $(ARTIFACTS)/$$GRAPH" ;\
	else \
		[ -e $$OTMPDIR/setup.stdout ] && ( cat $$OTMPDIR/setup.stdout ; rm -r $$OTMPDIR ) ;\
		$(ECHO) -e "  Dependency graph creation "$(FAILED) && exit 1 ;\
	fi

# Target to view dependency graph
view_dep: dep
	which $(DISPLAY) > /dev/null 2> /dev/null && [ -e $(ARTIFACTS)"/dep.png" ] && $(DISPLAY) $(ARTIFACTS)"/dep.png" || true

# Target to check for Python 3 compatibility of the environment
py3:
	OTMPDIR=$$(mktemp -d obj_XXXXXX) ;\
	RUNDIR=$(RUNS_DIR)/maclaurin_$${OTMPDIR//obj_/} ;\
	[ ! -e $(ARTIFACTS) ] && mkdir -p $(ARTIFACTS) ;\
	python3 $(PYTHON_DIR)/piernik_setup.py maclaurin -n -o $${OTMPDIR//obj_/} > $(ARTIFACTS)/py3.setup.stdout && \
		( $(cleanup_tmpdir) ; $(ECHO) -e "  Python 3 test "$(PASSED) ) || \
		( $(cleanup_tmpdir) ; $(ECHO) -e "  Python 3 test "$(FAILED) && exit 1 )

# Target to check compilation without HDF5 library
noHDF5:
	OTMPDIR=$$(mktemp -d obj_XXXXXX) ;\
	RUNDIR=$(RUNS_DIR)/crtest_$${OTMPDIR//obj_/} ;\
	[ ! -e $(ARTIFACTS) ] && mkdir -p $(ARTIFACTS) ;\
	$(SETUP) crtest --param problem.par.build -d I_KNOW_WHAT_I_AM_DOING -o $${OTMPDIR//obj_/} > $(ARTIFACTS)/noHDF5.setup.stdout && \
		( $(cleanup_tmpdir) ; $(ECHO) -e "  NoHDF5 test "$(PASSED) ) || \
		( $(cleanup_tmpdir) ; $(ECHO) -e "  NoHDF5 test "$(FAILED) && exit 1 )

# Target to test compilation with 64-bit integers
I64:
	OTMPDIR=$$(mktemp -d obj_XXXXXX) ;\
	RUNDIR=$(RUNS_DIR)/chimaera_$${OTMPDIR//obj_/} ;\
	[ ! -e $(ARTIFACTS) ] && mkdir -p $(ARTIFACTS) ;\
	$(SETUP) $(P) -o $${OTMPDIR//obj_/} --f90flags="-fdefault-integer-8 -Werror=conversion" > $(ARTIFACTS)/I64.setup.stdout && \
		( $(cleanup_tmpdir) ; $(ECHO) -e "  64-bit integer test "$(PASSED) ) || \
		( $(cleanup_tmpdir) ; $(ECHO) -e "  64-bit integer test "$(FAILED) && exit 1 )

# Target to run IO version 2 restart compatision test
IOv2:
	OTMPDIR=$$(mktemp -d obj_XXXXXX) ;\
	RUNDIR=$(RUNS_DIR)/advection_test_$${OTMPDIR//obj_/} ;\
	[ ! -e $(ARTIFACTS) ] && mkdir -p $(ARTIFACTS) ;\
	$(SETUP) advection_test -p problem.par.restart_test_v2 -o $${OTMPDIR//obj_/} > $(ARTIFACTS)/IOv2.setup.stdout && \
		( cd $${RUNDIR} ;\
			$(ECHO) -e "run_id = ts1\n  Start:   t = 0.0, nproc = 1" ;\
			$(MPIEXEC) -n 1 ./piernik > ts1.out ;\
			$(ECHO) -e "  Restart: t = 1.0, nproc = 1" ;\
			$(MPIEXEC) -n 1 ./piernik -n '&END_CONTROL tend = 2.0 /' >> ts1.out ;\
			$(ECHO) -e "  Finish:  t = 2.0\n\nrun_id = ts2\n  Start:   t = 0.0, nproc = 5" ;\
			$(MPIEXEC) -n 5 ./piernik -n '&OUTPUT_CONTROL run_id = "ts2" /' > ts2.out ;\
			$(ECHO) -e "  Restart: t = 1.0, nproc = 3" ;\
			$(MPIEXEC) -n 3 ./piernik -n '&END_CONTROL tend = 2.0 /' -n '&OUTPUT_CONTROL run_id = "ts2" /' >> ts2.out ;\
			$(ECHO) -e "  Finish:  t = 2.0" ;\
			../../$(BIN_DIR)/gdf_distance moving_pulse_ts{1,2}_0002.h5 2>&1 | tee compare.log ;\
		) >> $(ARTIFACTS)/IOv2.setup.stdout && \
		( [ $$( grep "^Total difference between" $${RUNDIR}/compare.log | awk '{print $$NF}' ) == 0 ] || exit 1 ) && \
		( $(cleanup_tmpdir) ; $(ECHO) -e "  IO v2 test "$(PASSED) ) || \
		( $(cleanup_tmpdir) ; $(ECHO) -e "  IO v2 test "$(FAILED) && exit 1 )

# Target to run Jeans instability test with accuracy check and restart reproducibility (restart done in a bit different way than it is done in IOv2)
Jeans:
	OTMPDIR=$$(mktemp -d obj_XXXXXX) ;\
	RUNDIR=$(RUNS_DIR)/jeans_$${OTMPDIR//obj_/} ;\
	[ ! -e $(ARTIFACTS) ] && mkdir -p $(ARTIFACTS) ;\
	$(SETUP) jeans -o $${OTMPDIR//obj_/} -o $${OTMPDIR//obj_/} > $(ARTIFACTS)/Jeans.setup.stdout && \
		( cd $${RUNDIR} ;\
			$(MPIEXEC) -n 1 ./piernik ;\
			gnuplot jeans.gnuplot 2>&1 || exit 1 ;\
			$(RM) *.res ;\
			$(MPIEXEC) -n 1 ./piernik -n '&END_CONTROL nend = 10/ &OUTPUT_CONTROL dt_hdf = 1.0 dt_res = 1.0 run_id = "rs1"/' ;\
			$(MPIEXEC) -n 1 ./piernik -n '&END_CONTROL nend = 20/ &OUTPUT_CONTROL dt_hdf = 1.0 dt_res = 1.0 run_id = "rs1"/' ;\
			$(MPIEXEC) -n 1 ./piernik -n '&END_CONTROL nend = 20/ &OUTPUT_CONTROL dt_hdf = 1.0 dt_res = 1.0 run_id = "rs2"/' ;\
			../../$(BIN_DIR)/gdf_distance jeans_rs{1,2}_0001.h5 2>&1 | tee compare.log ;\
		) >> $(ARTIFACTS)/Jeans.setup.stdout && \
		( [ $$( grep "^Total difference between" $${RUNDIR}/compare.log | awk '{print $$NF}' ) == 0 ] || exit 1 ) && \
		NORM=$$( awk '{exp_a = -0.0002475700224982; exp_p = 0.00574478339226114; if (NR==2) { printf("Amplitude error = %.16f, period error = %.16f", 1.*$$1, 1.*$$3); if ($$1 != exp_a || $$3 != exp_p) printf(" '$(DIFF)' (expected: %.16f, %.16f)", exp_a, exp_p); else printf(" '$(OK)'");} }' $${RUNDIR}/jeans.csv ) ;\
		cp $${RUNDIR}/jeans.{png,csv} $(ARTIFACTS) &&\
		( $(cleanup_tmpdir) ; $(ECHO) -e "  Jeans test "$(PASSED)", $${NORM}, ${ARTIFACTS}/jeans.png created" ) || \
		( $(cleanup_tmpdir) ; $(ECHO) -e "  Jeans test "$(FAILED) && exit 1 )

# Target to run Maclaurin test with accuracy check
Maclaurin:
	OTMPDIR=$$(mktemp -d obj_XXXXXX) ;\
	RUNDIR=$(RUNS_DIR)/maclaurin_$${OTMPDIR//obj_/} ;\
	[ ! -e $(ARTIFACTS) ] && mkdir -p $(ARTIFACTS) ;\
	$(SETUP) maclaurin -o $${OTMPDIR//obj_/} > $(ARTIFACTS)/Maclaurin.setup.stdout && \
		( cd $${RUNDIR} ;\
			$(MPIEXEC) -n 1 ./piernik ;\
			$(PYTHON) ../../$(PYTHON_DIR)/maclaurin.py ./maclaurin_sph_0000.h5 || exit 1 ;\
			$(ECHO) "L2 error norm,min. error,max. error" > maclaurin.csv ;\
			grep -a L2 *log | tail -n 1 | sed 's/.*= *\(.*\),.*= *\([^ ]*\) *\(.*\)$$/\1 , \2 , \3/' >> maclaurin.csv ;\
		) >> $(ARTIFACTS)/Maclaurin.setup.stdout && \
		NORM=$$( awk '{exp_n = 0.003359; if (NR==2) { printf("L2 error norm = %.6f", $$1); if ($$1 != exp_n) printf(" '$(DIFF)' (expected: %.6f)", exp_n); else printf(" '$(OK)'");} }' $${RUNDIR}/maclaurin.csv ) ;\
		cp $${RUNDIR}/maclaurin.{png,csv} $(ARTIFACTS) &&\
		( $(cleanup_tmpdir) ; $(ECHO) -e "  Maclaurin test "$(PASSED)", $${NORM}, ${ARTIFACTS}/maclaurin.png created." ) || \
		( $(cleanup_tmpdir) ; $(ECHO) -e "  Maclaurin test "$(FAILED) && exit 1 )

# Target to run 3-body test with accuracy and restart checks (here restart also checks particles)
3body:
	OTMPDIR=$$(mktemp -d obj_XXXXXX) ;\
	RUNDIR=$(RUNS_DIR)/3body_$${OTMPDIR//obj_/} ;\
	[ ! -e $(ARTIFACTS) ] && mkdir -p $(ARTIFACTS) ;\
	$(SETUP) 2body/3body -o $${OTMPDIR//obj_/} > $(ARTIFACTS)/3body.setup.stdout && \
		( cd $${RUNDIR} ;\
			( $(MPIEXEC) -n 3 ./piernik ;\
			../../$(PROBLEMS_DIR)/2body/3body/particle_error.py leapfrog_tst_0001.h5 | tee 3body.csv ) & pid_1=$$! ;\
			# Run the two restart groups in parallel (the rs1 group contains two sequential mpiexec calls) ;\
			RESTSTDIR=res ;\
			( $(ECHO) "Performing restart tests" ;\
				mkdir $${RESTSTDIR};\
				cd $${RESTSTDIR} ;\
				ln -s ../piernik ;\
				cp ../problem.par . ;\
				sed -i 's-problems-../problems-' problem.par ;\
				$(MPIEXEC) -n 3 ./piernik -n '&END_CONTROL nend = 30/' ; \
				$(MPIEXEC) -n 3 ./piernik ) >> stdout & pid_rs=$$! ;\
			# Wait for both jobs to finish before comparing ;\
			wait $$pid_1 $$pid_rs ;\
			cat stdout ;\
			../../$(BIN_DIR)/gdf_distance {,$${RESTSTDIR}/}leapfrog_tst_0001.h5 2>&1 | tee compare.log ;\
		) >> $(ARTIFACTS)/3body.setup.stdout && \
		[ $$( grep "^Total difference between" $${RUNDIR}/compare.log | awk '{print $$NF}' ) == 0 ] && \
		( NORM=$$( awk '{exp_p = 0.004358; exp_m = 3.05544e-08; exp_am = 0.00385792; if (NR==2) { printf("Period error = %.6f, momentum error = %.6g, angular momentum error = %.8f", 1.*$$1, 1.*$$3, 1*$$5); if ($$1 != exp_p || $$3 != exp_m || $$5 != exp_am) printf(" '$(DIFF)' (expected: %.6f, %.6g, %.8f)", exp_p, exp_m, exp_am); else printf(" '$(OK)'");} }' $${RUNDIR}/3body.csv ) ;\
			cp $${RUNDIR}/3body.csv $(ARTIFACTS) ;\
			$(cleanup_tmpdir) ;\
			$(ECHO) -e "  3-body test "$(PASSED)", $${NORM}" ) || \
		( $(cleanup_tmpdir) ; $(ECHO) -e "  3-body test "$(FAILED) && exit 1 )

# Target to run Propagation of Circularly polarized Alfvén Waves
CPAW:
	OTMPDIR=$$(mktemp -d obj_XXXXXX) ;\
	RUNDIR=$(RUNS_DIR)/cpaw_$${OTMPDIR//obj_/} ;\
	[ ! -e $(ARTIFACTS) ] && mkdir -p $(ARTIFACTS) ;\
	$(SETUP) cpaw -o $${OTMPDIR//obj_/} > $(ARTIFACTS)/CPAW.setup.stdout && \
		( cd $${RUNDIR} ;\
			$(MPIEXEC) -n 1 ./piernik -n '&BASE_DOMAIN n_d = 128, 64, 1 /';\
			$(ECHO) "L2 error norm" > CPAW.csv ;\
			grep -a L2 *log | tail -n 1 | sed 's/.*= *\(.*\).*$$/\1/' >> CPAW.csv ;\
		) >> $(ARTIFACTS)/CPAW.setup.stdout && \
		NORM=$$( awk '{exp_n = 0.00102726; if (NR==2) { printf("L2 error norm = %.8f", $$1); if ($$1 != exp_n) printf(" '$(DIFF)' (expected: %.8f)", exp_n); else printf(" '$(OK)'");} }' $${RUNDIR}/CPAW.csv ) ;\
		cp $${RUNDIR}/CPAW.csv $(ARTIFACTS) &&\
		( $(cleanup_tmpdir) ; $(ECHO) -e "  CPAW test "$(PASSED)", $${NORM}" ) || \
		( $(cleanup_tmpdir) ; $(ECHO) -e "  CPAW test "$(FAILED) && exit 1 )


# Target to run all CI artifact tests
artifacts:
	$(ECHO) -e $(BLUE)"Starting artifact tests ..."$(RESET)
	[ -e $(ARTIFACTS) ] && rm -rf $(ARTIFACTS) || true
	[ ! -e $(ARTIFACTS) ] && mkdir -p $(ARTIFACTS)
	$(MAKE) -k $(ARTIFACT_TESTS) || \
		( $(ECHO) -e $(RED)"Some artifact tests failed."$(RESET)" Details can be found in "$(ARTIFACTS)" directory." && exit 1 )
	$(ECHO) -e $(BLUE)"All artifact tests "$(PASSED)". Details can be found in "$(ARTIFACTS)" directory."

# Target to run all gold tests (old version)
oldgold:
	./jenkins/gold_test_list.sh

# Target to run gold tests in serial mode
gold-serial:
	SERIAL=1 ./jenkins/gold_test_list.sh

# Target to clean gold test files and CI artifacts
gold-clean:
	$(RM) -rf $(GOLDSPACE)/* $(ARTIFACTS)

# Function to run a single gold test
define run_gold_test
	$(GOLD_TEST_SCRIPT) $(CONFIG_DIR)/$(1).config > $(GOLDSPACE)/$(1).gold_stdout 2> $(GOLDSPACE)/$(1).gold_stderr && \
		( grep Warning $(GOLDSPACE)/$(1).gold_stderr; $(ECHO) -e "  $(1) test "$(PASSED) ) || \
		( $(ECHO) -e "  $(1) test "$(FAILED)" (more details in $(GOLDSPACE)/$(1).gold_std*)" && exit 1 )
endef

# Targets to run specific gold tests
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

# Target to run custom gold test, e.g. before making it official. Multiple tests will be executed serially here, so don't stack too many of them.
custom_gold:
	for i in jenkins/gold_configs/*.config ; do \
		T=$$( basename $$i ) ; \
		case $$T in \
			(mcrtest_CRESP.config|mcrtest_new.config|mcrwind.config|MHDsedovAMR.config|resist.config|streaming_instability.config) ;; \
			(*) $(call run_gold_test,$${T//.config/}) ;; \
		esac \
	done

# Target to run all gold tests.
# Call `make gold-clean` explicitly if a cleanup is required.
# Typically it shouldn't be required and existing files from previous gold runs can save some execution time here.
gold:
	$(ECHO) -e $(BLUE)"Starting gold tests ..."$(RESET)
	mkdir -p $(GOLDSPACE) # Create the directory if it doesn't exist
	$(MAKE) -k $(GOLD_TESTS) || \
		( $(ECHO) -e $(RED)"Some gold tests failed."$(RESET)" Details can be found in $(GOLDSPACE) directory." && exit 1 )
	$(ECHO) -e $(BLUE)"All gold tests "$(PASSED)". Details can be found in $(GOLDSPACE) directory."

# This set of tests is meant to be run locally, when Jenkins server is not available or one wants to test things before pushing to the repository.
# It is not meant to replace Jenkins, but to provide a way to run all tests locally.
# New Jenkins configuration should use them instead of custom scripts.
CI: QA
	$(MAKE) -k artifacts
	$(MAKE) -k gold
