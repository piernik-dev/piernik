# Running the gold tests by hand

You may also want GNU Parallel to be installed in your system for extra performance.
Make sure that current directory has at least 1 GB of free space.
Then call `make gold` and wait for results summary on `stdout` and `stderr`.
In case `make gold` fails due to limited resources you may want to use `make gold-serial`.
Each tested problem creates some more output files which may be inspected or removed by `make gold-clean`.

# The files

* `gold_test_list.sh` – The script called by `make gold`
* `gold_test.sh` – The script called for every gold test by `gold_test_list`.
* `gold_configs/` – The directory with configuration of the problems.
    Files placed here will be automatically used by `* gold CI` tasks on the Jenkins server.
    Each file there should assign values to the following variables:
    * `GOLD_COMMIT` – SHA1 of commit considered to be good. Update carefully, when everything else fails.
    * `PROBLEM_NAME` – Name of the problem in `problems/` directory.
    * `SETUP_PARAMS` – Options used when calling setup for current test.
    * `GOLD_PARAMS` – Options used when calling setup for gold test.
    * `OUTPUT` – Name of the output file to compare in both runs.
    * `NTHR` – Number of processes to be used in the run (optional, default = 1).
