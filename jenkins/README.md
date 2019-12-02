# Running the gold tests by hand

To run the gold tests locally, make sure that `yt` is available because `bin/gdf_distance` critically depends on it.
Make sure that `/tmp/` directory has about 500 MB of free space.
Then call `make gold` and wait for results summary on `stdout` and `stderr`.
Each tested problem creates some more output files which may be inspected or removed by `make gold-clean`.

# The files

* gold_test_list – The script called by `make gold` and by the `user_gold` test on jenkins
* gold_test.sh – The script called for every gold test by `gold_test_list`.
    It creates `/tmp/jenkins_gold/` directory and keeps there some non-critical data like results of gold runs that can be regenerated easily.
    Can be used to selectively run only some tests.
* gold_configs/ – The directory with configuration of the problems.
    Files placed here will be automatically used by `Piernik CI` task on Jenkins server.
    Each file should assign values to the following variables :
    * GOLD_COMMIT – SHA1 of commit considered to be good. Update czrefully, when everything else fails.
    * PROBLEM_NAME – name of the problem in `problems/` directory
    * SETUP_PARAMS – options used when calling setup for current test
    * GOLD_PARAMS – options used when calling setup for gold test
    * OUTPUT – name of the output file to compare in both runs
    * NTHR – (optional) number of processes to be used in the run (default = 1)