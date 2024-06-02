# Purpose
The benchmarking scripts are intended to run some Piernik setups on a **single machine** with various number of threads and collect some timings. The collected data may allow for:

* Identifying scaling problems of particular algorithms.
* Identifying weaknesses of particular architecture, such as poor memory bandwidth or throttling.
* Presenting performance and scalability improvements.

# Usage
* `benchmarking/standard_set_of_tests.sh results_directory` runs standard-sized set of tests 3 times and two enlarged tests (by 50% and by 100% in each direction).
* `benchmarking/piernik_bench.sh | tee results_file` runs a set of tests once.

By default, the tests are run starting from one thread and increasing the parallelism up to maximum number of hardware threads available (which may be more than the number of cores). One can override this by providing the selected number of threads to be tested as and argument to these scripts.

The results may be visualized by the `bench_plot.py` script from [Piernik benchmarks](https://github.com/gawrysz/Piernik-benchmarks). You can find there also a collection of performance characteristics took on various desktop, server and even some mobile platforms.

# Interpretation of results

## Scaling test types

* **Weak scaling** – the size of the problem is increased proportionally to the number of threads. Here we expect (*and desire*) flat characteristics for number of threads less or equal to number of cores. Deviations from flatness in the *weak scaling test* usually come from communication overhead.
* **Strong scaling** – the size of the problem remains constant as the number of threads increases. Here we expect the execution time to decrease with the number of cores in use, hence for the presentation purpose we multiply the recorded execution time by the number of threads in use. After this operation again we expect flat characteristics. Deviations from flatness in the *strong scaling test* may come from both communication overhead and difficulties in load balancing after partitioning the domain to too small pieces.
* **Flood scaling** – this is quite uncommon test, which runs a number of independent, single-threaded problems. Also here we expect flat characteristics. Deviations from flatness in the *flood scaling test* usually come from limited memory bandwidth or cache concurrency. One may take the *flood scaling* as the limiting case of *weak scaling* without the communication.
* As a bonus there are **compilation time** and **compilation CPU usage** graphs that may give some more hints about performance of the computer on a task that does not use much of FPU.

On many CPUs one can see a bump in the half of the graph and significantly higher execution times on the right side. This is due to HT or SMT technologies, which allow to run two threads on a single core. It can give some benefit as long as there are underutilized internal execution units in the CPU cores but in the case of *floating point*-dominated code, the typical outcome are the cores slowed down nearly by a half.

The current configuration is not well suited for CPUs which consist of *p-cores* and *e-cores*. It would require some extra effort in scheduling the processes and smart use of the `BALANCE` facility in Piernik to reach and show the peak performance available on such hybrid CPUs.

## Testing problems
The problems and their configuration was chosen to test various aspects of the algorithms used in the Piernik code and to expose the hardware bottlenecks.

* **`sedov`** – The most basic hydrodynamic problem run with the base domain of 64<sup>3</sup> and *cartesian decomposition*. It tends to scale well for small number of threads and exhibits peaks of poor performance when the number of threads is a prime number. This is because such a decomposition for a prime number results in creating slices in only one direction and relatively big communication overhead is combined with poor load balancing (imagine chopping 64<sup>3</sup> domain into 63 slices in _z_-direction ...).
* **`maclaurin`** – The multigrid self-gravitating setup with very demanding settings for the multipole solver (`lmax = 64`). It is using *block decomposition*. The blocks are 32<sup>3</sup> and the domain is by default 64<sup>3</sup> for the *weak* and *flood* tests and 128<sup>3</sup> for the *strong* test.
    * It is often scaling poorly in *flood test* because it is very sensitive to the memory bandwidth due to high memory usage in the multipole solver. This test shows best the benefits of having multiple memory channels in use.
    * It scales poorly in the *weak test* because, apart from the sensitivity to the memory bandwidth, the way that the domain is enlarged, forces increasing memory usage in the multipole solver proportionally to the number of threads involved.
    * It has problems with *strong scaling* for a high thread count because there are not enough blocks for efficient load balancing.
* **`crtest`** – The anisotropic multigrid diffusion problem run with the base domain of 32<sup>3</sup> by default and *noncartesian decomposition*. In contrast to `maclaurin`, this is a very CPU-bound problem, nearly independent of memory characteristics. Departures from flat characteristic here are often due to CPU lowering its operating frequency under load. The *strong* scaling is not impressive on the default, relatively small domain. It is mainly due to communication overhead, and it is improving on enlarged tests.

### To Do
* N-body problem with millions of particles for testing scalability of the N-body solver.
* CRESP problem with tens of CR fields stressing the previously unexplored regime of enormous `cg%u(:,:,:,:)` arrays.

# The past and the future
The `benchmarking` branch was started somewhere in 2016 (after PR [#188](https://github.com/piernik-dev/piernik/pull/188)). In order to maintain comparability of the results gathered in the `Piernik-benchmarks` directory, only the critical bugfixes were applied there. Unfortunately, even if the code does not change, the hardware changes, which demands use of new releases of the Linux kernels and operating systems. This has some impact on the measured performance, so we may decide to break the continuity at some point. It will be a good opportunity to update the set of problems and their parameters to be prepared for the future CPUs with hundreds of cores or hybrid CPUs with different populations of cores.
