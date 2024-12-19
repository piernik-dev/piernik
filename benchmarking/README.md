# Purpose
The benchmarking scripts are designed to run Piernik setups on a **single machine** with various numbers of threads and collect timing data. The collected data can help in:

* Identifying scaling issues of specific algorithms.
* Identifying weaknesses in particular architectures, such as poor memory bandwidth or throttling.
* Demonstrating performance and scalability improvements.

# Usage
* `benchmarking/standard_set_of_tests.sh results_directory` runs a standard-sized set of tests 3 times, and two enlarged tests (by 50% and by 100% in each direction).
* `benchmarking/piernik_bench.sh | tee results_file` runs a set of tests once.

By default, the tests start from one thread and increase the parallelism up to the maximum number of hardware threads available (which may exceed the number of cores). This behavior can be overridden by providing the desired number of threads as arguments to the `benchmarking/piernik_bench.sh` script.

The results can be visualized using the `bench_plot.py` script from [Piernik benchmarks](https://github.com/gawrysz/Piernik-benchmarks). This repository also contains a collection of performance characteristics for different setups.

# Interpretation of results

## Scaling test types

* **Weak scaling** – the size of the problem increases proportionally with the number of threads. Ideally, the execution time should remain constant for the number of threads less than or equal to the number of cores. Deviations from flatness in the *weak scaling test* usually come from communication overhead.
* **Strong scaling** – the size of the problem remains constant as the number of threads increases. We expect the execution time to decrease with the number of cores in use, hence for the presentation purpose we multiply the recorded execution time by the number of threads. Deviations from flatness in the *strong scaling test* may come from both communication overhead and difficulties in load balancing after partitioning the domain to too small pieces.
* **Flood scaling** – an uncommon test that runs multiple independent, single-threaded tasks. Flat characteristics are expected here as well. Deviations from flatness may indicate issues with memory bandwidth or other shared resources among threads. One may take the *flood scaling* as the limiting case of *weak scaling* without the communication.
* As a bonus, **compilation time** and **compilation CPU usage** graphs provide additional insights into the computer's performance on tasks that do not heavily utilize the FPU (Floating Point Unit).

On many CPUs, you may observe a bump in the middle of the graph and significantly higher execution times on the right side due to HT (Hyper-Threading) or SMT (Simultaneous Multi-Threading) technologies, which allow running two threads on a single core. When both threads are active, they share the same set of execution units, potentially leading to contention and increased execution time.

The current configuration is not well-suited for CPUs with a mix of performance cores (*p-cores*) and efficiency cores (*e-cores*). Properly utilizing such architectures would require additional effort in process scheduling and enabling features in the `BALANCE` facility in Piernik.

## Testing problems
The problems and their configurations are chosen to test various aspects of the algorithms used in the Piernik code and to expose hardware bottlenecks.

* **`sedov`** – The most basic hydrodynamic problem, run with a base domain of 64<sup>3</sup> and *cartesian decomposition*. It tends to scale well with a small number of threads but exhibits performance degradation peaks when the number of threads is a prime number. This effect comes from relatively big communication overhead combined with poor load balancing (imagine chopping 64<sup>3</sup> domain into 63 slices in _z_-direction).
* **`maclaurin`** – A multigrid self-gravitating setup with very demanding settings for the multipole solver (`lmax = 64`). It uses *block decomposition* with blocks of 32<sup>3</sup> and a domain size of 64<sup>3</sup> for the *weak* and *flood* tests and 128<sup>3</sup> for the *strong* test.
    * It often scales poorly in the *flood test* due to high sensitivity to memory bandwidth, as the multipole solver uses substantial amount of memory. This test best demonstrates the benefits of multiple memory channels and high-bandwidth memory.
    * It scales poorly in the *weak test* because, besides memory bandwidth sensitivity, the domain enlargement forces proportional increases in memory usage by the multipole solver.
    * It struggles with *strong scaling* at high thread counts due to insufficient number of blocks for efficient load balancing.
* **`crtest`** – An anisotropic multigrid diffusion problem run with a base domain of 32<sup>3</sup> by default and *non-cartesian decomposition*. Unlike `maclaurin`, this problem is very CPU-bound and scales well with the number of threads up to the number of physical cores. Departures from flat characteristic here are often due to CPU lowering its operating frequency under load. The *strong scaling* usually is improving on enlarged tests.

### To Do
* Implement an N-body problem with millions of particles to test the scalability of the N-body solver.
* Implement the CRESP problem with tens of CR fields to stress the previously unexplored regime of enormous `cg%u(:,:,:,:)` arrays.

# The past and the future
The `benchmarking` branch was initiated around 2016 (after PR [#188](https://github.com/piernik-dev/piernik/pull/188)). To maintain the comparability of results gathered in the `Piernik-benchmarking` repository, only the critical bugfixes were applied there.

Unfortunately, even if the code does not change, the hardware changes, which demands use of new releases of the Linux kernel and compilers. This has some impact on the measured performance, so we may decide to break the continuity at some point. It seems to be a good idea to update the set of tests and their parameters to be prepared for the future CPUs with hundreds of cores or hybrid CPUs with different populations of cores.