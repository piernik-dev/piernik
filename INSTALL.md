# Quick start with Piernik

For a typical desktop Linux installation you may need to add few packages to be able to compile and run Piernik.

## Ubuntu

For Ubuntu try:

    sudo apt install git make mpich libhdf5-mpich-dev libfftw3-dev pkg-config python pycodestyle python-numpy gfortran python3-requests

If the above complains that `E: Unable to locate package mpich` then do:

    sudo add-apt-repository universe

and try again to install the pacages mentioned above. You may use `libhdf5-openmpi-dev` and `openmpi-bin` if you prefer OpenMPI over MPICH.

## Fedora 32 and newer

For Fedora try:

    sudo dnf install make hdf5-openmpi-devel fftw-devel python environment-modules python3-pycodestyle python2-numpy python3-h5py python3-requests

You may use `hdf5-mpich-devel` if you prefer MPICH over OpenMPI.

* To enable MPI compiler wrappers in Fedora, use:

        module load mpi

    The above line is worth putting in your `.profile` so you won't have to remember it.

* Old approach to the Fortran MPI interface is now strongly disfavoured and Piernik does not compile by default. Use `-fallow-argument-mismatch` to turn MPI-related errors into warnings.
* The wrapper `h5pfc` from OpenMPI requires to set `-I/usr/lib64/gfortran/modules/openmpi/`.
* A compiler option `-frecursive` is suggested to be safer than the default in few places by the MPI wrappers.
* MPICH seems to have random bugs that manifest e.g by crashing at `MPI_Waitall` in `cg_list_bnd::internal_boundaries_MPI_1by1`.
* MPICH uses `-Werror=format-security` compiler flag which is no longer valid for `gfortran`. Use

      sudo sed -i 's/-Werror=format-security //' /usr/lib64/mpich/bin/mpif90

  to fix the warnings

* OpenMPI refuses to run on all threads of a CPU with SMT/HT. Override with `--use-hwthread-cpus` or `--oversubscribe` (only when necessary).
* If the wrapper `h5pfc` complains about "Nonexistent include directory" you may need to edit it (as a sudoer/root) and fix the `includedir` shell variable (seen in Fedora 37 and later).

### older Fedora releases (24 .. 31)

* Install `python-pycodestyle` instead of `python3-pycodestyle` if necessary for `make qa`.
* If you are affected by a bug in dependencies in Fedora packages that results in failure of MPI compilation due to missing `/usr/lib/rpm/redhat/redhat-hardened-cc1`, you can fix it by installing the missing package:

        sudo dnf install redhat-rpm-config

    You can also get rid of these *hardening* and *fortyfying* features (as these are unlikely to be important for a Piernik user and may negatively impact code performance) by stripping them off from the wrappers:

        sudo sed -i 's/-specs=[^ "]*//g' /usr/lib64/mpich/bin/h5pfc /usr/lib64/mpich/bin/mpif90

* If you see warnings saying that `-Werror=format-security` is not valid for Fortran, do:

        sudo sed -i 's/-Werror=format-security//' /usr/lib64/mpich/bin/h5pfc /usr/lib64/mpich/bin/mpif90

    Remember that doing so is a sort of hack, which may take revenge on you in a distant future.

If you have installed OpenMPI libraries, remember to replace `mpich` with `openmpi` in the `sed` calls above.

### other systems

On other systems you need to find your own way (and you may choose to describe it here). You will need:

* git
* decent Fortran 2008 compiler (for such things as polymorphism or modern MPI interface)
* MPI with Fortran support, preferably version that provides mpi_f08.mod. If it is limited to mpi.mod, it must provide interface to all MPI calls
* hdf5 1.8.8 or newer with high-level library, Fortran 2003 interface and MPI support (--enable-shared --enable-fortran --enable-fortran2003 --enable-parallel) – it allows you to use a wrapper ‘h5pfc’ as a compiler
* Python (unfortunately we still have some 2.7-based scripts), including packages such as h5py, numpy
* There are optional tasks that depend on gnuplot, parallel, graphviz and yt

## Piernik

Now make sure that you've forked Piernik project on GitHub and in a directory of your choice do:

    git clone https://github.com/<your_github_username>/piernik

Now you can `cd piernik` and test whether everything compiles well by calling:

    ./setup maclaurin

which should finish with a message:

    'maclaurin' ready in 'runs/maclaurin'.

The `./setup` command without any arguments will give you a list of recognized options. You can put some them in `.setuprc` file if you want to use them always.

Also look at the top of `Makefile` to find some useful tricks.

The file `bin/bash_completion.sh` may also make your life with Piernik a little bit easier by allowing to complete some setup options and problem names.

## Advanced features of Piernik

If you want to run so called `gold tests`, you will need Gnu Parallel (optional) and h5py:

    sudo dnf install parallel python-h5py

or

    apt install parallel python3-h5py

(on Fedora or Ubuntu, respectively).

Then you can execute so called "gold tests" locally by invoking `make gold`
or `make gold-serial` (when `make gold` requires too many resources, e.g. on
a laptop).

# Setting up Intel Fortran compiler

First, you have to obtain the software, as it most likely is not packaged in the regular repositories of your distro. You can obtain the packages [here](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html). Next, install `intel-hpckit` (a huge bunch of packages, IIRC ~15 GB) or at least `intel-oneapi-mpi-devel`, `intel-oneapi-compiler-fortran` `intel-oneapi-compiler-dpcpp-cpp-and-cpp-classic` and their dependencies. Make sure that you have working C++ compiler in your system, such as `gcc-c++` even when you don't need C++ interface for HDF5.

Then, you have to set up paths in your environment:

    source /opt/intel/oneapi/setvars.sh

which should set up the latest versions of everything.

For Piernik you need also to compile HDF5 as it is not bundled in oneAPI repository. You can find the HDF5 sources [here](https://www.hdfgroup.org/downloads/hdf5/source-code/), download them and unpack somewhere. Then configure and compile, eg.:

    ./configure --prefix=${HOME}/intel/HDF5 --enable-fortran --enable-shared --enable-parallel  --with-pic CC=mpiicc FC=mpiifort CXX=mpiicpc CFLAGS="-fPIC -O3 -xHost -ip -fno-alias -align" FFLAGS="-fPIC -O3 -xHost -ip -fno-alias -align" CXXFLAGS="-fPIC -O3 -xHost -ip -fno-alias -align" FFLAGS="-I/opt/intel/oneapi/mpi/latest/include -L/opt/intel/oneapi/mpi/latest/lib"
    make -j
    make install

* Please note that the classic compilers like `icc` will be soon deprecated so some modifications may be needed in order to use new family of compilers (`icx` and others).
* While it may be tempting to install HDF5 system-wide, in `/opt/intel/HDF5`, there are problems arising when you try to do `sudo make install` because some relinking fails as the subshell calls don't know proper paths.

Make sure that your shell knows about the new HDF5 scripts and libraries:

    export PATH=${HOME}/intel/HDF5/bin/:$PATH
    export LD_LIBRARY_PATH=${HOME}/intel/HDF5/lib/:$LD_LIBRARY_PATH

Last but not least, you have to prepare compiler configuration file. The absolute minimum is:

    PROG     = piernik
    F90      = h5pfc
    F90FLAGS = -r8

Of course there are may other useful options for optimization or checking which you may add. I do not recommend to use `-ipo` as I found some weird Piernik crashes when this option was used.

Let's hope that at this point your Intel compilers are able to produce correct code.
