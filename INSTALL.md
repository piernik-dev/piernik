# Quick start with Piernik

For a typical desktop Linux installation you may need to add few packages to be able to compile and run Piernik.

## Ubuntu

For Ubuntu try:

    sudo apt install git make mpich libhdf5-mpich-dev libfftw3-dev pkg-config python pycodestyle python-numpy gfortran

If the above complains that `E: Unable to locate package mpich` then do:

    sudo add-apt-repository universe

and try again to install the pacages mentioned above. You may use `libhdf5-openmpi-dev` and `openmpi-bin` if you prefer OpenMPI over MPICH.

## Fedora 32 and newer

For Fedora try:

    sudo dnf install make hdf5-openmpi-devel fftw-devel python environment-modules python3-pycodestyle python2-numpy python3-h5py

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

### older Fedora releases

* Install `python-pycodestyle` instead of `python3-pycodestyle` if necessary for `make qa`.
* If you are affected by a bug in dependencies in Fedora packages that results in failure of MPI compilation due to missing `/usr/lib/rpm/redhat/redhat-hardened-cc1`, you can fix it by installing the missing package:

        sudo dnf install redhat-rpm-config

    You can also get rid of these *hardening* and *fortyfying* features (as these are unlikely to be important for a Piernik user and may negatively impact code performance) by stripping them off from the wrappers:

        sudo sed -i 's/-specs=[^ "]*//g' /usr/lib64/mpich/bin/h5pfc /usr/lib64/mpich/bin/mpif90

* If you see warnings saying that `-Werror=format-security` is not valid for Fortran, do:

        sudo sed -i 's/-Werror=format-security//' /usr/lib64/mpich/bin/h5pfc /usr/lib64/mpich/bin/mpif90

    Remember that doing so is a sort of hack, which may take revenge on you in a distant future.

If you have installed OpenMPI libraries, remember to replace `mpich` with `openmpi` in the `sed` calls above.

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

    apt install parallel python-h5py

(on Fedora or Ubuntu, respectively).

Then you can execute so called "gold tests" locally by invoking `make gold`
or `make gold-serial` (when `make gold` requires too many resources, e.g. on
a laptop).
