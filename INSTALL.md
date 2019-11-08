# Quick start with Piernik

For typical desktop Linux installation you may need to add few packages to be able to compile and run Piernik.

## Ubuntu

For Ubuntu try:

    sudo apt install git make mpich libhdf5-mpich-dev libfftw3-dev\
    pkg-config python pep8 python-numpy gfortran

If the above complains that `E: Unable to locate package mpich` then do:

    sudo add-apt-repository universe

and try again to install the pacages mentioned above. You may use `libhdf5-openmpi-dev` and `openmpi-bin` if you prefer OpenMPI over MPICH.

## Fedora

For Fedora try:

    sudo dnf install make hdf5-mpich-devel fftw-devel python\
    environment-modules python-pep8 python2-numpy

You may use `hdf5-openmpi-devel` if you prefer OpenMPI over MPICH.

To enable MPI compiler wrappers in Fedora, use:

    module load mpi

The above line is worth putting in your `.profile` so you won't have to remember it.

Unfortunately there is a bug in dependencies in Fedora packages. Typical MPI compilation by default fails due to missing `/usr/lib/rpm/redhat/redhat-hardened-cc1`. You can fix it by installing the missing package:

    sudo dnf install redhat-rpm-config

You can also get rid of these *hardening* and *fortyfying* features (as these are unlikely to be important for a Piernik user and may negatively impact code performance) by stripping them off from the wrappers:

    sudo sed -i 's/-specs=[^ "]*//g' \
    /usr/lib64/mpich/bin/h5pfc /usr/lib64/mpich/bin/mpif90

If you see warnings saying that `-Werror=format-security` is not valid for Fortran, do:

    sudo sed -i 's/-Werror=format-security//' \
    /usr/lib64/mpich/bin/h5pfc /usr/lib64/mpich/bin/mpif90

Remember that doing so is a sort of hack, which may take revenge on you in a distant future.

If you installed OpenMPI libraries, remember to replace `mpich` with `openmpi` in the `sed` calls above.

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
