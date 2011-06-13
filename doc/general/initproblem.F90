!>
!! \page initproblem Problem description
!!
!! Each problem should provide a initproblem.F90 file with module initproblem and three subroutines:
!!
!! \par problem_pointers
!!
!! This routine iscalled very early, just after init_mpi and is supposed to set up procedure pointers, if necessary.
!! Such approach allows implementation of problem-specific hacks also during the initialization of the Piernik code.
!! In extreme cases it is possible to change behavior of the Piernik code quite radically, without precompiler directives.
!!
!! It assumes that all procedure pointers are initialized to Null().
!!
!! \par read_problem_par
!!
!! This routine reads the PROBLEM_CONTROL namelist by master process and communicates read parameters to other processes.
!! It can also calculate things that might be needed by the problem, except the initial conditions.
!!
!! This routine is called before it is decided whether to read restart file or start from scratch.
!!
!! \par init_prob
!!
!! This routine fills the grid with initial conditions, according to the parameters read in read_problem_par.
!!
!! It is not called when the code starts from the restart file.
!<
