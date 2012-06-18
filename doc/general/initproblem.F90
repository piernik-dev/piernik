!>
!! \page initproblem Problem description
!!
!! Each problem should provide a initproblem.F90 file with module initproblem and three subroutines:
!!
!! \par init_prob
!!
!! This routine fills the grid with initial conditions, according to the parameters read in read_problem_par.
!!
!! It is not called when the code starts from the restart file.
!!
!! \par read_problem_par
!!
!! This routine reads the PROBLEM_CONTROL namelist by master process and communicates read parameters to other processes.
!! It can also calculate things that might be needed by the problem, except the initial conditions.
!! \attention Pointers assignment dependent on parameters from problem.par file should be placed at the end of this routine.
!!
!! This routine is called before it is decided whether to read restart file or start from scratch.
!!
!! \par problem_pointers
!!
!! This routine is called very early, just after init_mpi and is supposed to set up procedure pointers, if necessary.
!! Such approach allows implementation of problem-specific hacks also during the initialization of the Piernik code.
!! In extreme cases it is possible to change behavior of the Piernik code quite radically, without precompiler directives.
!!
!! There are to sets of procedure pointers for users: @n @n
!! imported from user_hooks module: @n @n
!! \b problem_customize_solution - \copy user_hooks::problem_customize_solution @n @n
!! \b problem_grace_passed       - \copy user_hooks::problem_grace_passed @n @n
!! \b user_vars_arr_in_restart   - \copy user_hooks::user_vars_arr_in_restart @n @n
!! \b problem_post_restart       - \copy user_hooks::problem_post_restart @n @n
!! \b finalize_problem           - \copy user_hooks::finalize_problem @n @n
!! \b cleanup_problem            - \copy user_hooks::cleanup_problem @n @n
!! \b problem_refine_derefine    - \copy user_hooks::problem_refine_derefine @n @n
!! \b custom_emf_bnd             - \copy user_hooks::custom_emf_bnd @n @n
!! \b at_user_dims_settings      - \copy user_hooks::at_user_dims_settings @n @n
!! \b at_user_area_settings      - \copy user_hooks::at_user_area_settings @n @n
!! @n @n
!!
!! imported from dataio_user module: @n @n
!! \b user_attrs_pre       - \copy dataio_user::user_attrs_pre @n @n
!! \b user_attrs_rd        - \copy dataio_user::user_attrs_rd @n @n
!! \b user_attrs_wr        - \copy dataio_user::user_attrs_wr @n @n
!! \b user_plt_attrs       - \copy dataio_user::user_plt_attrs @n @n
!! \b user_reg_var_restart - \copy dataio_user::user_reg_var_restart @n @n
!! \b user_plt_hdf5        - \copy dataio_user::user_plt_hdf5 @n @n
!! \b user_vars_hdf5       - \copy dataio_user::user_vars_hdf5 @n @n
!! \b user_tsl             - \copy dataio_user::user_tsl @n @n
!! \b user_post_write_data - \copy dataio_user::user_post_write_data @n @n
!! @n @n
!!
!! It assumes that all procedure pointers are initialized to Null().
!!
!! \warning Pointers assignment may be dependent on parameters from problem.par file (unrecommended). In such a case they must
!! be assigned not before read_problem_par routine is called.
!<
