! $Id$
!
! PIERNIK Code Copyright (C) 2006 Michal Hanasz
!
!    This file is part of PIERNIK code.
!
!    PIERNIK is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    PIERNIK is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with PIERNIK.  If not, see <http://www.gnu.org/licenses/>.
!
!    Initial implemetation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.def"
!>
!! \brief (KK)
!!
!! In this module following namelists of parameters are specified:
!! \copydetails mpisetup::mpistart
!<
module mpisetup

   implicit none

   include 'mpif.h'

   integer :: nproc, proc, ierr , rc, info
   integer :: status(MPI_STATUS_SIZE,4)
   integer, dimension(4) :: req, err

   real, parameter       :: dt_default_grow = 2.

   real                  :: t, dt, dt_old, dtm
   integer               :: nstep

   integer, parameter    :: ndims = 3       ! 3D grid
   integer               :: comm, comm3d
   integer, dimension(3) :: psize, pcoords, coords
   logical               :: periods(3), reorder
   integer ::   procxl, procxr, procyl, procyr, proczl, proczr, procxyl, procyxl, procxyr, procyxr

   integer, parameter               :: buffer_dim=200
   integer, parameter               :: cbuff_len=32
   character(len=cbuff_len), dimension(buffer_dim) :: cbuff
   integer,   dimension(buffer_dim) :: ibuff
   real,      dimension(buffer_dim) :: rbuff
   logical,   dimension(buffer_dim) :: lbuff

   logical :: mpi_magic
   integer :: pxsize    !< number of MPI blocks in x-dimension
   integer :: pysize    !< number of MPI blocks in y-dimension
   integer :: pzsize    !< number of MPI blocks in z-dimension

   namelist /MPI_BLOCKS/ pxsize, pysize, pzsize, mpi_magic

   character(len=4) :: bnd_xl    !< type of boundary conditions for the left x-boundary
   character(len=4) :: bnd_xr    !< type of boundary conditions for right the x-boundary
   character(len=4) :: bnd_yl    !< type of boundary conditions for the left y-boundary
   character(len=4) :: bnd_yr    !< type of boundary conditions for the right y-boundary
   character(len=4) :: bnd_zl    !< type of boundary conditions for the left z-boundary
   character(len=4) :: bnd_zr    !< type of boundary conditions for the right z-boundary
   character(len=4) :: bnd_xl_dom, bnd_xr_dom, bnd_yl_dom, bnd_yr_dom, bnd_zl_dom, bnd_zr_dom !< computational domain boundaries

   namelist /BOUNDARIES/ bnd_xl, bnd_xr, bnd_yl, bnd_yr, bnd_zl, bnd_zr

   real    :: dt_initial
   real    :: dt_max_grow
   real    :: dt_min
   real    :: cfl
   real    :: smalld
   real    :: smallc
   real    :: smallei
   real    :: cfr_smooth
   integer :: integration_order

   namelist /NUMERICAL_SETUP/  cfl, smalld, smallei, integration_order, cfr_smooth, dt_initial, dt_max_grow, dt_min, smallc

   integer, dimension(3) :: domsize   !< local copy of nxd, nyd, nzd which can be used before init_grid()

   logical     :: mpi
   character(len=80)   :: cwd

   integer :: MPI_XZ_LEFT_BND=-1, MPI_XZ_RIGHT_BND=-1
   integer :: MPI_XZ_LEFT_DOM=-1, MPI_XZ_RIGHT_DOM=-1
   integer :: MPI_XY_LEFT_BND=-1, MPI_XY_RIGHT_BND=-1
   integer :: MPI_XY_LEFT_DOM=-1, MPI_XY_RIGHT_DOM=-1
   integer :: MPI_YZ_LEFT_BND=-1, MPI_YZ_RIGHT_BND=-1
   integer :: MPI_YZ_LEFT_DOM=-1, MPI_YZ_RIGHT_DOM=-1

   integer :: MAG_XZ_LEFT_BND=-1, MAG_XZ_RIGHT_BND=-1
   integer :: MAG_XZ_LEFT_DOM=-1, MAG_XZ_RIGHT_DOM=-1
   integer :: MAG_XY_LEFT_BND=-1, MAG_XY_RIGHT_BND=-1
   integer :: MAG_XY_LEFT_DOM=-1, MAG_XY_RIGHT_DOM=-1
   integer :: MAG_YZ_LEFT_BND=-1, MAG_YZ_RIGHT_BND=-1
   integer :: MAG_YZ_LEFT_DOM=-1, MAG_YZ_RIGHT_DOM=-1

   integer :: ARR_XZ_LEFT_BND=-1, ARR_XZ_RIGHT_BND=-1
   integer :: ARR_XZ_LEFT_DOM=-1, ARR_XZ_RIGHT_DOM=-1
   integer :: ARR_XY_LEFT_BND=-1, ARR_XY_RIGHT_BND=-1
   integer :: ARR_XY_LEFT_DOM=-1, ARR_XY_RIGHT_DOM=-1
   integer :: ARR_YZ_LEFT_BND=-1, ARR_YZ_RIGHT_BND=-1
   integer :: ARR_YZ_LEFT_DOM=-1, ARR_YZ_RIGHT_DOM=-1

   contains

!-----------------------------------------------------------------------------
!>
!! \brief Routine to start MPI
!!
!! \n \n
!! @b MPI_BLOCKS
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>pxsize</td><td>1</td><td>integer</td><td>\copydoc mpisetup::pxsize</td></tr>
!! <tr><td>pysize</td><td>1</td><td>integer</td><td>\copydoc mpisetup::pysize</td></tr>
!! <tr><td>pzsize</td><td>1</td><td>integer</td><td>\copydoc mpisetup::pzsize</td></tr>
!! </table>
!! \n \n
!! @b BOUNDARIES
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>bnd_xl</td><td>'per'</td><td>'per', 'ref', 'out', 'outd', 'outh', 'cor'</td><td>\copydoc mpisetup::bnd_xl</td></tr>
!! <tr><td>bnd_xr</td><td>'per'</td><td>'per', 'ref', 'out', 'outd', 'outh'       </td><td>\copydoc mpisetup::bnd_xr</td></tr>
!! <tr><td>bnd_yl</td><td>'per'</td><td>'per', 'ref', 'out', 'outd', 'outh', 'cor'</td><td>\copydoc mpisetup::bnd_yl</td></tr>
!! <tr><td>bnd_yr</td><td>'per'</td><td>'per', 'ref', 'out', 'outd', 'outh'       </td><td>\copydoc mpisetup::bnd_yr</td></tr>
!! <tr><td>bnd_zl</td><td>'per'</td><td>'per', 'ref', 'out', 'outd', 'outh'       </td><td>\copydoc mpisetup::bnd_zl</td></tr>
!! <tr><td>bnd_zr</td><td>'per'</td><td>'per', 'ref', 'out', 'outd', 'outh'       </td><td>\copydoc mpisetup::bnd_zr</td></tr>
!! </table>
!! \n \n
!! @b NUMERICAL_SETUP
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>cfl              </td><td>0.7   </td><td>real value between 0.0 and 1.0       </td><td>\copydoc mpisetup::cfl              </td></tr>
!! <tr><td>smalld           </td><td>1.e-10</td><td>real value                           </td><td>\copydoc mpisetup::smalld           </td></tr>
!! <tr><td>smallei          </td><td>1.e-10</td><td>real value                           </td><td>\copydoc mpisetup::smallei          </td></tr>
!! <tr><td>integration_order</td><td>2     </td><td>1 or 2 (or 3 - currently unavailable)</td><td>\copydoc mpisetup::integration_order</td></tr>
!! <tr><td>cfr_smooth       </td><td>0.0   </td><td>real value                           </td><td>\copydoc mpisetup::cfr_smooth       </td></tr>
!! <tr><td>dt_initial       </td><td>-1.   </td><td>positive real value or -1.           </td><td>\copydoc mpisetup::dt_initial       </td></tr>
!! <tr><td>dt_max_grow      </td><td>2.    </td><td>real value > 1.1                     </td><td>\copydoc mpisetup::dt_max_grow      </td></tr>
!! <tr><td>dt_min           </td><td>0.    </td><td>positive real value                  </td><td>\copydoc mpisetup::dt_min           </td></tr>
!! </table>
!! \n \n
!<
      subroutine mpistart

         use errh, only : die, namelist_errh

         implicit none

         namelist /DOMAIN_SIZES/ nxd, nyd, nzd, nb

         integer :: iproc, ierrh
         integer :: nxd, nyd, nzd, nb

         integer, parameter    :: cwdlen = 512 ! allow for moderately long CWD
         integer, parameter    :: hnlen = 32   ! hostname length limit
         character(LEN=cwdlen) :: cwd_proc
         character(LEN=hnlen)  :: host_proc
         integer               :: pid_proc

         character(LEN=cwdlen), allocatable, dimension(:) :: cwd_all
         character(LEN=hnlen) , allocatable, dimension(:) :: host_all
         integer              , allocatable, dimension(:) :: pid_all

         integer(kind=1)       :: getcwd, hostnm
         integer(kind=4)       :: getpid
         character(LEN=cwdlen) :: par_file, tmp_log_file
         integer :: cwd_status
         logical :: par_file_exist

         call MPI_INIT( ierr )
         call MPI_COMM_RANK(MPI_COMM_WORLD, proc, ierr)
         comm = MPI_COMM_WORLD
         info = MPI_INFO_NULL
         call MPI_COMM_SIZE(comm, nproc, ierr)

         if (allocated(cwd_all) .or. allocated(host_all) .or. allocated(pid_all)) call die("[mpisetup:mpistart] cwd_all, host_all or pid_all already allocated")

         ! BEWARE if proc /= 0 it is probably enough to allocate only one element or none at all (may depend on MPI implementation)
         allocate(cwd_all(0:nproc))
         allocate(host_all(0:nproc))
         allocate(pid_all(0:nproc))

         pid_proc   = getpid()
         status     = hostnm(host_proc)
         cwd_status = getcwd(cwd_proc)

         if (cwd_status /= 0) call die("[mpisetup:mpistart] problems accessing current working directory.")
#ifdef DEBUG
         write(*,'(3a,i6,3a)') 'mpisetup: host="',trim(host_proc),'", PID=',pid_proc,' CWD="',trim(cwd_proc),'"'
#endif /* DEBUG */

         if(proc == 0) then
            par_file = trim(cwd)//'/problem.par'
            inquire(file=par_file, exist=par_file_exist)
            if(.not. par_file_exist) call die('[mpisetup:mpistart] Cannot find "problem.par" in the working directory')
            tmp_log_file = trim(cwd)//'/tmp.log'
         endif

         call MPI_Gather(cwd_proc,  cwdlen, MPI_CHARACTER, cwd_all,  cwdlen, MPI_CHARACTER, 0, comm, err)
         call MPI_Gather(host_proc, hnlen,  MPI_CHARACTER, host_all, hnlen,  MPI_CHARACTER, 0, comm, err)
         call MPI_Gather(pid_proc,  1,      MPI_INTEGER,   pid_all,  1,      MPI_INTEGER,   0, comm, err)

         pxsize = 1
         pysize = 1
         pzsize = 1
         nxd    = 1
         nyd    = 1
         nzd    = 1
         nb     = 4

         mpi_magic = .true.

         bnd_xl = 'per'
         bnd_xr = 'per'
         bnd_yl = 'per'
         bnd_yr = 'per'
         bnd_zl = 'per'
         bnd_zr = 'per'

         cfl         = 0.7
         cfr_smooth  = 0.0
         smalld      = 1.e-10
         smallc      = 1.e-10
         smallei     = 1.e-10
         dt_initial  = -1.              !< negative value indicates automatic choice of initial timestep
         dt_max_grow = dt_default_grow !< for sensitive setups consider setting this as low as 1.1
         dt_min      = 0.

         integration_order  = 2

         if(proc == 0) then
            open(1,file=par_file)
               read(unit=1,nml=MPI_BLOCKS,iostat=ierrh)
               call namelist_errh(ierrh,'MPI_BLOCKS')
            close(1)
            open(1,file=par_file)
               read(unit=1,nml=BOUNDARIES,iostat=ierrh)
               call namelist_errh(ierrh,'BOUNDARIES')
            close(1)
            open(1,file=par_file)
               read(unit=1,nml=NUMERICAL_SETUP,iostat=ierrh)
               call namelist_errh(ierrh,'NUMERICAL_SETUP')
            close(1)
            open(1,file=par_file)
               read(unit=1,nml=DOMAIN_SIZES,iostat=ierrh)
               call namelist_errh(ierrh,'DOMAIN_SIZES')
            close(1)
         endif

         nxd = max(1, nxd)
         nyd = max(1, nyd)
         nzd = max(1, nzd)

         if(proc == 0) then

            cbuff(1) = bnd_xl
            cbuff(2) = bnd_xr
            cbuff(3) = bnd_yl
            cbuff(4) = bnd_yr
            cbuff(5) = bnd_zl
            cbuff(6) = bnd_zr

            ibuff(1) = pxsize
            ibuff(2) = pysize
            ibuff(3) = pzsize

            ibuff(4) = integration_order

            ibuff(5) = nxd
            ibuff(6) = nyd
            ibuff(7) = nzd
            ibuff(8) = nb

            rbuff(1) = smalld
            rbuff(10)= smallc
            rbuff(2) = smallei
            rbuff(3) = cfl
            rbuff(4) = cfr_smooth
            rbuff(5) = dt_initial
            rbuff(6) = dt_max_grow
            rbuff(7) = dt_min

            lbuff(1) = mpi_magic

            call MPI_BCAST(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
            call MPI_BCAST(ibuff,           buffer_dim, MPI_INTEGER,          0, comm, ierr)
            call MPI_BCAST(rbuff,           buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
            call MPI_BCAST(lbuff,           buffer_dim, MPI_LOGICAL,          0, comm, ierr)

         else

            call MPI_BCAST(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
            call MPI_BCAST(ibuff,           buffer_dim, MPI_INTEGER,          0, comm, ierr)
            call MPI_BCAST(rbuff,           buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
            call MPI_BCAST(lbuff,           buffer_dim, MPI_LOGICAL,          0, comm, ierr)

            mpi_magic   = lbuff(1)

            smalld      = rbuff(1)
            smallc      = rbuff(10)
            smallei     = rbuff(2)
            cfl         = rbuff(3)
            cfr_smooth  = rbuff(4)
            dt_initial  = rbuff(5)
            dt_max_grow = rbuff(6)
            dt_min      = rbuff(7)

            bnd_xl = cbuff(1)(1:4)
            bnd_xr = cbuff(2)(1:4)
            bnd_yl = cbuff(3)(1:4)
            bnd_yr = cbuff(4)(1:4)
            bnd_zl = cbuff(5)(1:4)
            bnd_zr = cbuff(6)(1:4)

            pxsize = ibuff(1)
            pysize = ibuff(2)
            pzsize = ibuff(3)

            integration_order = ibuff(4)

            nxd = ibuff(5)
            nyd = ibuff(6)
            nzd = ibuff(7)
            nb  = ibuff(8)

         endif

         domsize(:) = [nxd, nyd, nzd]

         bnd_xl_dom = bnd_xl
         bnd_xr_dom = bnd_xr
         bnd_yl_dom = bnd_yl
         bnd_yr_dom = bnd_yr
         bnd_zl_dom = bnd_zl
         bnd_zr_dom = bnd_zr

         psize(1)   = pxsize
         psize(2)   = pysize
         psize(3)   = pzsize

         if (pxsize*pysize*pzsize /= nproc) then
            if(mpi_magic) then
               call divide_domain_voodoo(nproc)
            else
               if(proc == 0)  write(*,*)'nproc =',nproc,' MUST BE EQUAL TO   pxsize*pysize*pzsize =',pxsize*pysize*pzsize
               call MPI_Barrier(MPI_COMM_WORLD, ierr)
               call MPI_Finalize(ierr)
               stop
            endif
         endif

         if (pxsize*pysize*pzsize /= 1) mpi = .true.

         periods(:) = .false.

         if (bnd_xl(1:3) == 'per' .or. bnd_xl(1:3) == 'she') then
            periods(1) = .true.  ! x periodic
            if (bnd_xr(1:3) /= bnd_xl(1:3)) call die("[mpisetup:mpistart] Periodic or shear BC do not match in X-direction")
         endif

         if (bnd_yl(1:3) == 'per') then
            periods(2) = .true.  ! y periodic
            if (bnd_yr(1:3) /= bnd_yl(1:3)) call die("[mpisetup:mpistart] Periodic BC do not match in Y-direction")
         endif

         if (bnd_zl(1:3) == 'per') then
            periods(3) = .true.  ! z periodic
            if (bnd_zr(1:3) /= bnd_zl(1:3)) call die("[mpisetup:mpistart] Periodic BC do not match in Z-direction")
         endif

         reorder = .false.     ! allows processes reordered for efficiency

         call MPI_CART_CREATE(comm, ndims, psize, periods, reorder, comm3d, ierr)
         call MPI_CART_COORDS(comm3d, proc, ndims, pcoords, ierr)

         if(proc == 0) then
            write(*, '("------------------------------------------------------------------------------------------------------")')
            open(3, file=tmp_log_file, status='unknown')
               write(3,'(a,/)')"###############     Initialization     ###############"
               write(3,"(a35,i2)") 'START OF MHD CODE,  No. of procs = ', nproc
               write(3,*)
               write(3,*) 'PROCESSES:'
               do iproc = 0, nproc-1
                  write(3,"(a6,i2,a7,i6,a1,a,a7,a)") " proc=",iproc, &
                           ", pid= ",pid_all(iproc), "@",trim(host_all(iproc)), &
                           ",  cwd=",trim(cwd)
               enddo
               write(3,*)

               write(unit=3,nml=MPI_BLOCKS)
               write(unit=3,nml=BOUNDARIES)
               write(unit=3,nml=NUMERICAL_SETUP)

            close(3)
            write(*,"(a35,i5)") 'START OF MHD CODE,  No. of procs = ', nproc
!            write(*,nml=MPI_BLOCKS)
         endif

         deallocate(host_all)
         deallocate(pid_all)
         deallocate(cwd_all)

! Compute neighbors

         call MPI_CART_SHIFT(comm3d,0,1,procxl,procxr,ierr)   ! x dim
         call MPI_CART_SHIFT(comm3d,1,1,procyl,procyr,ierr)   ! y dim
         call MPI_CART_SHIFT(comm3d,2,1,proczl,proczr,ierr)   ! z dim

         if(bnd_xl(1:3) == 'cor' .and. bnd_yl(1:3) == 'cor' ) then
            if(pcoords(1) == 0 .and. pcoords(2) > 0) then
               coords = (/pcoords(2),pcoords(1),pcoords(3)/)
               call MPI_Cart_rank(comm3d,coords,procxyl,ierr)
            else
               procxyl = MPI_PROC_NULL
            endif
            if(pcoords(2) == 0 .and. pcoords(1) > 0 ) then
               coords = (/pcoords(2),pcoords(1),pcoords(3)/)
               call MPI_Cart_rank(comm3d,coords,procyxl,ierr)
            else
               procyxl = MPI_PROC_NULL
            endif
         endif

         if(bnd_xr(1:3) == 'cor' .and. bnd_yr(1:3) == 'cor' ) then
            if(pcoords(1) == psize(1)-1 .and. pcoords(2) < psize(2)-1) then
               coords = (/pcoords(2),pcoords(1),pcoords(3)/)
               call MPI_Cart_rank(comm3d,coords,procxyr,ierr)
            else
               procxyr = MPI_PROC_NULL
            endif
            if(pcoords(2) == psize(2)-1 .and. pcoords(1) < psize(2)-1 ) then
               coords = (/pcoords(2),pcoords(1),pcoords(3)/)
               call MPI_Cart_rank(comm3d,coords,procyxr,ierr)
            else
               procyxr = MPI_PROC_NULL
            endif
         endif


#ifdef SHEAR_BND
         if(pysize > 1) stop 'Shear-pediodic boundary conditions do not permit pysize > 1'

#ifndef FFTW
         if(pcoords(1) == 0) then
            bnd_xl = 'she'
         else
            bnd_xl = 'mpi'
         endif

         if(pcoords(1) == pxsize-1) then
            bnd_xr = 'she'
         else
            bnd_xr = 'mpi'
         endif
#endif

#else /* SHEAR_BND */
         if(procxl /= MPI_PROC_NULL .and. procxl /= proc) bnd_xl = 'mpi'
         if(procxr /= MPI_PROC_NULL .and. procxr /= proc) bnd_xr = 'mpi'
#endif /* SHEAR_BND */

         if(procyl /= MPI_PROC_NULL .and. procyl /= proc) bnd_yl = 'mpi'
         if(procyr /= MPI_PROC_NULL .and. procyr /= proc) bnd_yr = 'mpi'

         if(proczl /= MPI_PROC_NULL .and. proczl /= proc) bnd_zl = 'mpi'
         if(proczr /= MPI_PROC_NULL .and. proczr /= proc) bnd_zr = 'mpi'

#ifdef DEBUG
         write(*,*) 'xdir: ',procxl, proc, procxr
         write(*,*) 'ydir: ',procyl, proc, procyr
         write(*,*) 'zdir: ',proczl, proc, proczr
         write(*,*)
#endif /* DEBUG */

         if(integration_order > 2) then
            stop 'For "ORIG" scheme integration_order must be 1 or 2'
         endif

         dt_old = -1.
         if (dt_max_grow < 1.01) then
            if (proc == 0) write(*,'(2(a,g10.3))')"[mpisetup:mpistart] dt_max_grow = ",dt_max_grow," is way too low. Resetting to ",dt_default_grow
            dt_max_grow = dt_default_grow
         end if

      end subroutine mpistart

!-----------------------------------------------------------------------------

      subroutine mpistop

         implicit none

         call MPI_Comm_free(comm3d, ierr)

         if (domsize(1) /= 1) then
            call MPI_Type_free(MPI_YZ_LEFT_BND, ierr)
            call MPI_Type_free(MPI_YZ_LEFT_DOM, ierr)
            call MPI_Type_free(MPI_YZ_RIGHT_DOM, ierr)
            call MPI_Type_free(MPI_YZ_RIGHT_BND, ierr)
            call MPI_Type_free(MAG_YZ_LEFT_BND, ierr)
            call MPI_Type_free(MAG_YZ_LEFT_DOM, ierr)
            call MPI_Type_free(MAG_YZ_RIGHT_DOM, ierr)
            call MPI_Type_free(MAG_YZ_RIGHT_BND, ierr)
         end if

         if (domsize(2) /= 1) then
            call MPI_Type_free(MPI_XZ_LEFT_BND, ierr)
            call MPI_Type_free(MPI_XZ_LEFT_DOM, ierr)
            call MPI_Type_free(MPI_XZ_RIGHT_DOM, ierr)
            call MPI_Type_free(MPI_XZ_RIGHT_BND, ierr)
            call MPI_Type_free(MAG_XZ_LEFT_BND, ierr)
            call MPI_Type_free(MAG_XZ_LEFT_DOM, ierr)
            call MPI_Type_free(MAG_XZ_RIGHT_DOM, ierr)
            call MPI_Type_free(MAG_XZ_RIGHT_BND, ierr)
         end if

         if (domsize(3) /= 1) then
            call MPI_Type_free(MPI_XY_LEFT_BND, ierr)
            call MPI_Type_free(MPI_XY_LEFT_DOM, ierr)
            call MPI_Type_free(MPI_XY_RIGHT_DOM, ierr)
            call MPI_Type_free(MPI_XY_RIGHT_BND, ierr)
            call MPI_Type_free(MAG_XY_LEFT_BND, ierr)
            call MPI_Type_free(MAG_XY_LEFT_DOM, ierr)
            call MPI_Type_free(MAG_XY_RIGHT_DOM, ierr)
            call MPI_Type_free(MAG_XY_RIGHT_BND, ierr)
         end if

         call MPI_BARRIER(comm,ierr)
         call sleep(5)
         call MPI_FINALIZE(ierr)

      end subroutine mpistop

!-----------------------------------------------------------------------------

      subroutine mpifind(var, what, loc_arr, loc_proc)

         implicit none

         character(len=3), intent(in) :: what
         real       :: var
         real, dimension(2)    :: rsend, rrecv
         integer, dimension(3) :: loc_arr
         integer               :: loc_proc

         rsend(1) = var
         rsend(2) = proc

         select case (what(1:3))
            case('min')
               CALL MPI_REDUCE(rsend, rrecv, 1, MPI_2DOUBLE_PRECISION, &
                                         MPI_MINLOC, 0, comm, ierr)
            case('max')
               CALL MPI_REDUCE(rsend, rrecv, 1, MPI_2DOUBLE_PRECISION, &
                                         MPI_MAXLOC, 0, comm, ierr)
            case default
               Write(*,*) 'mpifind: actual parameter "', what, '"is not allowed'
         end select

         if(proc == 0) then
            var = rrecv(1)
            loc_proc = rrecv(2)
         endif

         call MPI_BCAST(loc_proc, 1, MPI_INTEGER, 0, comm, ierr)

         if(loc_proc /= 0) then
            if(proc == loc_proc) then
               CALL MPI_SEND  (loc_arr, 3, MPI_INTEGER,     0, 11, comm, ierr)
            else if(proc == 0) then
               CALL MPI_RECV  (loc_arr, 3, MPI_INTEGER, loc_proc, 11, comm, status, ierr)
            endif
         endif

      end subroutine mpifind

!------------------------------------------------------------------------------------------
! Must be called by all procs to avoid communication and ensure that every proc has
! proper psize, pxsize, pysize, pzsize
!
! This routine tries to divide the computational domain into local domains.
! The goal is to minimize the ratio of longest to shortest edge to minimize the amount of inter-process communication.
! If the benchmarks show that some direction should be partitioned in more pieces than other directions, implement appropriate weighting in j1, j2 and j3 calculations.
!
! For some weird domains and PE counts this routine may find tiling that does not satisfy multigrid restrictions even if there is some. In such case divide domain manually.
!

      subroutine divide_domain_voodoo(np)

         use errh,      only: die
         use constants, only: some_primes

         implicit none

         integer, intent(in) :: np

         integer :: j1, j2, j3, jj, n, p
         integer, dimension(3) :: ldom

         ldom(1:3) = domsize(3:1:-1) ! Maxloc returns first occurrence of max, reversing direction order (to ZYX)  gives better cache utilisation.
         n = np
         psize(:) = 1

         do p = size(some_primes), 1, -1 ! start from largest defined primes, continue down to 2
            do while (mod(n, some_primes(p))==0)

               jj = 0
               j1 = sum(maxloc(ldom), 1) ! First try the longest edge; note the trick to make a scalar from 1-element vector without assigment to another variable
               if (mod(ldom(j1), some_primes(p))==0) then
                  jj = j1
               else
                  j2 = 1 + mod(j1 + 0, ndims)
                  j3 = 1 + mod(j1 + ndims -2, ndims)
                  if (ldom(j2) > ldom(j3) .and. mod(ldom(j2), some_primes(p))==0) jj = j2 ! middle edge ...
                  if (jj == 0 .and. mod(ldom(j3), some_primes(p))==0) jj = j3 ! try the shortest edge on last resort
               end if

               if (jj == 0) then
                  call die("[divide_domain_voodoo]: Can't find divisible edge")
               else
                  psize(jj) = psize(jj) * some_primes(p)
                  n         = n         / some_primes(p)
                  ldom(jj)  = ldom(jj)  / some_primes(p)
               end if

            end do
         end do

         if (n /= 1) call die("[divide_domain_voodoo]: I am not that intelligent") ! np has too big prime factors

         pxsize = psize(3) ! directions were reverted at ldom assigment
         pysize = psize(2)
         pzsize = psize(1)

         psize = [ pxsize, pysize, pzsize ]

         if (proc == 0 .and. np > 1) &
              write(*,'(a,3i4,a,3i6,a)')"[mpisetup:divide_domain_voodoo] Domain divided to [",psize(:)," ] pieces, each of [",ldom(3:1:-1)," ] cells."

      end subroutine divide_domain_voodoo

end module mpisetup
