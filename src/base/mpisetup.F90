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
module mpisetup
   implicit none
   include 'mpif.h'
   integer :: nproc, proc, ierr , rc, info
   integer :: status(MPI_STATUS_SIZE,4)
   integer, dimension(4) :: req, err

   real                  :: t,dt

   integer, parameter    :: ndims = 3       ! 3D grid
   integer               :: comm, comm3d
   integer, dimension(3) :: psize, pcoords, coords
   logical               :: periods(3), reorder
   integer ::   procxl, procxr, procyl, procyr, proczl, proczr, procxyl, procyxl, procxyr, procyxr
   integer ::   pxleft, pxright, pyleft, pyright, pzleft, pzright

   integer,   parameter             :: buffer_dim=200
   character, dimension(buffer_dim) :: cbuff*32
   integer,   dimension(buffer_dim) :: ibuff
   real,      dimension(buffer_dim) :: rbuff

   integer :: pxsize, pysize, pzsize
   namelist /MPI_BLOCKS/ pxsize, pysize, pzsize

   character(len=4) :: bnd_xl, bnd_xr, bnd_yl, bnd_yr, bnd_zl, bnd_zr
   character(len=4) :: bnd_xl_dom, bnd_xr_dom, bnd_yl_dom, bnd_yr_dom, bnd_zl_dom, bnd_zr_dom
   namelist /BOUNDARIES/ bnd_xl, bnd_xr, bnd_yl, bnd_yr, bnd_zl, bnd_zr

   logical     :: mpi
   character(len=80)   :: cwd

   integer :: MPI_XZ_LEFT_BND, MPI_XZ_RIGHT_BND
   integer :: MPI_XZ_LEFT_DOM, MPI_XZ_RIGHT_DOM
   integer :: MPI_XY_LEFT_BND, MPI_XY_RIGHT_BND
   integer :: MPI_XY_LEFT_DOM, MPI_XY_RIGHT_DOM
   integer :: MPI_YZ_LEFT_BND, MPI_YZ_RIGHT_BND
   integer :: MPI_YZ_LEFT_DOM, MPI_YZ_RIGHT_DOM

   integer :: MAG_XZ_LEFT_BND, MAG_XZ_RIGHT_BND
   integer :: MAG_XZ_LEFT_DOM, MAG_XZ_RIGHT_DOM
   integer :: MAG_XY_LEFT_BND, MAG_XY_RIGHT_BND
   integer :: MAG_XY_LEFT_DOM, MAG_XY_RIGHT_DOM
   integer :: MAG_YZ_LEFT_BND, MAG_YZ_RIGHT_BND
   integer :: MAG_YZ_LEFT_DOM, MAG_YZ_RIGHT_DOM

   contains

!-----------------------------------------------------------------------------
      subroutine mpistart
         use errh, only : namelist_errh
         implicit none
         integer :: iproc, ierrh

         character(LEN=80) :: cwd_proc,  cwd_all (0:4096)
         character(LEN=8)  :: host_proc, host_all(0:4096)
         integer           :: pid_proc,  pid_all (0:4096)

         integer(kind=1)    :: getcwd, hostnm
         integer(kind=4)    :: getpid
         character(LEN=100) :: par_file, tmp_log_file
         integer :: cwd_status
         logical :: par_file_exist

         call MPI_INIT( ierr )
         call MPI_COMM_RANK(MPI_COMM_WORLD, proc, ierr)
         comm = MPI_COMM_WORLD
         info = MPI_INFO_NULL
         call MPI_COMM_SIZE(comm, nproc, ierr)

         pid_proc = getpid()
         status = hostnm(host_proc)
         cwd_status =  getcwd(cwd_proc)
         if(cwd_status /= 0) stop 'mpisetup: problems accessing working directory'
#ifdef DEBUG
         write(*,*) 'pid  in mpisetup: ',pid_proc
         write(*,*) 'host in mpisetup: ',host_proc
         write(*,*) 'cwd  in mpisetup: ',cwd_proc
#endif /* DEBUG */
         if(proc == 0) then
            par_file = trim(cwd)//'/problem.par'
            inquire(file=par_file, exist=par_file_exist)
            if(.not. par_file_exist) stop '"problem.par" does not exist in the working directory'
            tmp_log_file = trim(cwd)//'/tmp.log'
         endif

         call MPI_GATHER ( cwd_proc, 80, MPI_CHARACTER, &
                           cwd_all,  80, MPI_CHARACTER, &
                           0, comm,err )

         call MPI_GATHER ( host_proc, 8, MPI_CHARACTER, &
                           host_all,  8, MPI_CHARACTER, &
                           0, comm,err )

         call MPI_GATHER ( pid_proc, 1, MPI_INTEGER, &
                           pid_all,  1, MPI_INTEGER, &
                           0, comm,err )

         pxsize = 1
         pysize = 1
         pzsize = 1

         bnd_xl = 'per'
         bnd_xr = 'per'
         bnd_yl = 'per'
         bnd_yr = 'per'
         bnd_zl = 'per'
         bnd_zr = 'per'

         if(proc == 0) then
            open(1,file=par_file)
               read(unit=1,nml=MPI_BLOCKS,iostat=ierrh)
               call namelist_errh(ierrh,'MPI_BLOCKS')
               read(unit=1,nml=BOUNDARIES,iostat=ierrh)
               call namelist_errh(ierrh,'BOUNDARIES')
            close(1)
         endif

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

            call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
            call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)

         else

            call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
            call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)

            bnd_xl = cbuff(1)(1:4)
            bnd_xr = cbuff(2)(1:4)
            bnd_yl = cbuff(3)(1:4)
            bnd_yr = cbuff(4)(1:4)
            bnd_zl = cbuff(5)(1:4)
            bnd_zr = cbuff(6)(1:4)

            pxsize = ibuff(1)
            pysize = ibuff(2)
            pzsize = ibuff(3)

         endif

         bnd_xl_dom = bnd_xl
         bnd_xr_dom = bnd_xr
         bnd_yl_dom = bnd_yl
         bnd_yr_dom = bnd_yr
         bnd_zl_dom = bnd_zl
         bnd_zr_dom = bnd_zr

         psize(1)   = pxsize
         psize(2)   = pysize
         psize(3)   = pzsize

         if(pxsize*pysize*pzsize /= nproc) then
            if(proc == 0)  write(*,*) &
               'nproc =',nproc,' MUST BE EQUAL TO   pxsize*pysize*pzsize =',pxsize*pysize*pzsize
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
            call MPI_Finalize(ierr)
            stop
         endif

         if(pxsize*pysize*pzsize /= 1) mpi = .true.


         if(bnd_xl(1:3) == 'per' .or. bnd_xl(1:3) == 'she') then
            periods(1) = .true.  ! x periodic
         else
            periods(1) = .false.
         endif
         if(bnd_yl(1:3) == 'per') then
            periods(2) = .true.  ! y periodic
         else
            periods(2) = .false.
         endif
         if(bnd_zl(1:3) == 'per') then
            periods(3) = .true.  ! z periodic
         else
            periods(3) = .false.
         endif

         reorder = .false.     ! allows processes reordered for efficiency

         call MPI_CART_CREATE(comm, ndims, psize, periods, reorder, comm3d, ierr)
         call MPI_CART_COORDS(comm3d, proc, ndims, pcoords, ierr)

         if(proc == 0) then
            open(3, file=tmp_log_file, status='unknown')
               write(3,"(a35,i2)") 'START OF MHD CODE,  No. of procs = ', nproc
               write(3,*)
               write(3,*) 'PROCESSES:'
               do iproc = 0, nproc-1
                  write(3,"(a6,i2,a7,i6,a1,a,a7,a)") ' proc=',iproc,',  &
                        pid=',pid_all(iproc), '@',trim(host_all(iproc)), ',  cwd=',trim(cwd)
               enddo
               write(3,*)

               write(unit=3,nml=MPI_BLOCKS)
               write(unit=3,nml=BOUNDARIES)

            close(3)
            write(*,*)
            write(*,"(a35,i2)") 'START OF MHD CODE,  No. of procs = ', nproc
            write(*,*)
            write(*,nml=MPI_BLOCKS)
            write(*,*)
         endif

! Compute neighbors

         pxleft  = pcoords(1) - 1
         pxright = pcoords(1) + 1
         pyleft  = pcoords(2) - 1
         pyright = pcoords(2) + 1
         pzleft  = pcoords(3) - 1
         pzright = pcoords(3) + 1

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


#ifdef SHEAR
         if(pysize > 1) stop 'Shear-pediodic boundary conditions do not permit pysize > 1'

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
#else /* SHEAR */
         if(procxl /= MPI_PROC_NULL .and. procxl /= proc) bnd_xl = 'mpi'
         if(procxr /= MPI_PROC_NULL .and. procxr /= proc) bnd_xr = 'mpi'
#endif /* SHEAR */

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

      end subroutine mpistart


!-----------------------------------------------------------------------------


      subroutine mpistop

         implicit none

         call MPI_FINALIZE()

      end subroutine mpistop

!-----------------------------------------------------------------------------


      subroutine mpifind(var, what, loc_arr, loc_proc)

         implicit none
         character what*(*)
         real       :: var
         real, dimension(2)    :: rsend, rrecv
         integer, dimension(3) :: loc_arr
         integer               :: loc_proc

         rsend(1) = var
         rsend(2) = proc

         select case (trim(what))
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

end module mpisetup


