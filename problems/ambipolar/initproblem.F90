#include "piernik.def"

module initproblem

! Initial condition for Parker instability in realistic galactic gravity
! Written by: M. Hanasz, February 2006
! Adopted for several fluids by: D. Woltanski, January 2008

   use arrays
   use start
   use grid
   use fluid_boundaries
   use hydrostatic

   real nbx0,nby0,nbz0, collf, n_x
   real x0,y0,z0,r0,bfaq
   character problem_name*32,run_id*3

   namelist /PROBLEM_CONTROL/  problem_name, run_id, &
                               bfaq, &
                               nbx0,nby0,nbz0, &
                               collf, n_x, &
                               x0,y0,z0, r0
   contains

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      implicit none

      problem_name = 'xxx'
      run_id  = 'aaa'
      bfaq    = 1.0
      nbx0    = 0.0
      nby0    = 1.0
      nbz0    = 0.0
      collf   = 1.0
      n_x     = 3.0
      x0      = 0.0
      y0      = 0.0
      z0      = 0.0

      if(proc .eq. 0) then
         open(1,file='problem.par')
            read(unit=1,nml=PROBLEM_CONTROL)
            write(*,nml=PROBLEM_CONTROL)
         close(1)
         open(3, file='tmp.log', position='append')
            write(3,nml=PROBLEM_CONTROL)
            write(3,*)
         close(3)
      endif

      if(proc .eq. 0) then

         cbuff(1) =  problem_name
         cbuff(2) =  run_id

         rbuff(1) = bfaq
         rbuff(2) = nbx0
         rbuff(3) = nby0
         rbuff(4) = nbz0
         rbuff(5) = collf
         rbuff(7) = n_x
         rbuff(8) = x0
         rbuff(9) = y0
         rbuff(10)= z0

         call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
         call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
         call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      else

         call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
         call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
         call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

         problem_name = cbuff(1)
         run_id       = cbuff(2)

         bfaq         = rbuff(1)
         nbx0         = rbuff(2)
         nby0         = rbuff(3)
         nbz0         = rbuff(4)
         collf        = rbuff(5)
         n_x          = rbuff(7)
         x0           = rbuff(8)
         y0           = rbuff(9)
         z0           = rbuff(10)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob

      use constants, only :cm, gram, sek, kpc, r_gc_sun, mp
      use arrays, only : nfluid
      use start,  only : r_gc
#ifdef GALACTIC_DISK
      use gravity_user, only : gpotdisk, gpothalo, gpotbulge
#endif /* GALACTIC_DISK */
      implicit none

      integer i,j,k
      real b0, vz
      real mass,emag,emaga
      real,dimension(nfluid) :: eint,einta,ekin,ekina
      real dcmol, dcneut, dcion, dchot, d0, zk

#ifdef GALACTIC_DISK
      allocate(gpotdisk(nx,ny,nz),gpothalo(nx,ny,nz),gpotbulge(nx,ny,nz))
#endif /* GALACTIC_DISK */

!    call read_problem_par

!   Secondary parameters

      allocate(dprof(nz))
      collfaq = collf*3.4e15*(8000.0/8000.0)**0.375*cm**3/gram/sek
      write(*,*) 'collfaq = ',collfaq,' [lengthunit^3/massunit/timeunit]'
      rc = r_gc
      dcmol = 2.6e20/(cm**2)*exp(-((rc - 4.5*kpc)**2-(r_gc_sun - 4.5*kpc)**2)/(2.9*kpc)**2)
      dcneut = 6.2e20/(cm**2)
      dcion =  1.46e20/(cm**2)*exp(-(rc**2 - r_gc_sun**2)/(37.0*kpc)**2) &
             + 1.20e18/(cm**2)*exp(-((rc - 4.0*kpc)**2 - (r_gc_sun - 4.0*kpc)**2)/(2.0*kpc)**2)
      dchot = 4.4e18/(cm**2)*(0.12*exp(-(rc-r_gc_sun)/4.9/kpc)+0.88*exp(-((rc-4.5*kpc)**2-(r_gc_sun-4.5*kpc)**2)/(2.9*kpc)**2))
      d0 = 1.36*mp*(dcmol+dcneut+dcion+dchot)
      b0 = bfaq*sqrt(2.*alpha*d0*c_si**2)


      einta = 0.0
      ekina = 0.0
      emaga = 0.0
      call hydrostatic_zeq(nx/2, ny/2, d0, dprof)
      do i = 1,nx
         do j = 1,ny
!        d0 = 1.36*mp*(dcmol+dcneut+dcion+dchot)
!        call hydrostatic_zeq(i, j, d0, dprof)
            do k = 1,nz
               zk=z(k)
               u(idna,i,j,k)   = 2.0*smalld+dprof(k)*0.8*(1.+0.25*abs(zk)/kpc)

               u(imxa,i,j,k) = 0.0
               u(imya,i,j,k) = 0.0
               u(imza,i,j,k) = 0.0


#ifndef ISO
               u(iena,i,j,k)   = c_si**2/(gamma-1.0) * u(idna,i,j,k) &
                               + 0.5*(u(imxa,i,j,k)**2+u(imya,i,j,k)**2+u(imza,i,j,k)**2)/u(idna,i,j,k)
               eint = c_si**2/(gamma-1.0)
               ekin = 0.5*(u(imxa,i,j,k)**2+u(imya,i,j,k)**2+u(imza,i,j,k)**2)/u(idna,i,j,k)
               einta = max(einta,eint)
               ekina = max(ekina,ekin)
#endif
            enddo
         enddo
      enddo

      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
               b(ibx,i,j,k)   = b0
               b(iby,i,j,k)   = 0.0
               b(ibz,i,j,k)   = 0.0
#ifndef ISO
               u(iena,i,j,k)   = u(iena,i,j,k) +0.5*sum(b(:,i,j,k)**2,1)
               emag = 0.5*sum(b(:,i,j,k)**2,1)
               emaga = max(emaga,emag)
#endif
            enddo
         enddo
      enddo
      write(*,*) 'eint = ',einta
      write(*,*) 'ekin = ',ekina
      write(*,*) 'emag = ',emaga

      if(allocated(dprof)) deallocate(dprof)


      return
   end subroutine init_prob

!-----------------------------------------------------------------------------



end module initproblem

