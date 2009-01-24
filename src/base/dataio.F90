! $Id$
#include "piernik.def"
!=====================================================================
!!
!!  dataio module responsible for data output
!!
!=====================================================================
!
module dataio

! Written by G. Kowal
! Modified for this code and extended by M.Hanasz
  use types
  use mpi_setup
#ifdef SN_SRC
  use sn_sources
#endif /* SN_SRC */

  implicit none

  integer, parameter    :: nvarsmx = 16 
  character(len=3)      :: new_id
  character(len=16)     :: restart, domain, mag_center
  integer               :: nrestart, resdel
  real                  :: dt_hdf, dt_res, dt_tsl, dt_log
  integer               :: min_disk_space_MB, sleep_minutes, sleep_seconds
  character(len=160)    :: user_message_file, system_message_file
  integer               :: ix,iy,iz, iv
  character(len=4), dimension(nvarsmx) :: vars

  integer               :: tsl_lun = 2, log_lun = 3
  integer               :: nhdf, nres, ntsl, nlog
  integer               :: step_hdf, step_res, nhdf_start, nres_start
  real                  :: t_start, last_hdf_time
  character(len=128)    :: log_file
  character(len=2)      :: pc1,pc2,pc3
  logical               :: tsl_firstcall

  logical               :: wait

  integer               :: nchar
  character             :: msg*16
  real                  :: msg_param

  character             :: hostfull*(80), host*8, fhost*10, fpid*10
  integer               :: pid, uid, ihost, scstatus
  real                  :: vx_max, vy_max, vz_max, va2max,va_max, &
                           cs2max, cs_max, &
                           dens_min, dens_max, pres_min, pres_max

#ifdef MAGNETIC
  real                  :: divb_max, b_min, b_max
#endif MAGNETIC
#ifndef ISO       
  real                  :: temp_min, temp_max, 
#endif /* ISO */     
#ifdef COSM_RAYS
  real                  :: encr_min, encr_max
#endif /* COSM_RAYS */

  namelist /RESTART_CONTROL/ restart, new_id, nrestart, resdel
  namelist /OUTPUT_CONTROL/ dt_hdf, dt_res, dt_tsl, dt_log, &
                            domain, vars, mag_center, ix, iy, iz,&
                            min_disk_space_MB, sleep_minutes, sleep_seconds, &
                            user_message_file, system_message_file


  contains

!---------------------------------------------------------------------
!
! inititalizes dataio parameters
!
!---------------------------------------------------------------------
!
  subroutine init_dataio
    implicit none
    integer(kind=1) :: getpid
    integer(kind=1) :: hostnm
    character(LEN=100) :: par_file, tmp_log_file



    restart = 'last'   ! 'last': autom. wybor ostatniego
                       ! niezaleznie od wartosci "nrestart"
                       ! cokolwiek innego: decyduje "nrestart"
    new_id  = ''
    nrestart=  3
    resdel  = 0

    dt_hdf = 0.0
    dt_res = 0.0
    dt_tsl = 0.0
    dt_log = 0.0
    domain = 'phys_domain'
    vars(:)   = '    '
    mag_center= 'no'
    min_disk_space_MB = 100
    sleep_minutes   = 0
    sleep_seconds   = 0
    user_message_file   = trim(cwd)//'/msg'
    system_message_file = '/tmp/piernik_msg'

    wait  = .false.
    tsl_firstcall = .true.

    nhdf  = 0
    ntsl  = 0
    nres  = 0
    nlog  = 0

    step_hdf  = -1
    step_res  = -1

    pc1 = '00'
    pc2 = '00'
    pc3 = '00'

    if(psize(1) .gt. 1) pc1 = '0x'
    if(psize(2) .gt. 1) pc2 = '0x'
    if(psize(3) .gt. 1) pc3 = '0x'

    pid = getpid()

    scstatus = hostnm(hostfull)
    ihost = index(hostfull,'.')
    if (ihost .eq. 0) ihost = index(hostfull,' ')
    host = hostfull(1:ihost-1)

   if(proc .eq. 0) then
      par_file = trim(cwd)//'/problem.par'
      tmp_log_file = trim(cwd)//'/tmp.log'
      open(1,file=par_file)
         read(unit=1,nml=OUTPUT_CONTROL)
      close(1)
      open(1,file=par_file)
         read(unit=1,nml=RESTART_CONTROL)
      close(1)
      open(3, file=tmp_log_file, position='append')
         write(unit=3,nml=OUTPUT_CONTROL)
         write(unit=3,nml=RESTART_CONTROL)
      close(3)

!  namelist /RESTART_CONTROL/ restart, new_id, nrestart, resdel

      cbuff(20) = restart
      cbuff(21) = new_id

      ibuff(20) = nrestart
      ibuff(21) = resdel

!  namelist /OUTPUT_CONTROL/ dt_hdf, dt_res, dt_tsl, domain, vars, mag_center, &
!                            min_disk_space_MB, sleep_minutes, ix, iy, iz&
!                            user_message_file, system_message_file

      ibuff(40) = min_disk_space_MB
      ibuff(41) = sleep_minutes
      ibuff(42) = sleep_seconds
      ibuff(43) = ix
      ibuff(44) = iy
      ibuff(45) = iz

      rbuff(40) = dt_hdf
      rbuff(41) = dt_res
      rbuff(42) = dt_tsl
      rbuff(43) = dt_log

      cbuff(40) = domain

      do iv = 1, nvarsmx
        cbuff(40+iv) = vars(iv)
      enddo

      cbuff(60) = mag_center
      cbuff(61) = user_message_file(1:32)
      cbuff(62) = system_message_file(1:32)

      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    else

      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

!  namelist /RESTART_CONTROL/ restart, new_id, nrestart, resdel

      restart             = trim(cbuff(20))
      new_id              = trim(cbuff(21))

      nrestart            = ibuff(20)
      resdel              = ibuff(21)

!  namelist /OUTPUT_CONTROL/ dt_hdf, dt_res, dt_tsl, domain, vars, mag_center, &
!                            min_disk_space_MB, sleep_minutes, ix, iy, iz &
!                            user_message_file, system_message_file

      min_disk_space_MB   = ibuff(40)
      sleep_minutes       = ibuff(41)
      sleep_seconds       = ibuff(42)
      ix                  = ibuff(43)
      iy                  = ibuff(44)
      iz                  = ibuff(45)


      dt_hdf              = rbuff(40)
      dt_res              = rbuff(41)
      dt_tsl              = rbuff(42)
      dt_log              = rbuff(43)

      domain              = trim(cbuff(40))
      do iv=1, nvarsmx
        vars(iv)          = trim(cbuff(40+iv))
      enddo

      mag_center          = trim(cbuff(60))

      user_message_file   = trim(cbuff(61))
      system_message_file = trim(cbuff(62))

   endif

   last_hdf_time = -dt_hdf


  end subroutine init_dataio

!---------------------------------------------------------------------
!
! controls data dumping
!
!---------------------------------------------------------------------
!
  subroutine write_data(output)
    use start, only:  t, dt, nstep
#ifdef USER_IO
    use init_problem, only : user_io_routine
#endif /* USER_IO */
    implicit none
    character  :: output*3

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if (output .eq. 'log' .or. output .eq. 'end') then
        call write_log
    endif

    if (output .eq. 'log' .or. output .eq. 'end') then
          call write_timeslice
    endif

!    CALL checkdf

#ifdef USER_IO
      call user_io_routine
#endif /* USER_IO */

#ifdef HDFSWEEP
    if (output .ne. 'gpt') then
      call write_log
      call write_timeslice
#else /* HDFSWEEP */
    if (dt_hdf .gt. 0.0 .and. nstep .gt. step_hdf .and. output .ne. 'gpt') then
#endif /* HDFSWEEP */

!      if ((nhdf-nhdf_start) .lt. (int((t-t_start) / dt_hdf) + 1) &
      if ((t-last_hdf_time) .ge. dt_hdf &
                .or. output .eq. 'hdf' .or. output .eq. 'end') then
        call write_hdf
        if((t-last_hdf_time) .ge. dt_hdf) last_hdf_time = last_hdf_time + dt_hdf
        if((t-last_hdf_time) .ge. dt_hdf) last_hdf_time = t ! dodatkowa regulacja w przypadku zmiany dt_hdf na mniejsze przez msg
        nhdf = nhdf + 1
        step_hdf = nstep
      endif
    endif



    if (dt_res .gt. 0.0 .and. nstep .gt. step_res) then
      if ((nres-nres_start) .lt. (int((t-t_start) / dt_res) + 1) &
                .or. output .eq. 'res' .or. output .eq. 'end') then
         if (nres > 0) then
           call write_restart
        endif
        nres = nres + 1
        step_res = nstep
      endif
    endif


  end subroutine write_data


!---------------------------------------------------------------------
!
! dumps data to hdf file
!
!---------------------------------------------------------------------
!
   subroutine write_hdf
      use mpi_setup, only: cwd
      use arrays, only : wa,outwa,outwb,outwc,b,u
      use start, only : nstep,t,dt
      use grid, only : dx,dy,dz,xmin,xmax,ymin,ymax,zmin,zmax,nxd,nyd,nzd,nb
      use grid, only : nx,ny,nz,nxb,nyb,nzb,x,y,z
      use init_problem, only : problem_name, run_id
      
      use fluidindex,  only : nfluid    
      use fluidindex,  only : ibx,iby,ibz
      use fluidindex,  only : nvar, iarr_all_dn,iarr_all_mx,iarr_all_my,iarr_all_mz

#ifndef STANDARD
      use constants, only : Gs
#endif /* STANDARD */

#ifndef ISO
      use fluidindex, only : iarr_all_en
#endif /* ISO */

#ifdef COSM_RAYS
      use initcosmicrays, only : iecr
#endif /* COSM_RAYS */

#ifdef GRAV
      use arrays, only : gp
#endif /* GRAV */

    implicit none

    character(len=128) :: file_name_hdf,file_name_disp
    character(LEN=4)   :: varname

    integer :: sd_id, sds_id, dim_id, iostatus
    integer :: rank, comp_type
    integer, dimension(3) :: dims, istart, stride
    integer, dimension(1) :: comp_prm
    integer :: sfstart, sfend, sfsnatt, sfcreate, sfwdata, sfscompress, sfendacc &
             , sfdimid, sfsdmname, sfsdscale, sfsdmstr

    integer :: iv, ifl, iw
    integer :: nxo = 1, nyo = 1, nzo = 1, &
               iso = 1, jso = 1, kso = 1, &
               ieo = 1, jeo = 1, keo = 1
	       
	       
    real(kind=4), dimension(:,:,:), allocatable :: tmparr

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if(domain .eq. 'full_domain') then
      nxo = nx
      nyo = ny
      nzo = nz
      iso = 1
      ieo = nx
      jso = 1
      jeo = ny
      if(nzd /= 1) then
        kso = 1
        keo = nz
      else
        kso = 1
        keo = 1
      endif

    else if(domain .eq. 'phys_domain') then
      nxo = nxb
      nyo = nyb
      nzo = nzb
      iso    = nb+1
      ieo    = nb+nxb
      jso    = nb+1
      jeo    = nb+nyb
      if(nzd /= 1) then
        kso    = nb+1
        keo    = nb+nzb
      else
        kso = 1
        keo = 1
      endif
    endif

    allocate(tmparr(iso:ieo,jso:jeo,kso:keo))

    rank = 3
    comp_type = 4
    comp_prm(1) = 6
    dims(1) = nxo
    dims(2) = nyo
    dims(3) = nzo
    istart(:) = 0
    stride(:) = 1

!  generate filename
!

    write (file_name_hdf,'(a,a1,a,a1,a3,a1,i2.2,a1,i2.2,a1,i2.2,a1,i4.4,a4)') &
              trim(cwd),'/',trim(problem_name),'_', run_id, &
              '_',pcoords(1),'_',pcoords(2),'_',pcoords(3),'_',nhdf, '.hdf'

    write (file_name_disp,'(a,a1,a3,a1,a2,a1,a2,a1,a2,a1,i4.4,a4)') &
              trim(problem_name),'_', run_id,'_',pc1,'_',pc2,'_',pc3,'_',nhdf, '.hdf'

    sd_id = sfstart(trim(file_name_hdf), 4)

! write attributes
!
    iostatus = sfsnatt( sd_id, 'problem' , 4, 32, problem_name)
    iostatus = sfsnatt( sd_id, 'run_id'  , 4, 32, run_id      )
    iostatus = sfsnatt( sd_id, 'domain'  , 4, 32, domain      )

    iostatus = sfsnatt( sd_id, 'psize1'  , 23, 1, psize(1) )
    iostatus = sfsnatt( sd_id, 'psize2'  , 23, 1, psize(2) )
    iostatus = sfsnatt( sd_id, 'psize3'  , 23, 1, psize(3) )

    iostatus = sfsnatt( sd_id, 'pcoords1', 23, 1, pcoords(1) )
    iostatus = sfsnatt( sd_id, 'pcoords2', 23, 1, pcoords(2) )
    iostatus = sfsnatt( sd_id, 'pcoords3', 23, 1, pcoords(3) )

    iostatus = sfsnatt( sd_id, 'dims1'   , 23, 1, dims(1)    )
    iostatus = sfsnatt( sd_id, 'dims2'   , 23, 1, dims(2)    )
    iostatus = sfsnatt( sd_id, 'dims3'   , 23, 1, dims(3)    )

    iostatus = sfsnatt( sd_id, 'nxd'     , 23, 1, nxd     )
    iostatus = sfsnatt( sd_id, 'nyd'     , 23, 1, nyd     )
    iostatus = sfsnatt( sd_id, 'nzd'     , 23, 1, nzd     )

    iostatus = sfsnatt( sd_id, 'nxb'     , 23, 1, nxb     )
    iostatus = sfsnatt( sd_id, 'nyb'     , 23, 1, nyb     )
    iostatus = sfsnatt( sd_id, 'nzb'     , 23, 1, nzb     )
    iostatus = sfsnatt( sd_id, 'nb'      , 23, 1, nb      )

    iostatus = sfsnatt( sd_id, 'xmin'    , 6,  1, xmin    )
    iostatus = sfsnatt( sd_id, 'xmax'    , 6,  1, xmax    )
    iostatus = sfsnatt( sd_id, 'ymin'    , 6,  1, ymin    )
    iostatus = sfsnatt( sd_id, 'ymax'    , 6,  1, ymax    )
    iostatus = sfsnatt( sd_id, 'zmin'    , 6,  1, zmin    )
    iostatus = sfsnatt( sd_id, 'zmax'    , 6,  1, zmax    )

    iostatus = sfsnatt( sd_id, 'nstep'   , 23, 1, nstep          )
    iostatus = sfsnatt( sd_id, 'time'    ,  6, 1, t              )
    iostatus = sfsnatt( sd_id, 'timestep',  6, 1, dt             )

!    do ifl=1,nfluid
!    write(gammaifl,'(a5,i1)') 'gamma',ifl
!    iostatus = sfsnatt( sd_id, gammaifl   ,  6, 1, gamma(ifl)     )
!    enddo


! write selected problem dependent parameters
#ifdef SN_SRC
    iostatus = sfsnatt( sd_id, 'nsn'      , 24, 1, nsn     )
#endif /* SN_SRC */

    iv = 1
    iw = 0
    ifl = 1

    do while (len_trim(vars(iv)) .ne. 0)

      select case(vars(iv))
      case ('dens')
        write(varname,'(a3,i1)') 'den',ifl
        wa(iso:ieo,jso:jeo,kso:keo) = u(iarr_all_dn(ifl),iso:ieo,jso:jeo,kso:keo)
        call next_fluid_or_var(ifl,iw,nfluid)

      case ('velx')
        write(varname,'(a3,i1)') 'vlx',ifl
        wa(iso:ieo,jso:jeo,kso:keo) = u(iarr_all_mx(ifl),iso:ieo,jso:jeo,kso:keo) / u(iarr_all_dn(ifl),iso:ieo,jso:jeo,kso:keo)
        call next_fluid_or_var(ifl,iw,nfluid)

      case ('vely')
        write(varname,'(a3,i1)') 'vly',ifl
        wa(iso:ieo,jso:jeo,kso:keo) = u(iarr_all_my(ifl),iso:ieo,jso:jeo,kso:keo) / u(iarr_all_dn(ifl),iso:ieo,jso:jeo,kso:keo)
        call next_fluid_or_var(ifl,iw,nfluid)

      case ('velz')
        write(varname,'(a3,i1)') 'vlz',ifl
        wa(iso:ieo,jso:jeo,kso:keo) = u(iarr_all_mz(ifl),iso:ieo,jso:jeo,kso:keo) / u(iarr_all_dn(ifl),iso:ieo,jso:jeo,kso:keo)
        call next_fluid_or_var(ifl,iw,nfluid)

#ifdef ISO
      case ('ener')
        write(varname,'(a3,i1)') 'ene',ifl
        wa(iso:ieo,jso:jeo,kso:keo) = 0.5*(u(iarr_all_mx(ifl),iso:ieo,jso:jeo,kso:keo)**2 &
                                    +u(iarr_all_my(ifl),iso:ieo,jso:jeo,kso:keo)**2 &
                                    +u(iarr_all_mz(ifl),iso:ieo,jso:jeo,kso:keo)**2)/u(iarr_all_dn(ifl),iso:ieo,jso:jeo,kso:keo)
        call next_fluid_or_var(ifl,iw,nfluid)

#else /* ISO */
      case ('ener')
        write(varname,'(a3,i1)') 'ene',ifl
        wa(iso:ieo,jso:jeo,kso:keo) = u(iarr_all_en(ifl),iso:ieo,jso:jeo,kso:keo)
        call next_fluid_or_var(ifl,iw,nfluid)
#endif /* ISO */

#ifdef COSM_RAYS
      case ('encr')
        varname = 'encr'            
        wa(iso:ieo,jso:jeo,kso:keo) = u(iecr, iso:ieo, jso:jeo, kso:keo)
#endif /* COSM_RAYS */

      case ('divb')
        wa(iso:ieo-1,jso:jeo-1,kso:keo-1) = &
           (b(ibx,iso+1:ieo,jso:jeo-1,kso:keo-1) - b(ibx,iso:ieo-1,jso:jeo-1,kso:keo-1))*dy*dz &
          +(b(iby,iso:ieo-1,jso+1:jeo,kso:keo-1) - b(iby,iso:ieo-1,jso:jeo-1,kso:keo-1))*dx*dz &
          +(b(ibz,iso:ieo-1,jso:jeo-1,kso+1:keo) - b(ibz,iso:ieo-1,jso:jeo-1,kso:keo-1))*dx*dy
        wa = abs(wa)
        wa(ieo,:,:) = wa(ieo-1,:,:)
        wa(:,jeo,:) = wa(:,jeo-1,:)
        wa(:,:,keo) = wa(:,:,keo-1)
#ifdef GRAV
      case ('gpot')
        varname = 'gpot'
        wa(iso:ieo,jso:jeo,kso:keo) = gp(iso:ieo,jso:jeo,kso:keo)
#endif /* GRAV */

      case ('magx')
        varname = 'magx'
        if(domain .eq. 'full_domain') then
          if(mag_center .eq. 'yes') then
            wa(:,:,:) = 0.5*b(ibx,:,:,:)
            wa(:,:,:) = wa(:,:,:)  + cshift(wa(:,:,:),shift=1,dim=1)
          else
            wa(:,:,:) = b(ibx,:,:,:)
          endif
        else if(domain .eq. 'phys_domain') then
          if(mag_center .eq. 'yes') then
            wa(iso:ieo,jso:jeo,kso:keo) = 0.5*(b(ibx,iso:ieo,jso:jeo,kso:keo))
            wa(iso:ieo,jso:jeo,kso:keo) = wa(iso:ieo,jso:jeo,kso:keo) + 0.5*(b(ibx,iso+1:ieo+1,jso:jeo,kso:keo))
          else
            wa(iso:ieo,jso:jeo,kso:keo) = b(ibx,iso:ieo,jso:jeo,kso:keo)
          endif
        endif

      case ('magy')
        varname = 'magy'
        if(domain .eq. 'full_domain') then
          if(mag_center .eq. 'yes') then
            wa(:,:,:) = 0.5*b(iby,:,:,:)
            wa(:,:,:) = wa(:,:,:)  + cshift(wa(:,:,:),shift=1,dim=2)
          else
            wa(:,:,:) = b(iby,:,:,:)
          endif
        else if(domain .eq. 'phys_domain') then
          if(mag_center .eq. 'yes') then
            wa(iso:ieo,jso:jeo,kso:keo) = 0.5*(b(iby,iso:ieo,jso:jeo,kso:keo))
            wa(iso:ieo,jso:jeo,kso:keo) = wa(iso:ieo,jso:jeo,kso:keo) + 0.5*(b(iby,iso:ieo,jso+1:jeo+1,kso:keo))
          else
            wa(iso:ieo,jso:jeo,kso:keo) = b(iby,iso:ieo,jso:jeo,kso:keo)
          endif
        endif

      case ('magz')
        varname = 'magz'
        if(domain .eq. 'full_domain') then
          if(mag_center .eq. 'yes') then
            wa(:,:,:) = 0.5*b(ibz,:,:,:)
            wa(:,:,:) = wa(:,:,:)  + cshift(wa(:,:,:),shift=1,dim=3)
          else
            wa(:,:,:) = b(ibz,:,:,:)
          endif
        else if(domain .eq. 'phys_domain') then
          if(mag_center .eq. 'yes') then
            wa(iso:ieo,jso:jeo,kso:keo) = 0.5*(b(ibz,iso:ieo,jso:jeo,kso:keo))
            wa(iso:ieo,jso:jeo,kso:keo) = wa(iso:ieo,jso:jeo,kso:keo) + 0.5*(b(ibz,iso:ieo,jso:jeo,kso+1:keo+1))
          else
            wa(iso:ieo,jso:jeo,kso:keo) = b(ibz,iso:ieo,jso:jeo,kso:keo)
          endif
        endif

      case default
        print *, 'Variable ', vars(iv), ' is not defined! Skipping.'
      end select

! write data
!
      sds_id = sfcreate(sd_id, varname, 5, rank, dims)
      iostatus = sfscompress(sds_id, comp_type, comp_prm)
      tmparr = real(wa(iso:ieo,jso:jeo,kso:keo),4)
      iostatus = sfwdata(sds_id, istart, stride, dims, tmparr)

! write coords
!
      dim_id = sfdimid( sds_id, 0 )
      iostatus = sfsdmname( dim_id, 'xc' )
      iostatus = sfsdscale( dim_id, dims(1), 5, real(x(iso:ieo),4))
      iostatus = sfsdmstr ( dim_id, 'X', 'pc', '' )

      dim_id = sfdimid( sds_id, 1 )
      iostatus = sfsdmname( dim_id, 'yc' )
      iostatus = sfsdscale( dim_id, dims(2), 5, real(y(jso:jeo),4))
      iostatus = sfsdmstr ( dim_id, 'Y', 'pc', '' )

      dim_id = sfdimid( sds_id, 2 )
      iostatus = sfsdmname( dim_id, 'zc' )
      iostatus = sfsdscale( dim_id, dims(3), 5, real(z(kso:keo),4))
      iostatus = sfsdmstr ( dim_id, 'Z', 'pc', '' )

      iostatus = sfendacc(sds_id)

      if(iw .eq. 0) iv = iv + 1
    end do

    iostatus = sfend(sd_id)

    if(proc.eq. 0) then
      open(log_lun, file=log_file, position='append')

      write(log_lun,*) 'Writing output   file: ', trim(file_name_disp)
      write(*,*)       'Writing output   file: ', trim(file_name_disp)

      close(log_lun)
    endif
    
    deallocate(tmparr)
    
  end subroutine write_hdf

  subroutine next_fluid_or_var(ifluid,ivar,nfluids)
    implicit none
    integer ifluid,ivar,nfluids
    if(ifluid .lt. nfluids) then
      ifluid=ifluid+1
      ivar=1
    else
      ivar=0
      ifluid=1
    endif
    return
  end subroutine next_fluid_or_var

!---------------------------------------------------------------------
!
! writes restart file
!
!---------------------------------------------------------------------
!
  subroutine write_restart

    use fluidindex,   only : nvar,nmag
    use arrays,       only : u,b
    use grid,         only : nxb,nyb,nzb,x,y,z,nx,ny,nz
    use grid,         only : xmin,xmax,ymin,ymax,zmin,zmax,nxd,nyd,nzd,nb
    use start,        only : t,dt,nstep
    use init_problem, only : problem_name, run_id
#ifdef GRAV
    use arrays, only : gp
#endif /* GRAV */

    implicit none

    character(len=128) :: file_name_res,file_name_disp,file_name_last
    character*160 syscom
    logical lastres_exist

    integer :: sd_id, sds_id, dim_id, scstatus
    integer :: iostatus, ranku, rankb, rank3d, comp_type
    integer, dimension(4) :: dimsu, dimsb, istart, stride
    integer, dimension(3) :: dims, dims3d
    integer, dimension(1) :: comp_prm
    integer :: sfstart, sfend, sfsnatt, sfcreate, sfwdata, sfscompress, sfendacc &
             , sfdimid, sfsdmname, sfsdscale

    integer(kind=1) :: system

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  prepare data dimensions
!
    comp_type   = 4
    comp_prm(1) = 6
    istart(:)   = 0
    stride(:)   = 1

    dims(1)     = nx
    dims(2)     = ny
    dims(3)     = nz


    ranku       = 4
    dimsu(1)    = nvar
    dimsu(2)    = dims(1)
    dimsu(3)    = dims(2)
    dimsu(4)    = dims(3)

    rankb       = 4
    dimsb(1)    = nmag
    dimsb(2)    = dims(1)
    dimsb(3)    = dims(2)
    dimsb(4)    = dims(3)

    rank3d      = 3
    dims3d(1)    = dims(1)
    dims3d(2)    = dims(2)
    dims3d(3)    = dims(3)

!  generate filename
!
    write (file_name_res,'(a,a1,a,a1,a3,a1,3(i2.2,a1),i3.3,a4)') &
              trim(cwd),'/',trim(problem_name),'_', run_id,  &
              '_',pcoords(1),'_',pcoords(2),'_',pcoords(3),'_',nres, '.res'

    write (file_name_disp,'(a,a1,a3,a1,3(a2,a1),i3.3,a4)') &
              trim(problem_name),'_', run_id, &
              '_',pc1,'_',pc2,'_',pc3,'_',nres, '.res'

    if((resdel .ne. 0) .and. (nres .gt. resdel)) then
    write (file_name_last,'(a,a1,a,a1,a3,a1,3(i2.2,a1),i3.3,a4)') &
              trim(cwd),'/',trim(problem_name),'_', run_id,  &
              '_',pcoords(1),'_',pcoords(2),'_',pcoords(3),'_',nres-resdel, '.res'
    endif

    sd_id = sfstart(file_name_res, 4)

!!  ATTRIBUTES
!!
! write config attributes
!
    iostatus = sfsnatt( sd_id, 'problem' , 4, 32, problem_name)
    iostatus = sfsnatt( sd_id, 'run_id'  , 4, 32, run_id      )
    iostatus = sfsnatt( sd_id, 'domain'  , 4, 32, domain      )

    iostatus = sfsnatt( sd_id, 'psize1'  , 23, 1, psize(1) )
    iostatus = sfsnatt( sd_id, 'psize2'  , 23, 1, psize(2) )
    iostatus = sfsnatt( sd_id, 'psize3'  , 23, 1, psize(3) )

    iostatus = sfsnatt( sd_id, 'pcoords1', 23, 1, pcoords(1) )
    iostatus = sfsnatt( sd_id, 'pcoords2', 23, 1, pcoords(2) )
    iostatus = sfsnatt( sd_id, 'pcoords3', 23, 1, pcoords(3) )

    iostatus = sfsnatt( sd_id, 'dimsu'   , 23, 1, nvar )
    iostatus = sfsnatt( sd_id, 'dims1'   , 23, 1, dims(1) )
    iostatus = sfsnatt( sd_id, 'dims2'   , 23, 1, dims(2) )
    iostatus = sfsnatt( sd_id, 'dims3'   , 23, 1, dims(3) )

    iostatus = sfsnatt( sd_id, 'nxd'     , 23,  1, nxd     )
    iostatus = sfsnatt( sd_id, 'nyd'     , 23,  1, nyd     )
    iostatus = sfsnatt( sd_id, 'nzd'     , 23,  1, nzd     )

    iostatus = sfsnatt( sd_id, 'nxb'     , 23,  1, nxb     )
    iostatus = sfsnatt( sd_id, 'nyb'     , 23,  1, nyb     )
    iostatus = sfsnatt( sd_id, 'nzb'     , 23,  1, nzb     )
    iostatus = sfsnatt( sd_id, 'nb'      , 23,  1, nb      )

    iostatus = sfsnatt( sd_id, 'xmin'    ,  6,  1, xmin   )
    iostatus = sfsnatt( sd_id, 'xmax'    ,  6,  1, xmax   )
    iostatus = sfsnatt( sd_id, 'ymin'    ,  6,  1, ymin   )
    iostatus = sfsnatt( sd_id, 'ymax'    ,  6,  1, ymax   )
    iostatus = sfsnatt( sd_id, 'zmin'    ,  6,  1, zmin   )
    iostatus = sfsnatt( sd_id, 'zmax'    ,  6,  1, zmax   )

! write evolution attributes
!
    iostatus = sfsnatt( sd_id, 'nstep'   , 24,  1, nstep   )
    iostatus = sfsnatt( sd_id, 'time'    ,  6,  1, t       )
    iostatus = sfsnatt( sd_id, 'timestep',  6,  1, dt      )

! write dataio attributes
!
    iostatus = sfsnatt( sd_id, 'nres'     , 24, 1, nres + 1 )
    iostatus = sfsnatt( sd_id, 'nhdf'     , 24, 1, nhdf     )
    iostatus = sfsnatt( sd_id, 'ntsl'     , 24, 1, ntsl     )
    iostatus = sfsnatt( sd_id, 'nlog'     , 24, 1, nlog     )
    iostatus = sfsnatt( sd_id, 'step_res' , 24, 1, nstep     )
    iostatus = sfsnatt( sd_id, 'step_hdf' , 24, 1, step_hdf )
    iostatus = sfsnatt( sd_id, 'last_hdf_time', 6, 1, last_hdf_time)

! write selected problem dependent parameters
#ifdef SN_SRC
    iostatus = sfsnatt( sd_id, 'nsn'      , 24, 1, nsn     )
#endif /* SN_SRC */

! write array of integer scalars
!
!    sds_id   = sfcreate(sd_id, 'intscal', 23, 1, nintscal)
!    iostatus = sfscompress(sds_id, comp_type, comp_prm)
!    iostatus = sfwdata(sds_id, istart, stride, nintscal, intscal)

! write array of real scalars
!
!    sds_id   = sfcreate(sd_id, 'rlscal', 6, 1, nrlscal)
!    iostatus = sfscompress(sds_id, comp_type, comp_prm)
!    iostatus = sfwdata(sds_id, istart, stride, nrlscal, rlscal)

! write initial vertical density profile
!
!    sds_id   = sfcreate(sd_id, 'dprof', 6, 1, nz)
!    iostatus = sfscompress(sds_id, comp_type, comp_prm)
!    iostatus = sfwdata(sds_id, istart, stride, nz, dprof)

! write fluid variables array
!
    sds_id = sfcreate(sd_id, 'fluid_vars', 6, ranku, dimsu)
    iostatus = sfscompress(sds_id, comp_type, comp_prm)
    iostatus = sfwdata(sds_id, istart, stride, dimsu, u)

! write magnetic field array
!
    sds_id = sfcreate(sd_id, 'mag_field', 6, rankb, dimsb)
    iostatus = sfscompress(sds_id, comp_type, comp_prm)
    iostatus = sfwdata(sds_id, istart, stride, dimsb, b)

#ifdef GRAV
! write gravitational potential
!
    sds_id = sfcreate(sd_id, 'grav_pot', 6, rank3d, dims3d)
    iostatus = sfscompress(sds_id, comp_type, comp_prm)
    iostatus = sfwdata(sds_id, istart, stride, dims3d, gp)
#endif /* GRAV */

! write coords
!
    dim_id = sfdimid( sds_id, 1 )
    iostatus = sfsdmname( dim_id, 'x' )
    iostatus = sfsdscale( dim_id, dimsu(2), 6, x)

    dim_id = sfdimid( sds_id, 2 )
    iostatus = sfsdmname( dim_id, 'y' )
    iostatus = sfsdscale( dim_id, dimsu(3), 6, y)

    dim_id = sfdimid( sds_id, 3 )
    iostatus = sfsdmname( dim_id, 'z' )
    iostatus = sfsdscale( dim_id, dimsu(4), 6, z)

    iostatus = sfendacc(sds_id)

    iostatus = sfend(sd_id)


    if(proc.eq. 0) then
      open(log_lun, file=log_file, position='append')

      write(log_lun,*) 'Writing restart  file: ', trim(file_name_disp)
      write(*,*)       'Writing restart  file: ', trim(file_name_disp)

      close(log_lun)
    endif

    if((resdel .ne. 0) .and. (nres .gt. resdel)) then
        inquire(file=file_name_last, exist=lastres_exist)
        if(lastres_exist) then
          syscom='rm -f'//file_name_last
          scstatus = system(syscom)
          write(*,*) trim(file_name_last),' removed'
        else
          write(*,*) trim(file_name_last),' does not exist'
        endif
    endif

  end subroutine write_restart

!---------------------------------------------------------------------
!
! read restart file
!
!---------------------------------------------------------------------
!
  subroutine read_restart !(all)
    use fluidindex, only : nvar, nmag
    use arrays, only : u,b
    use grid, only : nx,ny,nz
    use start, only  : t, dt, nstep
    use init_problem, only : problem_name, run_id
#ifdef GRAV
    use arrays, only : gp
#endif /* GRAV */

    implicit none

    character(len=128) :: file_name_res,file_name_disp
    logical file_exist, log_exist

    integer :: sd_id, sds_id, attr_index, sds_index
    integer :: iostatus, ranku, rankb, rank3d
    integer, dimension(4) :: dimsu, dimsb, istart, stride
    integer, dimension(3) :: dims3d
    integer :: sfstart, sfend, sffattr, sfrnatt, sfn2index, sfselect, sfrdata, sfendacc

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  prepare data dimensions
!
    istart(:) = 0
    stride(:) = 1

    ranku    = 4
    dimsu(1) = nvar
    dimsu(2) = nx
    dimsu(3) = ny
    dimsu(4) = nz

    rankb    = 4
    dimsb(1) = nmag
    dimsb(2) = nx
    dimsb(3) = ny
    dimsb(4) = nz

    rank3d   = 3
    dims3d(1) = nx
    dims3d(2) = ny
    dims3d(3) = nz


!  generate filename
!
    write (file_name_res,'(a,a1,a,a1,a3,a1,3(i2.2,a1),i3.3,a4)') &
              trim(cwd),'/',trim(problem_name),'_', run_id,  &
              '_',pcoords(1),'_',pcoords(2),'_',pcoords(3),'_',nres, '.res'

    write (file_name_disp,'(a,a1,a3,a1,3(a2,a1),i3.3,a4)') &
              trim(problem_name),'_', run_id, &
              '_',pc1,'_',pc2,'_',pc3,'_',nres, '.res'

    if(proc.eq. 0) then
        write(*,*)       'Reading restart  file: ', trim(file_name_disp)
    endif
  if(proc==0) then
    inquire(file =log_file , exist = log_exist)
    if(file_exist .eqv. .true.) then
      open(log_lun, file=log_file, position='append')
        write(log_lun,*) 'Reading restart  file: ', trim(file_name_disp)
      close(log_lun)
    endif
  endif

    inquire(file = file_name_res, exist = file_exist)
    if(file_exist .eqv. .false.) then
      if(log_exist) then
        open(log_lun, file=log_file, position='append')
          write(log_lun,*) 'Restart  file: ', trim(file_name_res), &
                         '               does not exist.  ABORTING !!! '
        close(log_lun)
      endif

      write(*,*)       'Restart  file: ', trim(file_name_res), &
                         '               does not exist.  ABORTING !!! '
      call MPI_BARRIER(comm,ierr)
      call mpistop
      stop
    endif


    sd_id = sfstart(file_name_res, 1)

! read evolution attributes
!
      attr_index = sffattr( sd_id, 'nstep'     )
      iostatus = sfrnatt( sd_id, attr_index, nstep   )
      attr_index = sffattr( sd_id, 'time'     )
      iostatus = sfrnatt( sd_id, attr_index, t   )
      attr_index = sffattr( sd_id, 'timestep' )
      iostatus = sfrnatt( sd_id, attr_index, dt     )

! read dataio attributes
!
      attr_index = sffattr( sd_id, 'nres'      )
      iostatus = sfrnatt( sd_id, attr_index, nres     )
      attr_index = sffattr( sd_id, 'nhdf'      )
      iostatus = sfrnatt( sd_id, attr_index, nhdf     )
      attr_index = sffattr( sd_id, 'ntsl'      )
      iostatus = sfrnatt( sd_id, attr_index, ntsl     )
      attr_index = sffattr( sd_id, 'nlog'      )
      iostatus = sfrnatt( sd_id, attr_index, nlog     )
      attr_index = sffattr( sd_id, 'step_res'  )
      iostatus = sfrnatt( sd_id, attr_index, step_res )
      attr_index = sffattr( sd_id, 'step_hdf'  )
      iostatus = sfrnatt( sd_id, attr_index, step_hdf )
      attr_index = sffattr( sd_id, 'last_hdf_time' )
      iostatus = sfrnatt( sd_id, attr_index, last_hdf_time )

! read selected problem dependent parameters
#ifdef SN_SRC
      attr_index = sffattr( sd_id, 'nsn'       )
      iostatus = sfrnatt( sd_id, attr_index, nsn_last    )
#endif /* SN_SRC */

! read variables array
!
      sds_index = sfn2index(sd_id, 'fluid_vars')
      sds_id = sfselect(sd_id, sds_index)
      iostatus = sfrdata(sds_id, istart, stride, dimsu, u)

! read magnetic field
!
      sds_index = sfn2index(sd_id, 'mag_field')
      sds_id = sfselect(sd_id, sds_index)
      iostatus = sfrdata(sds_id, istart, stride, dimsb, b)

#ifdef GRAV
! read gravitational potential
!
      sds_index = sfn2index(sd_id, 'grav_pot')
      sds_id = sfselect(sd_id, sds_index)
      iostatus = sfrdata(sds_id, istart, stride, dims3d, gp)
#endif /* GRAV */

      iostatus = sfendacc(sds_id)

    iostatus = sfend(sd_id)

  end subroutine read_restart

!------------------------------------------------------------------------

    subroutine find_last_restart(restart_number)

      use init_problem, only : problem_name, run_id

      implicit none

      character*120 file_name
      integer restart_number,nres
      logical exist
      character(len=128) :: file_name_base

      restart_number = 0

      write (file_name_base,'(a,a1,a3,a1)') trim(problem_name),'_',run_id,'_'

      call rm_file('restart_list.tmp')

      do nres =999,0,-1
        write (file_name,'(a,a1,a,a1,a3,a1,3(i2.2,a1),i3.3,a4)') &
               trim(cwd),'/',trim(problem_name),'_', run_id,'_',0,'_',0,'_',0,'_',nres,'.res'
        inquire(file = file_name, exist = exist)
        if(exist) then
           restart_number = nres
        return
        endif
      enddo

    end subroutine find_last_restart


!---------------------------------------------------------------------
!
! writes integrals to text file
!
!---------------------------------------------------------------------
!
  subroutine write_timeslice

    use fluidindex,   only : nfluid
    use fluidindex,   only : ibx,iby,ibz
    use fluidindex, only : nvar, iarr_all_dn,iarr_all_mx,iarr_all_my,iarr_all_mz
    use grid, only  : dvol,dx,dy,dz,is,ie,js,je,ks,ke,x,y,z,nxd,nyd,nzd
    use start, only : proc, dt, t, nstep
    use arrays, only : u,b,wa
    use init_problem, only : problem_name, run_id

#ifdef IONIZED
    use initionized, only : gamma_ion, cs_iso_ion,cs_iso_ion2
#endif /* IONIZED */

#ifdef NEUTRAL
    use initneutral, only : gamma_neu, cs_iso_neu,cs_iso_neu2
#endif /* NEUTRAL */

#ifndef ISO
    use fluidindex, only : iarr_all_en
#endif /* ISO */

#ifdef COSM_RAYS
    use initcosmicrays, only : iecr
#endif /* COSM_RAYS */

#ifndef STANDARD
    use constants, only : Gs
#endif /* STANDARD */

#ifdef GRAV
    use arrays, only : gp
#endif /* GRAV */

#ifdef ISO
    use start, only : csi2
#endif /* ISO */

#ifdef RESISTIVE
    use resistivity
#endif /* RESISTIVE */

#ifdef COSM_RAYS
    use initcosmicrays, only : iecr
#endif /* COSM_RAYS */

#ifdef SNE_DISTR
    use sn_distr, only : emagadd, tot_emagadd
#endif /* SNE_DISTR */

    implicit none

    character(len=128) :: tsl_file

    real :: mass = 0.0, momx = 0.0, momy = 0.0,  momz = 0.0, &
#ifndef ISO
            ener = 0.0, &
#endif /* ISO */
            ekin = 0.0,  emag = 0.0, &
            tot_mass = 0.0, tot_momx = 0.0, tot_momy = 0.0, tot_momz = 0.0, &
            tot_ener = 0.0, tot_eint = 0.0, tot_ekin = 0.0, tot_emag = 0.0, &
            tot_epot = 0.0, mflx = 0.0, mfly = 0.0, mflz = 0.0, &
            tot_mflx = 0.0, tot_mfly = 0.0, tot_mflz = 0.0
#ifdef GRAV
    integer :: i,j
    real :: epot =0.0, amomz = 0.0, tot_amomz = 0.0
#endif /* GRAV */
#ifdef COSM_RAYS
    real :: encr = 0.0, tot_encr = 0.0
#endif /* COSM_RAYS */
#ifdef SNE_DISTR
    real :: sum_emagadd = 0.0
#endif /* SNE_DISTR */
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
    if (proc .eq. 0) then
      write (tsl_file,'(a,a1,a,a1,a3,a1,i3.3,a4)') &
              trim(cwd),'/',trim(problem_name),'_', run_id,'_',nrestart,'.tsl'

      if (tsl_firstcall) then
        open(tsl_lun, file=tsl_file)

        write (tsl_lun, '(a1,a8,50a16)') '#','nstep', 'time', 'timestep', 'mass', &
                                           'momx', 'momy', 'momz', 'amomz', &
                                           'ener', 'epot', 'eint', 'ekin', &
#ifdef MAGNETIC					   
                                           'emag', 'mflx', 'mfly', 'mflz', &
#endif /* MAGNETIC */					   
#ifdef COSM_RAYS
                                           'encr_tot', 'encr_min',  'encr_max',&
#endif /* COSM_RAYS */
#ifdef RESISTIVE
                                           'eta_max', &
#endif /* RESISTIVE */
! some quantities computed in "write_log".One can add more, or change.
                                           'vx_max', 'vy_max', 'vz_max', 'va_max', 'cs_max', &
                                           'dens_min', 'dens_max', 'pres_min', 'pres_max', &
#ifndef ISO
                                           'temp_min', 'temp_max',  &
#endif /* ISO */
#ifdef SNE_DISTR
                                           'sum_emagadd', 'tot_emagadd', &
#endif /* SNE_DISTR */
#ifdef MAGNETIC                            'b_min', 'b_max', 'divb_max', &

#endif /* MAGNETIC */
                                           '     '

        write (tsl_lun, '(a1)') '#'
        tsl_firstcall = .false.
      else
        open(tsl_lun, file=tsl_file, position='append')
      endif
    endif

    mass = sum(u(iarr_all_dn,is:ie,js:je,ks:ke)) * dvol
    call mpi_allreduce(mass, tot_mass, 1, mpi_real8, mpi_sum, comm3d, ierr)

    momx = sum(u(iarr_all_mx,is:ie,js:je,ks:ke)) * dvol
    call mpi_allreduce(momx, tot_momx, 1, mpi_real8, mpi_sum, comm3d, ierr)

    momy = sum(u(iarr_all_my,is:ie,js:je,ks:ke)) * dvol
    call mpi_allreduce(momy, tot_momy, 1, mpi_real8, mpi_sum, comm3d, ierr)

    momz = sum(u(iarr_all_mz,is:ie,js:je,ks:ke)) * dvol
    call mpi_allreduce(momz, tot_momz, 1, mpi_real8, mpi_sum, comm3d, ierr)

#ifdef GRAV
!DW+
    amomz = 0.0
    do j=js,je
      do i=is,ie
         amomz = amomz + (x(i)*sum(u(iarr_all_my,i,j,ks:ke)) &
                         -y(j)*sum(u(iarr_all_mx,i,j,ks:ke)))*dvol
      enddo
    enddo
    call mpi_allreduce(amomz, tot_amomz, 1, mpi_real8, mpi_sum, comm3d, ierr)
!DW-
    epot = sum(u(iarr_all_dn(1),is:ie,js:je,ks:ke) *gp(is:ie,js:je,ks:ke)) * dvol
    call mpi_allreduce(epot, tot_epot, 1, mpi_real8, mpi_sum, comm3d, ierr)
#endif /* GRAV */

    wa(is:ie,js:je,ks:ke) &
        = 0.5 * (u(iarr_all_mx(1),is:ie,js:je,ks:ke)**2   &
               + u(iarr_all_my(1),is:ie,js:je,ks:ke)**2  &
               + u(iarr_all_mz(1),is:ie,js:je,ks:ke)**2)/ &
                 u(iarr_all_dn(1),is:ie,js:je,ks:ke)
    ekin = sum(wa(is:ie,js:je,ks:ke)) * dvol
    call mpi_allreduce(ekin, tot_ekin, 1, mpi_real8, mpi_sum, comm3d, ierr)

    wa(is:ie,js:je,ks:ke) &
       = 0.5 * (b(ibx,is:ie,js:je,ks:ke)**2 + &
                b(iby,is:ie,js:je,ks:ke)**2 + &
                b(ibz,is:ie,js:je,ks:ke)**2)
    emag = sum(wa(is:ie,js:je,ks:ke)) * dvol
    call mpi_allreduce(emag, tot_emag, 1, mpi_real8, mpi_sum, comm3d, ierr)

    wa(is:ie,js:je,ks:ke) = b(ibx,is:ie,js:je,ks:ke)
    mflx = sum(wa(is:ie,js:je,ks:ke)) * dy*dz/nxd
    call mpi_allreduce(mflx, tot_mflx, 1, mpi_real8, mpi_sum, comm3d, ierr)

    wa(is:ie,js:je,ks:ke) = b(iby,is:ie,js:je,ks:ke)
    mfly = sum(wa(is:ie,js:je,ks:ke)) * dx*dz/nyd
    call mpi_allreduce(mfly, tot_mfly, 1, mpi_real8, mpi_sum, comm3d, ierr)

    wa(is:ie,js:je,ks:ke) = b(ibz,is:ie,js:je,ks:ke)
    mflz = sum(wa(is:ie,js:je,ks:ke)) * dx*dy/nzd
    call mpi_allreduce(mflz, tot_mflz, 1, mpi_real8, mpi_sum, comm3d, ierr)



#ifdef ISO
    tot_eint = csi2*tot_mass
    tot_ener = tot_eint+tot_ekin+tot_emag
#else /* ISO */
    ener = sum(u(iarr_all_en,is:ie,js:je,ks:ke)) * dvol
    call mpi_allreduce(ener, tot_ener, 1, mpi_real8, mpi_sum, comm3d, ierr)
    tot_eint = tot_ener - tot_ekin - tot_emag
#endif /* ISO */
#ifdef GRAV
    tot_ener = tot_ener + tot_epot
#endif /* GRAV */

#ifdef COSM_RAYS
    encr = sum(u(iecr,is:ie,js:je,ks:ke)) * dvol
    call mpi_allreduce(encr, tot_encr, 1, mpi_real8, mpi_sum, comm3d, ierr)
#endif /* COSM_RAYS */
#ifdef SNE_DISTR
    call mpi_allreduce(emagadd, sum_emagadd, 1, mpi_real8, mpi_sum, comm3d, ierr)
    tot_emagadd = tot_emagadd + sum_emagadd
#endif /* SNE_DISTR */

    if (proc .eq. 0) then
      write (tsl_lun, '(1x,i8,50(1x,1pe15.8))') &
                      nstep, &
                      t, dt, tot_mass, &
                      tot_momx, tot_momy, tot_momz, &
#ifdef GRAV
                      tot_amomz, &
#endif /* GRAV */
                      tot_ener, tot_epot, tot_eint, tot_ekin, tot_emag, &
                      tot_mflx, tot_mfly, tot_mflz, &
#ifdef COSM_RAYS
                      tot_encr, encr_min, encr_max, &
#endif /* COSM_RAYS */
#ifdef MAGNETIC
                      b_min, b_max, divb_max &

#ifdef RESISTIVE
                      eta_max, &
#endif /* RESISTIVE */
#endif /* MAGNETIC */

! some quantities computed in "write_log".One can add more, or change.
                      vx_max, vy_max, vz_max, va_max, cs_max, &
                      dens_min, dens_max, pres_min, pres_max, &
#ifndef ISO
                      temp_min, temp_max, &
#endif /* ISO */
                      0.0
		                           
      close(tsl_lun)
    endif


  end subroutine write_timeslice

!---------------------------------------------------------------------
!
! writes timestep diagnostics to the logfile
!
!---------------------------------------------------------------------
!
  subroutine  write_log

    use fluidindex, only : ibx, iby, ibz, nfluid
    use arrays, only : wa,u,b
    use grid, only   : dx,dy,dz,dxmn,nb,is,ie,js,je,ks,ke,nx,ny,nz
    use constants, only : small, hydro_mass, k_B
    use start, only : t,dt,nstep,smallei,cfl

#ifdef IONIZED
    use initionized, only : gamma_ion, cs_iso_ion,cs_iso_ion2
    use initionized, only : idni,imxi,imyi,imzi
#ifndef ISO
    use initionized, only : ieni
#endif /* ISO */
#endif /* IONIZED */

#ifdef NEUTRAL
    use initneutral, only : gamma_neu, cs_iso_neu,cs_iso_neu2
    use initneutral, only : idnn,imxn,imyn,imzn
#ifndef ISO
    use initneutral, only : ienn
#endif /* ISO */
#endif /* NEUTRAL */

#ifdef DUST
    use initdust, only : idnd,imxd,imyd,imzd
#endif /* DUST */

#ifdef COSM_RAYS
    use timestepcosmicrays, only  : dt_crs
    use initcosmicrays, only : iecr
#endif /* COSM_RAYS */

#ifdef RESISTIVE
    use resistivity
#endif /* RESISTIVE */

    implicit none

#ifdef MAGNETIC
    integer, dimension(3) :: loc_b_min, loc_b_max, loc_divb_max
    integer               :: proc_b_min, proc_b_max, proc_divb_max    
    real                  :: b_min, b_max    
#endif /* MAGNETIC */ 

#ifdef IONIZED 
    integer, dimension(3) :: loc_vxi_max, loc_vyi_max, loc_vzi_max, &
                             loc_csi_max, loc_deni_min, loc_deni_max, &
                             loc_prei_min, loc_prei_max, loc_temi_min, loc_temi_max, &
                             loc_vai_max 
    integer               :: proc_vxi_max, proc_vyi_max, proc_vzi_max, &
                             proc_csi_max, proc_deni_min, proc_deni_max, &
                             proc_prei_min, proc_prei_max, proc_temi_min, proc_temi_max, &
                             proc_vai_max, 
    real                  :: deni_min, deni_max, vxi_max, vyi_max, vzi_max, &
                             prei_min, prei_max, temi_min, temi_max, vai_max, csi_max                            
#endif /* IONIZED */

#ifdef NEUTRAL 
    integer, dimension(3) :: loc_vxn_max, loc_vyn_max, loc_vzn_max, &
                             loc_csn_max, loc_denn_min, loc_denn_max, &
                             loc_pren_min, loc_pren_max, loc_temn_min, loc_temn_max
    integer               :: proc_vxn_max, proc_vyn_max, proc_vzn_max, &
                             proc_csn_max, proc_denn_min, proc_denn_max, &
                             proc_pren_min, proc_pren_max, proc_temn_min, proc_temn_max
    real                  :: denn_min, denn_max, vxn_max, vyn_max, vzn_max, &
                             pren_min, pren_max, temn_min, temn_max, csn_max
#endif /* NEUTRAL */

#ifdef DUST
    integer, dimension(3) :: loc_dend_min,  loc_dend_max, &
                             loc_vxd_max, loc_vyd_max, loc_vzd_max

    integer               :: proc_dend_min,  proc_dend_max,  &
                             proc_vxd_max,   proc_vyd_max,   proc_vzd_max
    real                  :: dend_min, dend_max, vxd_max, vyd_max, vzd_max
#endif /* DUST */

#ifdef COSM_RAYS
    integer, dimension(3) :: loc_encr_min, loc_encr_max
    integer               :: proc_encr_min, proc_encr_max
#endif /* COSM_RAYS */

#ifdef RESISTIVE
    integer               :: proc_eta_max
#endif /* RESISTIVE */

! Timestep diagnostics
#ifdef NEUTRAL
    wa            = u(idnn,:,:,:)
    denn_min      = minval(wa(is:ie,js:je,ks:ke))
    loc_denn_min  = minloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
    call mpifind(denn_min, 'min', loc_denn_min, proc_denn_min)

    denn_max      = maxval(wa(is:ie,js:je,ks:ke))
    loc_denn_max  = maxloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
    call mpifind(denn_max, 'max', loc_denn_max, proc_denn_max)

    wa          = abs(u(imxn,:,:,:)/u(idnn,:,:,:))
    vxn_max     = maxval(wa(is:ie,js:je,ks:ke))
    loc_vxn_max = maxloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
    call mpifind(vxn_max, 'max', loc_vxn_max, proc_vxn_max)

    wa          = abs(u(imyn,:,:,:)/u(idnn,:,:,:))
    vyn_max     = maxval(wa(is:ie,js:je,ks:ke))
    loc_vyn_max = maxloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
    call mpifind(vyn_max, 'max', loc_vyn_max, proc_vyn_max)

    wa           = abs(u(imzn,:,:,:)/u(idnn,:,:,:))
    vzn_max      = maxval(wa(is:ie,js:je,ks:ke))
    loc_vzn_max  = maxloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
    call mpifind(vzn_max, 'max', loc_vzn_max, proc_vzn_max)
#ifdef ISO
    pren_min      = cs_iso_neu2*denn_min
    loc_pren_min  = loc_denn_min
    proc_pren_min = proc_denn_min
    pren_max      = cs_iso_neu2*denn_max
    loc_pren_max  = loc_denn_max
    proc_pren_max = proc_denn_max
    csn_max       = cs_iso_neu
    loc_csn_max   = 0
    proc_csn_max  = 0
    temn_min      = hydro_mass / k_B * cs_iso_neu2
    loc_temn_min  = 0
    proc_temn_min = 0
    temn_max      = hydro_mass / k_B * cs_iso_neu2
    loc_temn_max  = 0
    proc_temn_max = 0
#else /* ISO */
    wa(:,:,:) = (u(ienn,:,:,:) &                ! eint
                - 0.5*((u(imxn,:,:,:)**2 +u(imyn,:,:,:)**2 &
                  + u(imzn,:,:,:)**2)/u(idnn,:,:,:)))
    wa(:,:,:) = max(wa(:,:,:),smallei)
    wa(:,:,:) = (gamma_neu-1.0)*wa(:,:,:)           ! pres

    pren_min      = minval(wa(is:ie,js:je,ks:ke))
    loc_pren_min  = minloc(wa(is:ie,js:je,ks:ke)) + (/nb,nb,nb/)
    call mpifind(pren_min, 'min', loc_pren_min, proc_pren_min)

    pren_max      = maxval(wa(is:ie,js:je,ks:ke))
    loc_pren_max  = maxloc(wa(is:ie,js:je,ks:ke)) + (/nb,nb,nb/)
    call mpifind(pren_max, 'max', loc_pren_max, proc_pren_max)

    temn_max      = maxval( hydro_mass / k_B * wa(is:ie,js:je,ks:ke) &
                                             /u(idnn,is:ie,js:je,ks:ke))
    loc_temn_max  = maxloc(wa(is:ie,js:je,ks:ke)    &
                         /u(idnn,is:ie,js:je,ks:ke)  ) + (/nb,nb,nb/)
    call mpifind(temn_max, 'max', loc_temn_max, proc_temn_max)

    temn_min      = minval( hydro_mass / k_B * wa(is:ie,js:je,ks:ke) &
                                             /u(idnn,is:ie,js:je,ks:ke))
    loc_temn_min  = minloc(wa(is:ie,js:je,ks:ke) &
                         /u(idnn,is:ie,js:je,ks:ke)  ) &
                     + (/nb,nb,nb/)
    call mpifind(temn_min, 'min', loc_temn_min, proc_temn_min)

    wa(:,:,:) = gamma_neu*wa(:,:,:)
    csn_max        = sqrt(maxval(wa(is:ie,js:je,ks:ke) &
                            /u(idnn,is:ie,js:je,ks:ke)))
    loc_csn_max    = maxloc(wa(is:ie,js:je,ks:ke) &
                            /u(idnn,is:ie,js:je,ks:ke)) &
                     + (/nb,nb,nb/)
    call mpifind(csn_max, 'max', loc_csn_max, proc_csn_max)
#endif /* ISO */

#endif /* NEUTRAL */


#ifdef IONIZED
    wa            = u(idni,:,:,:)
    deni_min      = minval(wa(is:ie,js:je,ks:ke))
    loc_deni_min  = minloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
    call mpifind(deni_min, 'min', loc_deni_min, proc_deni_min)

    deni_max      = maxval(wa(is:ie,js:je,ks:ke))
    loc_deni_max  = maxloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
    call mpifind(deni_max, 'max', loc_deni_max, proc_deni_max)

    wa          = abs(u(imxi,:,:,:)/u(idni,:,:,:))
    vxi_max     = maxval(wa(is:ie,js:je,ks:ke))
    loc_vxi_max = maxloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
    call mpifind(vxi_max, 'max', loc_vxi_max, proc_vxi_max)

    wa          = abs(u(imyi,:,:,:)/u(idni,:,:,:))
    vyi_max     = maxval(wa(is:ie,js:je,ks:ke))
    loc_vyi_max = maxloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
    call mpifind(vyi_max, 'max', loc_vyi_max, proc_vyi_max)

    wa           = abs(u(imzi,:,:,:)/u(idni,:,:,:))
    vzi_max      = maxval(wa(is:ie,js:je,ks:ke))
    loc_vzi_max  = maxloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
    call mpifind(vzi_max, 'max', loc_vzi_max, proc_vzi_max)

#ifdef MAGNETIC
    wa(:,:,:)  = b(1,:,:,:)*b(1,:,:,:) + b(2,:,:,:)*b(2,:,:,:) + &
                 b(3,:,:,:)*b(3,:,:,:)
    b_min      = sqrt(minval(wa(is:ie,js:je,ks:ke)))
    loc_b_min  = minloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
    call mpifind(b_min, 'min', loc_b_min, proc_b_min)

    b_max      = sqrt(maxval(wa(is:ie,js:je,ks:ke)))
    loc_b_max  = maxloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
    call mpifind(b_max, 'max', loc_b_max, proc_b_max)

    vai_max     = sqrt(maxval(wa(is:ie,js:je,ks:ke) &
                       /u(idni,is:ie,js:je,ks:ke)))
    loc_vai_max = maxloc(wa(is:ie,js:je,ks:ke)     &
                       /u(idni,is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
    call mpifind(va_max, 'max', loc_vai_max, proc_vai_max)
#endif /* MAGNETIC */

#ifdef ISO
    prei_min      = cs_iso_ion2*deni_min
    loc_prei_min  = loc_deni_min
    proc_prei_min = proc_deni_min
    prei_max      = cs_iso_ion2*deni_max
    loc_prei_max  = loc_deni_max
    proc_prei_max = proc_deni_max
    csi_max       = cs_iso_ion
    loc_csi_max   = 0
    proc_csi_max  = 0
    temi_min      = hydro_mass / k_B * cs_iso_ion2
    loc_temi_min  = 0
    proc_temi_min = 0
    temi_max      = hydro_mass / k_B * cs_iso_ion2
    loc_temi_max  = 0
    proc_temi_max = 0
#else /* ISO */
    wa(:,:,:) = (u(ieni,:,:,:) &                ! eint
                - 0.5*((u(imxi,:,:,:)**2 +u(imyi,:,:,:)**2 &
                  + u(imzi,:,:,:)**2)/u(idni,:,:,:)))
#ifdef MAGNETIC
    wa(:,:,:) = wa(:,:,:) - 0.5*(b(ibx,:,:,:)**2 + b(iby,:,:,:)**2 + &
                 b(ibz,:,:,:)**2)
#endif /* MAGNETIC */
    wa(:,:,:) = max(wa(:,:,:),smallei)
    wa(:,:,:) = (gamma_ion-1.0)*wa(:,:,:)           ! pres

    prei_min      = minval(wa(is:ie,js:je,ks:ke))
    loc_prei_min  = minloc(wa(is:ie,js:je,ks:ke)) + (/nb,nb,nb/)
    call mpifind(prei_min, 'min', loc_prei_min, proc_prei_min)

    prei_max      = maxval(wa(is:ie,js:je,ks:ke))
    loc_prei_max  = maxloc(wa(is:ie,js:je,ks:ke)) + (/nb,nb,nb/)
    call mpifind(prei_max, 'max', loc_prei_max, proc_prei_max)

    temi_max      = maxval( hydro_mass / k_B * wa(is:ie,js:je,ks:ke) &
                                             /u(idni,is:ie,js:je,ks:ke))
    loc_temi_max  = maxloc(wa(is:ie,js:je,ks:ke)    &
                         /u(idni,is:ie,js:je,ks:ke)  ) + (/nb,nb,nb/)
    call mpifind(temi_max, 'max', loc_temi_max, proc_temi_max)

    temi_min      = minval( hydro_mass / k_B * wa(is:ie,js:je,ks:ke) &
                                             /u(idni,is:ie,js:je,ks:ke))
    loc_temi_min  = minloc(wa(is:ie,js:je,ks:ke) &
                         /u(idni,is:ie,js:je,ks:ke)  ) &
                     + (/nb,nb,nb/)
    call mpifind(temi_min, 'min', loc_temi_min, proc_temi_min)

    wa(:,:,:) = gamma_ion*wa(:,:,:)
    csi_max        = sqrt(maxval(wa(is:ie,js:je,ks:ke) &
                            /u(idni,is:ie,js:je,ks:ke)))
    loc_csi_max    = maxloc(wa(is:ie,js:je,ks:ke) &
                            /u(idni,is:ie,js:je,ks:ke)) &
                     + (/nb,nb,nb/)
    call mpifind(csi_max, 'max', loc_csi_max, proc_csi_max)
#endif /* ISO */

#endif /* IONIZED */

#ifdef DUST
    wa            = u(idnd,:,:,:)
    dend_min      = minval(wa(is:ie,js:je,ks:ke))
    loc_dend_min  = minloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
    call mpifind(dend_min, 'min', loc_dend_min, proc_dend_min)

    dend_max      = maxval(wa(is:ie,js:je,ks:ke))
    loc_dend_max  = maxloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
    call mpifind(dend_max, 'max', loc_dend_max, proc_dend_max)

    wa          = abs(u(imxd,:,:,:)/u(idnd,:,:,:))
    vxd_max     = maxval(wa(is:ie,js:je,ks:ke))
    loc_vxd_max = maxloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
    call mpifind(vxd_max, 'max', loc_vxd_max, proc_vxd_max)

    wa          = abs(u(imyd,:,:,:)/u(idnd,:,:,:))
    vyd_max     = maxval(wa(is:ie,js:je,ks:ke))
    loc_vyd_max = maxloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
    call mpifind(vyd_max, 'max', loc_vyd_max, proc_vyd_max)

    wa           = abs(u(imzd,:,:,:)/u(idnd,:,:,:))
    vzd_max      = maxval(wa(is:ie,js:je,ks:ke))
    loc_vzd_max  = maxloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
    call mpifind(vzd_max, 'max', loc_vzd_max, proc_vzd_max)
#endif /* DUST */


#ifdef RESISTIVE
      call mpifind(eta_max, 'max', loc_eta_max, proc_eta_max)
#endif /* RESISTIVE */


#ifdef MAGNETIC
    wa(1:nx-1,1:ny-1,1:max(nz-1,1)) = &
                 (b(ibx,2:nx,1:ny-1,1:max(nz-1,1)) - b(ibx,1:nx-1,1:ny-1,1:max(nz-1,1)))*dy*dz &
                +(b(iby,1:nx-1,2:ny,1:max(nz-1,1)) - b(iby,1:nx-1,1:ny-1,1:max(nz-1,1)))*dx*dz &
                +(b(ibz,1:nx-1,1:ny-1,min(2,nz):nz) - b(ibz,1:nx-1,1:ny-1,1:max(nz-1,1)))*dx*dy
    wa = abs(wa)

    wa(ie,:,:) = wa(ie-1,:,:)
    wa(:,je,:) = wa(:,je-1,:)
    wa(:,:,ke) = wa(:,:,max(ke-1,1))

    divb_max      = maxval(wa(is:ie,js:je,ks:ke))
    loc_divb_max  = maxloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
    call mpifind(divb_max, 'max', loc_divb_max, proc_divb_max)
#endif /* MAGNETIC */

#ifdef COSM_RAYS
    wa            = u(iecr,:,:,:)
    encr_min      = minval(wa(is:ie,js:je,ks:ke))
    loc_encr_min  = minloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
    call mpifind(encr_min, 'min', loc_encr_min, proc_encr_min)

    wa            =  u(iecr,:,:,:)
    encr_max      = maxval(wa(is:ie,js:je,ks:ke))
    loc_encr_max  = maxloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
    call mpifind(encr_max, 'max', loc_encr_max, proc_encr_max)
#endif /* COSM_RAYS */

    if(proc .eq. 0)  then

      open(log_lun, file=log_file, position='append')
#ifdef IONIZED
        write(log_lun,771) 'min(dens)   ION  =', deni_min,  proc_deni_min,  loc_deni_min
        write(log_lun,771) 'max(dens)   ION  =', deni_max,  proc_deni_max,  loc_deni_max
#ifndef ISO
        write(log_lun,771) 'min(temp)   ION  =', temi_min,  proc_temi_min,  loc_temi_min
        write(log_lun,771) 'max(temp)   ION  =', temi_max,  proc_temi_max,  loc_temi_max
#endif /* ISO */
        write(log_lun,771) 'min(pres)   ION  =', prei_min,  proc_prei_min,  loc_prei_min
        write(log_lun,771) 'max(pres)   ION  =', prei_max,  proc_prei_max,  loc_prei_max
        write(log_lun,777) 'max(|vx|)   ION  =', vxi_max, 'dt=',cfl*dx/(vxi_max+small),   proc_vxi_max, loc_vxi_max
        write(log_lun,777) 'max(|vy|)   ION  =', vyi_max, 'dt=',cfl*dy/(vyi_max+small),   proc_vyi_max, loc_vyi_max
        write(log_lun,777) 'max(|vz|)   ION  =', vzi_max, 'dt=',cfl*dz/(vzi_max+small),   proc_vzi_max, loc_vzi_max
        write(log_lun,777) 'max(c_s )   ION  =', csi_max, 'dt=',cfl*dxmn/(csi_max+small), proc_csi_max, loc_csi_max
#ifdef MAGNETIC
        write(log_lun,777) 'max(c_f)    ION  =', sqrt(csi_max**2+vai_max**2), 'dt=',cfl*dxmn/sqrt(csi_max**2+vai_max**2)
        write(log_lun,777) 'max(v_a)    ION  =', vai_max, 'dt=',cfl*dxmn/(vai_max+small), proc_vai_max, loc_vai_max
        write(log_lun,770) 'min(|b|)    MAG  =', b_min,     proc_b_min,     loc_b_min
        write(log_lun,770) 'max(|b|)    MAG  =', b_max,     proc_b_max,     loc_b_max
        write(log_lun,770) 'max(|divb|) MAG  =', divb_max,  proc_divb_max,  loc_divb_max
#else /* MAGNETIC */
        write(log_lun,777) 'max(c_s)    ION  =', sqrt(csi_max**2), 'dt=',cfl*dxmn/sqrt(csi_max**2)
#endif /* MAGNETIC */
#endif /* IONIZED */
#ifdef NEUTRAL
        write(log_lun,771) 'min(dens)   NEU  =', denn_min,  proc_denn_min,  loc_denn_min
        write(log_lun,771) 'max(dens)   NEU  =', denn_max,  proc_denn_max,  loc_denn_max
#ifndef ISO
        write(log_lun,771) 'min(temp)   NEU  =', temn_min,  proc_temn_min,  loc_temn_min
        write(log_lun,771) 'max(temp)   NEU  =', temn_max,  proc_temn_max,  loc_temn_max
#endif /* ISO */
        write(log_lun,771) 'min(pres)   NEU  =', pren_min,  proc_pren_min,  loc_pren_min
        write(log_lun,771) 'max(pres)   NEU  =', pren_max,  proc_pren_max,  loc_pren_max
        write(log_lun,777) 'max(|vx|)   NEU  =', vxn_max, 'dt=',cfl*dx/(vxn_max+small),   proc_vxn_max, loc_vxn_max
        write(log_lun,777) 'max(|vy|)   NEU  =', vyn_max, 'dt=',cfl*dy/(vyn_max+small),   proc_vyn_max, loc_vyn_max
        write(log_lun,777) 'max(|vz|)   NEU  =', vzn_max, 'dt=',cfl*dz/(vzn_max+small),   proc_vzn_max, loc_vzn_max
        write(log_lun,777) 'max(c_s )   NEU  =', csn_max, 'dt=',cfl*dxmn/(csn_max+small), proc_csn_max, loc_csn_max
#endif /* NEUTRAL */
#ifdef DUST
        write(log_lun,771) 'min(dens)   DST  =', dend_min,  proc_dend_min,  loc_dend_min
        write(log_lun,771) 'max(dens)   DST  =', dend_max,  proc_dend_max,  loc_dend_max
        write(log_lun,777) 'max(|vx|)   DST  =', vxd_max, 'dt=',cfl*dx/(vxd_max+small),   proc_vxd_max, loc_vxd_max
        write(log_lun,777) 'max(|vy|)   DST  =', vyd_max, 'dt=',cfl*dy/(vyd_max+small),   proc_vyd_max, loc_vyd_max
        write(log_lun,777) 'max(|vz|)   DST  =', vzd_max, 'dt=',cfl*dz/(vzd_max+small),   proc_vzd_max, loc_vzd_max
#endif /* DUST */
#ifdef COSM_RAYS
        write(log_lun,771) 'min(encr)   CRS  =', encr_min,        proc_encr_min, loc_encr_min
        write(log_lun,777) 'max(encr)   CRS  =', encr_max,      'dt=',dt_crs,     proc_encr_max, loc_encr_max
#endif /* COSM_RAYS */
#ifdef RESISTIVE
        write(log_lun,776) 'max(eta)    RES  =', eta_max ,      'dt=',dt_resist, proc_eta_max,  loc_eta_max
#endif /* RESISTIVE */

      close(log_lun)

    endif

    if(proc .eq. 0)  then
      open(log_lun, file=log_file, position='append')
        write(log_lun,'(a80)') '================================================================================'
        write(log_lun,900) nstep,dt,t
        write(log_lun,'(a80)') '================================================================================'
      close(log_lun)
    endif



770 format(5x,a18,(1x,e15.9),4(1x,i4))

771 format(5x,a18,(1x,e15.9),16x,5(1x,i4))
777 format(5x,a18,(1x,e15.9),2x,a3,(1x,e10.4),5(1x,i4))

900 format('   nstep = ',i7,'   dt = ',f22.16,'   t = ',f22.16,2(1x,i4))
#ifdef RESISTIVE
776 format(5x,a18,(1x,e10.4),2x,a3,(1x,e10.4),4(1x,i4))
#endif /* RESISTIVE */

!#endif /* NOT_WORKING */

  end subroutine write_log


!
!=======================================================================
!
!                    B E G I N   S U B R O U T I N E
!                             R E A D F M S G
!
!=======================================================================
!
    subroutine read_file_msg


!     written by: michal hanasz
!     date:       26. june 2003
!
!-------------------------------------------------------------------------
!     configurable parameters: problem.par
!-------------------------------------------------------------------------
!      user_message_file           ! 1st (user) message file (eg.'./msg')
!      system_message_file         ! 2nd (ups)  message file (eg.'/etc/ups/user/msg')
!-------------------------------------------------------------------------
      implicit none
      character user_last_msg_file*80
      character system_last_msg_file*80


      character user_msg_time(10)*80,system_msg_time(10)*80

      character, save :: user_msg_time_old*80,system_msg_time_old*80
      character*160 syscom
      integer i
      integer(kind=1) :: system

      msg=''
      user_msg_time(9)=''
      nchar=0

      open(91,file=user_message_file,status='old',err=224)
      read(91,fmt=*,err=224,end=224) msg, msg_param

      close(91)
      goto 225
224   continue
      close(91)

      open(92,file=user_message_file,status='old',err=333)
      read(92,fmt=*,err=333,end=333) msg
      close(92)
      goto 225
333   continue
      close(92)
      goto 888

225   continue
      nchar=len_trim(msg)

      user_last_msg_file='./user_last_msg.tmp'
      call rm_file(user_last_msg_file)

      syscom='ls -l --full-time msg >'//user_last_msg_file
      scstatus = system(syscom)
      open(93,file=user_last_msg_file,status='old',err=888)
        read(93,fmt=*,err=888,end=888) (user_msg_time(i), i=1,10)
      close(93)

!---  do the requested action only once for a given user message file

      if(user_msg_time(7) .eq. user_msg_time_old) then
        msg=''
        nchar=0
      else
        system_msg_time_old= user_msg_time(7)
        return
      endif

888   continue

      call rm_file(user_message_file)

!------------------------------------------------------------------------


       open(96,file=system_message_file,status='old',err=424)
       read(96,fmt=*,err=424,end=424) msg, msg_param
       close(96)
       goto 425
424    continue
       close(96)

       open(97,file=system_message_file,status='old',err=533)
       read(97,fmt=*,err=533,end=533) msg
       close(97)
       goto 425
533    continue
       close(97)
       goto 999

425    continue
       nchar=len_trim(msg)

       system_last_msg_file='./system_last_msg.tmp'
       call rm_file(system_last_msg_file)

       syscom='ls -l --full-time '//system_message_file//' > '//system_last_msg_file
       scstatus = system(syscom)

       open(98,file=system_last_msg_file,status='old',err=999)
         read(98,fmt=*,err=999,end=999) (system_msg_time(i), i=1,10)
       close(98)

!---  do the requested action only once for a given system message file

       if(system_msg_time(7) .eq. system_msg_time_old) then
         msg=''
         nchar=0
       else
         system_msg_time_old= system_msg_time(7)
         return
       endif

999    continue


       return
    end subroutine read_file_msg

!------------------------------------------------------------------------

    subroutine rm_file(file_name)
      implicit none
      character*(*) file_name
      character*160 syscom
      logical exist
      integer(kind=1) :: system

!     delete file if exists
      inquire(file = file_name, exist = exist)
      if(exist) then
        syscom='rm -f '//trim(file_name)
        scstatus = system(syscom)
      endif

    end subroutine rm_file

end module dataio
