
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
! modifications and extensions by D.Woltanski in lines: 23, 65, 132, 133, 136-137,
! 474-476, 478, 527-532, 583, 585, 658-669, 829-830, 932, 936, 967, 971-975
  use types
  use mpi_setup
!  use utils
#ifdef RESISTIVITY
!  use resistivity
#endif /* RESISTIVITY */
#ifdef SN_SRC
  use sn_sources
#endif /* SN_SRC */

  implicit none

  integer :: tsl_lun = 2, log_lun = 3
  integer :: nhdf, nres, ntsl, nlog
  integer :: step_hdf, step_res, nhdf_start, nres_start
  real    :: t_start, last_hdf_time
  character(len=128) :: log_file
  character:: pc1*2,pc2*2,pc3*2
  logical tsl_firstcall

  logical wait

  type(hdf) :: chdf


  integer nchar
  character msg*16
  real msg_param

!  integer ctoi
!  external ctoi
  character hostfull*(80), host*8, fhost*10, fpid*10
  integer :: pid, uid, ihost, scstatus
  real vx_max, vy_max, vz_max, va2max,va_max, cs2max, cs_max, &
       dens_min, dens_max, pres_min, pres_max, b_min, b_max,  &
       temp_min, temp_max, divb_max
#ifdef COSM_RAYS
  real encr_min, encr_max
#endif /* COSM_RAYS */



  contains

     subroutine set_container(chdf)
       use types
       implicit none
       type(hdf), intent(out) :: chdf

       chdf%nhdf = nhdf
       chdf%ntsl = ntsl
       chdf%nres = nres
       chdf%nlog = nlog
       chdf%step_hdf = step_hdf
       chdf%log_lun = log_lun
       chdf%last_hdf_time = last_hdf_time
       chdf%log_file = log_file

     end subroutine set_container

     subroutine get_container(chdf)
       use types
       implicit none
       type(hdf), intent(in) :: chdf

       nhdf = chdf%nhdf
       ntsl = chdf%ntsl
       nres = chdf%nres
       nlog = chdf%nlog
       step_hdf =  chdf%step_hdf
       log_lun = chdf%log_lun
       last_hdf_time = chdf%last_hdf_time
       log_file = chdf%log_file

     end subroutine get_container


!=======================================================================
!
!    \\\\\\\\\\        B E G I N   F U N C T I O N        //////////
!    //////////                  C T O I                  \\\\\\\\\\
!
!=======================================================================
!
       function ctoi ( string, istrt, ifin )
!
!    jms:zeus3d.ctoi <------------- converts CHAR*8 string to an integer
!    from jms:zeus2d.strtoi                              ?????????, 19??
!
!    written by: Jim Stone
!    modified 1: January, 1990 by David Clarke; incorporated into ZEUS3D
!
!  PURPOSE:  Converts a segment of a character*8 string into an integer.
!
!  EXTERNALS: [NONE]
!
!-----------------------------------------------------------------------
!
!*ca imp
       implicit none

       character*10   string
       integer       istrt   , ifin    , ishift  , ival    , i
       integer       ctoi

       integer       asciic  (48:57)
!
!      Data Statements
!
       data          asciic  / 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 /
       save          asciic
!
!-----------------------------------------------------------------------
!
       ishift = 1
       ival   = 0
       do 10 i=ifin,istrt,-1
         ival   = ival + ishift * asciic ( ichar(string(i:i)) )
         ishift = ishift * 10
10     continue
       ctoi = ival

       return
       end function ctoi
!
!=======================================================================
!
!    \\\\\\\\\\          E N D   F U N C T I O N          //////////
!    //////////                  C T O I                  \\\\\\\\\\
!
!=======================================================================
!---------------------------------------------------------------------
!
! inititalizes dataio parameters
!
!---------------------------------------------------------------------
!
  subroutine init_dataio
    use start, only : dt_hdf
    implicit none
    integer(kind=1) :: getpid
    integer(kind=1) :: hostnm
    wait  = .false.
    tsl_firstcall = .true.

    nhdf  = 0
    ntsl  = 0
    nres  = 0
    nlog  = 0

    step_hdf  = -1
    step_res  = -1
    last_hdf_time = -dt_hdf

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

  end subroutine init_dataio

!---------------------------------------------------------------------
!
! controls data dumping
!
!---------------------------------------------------------------------
!
  subroutine write_data(output)
    use start, only:  dt_hdf, dt_res, dt_tsl, dt_log, t, dt, nstep
#ifdef HDF5
    use init_problem, only:  problem_name, run_id
    use dataio_hdf5, only: write_hdf5, write_restart_hdf5
#endif /* HDF5 */
#ifdef USER_IO
    use init_problem, only : user_io_routine
#endif /* USER_IO */
    implicit none
    character  :: output*3
#ifdef HDF5
    character :: filename*128
#endif /* HDF5 */


!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if (output .eq. 'log' .or. output .eq. 'end') then
        call write_log
    endif

    if (output .eq. 'log' .or. output .eq. 'end') then
          call write_timeslice
    endif

!    CALL checkdf

#ifdef GALACTIC_DISK
    if(output .eq. 'gpt') then
!      call write_gpot
    endif
#endif /* GALACTIC_DISK */
#ifdef BLOB
      call write_blob
#endif /* BLOB */
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
#ifndef HDF5
        call write_hdf
#else /* HDF5 */
     call set_container(chdf)
     call write_hdf5(chdf)
#endif /* HDF5 */
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
#ifdef HDF5
           if(proc==0) then
              write (filename,'(a,a1,a3,a1,i4.4,a4)') &
                trim(problem_name),'_', run_id,'_',nres,'.res'
           endif
           call MPI_BCAST(filename, 128, MPI_CHARACTER, 0, comm, ierr)
           call set_container(chdf); chdf%nstep = nstep
           call write_restart_hdf5(filename,chdf)
#else /* HDF5 */
           call write_restart
#endif /* HDF5 */
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
    use start, only : vars,xmin,xmax,ymin,ymax,zmin,zmax,nstep,t, &
         domain, nxd,nyd,nzd,nb,mag_center, gamma, dimensions, dt
    use grid, only : dx,dy,dz
    use arrays, only : nx,ny,nz,nxb,nyb,nzb,x,y,z,wa,outwa,outwb,outwc,b,u, &
         idna,imxa,imya,imza,ibx,iby,ibz
#ifdef COSM_RAYS
    use arrays, only : iecr
#endif /* COSM_RAYS */
    use init_problem, only : problem_name, run_id
#ifndef ISO
    use arrays, only : iena
    use start, only  : gamma
#endif /* ISO */
#ifndef SPLIT
    use arrays, only : Lu
#endif /* SPLIT */
#ifdef GRAV
    use arrays, only : gp
#endif /* GRAV */
    implicit none

    character(len=128) :: file_name_hdf,file_name_disp

    integer :: sd_id, sds_id, dim_id, iostatus
    integer :: rank, comp_type
    integer, dimension(3) :: dims, istart, stride
    integer, dimension(1) :: comp_prm
    integer :: sfstart, sfend, sfsnatt, sfcreate, sfwdata, sfscompress, sfendacc &
             , sfdimid, sfsdmname, sfsdscale, sfsdmstr

    integer :: iv, ibe, jbe
    integer :: nxo, nyo, nzo, iso,ieo,jso,jeo,kso,keo
    real(kind=4), dimension(:,:,:), allocatable :: temp

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if(domain .eq. 'full_domain') then
      nxo = nx
      nyo = ny
      nzo = nz
      iso = 1
      ieo = nx
      jso = 1
      jeo = ny
      if(dimensions .eq. '3d') then
        kso = 1
        keo = nz
      elseif(dimensions .eq. '2dxy') then
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
      if(dimensions .eq. '3d') then
        kso    = nb+1
        keo    = nb+nzb
      elseif(dimensions .eq. '2dxy') then
        kso = 1
        keo = 1
      endif
    endif

    allocate(temp(iso:ieo,jso:jeo,kso:keo))

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

    iostatus = sfsnatt( sd_id, 'gamma'   ,  6, 1, gamma          )


! write selected problem dependent parameters
#ifdef SN_SRC
    iostatus = sfsnatt( sd_id, 'nsn'      , 24, 1, nsn     )
#endif /* SN_SRC */

    iv = 1

    do while (len_trim(vars(iv)) .ne. 0)

      select case(vars(iv))
      case ('dens')
        wa(iso:ieo,jso:jeo,kso:keo) = u(idna,iso:ieo,jso:jeo,kso:keo)
      case ('velx')
        wa(iso:ieo,jso:jeo,kso:keo) = u(imxa,iso:ieo,jso:jeo,kso:keo) / u(idna,iso:ieo,jso:jeo,kso:keo)

      case ('vely')
        wa(iso:ieo,jso:jeo,kso:keo) = u(imya,iso:ieo,jso:jeo,kso:keo) / u(idna,iso:ieo,jso:jeo,kso:keo)

      case ('velz')
        wa(iso:ieo,jso:jeo,kso:keo) = u(imza,iso:ieo,jso:jeo,kso:keo) / u(idna,iso:ieo,jso:jeo,kso:keo)

#ifdef ISO
      case ('ener')
        wa(iso:ieo,jso:jeo,kso:keo) = 0.5*(u(imxa,iso:ieo,jso:jeo,kso:keo)**2 &
                                    +u(imya,iso:ieo,jso:jeo,kso:keo)**2 &
                                    +u(imza,iso:ieo,jso:jeo,kso:keo)**2)/u(idna,iso:ieo,jso:jeo,kso:keo)
#else /* ISO */
      case ('ener')
        wa(iso:ieo,jso:jeo,kso:keo) = u(iena,iso:ieo,jso:jeo,kso:keo)
      case ('eint')
        wa(iso:ieo,jso:jeo,kso:keo) = u(iena,iso:ieo,jso:jeo,kso:keo) &
                                    - 0.5*(u(imxa,iso:ieo,jso:jeo,kso:keo)**2 &
                                          +u(imya,iso:ieo,jso:jeo,kso:keo)**2 &
                                          +u(imza,iso:ieo,jso:jeo,kso:keo)**2)/u(idna,iso:ieo,jso:jeo,kso:keo)
#endif /* ISO */
#ifdef COSM_RAYS
      case ('encr')
        wa(iso:ieo,jso:jeo,kso:keo) = u(iecr,iso:ieo,jso:jeo,kso:keo)
#endif /* COSM_RAYS */

      case ('divb')
        wa(iso:ieo-1,jso:jeo-1,kso:keo-1) = &
           (b(1,iso+1:ieo, jso:jeo-1, kso:keo-1) - b(1,iso:ieo-1, jso:jeo-1, kso:keo-1))*dy*dz &
          +(b(2,iso:ieo-1, jso+1:jeo, kso:keo-1) - b(2,iso:ieo-1, jso:jeo-1, kso:keo-1))*dx*dz &
          +(b(3,iso:ieo-1, jso:jeo-1, kso+1:keo) - b(3,iso:ieo-1, jso:jeo-1, kso:keo-1))*dx*dy
        wa = abs(wa)
        wa(ieo,:,:) = wa(ieo-1,:,:)
        wa(:,jeo,:) = wa(:,jeo-1,:)
        wa(:,:,keo) = wa(:,:,keo-1)

      case ('omga')
        do ibe=iso,ieo
          do jbe=jso,jeo
            wa(ibe,jbe,kso:keo) = (u(imya,ibe,jbe,kso:keo) / u(idna,ibe,jbe,kso:keo) * x(ibe) &
                                 - u(imxa,ibe,jbe,kso:keo) / u(idna,ibe,jbe,kso:keo) * y(jbe))/(x(ibe)**2+y(jbe)**2)
          enddo
        enddo

      case ('vrot')
        do ibe=iso,ieo
          do jbe=jso,jeo
            wa(ibe,jbe,kso:keo) = (u(imya,ibe,jbe,kso:keo) / u(idna,ibe,jbe,kso:keo) * x(ibe) &
                                 - u(imxa,ibe,jbe,kso:keo) / u(idna,ibe,jbe,kso:keo) * y(jbe))/sqrt(x(ibe)**2+y(jbe)**2)
          enddo
        enddo

      case ('vout')
        do ibe=iso,ieo
          do jbe=jso,jeo
            wa(ibe,jbe,kso:keo) = (u(imxa,ibe,jbe,kso:keo) / u(idna,ibe,jbe,kso:keo) * x(ibe) &
                                 + u(imya,ibe,jbe,kso:keo) / u(idna,ibe,jbe,kso:keo) * y(jbe))/sqrt(x(ibe)**2+y(jbe)**2)
          enddo
        enddo

      case ('dcol')
        do ibe=iso,ieo
          do jbe=jso,jeo
            wa(ibe,jbe,kso:keo) = sum(u(idna,ibe,jbe,kso:keo))*(z(kso)-z(kso-1))
          enddo
        enddo

#ifndef ISO
      case ('csnd')
        wa(iso:ieo,jso:jeo,kso:keo) = sqrt((u(iena,iso:ieo,jso:jeo,kso:keo) &
       -0.5*(u(imxa,iso:ieo,jso:jeo,kso:keo)**2+u(imya,iso:ieo,jso:jeo,kso:keo)**2+u(imza,iso:ieo,jso:jeo,kso:keo)**2) &
        / u(idna,iso:ieo,jso:jeo,kso:keo))*(gamma-1.)*gamma/u(idna,iso:ieo,jso:jeo,kso:keo))
#endif /* ISO */
#ifdef GRAV
      case ('gpot')
        wa(iso:ieo,jso:jeo,kso:keo) = gp(iso:ieo,jso:jeo,kso:keo)
#endif /* GRAV */
      case ('magx')
        if(domain .eq. 'full_domain') then
          if(mag_center .eq. 'yes') then
            wa(:,:,:) = 0.5*(b(ibx,:,:,:))
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
        if(domain .eq. 'full_domain') then
          if(mag_center .eq. 'yes') then
            wa(:,:,:) = 0.5*(b(iby,:,:,:))
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
        if(domain .eq. 'full_domain') then
          if(mag_center .eq. 'yes') then
            wa(:,:,:) = 0.5*(b(ibz,:,:,:))
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

      case ('curx')
       wa(:,:,:) = (b(ibz,:,:,:)-cshift(b(ibz,:,:,:),shift=-1,dim=2))/dy &
                 - (b(iby,:,:,:)-cshift(b(iby,:,:,:),shift=-1,dim=3))/dz
       if(domain .eq. 'full_domain') then
         if(mag_center .eq. 'yes') then
           wa(:,:,:) = 0.5*(wa(:,:,:) + cshift(wa(:,:,:),shift=1,dim=1))
         endif
       else if(domain .eq. 'phys_domain') then
         if(mag_center .eq. 'yes') then
           wa(iso:ieo,jso:jeo,kso:keo) = 0.5*(wa(iso:ieo,jso:jeo,kso:keo)+wa(iso+1:ieo+1,jso:jeo,kso:keo))
         endif
       endif

      case ('cury')
       wa(:,:,:) = (b(ibx,:,:,:)-cshift(b(ibx,:,:,:),shift=-1,dim=3))/dz &
                 - (b(ibz,:,:,:)-cshift(b(ibz,:,:,:),shift=-1,dim=1))/dx
       if(domain .eq. 'full_domain') then
         if(mag_center .eq. 'yes') then
           wa(:,:,:) = 0.5*(wa(:,:,:) + cshift(wa(:,:,:),shift=1,dim=2))
         endif
       else if(domain .eq. 'phys_domain') then
         if(mag_center .eq. 'yes') then
           wa(iso:ieo,jso:jeo,kso:keo) = 0.5*(wa(iso:ieo,jso:jeo,kso:keo)+wa(iso:ieo,jso+1:jeo+1,kso:keo))
         endif
       endif

      case ('curz')
       wa(:,:,:) = (b(iby,:,:,:)-cshift(b(iby,:,:,:),shift=-1,dim=1))/dx &
                 - (b(ibx,:,:,:)-cshift(b(ibx,:,:,:),shift=-1,dim=2))/dy
       if(domain .eq. 'full_domain') then
         if(mag_center .eq. 'yes') then
           wa(:,:,:) = 0.5*(wa(:,:,:)+ cshift(wa(:,:,:),shift=1,dim=3))
         endif
       else if(domain .eq. 'phys_domain') then
         if(mag_center .eq. 'yes') then
           wa(iso:ieo,jso:jeo,kso:keo) = 0.5*(wa(iso:ieo,jso:jeo,kso:keo) + wa(iso:ieo,jso:jeo,kso+1:keo+1))
         endif
       endif
#ifndef SPLIT
      case ('flx1')
        wa(:,:,:) = Lu(1,:,:,:)
      case ('flx2')
        wa(:,:,:) = Lu(2,:,:,:)
      case ('flx3')
        wa(:,:,:) = Lu(3,:,:,:)
      case ('flx4')
        wa(:,:,:) = Lu(4,:,:,:)
#ifndef ISO
      case ('flx5')
        wa(:,:,:) = Lu(5,:,:,:)
#endif /* ~ISO */
#endif /* SPLIT */
!DW+
      case ('grvx')
        wa(iso:ieo,jso:jeo,kso:keo) = outwa(iso:ieo,jso:jeo,kso:keo)

      case ('grvy')
        wa(iso:ieo,jso:jeo,kso:keo) = outwb(iso:ieo,jso:jeo,kso:keo)

      case ('grvz')
        wa(iso:ieo,jso:jeo,kso:keo) = outwc(iso:ieo,jso:jeo,kso:keo)
!DW-
      case default
        print *, 'Variable ', vars(iv), ' is not defined! Skipping.'
      end select

! write data
!
      sds_id = sfcreate(sd_id, vars(iv), 5, rank, dims)
      iostatus = sfscompress(sds_id, comp_type, comp_prm)
      temp = real(wa(iso:ieo,jso:jeo,kso:keo),4)
      iostatus = sfwdata(sds_id, istart, stride, dims, temp)

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

      iv = iv + 1
    end do

    iostatus = sfend(sd_id)


    if(proc.eq. 0) then
      open(log_lun, file=log_file, position='append')

      write(log_lun,*) 'Writing output   file: ', trim(file_name_disp)
      write(*,*)       'Writing output   file: ', trim(file_name_disp)

      close(log_lun)
    endif
    deallocate(temp)


  end subroutine write_hdf

#ifdef GALACTIC_DISK
  subroutine write_gpot

    use arrays, only       : wa, gp,nxb,nyb,nzb,nx,ny,nz,x,y,z
    use init_problem, only : problem_name, run_id
    use gravity, only      : gpotdisk,gpothalo,gpotbulge
    use grid, only         : dx,dy,dz
    use start, only        : nstep,t,dt,nxd,nyd,nzd,nb,xmin,xmax,ymin,ymax,zmin,zmax,dimensions,gamma,domain


    implicit none

    character(len=128) :: file_name_hdf,file_name_disp
    character, dimension(24) :: gvars*4

    integer :: sd_id, sds_id, dim_id, iostatus
    integer :: rank, comp_type
    integer, dimension(3) :: dims, istart, stride
    integer, dimension(1) :: comp_prm
    integer :: sfstart, sfend, sfsnatt, sfcreate, sfwdata, sfscompress, sfendacc &
             , sfdimid, sfsdmname, sfsdscale, sfsdmstr

    integer :: iv, i, j, ibe, jbe
    integer :: nxo, nyo, nzo, iso,ieo,jso,jeo,kso,keo

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if(domain .eq. 'full_domain') then
      nxo   = nx
      nyo   = ny
      nzo   = nz
      iso   = 1
      ieo   = nx
      jso   = 1
      jeo   = ny
      if(dimensions .eq. '3d') then
      kso   = 1
      keo   = nz
      elseif(dimensions .eq. '2dxy') then
        kso = 1
        keo = 1
      endif

    else if(domain .eq. 'phys_domain') then
      nxo   = nxb
      nyo   = nyb
      nzo   = nzb
      iso    = nb+1
      ieo    = nb+nxb
      jso    = nb+1
      jeo    = nb+nyb
      if(dimensions .eq. '3d') then
        kso    = nb+1
        keo    = nb+nzb
      elseif(dimensions .eq. '2dxy') then
        kso = 1
        keo = 1
      endif
    endif


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

    write (file_name_hdf,'(a,a1,a,a5,i2.2,a1,i2.2,a1,i2.2,a1,a8)') &
              trim(cwd),'/',trim(problem_name),'_gpt_', &
	      pcoords(1),'_',pcoords(2),'_',pcoords(3),'_', '0000.hdf'

    write (file_name_disp,'(a,a1,a,a5,a2,a1,a2,a1,a2,a1,a8)') &
              trim(cwd),'/',trim(problem_name),'_gpt_',pc1,'_',pc2,'_',pc3,'_', '0000.hdf'


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

    iostatus = sfsnatt( sd_id, 'gamma'   ,  6, 1, gamma          )

    do iv=1,24

      select case(iv)

      case (1)
        gvars(1)='gpot'
        wa(iso:ieo,jso:jeo,kso:keo) = gp(iso:ieo,jso:jeo,kso:keo)

      case (2)
        gvars(2)='gdsk'
        wa(iso:ieo,jso:jeo,kso:keo) = gpotdisk(iso:ieo,jso:jeo,kso:keo)

      case (3)
        gvars(3)='ghal'
        wa(iso:ieo,jso:jeo,kso:keo) = gpothalo(iso:ieo,jso:jeo,kso:keo)

      case (4)
        gvars(4)='gblg'
        wa(iso:ieo,jso:jeo,kso:keo) = gpotbulge(iso:ieo,jso:jeo,kso:keo)

      case (5)
        gvars(5)='grxl'
        wa(iso:ieo,jso:jeo,kso:keo) = (gp(iso:ieo,jso:jeo,kso:keo)-gp(iso-1:ieo-1,jso:jeo,kso:keo))/dx

      case (6)
        gvars(6)='gryl'
        wa(iso:ieo,jso:jeo,kso:keo) = (gp(iso:ieo,jso:jeo,kso:keo)-gp(iso:ieo,jso-1:jeo-1,kso:keo))/dy

      case (7)
        gvars(7)='grzl'
        wa(iso:ieo,jso:jeo,kso:keo) = (gp(iso:ieo,jso:jeo,kso:keo)-gp(iso:ieo,jso:jeo,kso-1:keo-1))/dz

      case (8)
        gvars(8)='gdxl'
        wa(iso:ieo,jso:jeo,kso:keo) = (gpotdisk(iso:ieo,jso:jeo,kso:keo)-gpotdisk(iso-1:ieo-1,jso:jeo,kso:keo))/dx

      case (9)
        gvars(9)='gdyl'
        wa(iso:ieo,jso:jeo,kso:keo) = (gpotdisk(iso:ieo,jso:jeo,kso:keo)-gpotdisk(iso:ieo,jso-1:jeo-1,kso:keo))/dy

      case (10)
        gvars(10)='gdzl'
        wa(iso:ieo,jso:jeo,kso:keo) = (gpotdisk(iso:ieo,jso:jeo,kso:keo)-gpotdisk(iso:ieo,jso:jeo,kso-1:keo-1))/dz

      case (11)
        gvars(11)='ghxl'
        wa(iso:ieo,jso:jeo,kso:keo) = (gpothalo(iso:ieo,jso:jeo,kso:keo)-gpothalo(iso-1:ieo-1,jso:jeo,kso:keo))/dx

      case (12)
        gvars(12)='ghyl'
        wa(iso:ieo,jso:jeo,kso:keo) = (gpothalo(iso:ieo,jso:jeo,kso:keo)-gpothalo(iso:ieo,jso-1:jeo-1,kso:keo))/dy

      case (13)
        gvars(13)='ghzl'
        wa(iso:ieo,jso:jeo,kso:keo) = (gpothalo(iso:ieo,jso:jeo,kso:keo)-gpothalo(iso:ieo,jso:jeo,kso-1:keo-1))/dz

      case (14)
        gvars(14)='gbxl'
        wa(iso:ieo,jso:jeo,kso:keo) = (gpotbulge(iso:ieo,jso:jeo,kso:keo)-gpotbulge(iso-1:ieo-1,jso:jeo,kso:keo))/dx

      case (15)
        gvars(15)='gbyl'
        wa(iso:ieo,jso:jeo,kso:keo) = (gpotbulge(iso:ieo,jso:jeo,kso:keo)-gpotbulge(iso:ieo,jso-1:jeo-1,kso:keo))/dy

      case (16)
        gvars(16)='gbzl'
        wa(iso:ieo,jso:jeo,kso:keo) = (gpotbulge(iso:ieo,jso:jeo,kso:keo)-gpotbulge(iso:ieo,jso:jeo,kso-1:keo-1))/dz

      case (17)
        gvars(17)='vrgx'
	do ibe = iso,ieo
	do jbe = jso,jeo
        wa(ibe,jbe,kso:keo) = (gp(ibe+1,jbe,kso:keo)-gp(ibe-1,jbe,kso:keo))/dx/2.0
	wa(ibe,jbe,kso:keo) = sqrt(abs(wa(ibe,jbe,kso:keo)/x(ibe)))*sqrt(x(ibe)**2+y(jbe)**2)
	enddo
	enddo

      case (18)
        gvars(18)='vrgy'
	do ibe = iso,ieo
	do jbe = jso,jeo
        wa(ibe,jbe,kso:keo) = (gp(ibe,jbe+1,kso:keo)-gp(ibe,jbe-1,kso:keo))/dy/2.0
	wa(ibe,jbe,kso:keo) = sqrt(abs(wa(ibe,jbe,kso:keo)/y(jbe)))*sqrt(x(ibe)**2+y(jbe)**2)
	enddo
	enddo

      case (19)
        gvars(19)='vrdx'
	do ibe = iso,ieo
	do jbe = jso,jeo
        wa(ibe,jbe,kso:keo) = (gpotdisk(ibe+1,jbe,kso:keo)-gpotdisk(ibe-1,jbe,kso:keo))/dx/2.0
	wa(ibe,jbe,kso:keo) = sqrt(abs(wa(ibe,jbe,kso:keo)/x(ibe)))*sqrt(x(ibe)**2+y(jbe)**2)
	enddo
	enddo

      case (20)
        gvars(20)='vrdy'
	do ibe = iso,ieo
	do jbe = jso,jeo
        wa(ibe,jbe,kso:keo) = (gpotdisk(ibe,jbe+1,kso:keo)-gpotdisk(ibe,jbe-1,kso:keo))/dy/2.0
	wa(ibe,jbe,kso:keo) = sqrt(abs(wa(ibe,jbe,kso:keo)/y(jbe)))*sqrt(x(ibe)**2+y(jbe)**2)
	enddo
	enddo

      case (21)
        gvars(21)='vrhx'
	do ibe = iso,ieo
	do jbe = jso,jeo
        wa(ibe,jbe,kso:keo) = (gpothalo(ibe+1,jbe,kso:keo)-gpothalo(ibe-1,jbe,kso:keo))/dx/2.0
	wa(ibe,jbe,kso:keo) = sqrt(abs(wa(ibe,jbe,kso:keo)/x(ibe)))*sqrt(x(ibe)**2+y(jbe)**2)
	enddo
	enddo

      case (22)
        gvars(22)='vrhy'
	do ibe = iso,ieo
	do jbe = jso,jeo
        wa(ibe,jbe,kso:keo) = (gpothalo(ibe,jbe+1,kso:keo)-gpothalo(ibe,jbe-1,kso:keo))/dy/2.0
	wa(ibe,jbe,kso:keo) = sqrt(abs(wa(ibe,jbe,kso:keo)/y(jbe)))*sqrt(x(ibe)**2+y(jbe)**2)
	enddo
	enddo

      case (23)
        gvars(23)='vrbx'
	do ibe = iso,ieo
	do jbe = jso,jeo
        wa(ibe,jbe,kso:keo) = (gpotbulge(ibe+1,jbe,kso:keo)-gpotbulge(ibe-1,jbe,kso:keo))/dx/2.0
	wa(ibe,jbe,kso:keo) = sqrt(abs(wa(ibe,jbe,kso:keo)/x(ibe)))*sqrt(x(ibe)**2+y(jbe)**2)
	enddo
	enddo

      case (24)
        gvars(24)='vrby'
	do ibe = iso,ieo
	do jbe = jso,jeo
        wa(ibe,jbe,kso:keo) = (gpotbulge(ibe,jbe+1,kso:keo)-gpotbulge(ibe,jbe-1,kso:keo))/dy/2.0
	wa(ibe,jbe,kso:keo) = sqrt(abs(wa(ibe,jbe,kso:keo)/y(jbe)))*sqrt(x(ibe)**2+y(jbe)**2)
	enddo
	enddo

      case default
        print *, 'Variable ', gvars(iv), ' is not defined! Skipping.'
      end select
!       print *, gvars(iv), minval(wa), maxval(wa)

! write data
!
      sds_id = sfcreate(sd_id, gvars(iv), 5, rank, dims)
      iostatus = sfscompress(sds_id, comp_type, comp_prm)
      iostatus = sfwdata(sds_id, istart, stride, dims, real(wa(iso:ieo,jso:jeo,kso:keo),4))

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

    end do

    iostatus = sfend(sd_id)


    if(proc.eq. 0) then
!      open(log_lun, file=log_file, position='append')

!      write(log_lun,*) 'Writing output   file: ', trim(file_name_disp)
      write(*,*)       'Writing output   file: ', trim(file_name_disp)

!      close(log_lun)
    endif


  end subroutine write_gpot
#endif /* GALACTIC_DISK */

!---------------------------------------------------------------------
!
! writes restart file
!
!---------------------------------------------------------------------
!
  subroutine write_restart

    use arrays, only : nxb,nyb,nzb,x,y,z,u,b,nx,ny,nz,nm,nu
    use start, only  : xmin,xmax,ymin,ymax,zmin,zmax,nxd,nyd,nzd,t,dt, &
         resdel,nb,nstep,domain
    use init_problem, only : problem_name, run_id
#ifdef GRAV
    use arrays, only : gp
#endif /* GRAV */
#ifdef MASS_COMPENS
    use arrays, only : dinit
#endif /*  MASS_COMPENS */

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
    dimsu(1)    = nu
    dimsu(2)    = dims(1)
    dimsu(3)    = dims(2)
    dimsu(4)    = dims(3)

    rankb       = 4
    dimsb(1)    = nm
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

    iostatus = sfsnatt( sd_id, 'dimsu'   , 23, 1, nu )
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
#ifdef MASS_COMPENS
! write initial density distribution
!
    sds_id = sfcreate(sd_id, 'dinit', 6, rank3d, dims3d)
    iostatus = sfscompress(sds_id, comp_type, comp_prm)
    iostatus = sfwdata(sds_id, istart, stride, dims3d, dinit)
#endif /* MASS_COMPENS */

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

    use arrays, only : u,b,nx,ny,nz,nu,nm
    use start, only  : t, dt, nstep
    use init_problem, only : problem_name, run_id
#ifdef GRAV
    use arrays, only : gp
#endif /* GRAV */
#ifdef MASS_COMPENS
    use arrays, only : dinit
#endif /*  MASS_COMPENS */

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
    dimsu(1) = nu
    dimsu(2) = nx
    dimsu(3) = ny
    dimsu(4) = nz

    rankb    = 4
    dimsb(1) = nm
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
#ifdef MASS_COMPENS
! read initial density distribution
!
      sds_index = sfn2index(sd_id, 'dinit')
      sds_id = sfselect(sd_id, sds_index)
      iostatus = sfrdata(sds_id, istart, stride, dims3d, dinit)
#endif /* MASS_COMPENS */

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
#ifdef HDF5
        write (file_name,'(a,a1,a,a1,a3,a1,i4.4,a4)') &
               trim(cwd),'/',trim(problem_name),'_', run_id,'_',nres,'.res'
#else /* HDF5 */
        write (file_name,'(a,a1,a,a1,a3,a1,3(i2.2,a1),i3.3,a4)') &
               trim(cwd),'/',trim(problem_name),'_', run_id,'_',0,'_',0,'_',0,'_',nres,'.res'
#endif /* HDF5 */
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

    use arrays, only : is,ie,js,je,ks,ke,u,b,idna,imxa,imya,imza, wa, &
         ibx,iby,ibz,x,y,z,wa
#ifdef COSM_RAYS
    use arrays, only : iecr
#endif /* COSM_RAYS */
    use grid, only  : dvol,dx,dy,dz
    use start, only : proc, dt, t, nstep, nrestart,nxd,nyd,nzd
#ifdef GRAV
    use arrays, only : gp
#endif /* GRAV */
#ifdef ISO
    use start, only : csi2
#endif /* ISO */
    use init_problem, only : problem_name, run_id
#ifdef COOL_HEAT
    use thermal
#endif /* COOL_HEAT */
#ifdef RESIST
    use resistivity
#endif /* RESIST */
#ifndef ISO
    use arrays, only : iena
#endif /* ISO */
#ifdef COSM_RAYS
    use arrays, only : iecr
#endif /* COSM_RAYS */
#ifdef SNE_DISTR
    use sn_distr, only : emagadd, tot_emagadd
#endif /* SNE_DISTR */

    implicit none
    integer i,j

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
                                           'ener', 'epot', 'eint', 'ekin', 'emag', &
                                           'mflx', 'mfly', 'mflz', &
#ifdef COSM_RAYS
                                           'encr_tot', 'encr_min',  'encr_max',&
#endif /* COSM_RAYS */
#ifdef RESIST
                                           'eta_max', &
#endif /* RESIST */
! some quantities computed in "write_log".One can add more, or change.
                                           'vx_max', 'vy_max', 'vz_max', 'va_max', 'cs_max', &
                                           'dens_min', 'dens_max', 'pres_min', 'pres_max', &
#ifndef ISO
                                           'temp_min', 'temp_max',  &
#endif /* ISO */
#ifdef SNE_DISTR
                                           'sum_emagadd', 'tot_emagadd', &
#endif /* SNE_DISTR */
                                           'b_min', 'b_max', 'divb_max'


        write (tsl_lun, '(a1)') '#'
        tsl_firstcall = .false.
      else
        open(tsl_lun, file=tsl_file, position='append')
      endif
    endif

    mass = sum(u(idna,is:ie,js:je,ks:ke)) * dvol
    call mpi_allreduce(mass, tot_mass, 1, mpi_real8, mpi_sum, comm3d, ierr)

    momx = sum(u(imxa,is:ie,js:je,ks:ke)) * dvol
    call mpi_allreduce(momx, tot_momx, 1, mpi_real8, mpi_sum, comm3d, ierr)

    momy = sum(u(imya,is:ie,js:je,ks:ke)) * dvol
    call mpi_allreduce(momy, tot_momy, 1, mpi_real8, mpi_sum, comm3d, ierr)

    momz = sum(u(imza,is:ie,js:je,ks:ke)) * dvol
    call mpi_allreduce(momz, tot_momz, 1, mpi_real8, mpi_sum, comm3d, ierr)

#ifdef GRAV
!DW+
    amomz = 0.0
    do j=js,je
      do i=is,ie
         amomz = amomz + (x(i)*sum(u(imya,i,j,ks:ke)) &
                         -y(j)*sum(u(imxa,i,j,ks:ke)))*dvol
      enddo
    enddo
    call mpi_allreduce(amomz, tot_amomz, 1, mpi_real8, mpi_sum, comm3d, ierr)
!DW-
    epot = sum(u(idna,is:ie,js:je,ks:ke) *gp(is:ie,js:je,ks:ke)) * dvol
    call mpi_allreduce(epot, tot_epot, 1, mpi_real8, mpi_sum, comm3d, ierr)
#endif /* GRAV */

    wa(is:ie,js:je,ks:ke) &
        = 0.5 * (u(imxa,is:ie,js:je,ks:ke)**2   &
               + u(imya,is:ie,js:je,ks:ke)**2  &
               + u(imza,is:ie,js:je,ks:ke)**2)/ &
                 u(idna,is:ie,js:je,ks:ke)
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
    ener = sum(u(iena,is:ie,js:je,ks:ke)) * dvol
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
#ifdef RESIST
                      eta_max, &
#endif /* RESIST */
! some quantities computed in "write_log".One can add more, or change.
                      vx_max, vy_max, vz_max, va_max, cs_max, &
                      dens_min, dens_max, pres_min, pres_max, &
#ifndef ISO
                      temp_min, temp_max,  &
#endif /* ISO */
#ifdef SNE_DISTR
                      sum_emagadd, tot_emagadd, &
#endif /* SNE_DISTR */
                      b_min, b_max, divb_max
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

    use arrays, only : wa,is,ie,js,je,ks,ke,idna,imxa,imya,imza,u,b,nx,ny,nz
    use grid, only   : dx,dy,dz,dxmn
    use constants, only : small, hydro_mass, k_B
    use start, only : t,dt,nstep,sleep_minutes,sleep_seconds, smallei,nb, &
         gamma,cfl
#ifdef COSM_RAYS
    use start, only  : dt_cr
    use arrays, only : iecr
#endif /* COSM_RAYS */
#ifdef ISO
    use start, only : csi2,c_si
#endif /* ISO */
!    use init_problem
!   use thermal
#ifndef ISO
    use arrays, only : iena
#endif /* ISO */
#ifdef RESIST
    use resistivity
#endif /* RESIST */

    implicit none

    integer, dimension(3) :: loc_vx_max, loc_vy_max, loc_vz_max, loc_va_max, &
                             loc_cs_max, loc_dens_min, loc_dens_max, loc_pres_min, &
                             loc_pres_max, loc_b_min, loc_b_max, &
                             loc_temp_min, loc_temp_max, loc_divb_max
#ifdef COOL_HEAT
    integer, dimension(3) :: loc_dt_cool, loc_dt_heat
#endif /* COOL_HEAT */
#ifdef COSM_RAYS
    integer, dimension(3) :: loc_encr_min, loc_encr_max
#endif /* COSM_RAYS */

    integer               :: proc_vx_max, proc_vy_max, proc_vz_max, proc_va_max, &
                             proc_cs_max, proc_dens_min, proc_dens_max, proc_pres_min, &
                             proc_pres_max, proc_b_min, proc_b_max, &
                             proc_temp_min, proc_temp_max, proc_divb_max
#ifdef COOL_HEAT
    integer               :: proc_dt_cool, proc_dt_heat
#endif /* COOL_HEAT */
#ifdef RESIST
    integer               :: proc_eta_max
#endif /* RESIST */
#ifdef COSM_RAYS
    integer               :: proc_encr_min, proc_encr_max
#endif /* COSM_RAYS */

! Timestep diagnostics

    wa         =  u(idna,:,:,:)
    dens_min      = minval(wa(is:ie,js:je,ks:ke))
    loc_dens_min  = minloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
    call mpifind(dens_min, 'min', loc_dens_min, proc_dens_min)

    wa         =  u(idna,:,:,:)
    dens_max      = maxval(wa(is:ie,js:je,ks:ke))
    loc_dens_max  = maxloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
    call mpifind(dens_max, 'max', loc_dens_max, proc_dens_max)


    wa         = abs(u(imxa,:,:,:)/u(idna,:,:,:))
    vx_max      = maxval(wa(is:ie,js:je,ks:ke))
    loc_vx_max  = maxloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
    call mpifind(vx_max, 'max', loc_vx_max, proc_vx_max)

    wa         = abs(u(imya,:,:,:)/u(idna,:,:,:))
    vy_max      = maxval(wa(is:ie,js:je,ks:ke))
    loc_vy_max  = maxloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
    call mpifind(vy_max, 'max', loc_vy_max, proc_vy_max)

    wa         = abs(u(imza,:,:,:)/u(idna,:,:,:))
    vz_max      = maxval(wa(is:ie,js:je,ks:ke))
    loc_vz_max  = maxloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
    call mpifind(vz_max, 'max', loc_vz_max, proc_vz_max)

!    wa         = sum(b(:,:,:,:)**2,1)
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

    va_max      = sqrt(maxval(wa(is:ie,js:je,ks:ke) &
                       /u(idna,is:ie,js:je,ks:ke)))
    loc_va_max  = maxloc(wa(is:ie,js:je,ks:ke) &
                       /u(idna,is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
    call mpifind(va_max, 'max', loc_va_max, proc_va_max)

    wa(1:nx-1,1:ny-1,1:nz-1) = &
          (b(1,2:nx,1:ny-1,1:nz-1) - b(1,1:nx-1,1:ny-1,1:nz-1))*dy*dz &
         +(b(2,1:nx-1,2:ny,1:nz-1) - b(2,1:nx-1,1:ny-1,1:nz-1))*dx*dz &
         +(b(3,1:nx-1,1:ny-1,2:nz) - b(3,1:nx-1,1:ny-1,1:nz-1))*dx*dy

    wa(nx,:,:) = wa(nx-1,:,:)
    wa(:,ny,:) = wa(:,ny-1,:)
    wa(:,:,nz) = wa(:,:,nz-1)
    wa = abs(wa)

    divb_max      = maxval(wa(is:ie,js:je,ks:ke))
    loc_divb_max  = maxloc(wa(is:ie,js:je,ks:ke)) + (/nb,nb,nb/)
    call mpifind(divb_max, 'max', loc_divb_max, proc_divb_max)

#ifdef ISO
    pres_min     = csi2*dens_min
    loc_pres_min  = loc_dens_min
    proc_pres_min = proc_dens_min
    pres_max     = csi2*dens_max
    loc_pres_max  = loc_dens_max
    proc_pres_max = proc_dens_max
    cs_max        = c_si
    loc_cs_max    = 0
    proc_cs_max   = 0
    temp_min      = hydro_mass / k_B * csi2
    loc_temp_min  = 0
    proc_temp_min = 0
    temp_max      = hydro_mass / k_B * csi2
    loc_temp_max  = 0
    proc_temp_max = 0
#else /* ISO */
    wa            = u(iena,:,:,:)
!    wa            = wa - 0.5*(sum(u(imxa:imza,:,:,:)**2,1) /u(idna,:,:,:))
    wa            = wa - 0.5 * u(imxa,:,:,:)**2 / u(idna,:,:,:)
    wa            = wa - 0.5 * u(imya,:,:,:)**2 / u(idna,:,:,:)
    wa            = wa - 0.5 * u(imza,:,:,:)**2 / u(idna,:,:,:)
    wa            = wa - b(1,:,:,:)**2
    wa            = wa - b(2,:,:,:)**2
    wa            = wa - b(3,:,:,:)**2
    wa            = max(wa,smallei)
    wa            = (gamma-1)*wa                    ! pres

    pres_min      = minval(wa(is:ie,js:je,ks:ke))
    loc_pres_min  = minloc(wa(is:ie,js:je,ks:ke)) &
                     + (/nb,nb,nb/)
    call mpifind(pres_min, 'min', loc_pres_min, proc_pres_min)

    pres_max      = maxval(wa(is:ie,js:je,ks:ke))
    loc_pres_max  = maxloc(wa(is:ie,js:je,ks:ke)) &
                     + (/nb,nb,nb/)
    call mpifind(pres_max, 'max', loc_pres_max, proc_pres_max)

    temp_max      = maxval( hydro_mass / k_B * wa(is:ie,js:je,ks:ke) &
                                             /u(idna,is:ie,js:je,ks:ke))
    loc_temp_max  = maxloc(wa(is:ie,js:je,ks:ke)    &
                         /u(idna,is:ie,js:je,ks:ke)  ) &
                     + (/nb,nb,nb/)
    call mpifind(temp_max, 'max', loc_temp_max, proc_temp_max)


    temp_min      = minval( hydro_mass / k_B * wa(is:ie,js:je,ks:ke) &
                                             /u(idna,is:ie,js:je,ks:ke))
    loc_temp_min  = minloc(wa(is:ie,js:je,ks:ke)    &
                         /u(idna,is:ie,js:je,ks:ke)  ) &
                     + (/nb,nb,nb/)
    call mpifind(temp_min, 'min', loc_temp_min, proc_temp_min)

    cs_max        = sqrt(maxval(gamma*wa(is:ie,js:je,ks:ke) &
                            /u(idna,is:ie,js:je,ks:ke)))
    loc_cs_max    = maxloc(gamma*wa(is:ie,js:je,ks:ke) &
                            /u(idna,is:ie,js:je,ks:ke)) &
                     + (/nb,nb,nb/)
    call mpifind(cs_max, 'max', loc_cs_max, proc_cs_max)
#endif /* ISO */

#ifdef COOL_HEAT
      call mpifind(eint_src_min, 'min', loc_dt_cool, proc_dt_cool)
      call mpifind(eint_src_max, 'max', loc_dt_heat, proc_dt_heat)
      call mpifind(dt_cool,      'min', loc_dt_cool, proc_dt_cool)
      call mpifind(dt_heat,      'min', loc_dt_heat, proc_dt_heat)
#endif /* COOL_HEAT */
#ifdef RESIST
! Tu trzba sprawdzic czy poprawnie znajdowane jest max i loc dla wielu procesow
      call mpifind(eta_max,      'max', loc_eta_max, proc_eta_max)
#endif /* RESIST */
#ifdef COSM_RAYS
    wa         =  u(iecr,:,:,:)
    encr_min      = minval(wa(is:ie,js:je,ks:ke))
    loc_encr_min  = minloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
    call mpifind(encr_min, 'min', loc_encr_min, proc_encr_min)

    wa         =  u(iecr,:,:,:)
    encr_max      = maxval(wa(is:ie,js:je,ks:ke))
    loc_encr_max  = maxloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
    call mpifind(encr_max, 'max', loc_encr_max, proc_encr_max)
#endif /* COSM_RAYS */

    if(proc .eq. 0)  then

      open(log_lun, file=log_file, position='append')

        write(log_lun,770) 'min(dens)      =', dens_min,  proc_dens_min,  loc_dens_min
        write(log_lun,770) 'max(dens)      =', dens_max,  proc_dens_max,  loc_dens_max
        write(log_lun,770) 'min(temp)      =', temp_min,  proc_temp_min,  loc_temp_min
        write(log_lun,770) 'max(temp)      =', temp_max,  proc_temp_max,  loc_temp_max
        write(log_lun,770) 'min(pres)      =', pres_min,  proc_pres_min,  loc_pres_min
        write(log_lun,770) 'max(pres)      =', pres_max,  proc_pres_max,  loc_pres_max
        write(log_lun,770) 'min(|b|)       =', b_min,     proc_b_min,     loc_b_min
        write(log_lun,770) 'max(|b|)       =', b_max,     proc_b_max,     loc_b_max
        write(log_lun,770) 'max(|divb|)    =', divb_max,  proc_divb_max,  loc_divb_max

        write(log_lun,777) 'max(|velx|)    =', vx_max, 'dt=',cfl*dx/(vx_max+small),   proc_vx_max, loc_vx_max
        write(log_lun,777) 'max(|vely|)    =', vy_max, 'dt=',cfl*dy/(vy_max+small),   proc_vy_max, loc_vy_max
        write(log_lun,777) 'max(|velz|)    =', vz_max, 'dt=',cfl*dz/(vz_max+small),   proc_vz_max, loc_vz_max
        write(log_lun,777) 'max(v_alfven)  =', va_max, 'dt=',cfl*dxmn/(va_max+small), proc_va_max, loc_va_max
        write(log_lun,777) 'max(c_sound)   =', cs_max, 'dt=',cfl*dxmn/(cs_max+small), proc_cs_max, loc_cs_max
        write(log_lun,777) 'max(c_fast)    =', sqrt(cs_max**2+va_max**2), 'dt=',cfl*dxmn/sqrt(cs_max**2+va_max**2)
#ifdef COOL_HEAT
        write(log_lun,777) 'min(esrc/eint) =', eint_src_min , 'dt=',dt_cool,   proc_dt_cool,  loc_dt_cool
        write(log_lun,777) 'max(esrc/eint) =', eint_src_max , 'dt=',dt_heat,   proc_dt_heat,  loc_dt_heat
#endif /* COOL_HEAT */
#ifdef RESIST
        write(log_lun,777) 'max(eta)       =', eta_max ,      'dt=',dt_resist, proc_eta_max,  loc_eta_max
#endif /* RESIST */
#ifdef COSM_RAYS
        write(log_lun,777) 'min(encr)      =', encr_min,         '',  0.0,     proc_encr_min, loc_encr_min
        write(log_lun,777) 'max(encr)      =', encr_max,      'dt=',dt_cr,     proc_encr_max, loc_encr_max
#endif /* COSM_RAYS */

        write(log_lun,900) nstep,dt,t

      close(log_lun)

    endif


770 format(5x,a16,(1x,e15.9),4(1x,i4))

777 format(5x,a16,(1x,e10.4),2x,a3,(1x,e10.4),4(1x,i4))
900 format('   nstep = ',i7,'   dt = ',f22.16,'   t = ',f22.16,2(1x,i4))

  end subroutine write_log


!
!=======================================================================
!
!                    B E G I N   S U B R O U T I N E
!                           C H E C K _ D I S K
!
!=======================================================================
!
  subroutine check_disk

!     written by: michal hanasz
!     date:       26. june 2003
      use start, only : sleep_minutes,sleep_seconds, min_disk_space_MB
      implicit none

      integer scstatus,ispace,islash,i,i1,ldfout,lcwd
      logical diskfree,diskaccess
      integer min_disk_space, darray(3),tarray(3)
      character df_file*32,cwd_file*32,fsys*30,dfout*80,cspace*10,cpid*10,cproc*2
      character host_small_space*8, dfout_small_space*80
      character*160 syscom
      integer date_time (8) , sleep_time
      character (len = 10) big_ben (3)

      integer isend(2), irecv(2), loc_proc
      integer tsleep

      integer(kind=1) :: system

!-------------------------------------------------------------------------
!     Configurable parameters in problem.par
!-------------------------------------------------------------------------
!      min_disk_space_MB     ! Minimum free space on the working directory
!      sleep_minutes         ! Waiting time in minutes
!      sleep_seconds         ! Waiting time in minutes
!-------------------------------------------------------------------------
      min_disk_space = 1000*min_disk_space_MB

      write (cpid,'(i10)') pid
      write(cproc,'(i2.2)') proc
      df_file = trim(cpid)//'_'//trim(host)//'_'//cproc//'_'//'df.tmp'
      cwd_file= trim(cpid)//'_'//trim(host)//'_'//cproc//'_'//'cwd.tmp'

      diskfree   = .true.
      diskaccess = .true.
      sleep_time = 60*sleep_minutes+sleep_seconds

111   continue

      lcwd=LEN_TRIM(cwd)
      islash = INDEX(cwd(2:lcwd),'/')
      fsys=cwd(2:islash)

!---  check disk access - writing and deleting cwd_file

      call rm_file(cwd_file)
      open(99,file=cwd_file,status='new',err=221)
      close(99)
      call rm_file(cwd_file)

      goto 222
221   continue
        diskaccess=.false.
222   continue

      if(diskaccess) then
        tsleep = 0
      else
        tsleep = sleep_time
      endif

!---  broadcasting a message (tsleep) on disk inaccessibility

      isend(1) = tsleep
      isend(2) = proc

      call MPI_REDUCE(isend, irecv, 1, MPI_2INTEGER, &
                                       MPI_MINLOC, 0, comm, ierr)
      if(proc .eq.0) then
        tsleep   = irecv(1)
        loc_proc = irecv(2)
      endif

      call MPI_BCAST(loc_proc, 1, MPI_INTEGER, 0, comm, ierr)
      call MPI_BCAST(tsleep,   1, MPI_INTEGER, 0, comm, ierr)

      if(tsleep .gt. 0) then
        if(proc .eq. 0) then
          write(*,*)
          write(*,18) (darray(i),i=1,3),(tarray(i),i=1,3)
          write(*,*) 'WORKING DIRECTORY UNAVAILABLE, proc =', loc_proc
          write(*,*) 'WAITING ',sleep_minutes,' min', sleep_minutes,' sec'
        endif
        call sleep(tsleep)
        goto 111

      endif

!---  disk access restored
!------------------------------------------------------------------------


!--- check disk space on working partition

!     delete df_file if exists
      call rm_file(df_file)

      syscom='df /'//fsys// '| grep '//fsys// '> '//df_file
      scstatus = system(syscom)


      open(81,file=df_file,status='old',err=888)
        read(81,15,err=888,end=888) dfout
      close(81)

!     delete df_file
      call rm_file(df_file)

      ldfout = LEN_TRIM(dfout)
      cspace = dfout(41:50)

      do i=1,10
        i1 = index(cspace(i:10),' ')
        if(i1 .eq. 0) then
          ispace = ctoi(cspace,i,10)
          goto 777
        endif
      enddo
777   continue

!     find a minimum of disk space, process number, machine and broadcast

      isend(1) = ispace
      isend(2) = proc

      call MPI_REDUCE(isend, irecv, 1, MPI_2INTEGER, &
                                       MPI_MINLOC, 0, comm, ierr)
      if(proc .eq.0) then
        ispace   = irecv(1)
        loc_proc = irecv(2)
      endif

      call MPI_BCAST(loc_proc, 1, MPI_INTEGER, 0, comm, ierr)
      call MPI_BCAST(ispace, 1, MPI_INTEGER, 0, comm, ierr)

      if(ispace .lt. min_disk_space) then

        diskfree = .false.
        wait = .true.

        if(loc_proc .ne. 0) then
          if(proc .eq. loc_proc) then
            call MPI_SEND  (host,   8, MPI_CHARACTER,        0, 15, comm, ierr)
          else if(proc .eq. 0) then
            call MPI_RECV  ( host_small_space,  8, MPI_CHARACTER, loc_proc, 15, comm, status, ierr)
          endif

          if(proc .eq. loc_proc) then
            call MPI_SEND  (dfout, 80, MPI_CHARACTER,        0, 16, comm, ierr)
          else if(proc .eq. 0) then
            call MPI_RECV  (dfout_small_space, 80, MPI_CHARACTER, loc_proc, 16, comm, status, ierr)
          endif
        else
          host_small_space = host
          dfout_small_space = dfout
        endif

        if(proc .eq. 0) then

!--- inform user and write to logfile

          call date_and_time (big_ben (1), big_ben (2), big_ben (3), date_time)
          darray(1) = date_time(3)
          darray(2) = date_time(2)
          darray(3) = date_time(1)
          tarray(1) = date_time(5)
          tarray(2) = date_time(6)
          tarray(3) = date_time(7)

          write(*,*)
          write(*,18) (darray(i),i=1,3),(tarray(i),i=1,3)
          write(*,*)'!DISK SPACE TOO SMALL:   ', host_small_space
          write(*,17) dfout_small_space(1:79)
          write(*,13)
          write(*,14) min_disk_space
          write(*,*) 'WAITING ',sleep_minutes,' min', sleep_seconds,' sec'

          open(log_lun, file=log_file, position='append')
          write(log_lun,*)
          write(log_lun,18) (darray(i),i=1,3),(tarray(i),i=1,3)
          write(log_lun,*) '!DISK SPACE TOO SMALL:   ', host_small_space
          write(log_lun,17) dfout_small_space(1:79)
          write(log_lun,13)
          write(log_lun,14) min_disk_space
          write(log_lun,*)
          close(log_lun)
        endif !(proc .eq. 0)

          call sleep(sleep_time)

!---      go back and check again
!          goto 111

      else !(ispace .lt. min_disk_space)
        wait = .false.

! inform if disk space is restored

        if(diskfree .eqv. .false.) then
          if(proc .eq. 0) then
            write(*,*)
            write(*,18) (darray(i),i=1,3),(tarray(i),i=1,3)
            write(*,*)'!SUFFICIENT DISK SPACE: ', host_small_space
            write(*,17) dfout(1:79)
            write(*,13)
            write(*,*) 'CONTINUING '

            open(log_lun, file=log_file, position='append')
            write(log_lun,*)
            write(log_lun,18) (darray(i),i=1,3),(tarray(i),i=1,3)
            write(log_lun,*)'!SUFFICIENT DISK SPACE: ', host_small_space
            write(log_lun,17) dfout(1:79)
            write(log_lun,13)
            write(log_lun,*) 'CONTINUING '
            write(log_lun,*)
            close(log_lun)
          endif
          diskfree = .true.
        endif
      endif


888   continue

13    format(42x,'----------')
14    format(23x,'SHOULD BE AT LEAST:',i9)
15    format(a80)
17    format(1x,a79)
18    format(i2.2,'.',i2.2,'.',i4.4,', ',i2.2,':',i2.2,':',i2.2)



    end subroutine check_disk
!=======================================================================
!
!                     E N D   S U B R O U T I N E
!                             C H E C K D F
!
!=======================================================================
!
!=======================================================================
!
!                    B E G I N   S U B R O U T I N E
!                             R E A D F M S G
!
!=======================================================================
!
    subroutine read_file_msg

      use start, only : user_message_file, system_message_file


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

!
!=======================================================================
!
!     \\\\\\\\\\     B E G I N   S U B R O U T I N E     //////////
!     //////////                 C T O R                 \\\\\\\\\\
!
!=======================================================================
!
       subroutine ctor ( cstr, nch, rval, ierr )
!
!    dac:zeus3d.ctor <-------------- translates character string to real
!                                                          january, 1992
!
!      Written by: David Clarke
!      Modified 1:
!
!  PURPOSE:  This subroutine translates a character string to a real
!  constant.
!
!  INPUT VARIABLES:
!    cstr      character representation of real number
!    nch       number of characters in "cstr"
!
!  OUTPUT VARIABLES
!    rval      floating representation of real number
!    ierr      0 => floating conversion successful
!              1 => non-translatable character found in "cstr"
!
!  EXTERNALS: [NONE]
!
!-----------------------------------------------------------------------
!
!*ca imp
!
       implicit none
       character*1   char1
       character*(*) cstr
       integer       nch     , ierr    , i       , esign   , ipnt &
                   , npnt    , more    , iman    , iexp    , j
       REAL        rval    , msign   , term    , fact

       character*1   ch      (  10)
       integer       exponent(   4), mantissa(  16)
!
!      Data statements
!
       data          ch       /'0','1','2','3','4','5','6','7','8','9'/
!
!-----------------------------------------------------------------------
!
!      Initialise parameters.
!
       do 10 i=1,4
         exponent(i) = 0
10     continue
       do 20 i=1,16
         mantissa(i) = 0
20     continue
       rval  = 0.0
       msign = 1.0
       esign = 1
       ipnt  = 0
       npnt  = 0
       more  = 1
       iman  = 0
       iexp  = 0
       ierr  = 0
!
!      Scan "cstr".
!
       do 40 i=1,nch
         char1 = cstr(i:i)
         if (char1 .eq. ' ') go to 40
         if (char1 .eq. '-') then
           if (more .eq. 1) msign = -1.0
           if (more .eq. 2) esign = -1
           go to 40
         endif
         if (char1 .eq. '+') go to 40
         if (char1 .eq. '.') then
           npnt = 1
           ipnt = i - 2
           if (msign .eq. -1.0) ipnt = ipnt - 1
           go to 40
         endif
         if ( (char1 .eq. 'e') .or. (char1 .eq. 'd') ) then
           more = 2
           if (npnt .eq. 0) then
             npnt = 1
             ipnt = i - 2
             if (msign .eq. -1.0) ipnt = ipnt - 1
             if (ipnt .lt. 0) ipnt = 0
           endif
           go to 40
         endif
         do 30 j=1,10
           if (char1 .eq. ch(j)) then
             if (more .eq. 1) then
               iman = iman + 1
               mantissa(iman) = j - 1
             endif
             if (more .eq. 2) then
               iexp = iexp + 1
               exponent(iexp) = j - 1
             endif
             go to 40
           endif
30       continue
!
!      Illegal character found.  Set "ierr" to 1 and return.
!
         ierr = 1
         return
40     continue
       if (npnt .eq. 0) ipnt = iman - 1
!
!      Evaluate the mantissa.
!
       if (iman .gt. 0) then
         do 50 i=1,iman
           term =  real(mantissa(i)) * 10.0**ipnt
           rval = rval + term
           ipnt = ipnt - 1
50       continue
       else
         rval = 1.0
       endif
       rval = rval * msign
!
!      Evaluate the exponential factor.
!
       fact = 1.0
       if (iexp .gt. 0) then
         do 60 i=1,iexp
           ipnt = iexp - i
           fact = fact * 10.0 ** (exponent(i) * 10.0**ipnt)
60       continue
       endif
       if (esign .eq. -1) fact = 1.0 / fact
!
!      Return real representation of character string.
!
       rval = rval * fact

       return
       end subroutine ctor

!=======================================================================
!
!     \\\\\\\\\\       E N D   S U B R O U T I N E       //////////
!     //////////                 C T O R                 \\\\\\\\\\
!
!=======================================================================
!
#ifdef BLOB
!-----------------------------------------------------------------------
! estimation of the fraction of cloud matter is not destructed
!-----------------------------------------------------------------------

  subroutine write_blob

    use arrays, only : nxb,nyb,nzb,u
    use constants, only : pi
    use grid, only  : dx,dy,dz
    use start, only : gamma,t,nstep
    use init_problem, only : chi,denv,tkh,Mext,rblob

    implicit none
    integer i,j,k
    real mass, eint, csirel, masssum


   if(proc .eq. 0) open(5,file='cloudfrac.out',position='append')
   mass = 0.0
   do i=1,nxb
   do j=1,nyb
   do k=1,nzb
     if(u(1,i,j,k) .gt. 0.64*chi*denv) then
       eint=u(5,i,j,k)-0.5*(u(2,i,j,k)**2+u(3,i,j,k)**2+u(4,i,j,k)**2)/u(1,i,j,k)
       csirel =sqrt(eint*(gamma-1)*tkh*Mext*gamma/3.2/rblob/sqrt(chi)/u(1,i,j,k))
       if(csirel .lt. 0.9) then
         mass = mass + u(1,i,j,k)*dx*dy*dz
       endif
     endif
   enddo
   enddo
   enddo

   call MPI_REDUCE(mass,masssum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,comm,ierr)
   if(proc .eq. 0) write(5,555) nstep,t,t/tkh,masssum,masssum/(4./3.*pi*rblob**3*chi*denv)
555 format(i5,4(1x,e20.10))
   if(proc .eq. 0) close(5)

  end subroutine write_blob
#endif /* BLOB */


end module dataio




