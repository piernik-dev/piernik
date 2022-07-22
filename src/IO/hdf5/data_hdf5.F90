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
!    Initial implementation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!

#include "piernik.h"

!>
!! \brief Module that contains HDF5 I/O routines for writing single-precision data dumps
!<
module data_hdf5
! pulled by HDF5

   implicit none

   private
   public :: init_data, write_hdf5, gdf_translate

   interface
      subroutine h5_write(sequential)
         implicit none
         logical, intent(in) :: sequential
      end subroutine h5_write
   end interface

   procedure(h5_write), pointer :: write_hdf5 => NULL() !h5_write_to_single_file

contains

   subroutine init_data

      use dataio_pub, only: multiple_h5files, use_v2_io, warn
      use mpisetup,   only: master

      implicit none

      if (multiple_h5files) then
         write_hdf5 => h5_write_to_multiple_files
         if (use_v2_io .and. master) call warn('[data_hdf5:init_data] v2 I/O format for multiple h5 files not available')
      else
         write_hdf5 => h5_write_to_single_file
      endif

   end subroutine init_data

!>
!! \brief Set up unit labels and cgs coefficients for standard fields.
!!
!! \todo Provide user hook for defining unit labels and cgs coefficients.
!<

   function datafields_descr(var) result(f)

      use constants, only: fpi
      use gdf,       only: gdf_field_type
      use units,     only: cm, erg, gram, sek, miu0

      implicit none

      character(len=*), intent(in) :: var
      type(gdf_field_type)         :: f

      f%f2cgs = 1.0
      f%stag  = 0
      f%fn    = trim(var)
      f%fu    = 'fixme'
      select case (trim(var))
         case ("dend", "deni", "denn")
            f%fu = "\rm{g}/\rm{cm}^3"
            f%f2cgs = 1.0 / (gram/cm**3)
         case ("vlxd", "vlxn", "vlxi", "vlyd", "vlyn", "vlyi", "vlzd", "vlzn", "vlzi", "v", "c_s", "cs")
            f%fu = "\rm{cm}/\rm{s}"
            f%f2cgs = 1.0 / (cm/sek)
         case ("enen", "enei")
            f%fu = "\rm{erg}/\rm{cm}^3"
            f%f2cgs = 1.0 / (erg/cm**3)
         case ("ethn", "ethi")
            f%fu = "\rm{erg}/\rm{g}"
            f%f2cgs = 1.0 / (erg/gram)
         case ("pren", "prei")
            f%fu = "\rm{g}/\rm{cm}/\rm{s}^2"
            f%f2cgs = 1.0 / (gram/cm/sek**2)
         case ("temn", "temi")
            f%fu = "\rm{K}"
            f%f2cgs = 1.0
         case ("magx", "magy", "magz", "magB")
            f%fu = "\rm{Gs}"
            f%f2cgs = 1.0 / (fpi * sqrt(cm / (miu0 * gram)) * sek)
            f%stag = 1
         case ("divbc", "divbf", "divbc4", "divbf4", "divbc6", "divbf6", "divbc8", "divbf8")
            f%fu= "\rm{Gs}/\rm{cm}" ! I'm not sure if it is a best description
            f%f2cgs = 1.0 / (fpi * sqrt(cm / (miu0 * gram)) * sek * cm)
         case ("magdir")
            f%fu = "\rm{radians}"
#ifdef COSM_RAYS
         case ("cr01" : "cr99")
            f%fu = "\rm{erg}/\rm{cm}^3"
            f%f2cgs = 1.0 / (erg/cm**3)
#endif /* COSM_RAYS */
#ifdef CRESP
         case ("cr_e-n01" : "cr_e-n99")
             f%fu = "1/\rm{cm}^3"
             f%f2cgs = 1.0 / (1.0/cm**3) ! number density
         case ("cr_e-e01" : "cr_e-e99")
             f%fu = "\rm{erg}/\rm{cm}^3"
             f%f2cgs = 1.0 / (erg/cm**3)

         case ("cr_p+n01" : "cr_p+n99")
            f%fu = "1/\rm{cm}^3"
            f%f2cgs = 1.0 / (1.0/cm**3) ! number density
         case ("cr_p+e01" : "cr_p+e99")
            f%fu = "\rm{erg}/\rm{cm}^3"
            f%f2cgs = 1.0 / (erg/cm**3)

         case ("cr_C12n01" : "cr_C12n99")
            f%fu = "1/\rm{cm}^3"
            f%f2cgs = 1.0 / (1.0/cm**3) ! number density
         case ("cr_C12e01" : "cr_C12e99")
            f%fu = "\rm{erg}/\rm{cm}^3"
            f%f2cgs = 1.0 / (erg/cm**3)

         case ("cr_N14n01" : "cr_N14n99")
            f%fu = "1/\rm{cm}^3"
            f%f2cgs = 1.0 / (1.0/cm**3) ! number density
         case ("cr_N14e01" : "cr_N14e99")
            f%fu = "\rm{erg}/\rm{cm}^3"
            f%f2cgs = 1.0 / (erg/cm**3)

         case ("cr_O16n01" : "cr_O16n99")
            f%fu = "1/\rm{cm}^3"
            f%f2cgs = 1.0 / (1.0/cm**3) ! number density
         case ("cr_O16e01" : "cr_O16e99")
            f%fu = "\rm{erg}/\rm{cm}^3"
            f%f2cgs = 1.0 / (erg/cm**3)

         case ("cr_Li7n01" : "cr_Li7n99")
            f%fu = "1/\rm{cm}^3"
            f%f2cgs = 1.0 / (1.0/cm**3) ! number density
         case ("cr_Li7e01" : "cr_Li7e99")
            f%fu = "\rm{erg}/\rm{cm}^3"
            f%f2cgs = 1.0 / (erg/cm**3)

         case ("cr_Be9n01" : "cr_Be9n99")
            f%fu = "1/\rm{cm}^3"
            f%f2cgs = 1.0 / (1.0/cm**3) ! number density
         case ("cr_Be9e01" : "cr_Be9e99")
            f%fu = "\rm{erg}/\rm{cm}^3"
            f%f2cgs = 1.0 / (erg/cm**3)

         case ("cr_Be10n01" : "cr_Be10n99")
            f%fu = "1/\rm{cm}^3"
            f%f2cgs = 1.0 / (1.0/cm**3) ! number density
         case ("cr_Be10e01" : "cr_Be10e99")
            f%fu = "\rm{erg}/\rm{cm}^3"
            f%f2cgs = 1.0 / (erg/cm**3)





         case ("cref01" : "cref99")
            f%fu = "\rm{s}^3/\rm{g}^2\rm{cm}^6"
            f%f2cgs = sek**3 / gram**2 / cm**6
         case ("crep01" : "crep02")     ! dimensionless, p treated as Lorentz's gamma
            f%fu = ""
            f%f2cgs = 1.0
         case ("creq01" : "creq99")
            f%fu = ""                   ! dimensionless q
            f%f2cgs = 1.0
#endif /* CRESP */
         case ("gpot", "sgpt")
            f%fu = "\rm{cm}^2 / \rm{s}^2"
            f%f2cgs = 1.0 / (cm**2 / sek**2)
         case ("trcr")
      end select
   end function datafields_descr

   elemental function gdf_translate(var) result(newname)

      use constants,  only: dsetnamelen
      use dataio_pub, only: gdf_strict

      implicit none

      character(len=*), intent(in) :: var
      character(len=dsetnamelen)   :: newname

      if (gdf_strict) then
         select case (trim(var))
            case ("dend", "deni", "denn")
               newname = "density"
            case ("vlxd", "vlxn", "vlxi", "vlyd", "vlyn", "vlyi", "vlzd", "vlzn", "vlzi")
               write(newname, '("velocity_",A1)') var(3:3)
            case ("momxd", "momxn", "momxi", "momyd", "momyn", "momyi", "momzd", "momzn", "momzi")
               write(newname, '("momentum_",A1)') var(4:4)
            case ("enen", "enei")
               newname = "energy_density"
            case ("ethn", "ethi")
               newname = "specific_energy"
            case ("pren", "prei")
               newname = "pressure"
            case ("temn", "temi")
               newname = "temperature"
            case ("magx", "magy", "magz")
               write(newname, '("mag_field_",A1)') var(4:4)
            case ("divbc", "divbf")
               write(newname, '("magnetic_field_divergence_",A1)') var(5:5)
            case ("divbc4", "divbf4")
               write(newname, '("magnetic_field_divergence_",A1,"_O(4)")') var(5:5)
            case ("divbc6", "divbf6")
               write(newname, '("magnetic_field_divergence_",A1,"_O(6)")') var(5:5)
            case ("divbc8", "divbf8")
               write(newname, '("magnetic_field_divergence_",A1,"_O(8)")') var(5:5)
            case ("pmag%")
               newname = "p_mag_to_p_tot_ratio"
            case ("magB")
               newname = "magnetic_field_magnitude"
            case ("magdir")
               newname = "magnetic_field_direction"
            case ("v")
               newname = "total_velocity"
            case ("cs", "c_s")
#ifdef MAGNETIC
               newname = "fast_magnetosonic_speed"
#else /* !MAGNETIC */
               newname = "sound_speed"
#endif /* MAGNETIC */
            case ("Mach", "mach")
#ifdef MAGNETIC
               newname = "Mach_number_(fast)"
#else /* !MAGNETIC */
               newname = "Mach_number"
#endif /* MAGNETIC */
#ifdef NBODY
            case ("ppid")
               newname="id"
            case ("ener")
               newname="energy"
            case ("posx")
               newname="position_x"
            case ("posy")
               newname="position_y"
            case ("posz")
               newname="position_z"
            case ("velx")
               newname="velocity_x"
            case ("vely")
               newname="velocity_y"
            case ("velz")
               newname="velocity_z"
            case ("accx")
               newname="acceleration_x"
            case ("accy")
               newname="acceleration_y"
            case ("accz")
               newname="acceleration_z"
#endif /* NBODY */
            case default
               write(newname, '(A)') trim(var)
         end select
      else
         write(newname, '(A)') trim(var)
      endif
   end function gdf_translate

   subroutine create_units_description(gid)

      use common_hdf5,  only: hdf_vars, hdf_vars_avail
      use constants,    only: units_len, cbuff_len, I_FIVE, I_ONE
      use hdf5,         only: HID_T, h5dopen_f, h5dclose_f
      use helpers_hdf5, only: create_dataset, create_attribute
      use units,        only: lmtvB, s_lmtvB, get_unit

      implicit none

      integer(HID_T), intent(in)             :: gid
      integer(HID_T)                         :: dset_id
      integer(kind=4)                        :: error, i, ip  !< error perhaps should be of type integer(HID_T)
      character(len=cbuff_len), pointer      :: ssbuf
      character(len=units_len), pointer      :: sbuf
      character(len=units_len), target       :: s_unit
      real                                   :: val_unit

      character(len=cbuff_len), dimension(I_FIVE), parameter :: base_dsets = &
         &  ["length_unit  ", "mass_unit    ", "time_unit    ",  "velocity_unit", "magnetic_unit"]

      do i = lbound(base_dsets, 1), ubound(base_dsets, 1)
         call create_dataset(gid, base_dsets(i), lmtvB(i))
         call h5dopen_f(gid, base_dsets(i), dset_id, error)
         ssbuf => s_lmtvB(i)
         call create_attribute(dset_id, "unit", ssbuf)
         call h5dclose_f(dset_id, error)
      enddo
      ip = lbound(base_dsets, 1) - 1
      do i = lbound(hdf_vars, 1, kind=4), ubound(hdf_vars, 1, kind=4)
         if (.not.hdf_vars_avail(i)) cycle
         ip = ip + I_ONE
         call get_unit(gdf_translate(hdf_vars(i)), val_unit, s_unit)
         call create_dataset(gid, gdf_translate(hdf_vars(i)), val_unit)
         call h5dopen_f(gid, gdf_translate(hdf_vars(i)), dset_id, error)
         sbuf => s_unit
         call create_attribute(dset_id, "unit", sbuf)
         call h5dclose_f(dset_id, error)
      enddo

   end subroutine create_units_description

   subroutine create_datafields_descrs(place)

      use common_hdf5,  only: hdf_vars, hdf_vars_avail
      use gdf,          only: gdf_field_type, fmax
      use hdf5,         only: HID_T, h5gcreate_f, h5gclose_f
      use helpers_hdf5, only: create_attribute

      implicit none

      integer(HID_T), intent(in)             :: place

      integer                                :: i, ip
      integer(kind=4)                        :: error  !< error perhaps should be of type integer(HID_T)
      integer(HID_T)                         :: g_id
      type(gdf_field_type), target           :: f
      integer(kind=8), pointer, dimension(:) :: ibuf
      character(len=fmax), pointer           :: sbuf

      allocate(ibuf(1))
      ip = lbound(hdf_vars,1) - 1
      do i = lbound(hdf_vars,1), ubound(hdf_vars,1)
         if (.not.hdf_vars_avail(i)) cycle
         ip = ip + 1
         f = datafields_descr(hdf_vars(i))
         call h5gcreate_f(place, gdf_translate(hdf_vars(i)), g_id, error)
         call create_attribute(g_id, 'field_to_cgs', [f%f2cgs])
         ibuf = f%stag
         call create_attribute(g_id, 'staggering',   ibuf)
         sbuf => f%fu
         call create_attribute(g_id, 'field_units',  sbuf)
         sbuf => f%fn
         call create_attribute(g_id, 'field_name',   sbuf)
         call h5gclose_f(g_id, error)
      enddo
      deallocate(ibuf)
   end subroutine create_datafields_descrs
!>
!! \brief Routine calculating quantities for .hdf files
!<
   subroutine datafields_hdf5(var, tab, ierrh, cg)

      use common_hdf5,      only: common_shortcuts, hdf_vars
      use constants,        only: dsetnamelen, I_ONE
      use fluids_pub,       only: has_ion, has_neu, has_dst
      use fluidindex,       only: flind
      use fluidtypes,       only: component_fluid
      use func,             only: ekin, emag, sq_sum3
      use grid_cont,        only: grid_container
      use mpisetup,         only: proc
#ifdef MAGNETIC
      use constants,        only: xdim, ydim, zdim, half, two, I_TWO, I_FOUR, I_SIX, I_EIGHT
      use div_B,            only: divB_c_IO
      use domain,           only: dom
      use global,           only: cc_mag

#endif /* MAGNETIC */
#ifdef CRESP
      use initcrspectrum,   only: dfpq
      use named_array_list, only: wna
#endif /* CRESP */
#ifdef COSM_RAYS
      use cr_data,          only: cr_names, cr_spectral
#endif /* COSM_RAYS */
#ifndef ISO
      use units,            only: kboltz, mH
#endif /* !ISO */

      implicit none

      character(len=dsetnamelen),      intent(in)    :: var
      real, dimension(:,:,:), pointer, intent(inout) :: tab
      integer,                         intent(out)   :: ierrh
      type(grid_container),   pointer, intent(in)    :: cg

      class(component_fluid), pointer                :: fl_dni, fl_mach
      integer(kind=4)                                :: i_xyz, clast
      integer                                        :: ii, jj, kk, icr
#ifdef COSM_RAYS
      integer                                        :: i, ibin
      integer, parameter                             :: auxlen = dsetnamelen - 1
      character(len=auxlen)                          :: aux
      character(len=2)                               :: varn2
      !character(len=*)                               :: vname
#endif /* COSM_RAYS */

      call common_shortcuts(var, fl_dni, i_xyz)
      if (.not. associated(fl_dni)) tab = -huge(1.)
      ierrh = 0
      tab = 0.0

#ifdef MAGNETIC
      associate(emag_c => merge(emag(cg%b(xdim, RNG), cg%b(ydim, RNG),  cg%b(zdim, RNG)), &
           &                    emag(half*(cg%b(xdim, RNG) + cg%b(xdim, cg%is+dom%D_x:cg%ie+dom%D_x, cg%js        :cg%je,         cg%ks        :cg%ke        )), &
           &                         half*(cg%b(ydim, RNG) + cg%b(ydim, cg%is        :cg%ie,         cg%js+dom%D_y:cg%je+dom%D_y, cg%ks        :cg%ke        )), &
           &                         half*(cg%b(zdim, RNG) + cg%b(zdim, cg%is        :cg%ie,         cg%js        :cg%je,         cg%ks+dom%D_z:cg%ke+dom%D_z))), &
           &                    cc_mag))  ! fortran way of constructing ternary operators
#else /* !MAGNETIC */
      associate(emag_c => 0.)
#endif /* !MAGNETIC */
      !print *, lbound(hdf_vars,1, kind=4), ubound(hdf_vars,1, kind=4)
      select case (var)
#ifdef COSM_RAYS
         case ("cr01" : "cr99")
            read(var,'(A2,I2.2)') aux, i !> \deprecated BEWARE 0 <= i <= 99, no other indices can be dumped to hdf file
            tab(:,:,:) = cg%u(flind%crn%beg+i-1, RNG)
#endif /* COSM_RAYS */
#ifdef CRESP
         case ('cr_A000' : 'cr_zz99')
            clast = len(trim(var))
            varn2 = var(clast - 1:clast)
            if (var(clast - 2:clast - 2) == 'e') then

            !part of the code for spectrally resolved species : energy density

               read (varn2,'(I2.2)') ibin
               do i = 1, size(cr_names)
                  if (cr_names(i).eq.var(4:clast-3)) icr = i
               enddo
               tab(:,:,:) = 1.0 !cg%u(flind%crspcs(icr)%ebeg+ibin-1, RNG)
               !print *, 'icr : ', icr, ' max val : ', maxval(tab(:,:,:)), 'bin : ', ibin

            else if (var(clast - 2:clast - 2) == 'n') then

            !part of the code for spectrally resolved species : number density

               read (varn2,'(I2.2)') ibin
               do i = 1, size(cr_names)
                  if (cr_names(i).eq.var(4:clast-3)) icr = i
               enddo
               tab(:,:,:) = cg%u(flind%crspcs(icr)%nbeg+ibin-1, RNG)
               !print *, 'icr : ', icr, ' cgu : ', maxval(cg%u(flind%crspcs(icr)%nbeg+ibin-1, RNG)), 'bin : ', ibin

            else
               do i = 1, size(cr_names)
                  if (var == trim('cr_' // cr_names(i))) exit
               enddo
               tab(:,:,:) = cg%u(flind%crn%beg+ibin-1-count(cr_spectral), RNG)
            endif
  !        print *, var, aux
  !        print *, flind%crn%beg+i-1-count(cr_spectral)
         case ("cren01" : "cren99")
            read(var,'(A4,I2.2)') aux, i !> \deprecated BEWARE 0 <= i <= 99, no other indices can be dumped to hdf file
            tab(:,:,:) = cg%u(flind%crspc%nbeg+i-1, RNG)
         case ("cree01" : "cree99")
            read(var,'(A4,I2.2)') aux, i !> \deprecated BEWARE 0 <= i <= 99, no other indices can be dumped to hdf file
            tab(:,:,:) = cg%u(flind%crspc%ebeg+i-1, RNG)
         case ("cref01" : "cref99")
            read(var,'(A4,I2.2)') aux, i !> \deprecated BEWARE 0 <= i <= 99, no other indices can be dumped to hdf file
            tab(:,:,:) = cg%w(wna%ind(dfpq%f_nam))%arr(i,RNG)  !flind%crspc%fbeg+i-1, RNG)
         case ("crep01" : "crep02")
            read(var,'(A4,I2.2)') aux, i !> \deprecated BEWARE 0 <= i <= 99, no other indices can be dumped to hdf file
            tab(:,:,:) = cg%w(wna%ind(dfpq%p_nam))%arr(i,RNG)  !flind%crspc%fbeg+i-1, RNG)
         case ("creq01" : "creq99")
            read(var,'(A4,I2.2)') aux, i !> \deprecated BEWARE 0 <= i <= 99, no other indices can be dumped to hdf file
            tab(:,:,:) = cg%w(wna%ind(dfpq%q_nam))%arr(i,RNG)  !flind%crspc%fbeg+i-1, RNG)
#endif /* CRESP */
#ifdef TRACER
         case ("trcr")
            tab(:,:,:) = cg%u(flind%trc%beg, RNG)
#endif /* TRACER */
         case ("dend", "deni", "denn")
            if (associated(fl_dni)) tab(:,:,:) = cg%u(fl_dni%idn, RNG)
         case ("vlxd", "vlxn", "vlxi", "vlyd", "vlyn", "vlyi", "vlzd", "vlzn", "vlzi")
            if (associated(fl_dni)) tab(:,:,:) = cg%u(fl_dni%imx + i_xyz, RNG) / cg%u(fl_dni%idn, RNG)
         case ("momxd", "momxn", "momxi", "momyd", "momyn", "momyi", "momzd", "momzn", "momzi")
            if (associated(fl_dni)) tab(:,:,:) = cg%u(fl_dni%imx + i_xyz, RNG)
         case ("enen", "enei")
#ifdef ISO
            if (associated(fl_dni)) tab(:,:,:) = ekin(cg%u(fl_dni%imx, RNG), cg%u(fl_dni%imy, RNG), cg%u(fl_dni%imz, RNG), cg%u(fl_dni%idn, RNG))
#else /* !ISO */
            if (associated(fl_dni)) tab(:,:,:) = cg%u(fl_dni%ien, RNG)
#endif /* !ISO */
         case ("pren")
#ifndef ISO
            tab(:,:,:) = flind%neu%gam_1 * (cg%u(flind%neu%ien, RNG) - ekin(cg%u(flind%neu%imx, RNG), cg%u(flind%neu%imy, RNG), cg%u(flind%neu%imz, RNG), cg%u(flind%neu%idn, RNG)))
#endif /* !ISO */
         case ("prei")
#ifndef ISO
            tab(:,:,:) = flind%ion%gam_1 * (cg%u(flind%ion%ien, RNG) - ekin(cg%u(flind%ion%imx, RNG), cg%u(flind%ion%imy, RNG), cg%u(flind%ion%imz, RNG), cg%u(flind%ion%idn, RNG)) - emag_c)
#endif /* !ISO */
         case ("pmag%")
#ifndef ISO
#ifdef IONIZED
            tab(:,:,:) = emag_c / &
                 &      (flind%ion%gam_1 * (cg%u(flind%ion%ien, RNG) - ekin(cg%u(flind%ion%imx, RNG), cg%u(flind%ion%imy, RNG), cg%u(flind%ion%imz, RNG), cg%u(flind%ion%idn, RNG)) - emag_c) + &
                 &       emag_c)
#endif /* IONIZED */
#endif /* !ISO */
         case ("ethn")
#ifndef ISO
            tab(:,:,:) = (cg%u(flind%neu%ien, RNG) - &
                 &       ekin(cg%u(flind%neu%imx, RNG), cg%u(flind%neu%imy, RNG), cg%u(flind%neu%imz, RNG), cg%u(flind%neu%idn, RNG))) /         &
                 &       cg%u(flind%neu%idn, RNG)
#endif /* !ISO */
         case ("ethi")
#ifndef ISO
            tab(:,:,:) = (cg%u(flind%ion%ien, RNG) - &
                 &       ekin(cg%u(flind%ion%imx, RNG), cg%u(flind%ion%imy, RNG), cg%u(flind%ion%imz, RNG), cg%u(flind%ion%idn, RNG)) -          &
                 &       emag_c) / cg%u(flind%ion%idn, RNG)
#endif /* !ISO */
         case ("temn")
#ifndef ISO
            tab(:,:,:) = flind%neu%gam_1 * mH / kboltz * (cg%u(flind%neu%ien, RNG) - &
                 &       ekin(cg%u(flind%neu%imx, RNG), cg%u(flind%neu%imy, RNG), cg%u(flind%neu%imz, RNG), cg%u(flind%neu%idn, RNG))) /         &
                 &       cg%u(flind%neu%idn, RNG)
#endif /* !ISO */
         case ("temi")
#ifndef ISO
            tab(:,:,:) = flind%ion%gam_1 * mH / kboltz * (cg%u(flind%ion%ien, RNG) - &
                 &       ekin(cg%u(flind%ion%imx, RNG), cg%u(flind%ion%imy, RNG), cg%u(flind%ion%imz, RNG), cg%u(flind%ion%idn, RNG)) -          &
                 &       emag_c) / cg%u(flind%ion%idn, RNG)
#endif /* !ISO */
#ifdef MAGNETIC
         case ("magx", "magy", "magz")
            tab(:,:,:) = cg%b(xdim + i_xyz, RNG) ! beware: these are "raw", face-centered. Use them with care when you process plotfiles
         case ("magB")
            tab(:,:,:) = sqrt(two * emag_c)
         case ("magdir")
            tab(:,:,:) =  merge(atan2(cg%b(ydim, RNG), cg%b(xdim, RNG)), &
                 &              atan2(cg%b(ydim, RNG) + cg%b(ydim, cg%is        :cg%ie,         cg%js+dom%D_y:cg%je+dom%D_y, cg%ks        :cg%ke        ), &
                 &                    cg%b(xdim, RNG) + cg%b(xdim, cg%is+dom%D_x:cg%ie+dom%D_x, cg%js        :cg%je,         cg%ks        :cg%ke        )),  &
                 &              cc_mag)
            ! ToDo: magi - inclination
            ! ToDo: curlb - nabla x B
!! ToDo: autodetect centering, add option for dumping both just in case
!! face-centered div(B): RTVD and RIEMANN, both with constrained transport
         case ("divbf")
            tab(:,:,:) = divB_c_IO(cg, I_TWO,  .false.)
         case ("divbf4")
            tab(:,:,:) = divB_c_IO(cg, I_FOUR, .false.)
         case ("divbf6")
            tab(:,:,:) = divB_c_IO(cg, I_SIX,  .false.)
         case ("divbf8")
            tab(:,:,:) = divB_c_IO(cg, I_EIGHT,.false.)
!! cell-centered div(B): RIEMANN with divergence cleaning
         case ("divbc")
            tab(:,:,:) = divB_c_IO(cg, I_TWO,  .true.)
         case ("divbc4")
            tab(:,:,:) = divB_c_IO(cg, I_FOUR, .true.)
         case ("divbc6")
            tab(:,:,:) = divB_c_IO(cg, I_SIX,  .true.)
         case ("divbc8")
            tab(:,:,:) = divB_c_IO(cg, I_EIGHT,.true.)
#endif /* MAGNETIC */
         case ("v") ! perhaps this should be expanded to vi, vn or vd, depending on fluids present
            nullify(fl_mach)
            if (has_ion) then
               fl_mach => flind%ion
            else if (has_neu) then
               fl_mach => flind%neu
            else if (has_dst) then
               fl_mach => flind%dst
            endif
            tab(:,:,:) = sqrt(sq_sum3(cg%u(fl_mach%imx, RNG), cg%u(fl_mach%imy, RNG), cg%u(fl_mach%imz, RNG))) / cg%u(fl_mach%idn, RNG)
         case ("cs", "c_s")
            nullify(fl_mach)
            if (has_ion) then
               fl_mach => flind%ion
            else if (has_neu) then
               fl_mach => flind%neu
            endif
            if (associated(fl_mach)) then
               do kk = cg%ks, cg%ke
                  do jj = cg%js, cg%je
                     do ii = cg%is, cg%ie
                        tab(ii - cg%is + I_ONE, jj - cg%js + I_ONE, kk - cg%ks + I_ONE) = fl_mach%get_cs(ii, jj, kk, cg%u, cg%b, cg%cs_iso2)
                     enddo
                  enddo
               enddo
            else
               tab(:,:,:) = 0.
            endif
         case ("mach", "Mach")
            nullify(fl_mach)
            if (has_ion) then
               fl_mach => flind%ion
            else if (has_neu) then
               fl_mach => flind%neu
            endif
            if (associated(fl_mach)) then
               do kk = cg%ks, cg%ke
                  do jj = cg%js, cg%je
                     do ii = cg%is, cg%ie
                        tab(ii - cg%is + I_ONE, jj - cg%js + I_ONE, kk - cg%ks + I_ONE) = fl_mach%get_mach(ii, jj, kk, cg%u, cg%b, cg%cs_iso2)
                     enddo
                  enddo
               enddo
            else
               tab(:,:,:) = 0.
            endif
         case ("gpot")
            if (associated(cg%gpot)) tab(:,:,:) = cg%gpot(RNG)
         case ("sgpt")
            if (associated(cg%sgp)) tab(:,:,:) = cg%sgp(RNG)
         case ("level")
            tab(:,:,:) = cg%l%id
         case ("grid_id")
            tab(:,:,:) = cg%grid_id
         case ("proc")
            tab(:,:,:) = proc
         case default
            ierrh = -1
      end select
      end associate

   end subroutine datafields_hdf5

!
! ------------------------------------------------------------------------------------
!
   subroutine h5_write_to_single_file(sequential)

      use common_hdf5,     only: dump_announcement, dump_announce_time, set_common_attributes
      use constants,       only: cwdlen, PPP_IO, HDF
      use dataio_pub,      only: nhdf, use_v2_io, last_hdf_time
      use mpisetup,        only: report_to_master, report_string_to_master
      use piernik_mpi_sig, only: sig
      use ppp,             only: ppp_main
#if defined(MULTIGRID) && defined(SELF_GRAV)
      use multigrid_gravity, only: unmark_oldsoln
#endif /* MULTIGRID && SELF_GRAV */
#ifdef NBODY_1FILE
!      use particles_io_hdf5, only: write_nbody_hdf5
#endif /* NBODY_1FILE */

      implicit none

      logical,         intent(in) :: sequential
      character(len=cwdlen)       :: fname
      character(len=*), parameter :: wrd_label = "IO_write_datafile_v1"

      call ppp_main%start(wrd_label, PPP_IO)
      ! Initialize HDF5 library and Fortran interfaces.
      !
      call dump_announcement(HDF, nhdf, fname, last_hdf_time, sequential)

      call set_common_attributes(fname)
      if (use_v2_io) then
         call h5_write_to_single_file_v2(fname)
      else
         call h5_write_to_single_file_v1(fname)
      endif
#if defined(MULTIGRID) && defined(SELF_GRAV)
      call unmark_oldsoln
#endif /* MULTIGRID && SELF_GRAV */

      call dump_announce_time
      call report_to_master(sig%hdf_written, only_master=.True.)
      call report_string_to_master(fname, only_master=.True.)
#ifdef NBODY_1FILE
!      call write_nbody_hdf5(fname)
#endif /* NBODY_1FILE */

      call ppp_main%stop(wrd_label, PPP_IO)

   end subroutine h5_write_to_single_file

   subroutine h5_write_to_single_file_v2(fname)

      use constants,   only: PPP_IO
      use common_hdf5, only: write_to_hdf5_v2, O_OUT
      use gdf,         only: gdf_create_root_group
      use mpisetup,    only: master, piernik_MPI_Barrier
      use ppp,         only: ppp_main

      implicit none

      character(len=*), intent(in) :: fname
      character(len=*), parameter :: wrd_label = "IO_write_datafile_v2"

      call ppp_main%start(wrd_label, PPP_IO)

      call write_to_hdf5_v2(fname, O_OUT, create_empty_cg_datasets_in_output, write_cg_to_output)

      if (master) then
         call gdf_create_root_group(fname, 'field_types', create_datafields_descrs)
         call gdf_create_root_group(fname, 'dataset_units', create_units_description)
      endif
      call piernik_MPI_Barrier

      call ppp_main%stop(wrd_label, PPP_IO)

   end subroutine h5_write_to_single_file_v2

!> \brief Write all grid containers to the file

   subroutine write_cg_to_output(cgl_g_id, cg_n, cg_all_n_b, cg_all_n_o)

      use cg_leaves,   only: leaves
      use cg_list,     only: cg_list_element
      use common_hdf5, only: get_nth_cg, hdf_vars, cg_output, hdf_vars, hdf_vars_avail, enable_all_hdf_var
      use constants,   only: xdim, ydim, zdim, ndims, FP_REAL, PPP_IO, PPP_CG, I_ONE
      use dataio_pub,  only: die, nproc_io, can_i_write, h5_64bit, nstep_start
      use global,      only: nstep
      use grid_cont,   only: grid_container
      use hdf5,        only: HID_T, HSIZE_T, H5T_NATIVE_REAL, H5T_NATIVE_DOUBLE, h5sclose_f, h5dwrite_f, h5sselect_none_f, h5screate_simple_f
      use MPIF,        only: MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, MPI_COMM_WORLD
      use MPIFUN,      only: MPI_Recv, MPI_Send
      use mpisetup,    only: master, FIRST, proc, err_mpi, tag_ub
      use ppp,         only: ppp_main
#ifdef NBODY_1FILE
      use cg_particles_io, only: pdsets, nbody_datafields
      use domain,          only: is_multicg
      use MPIF,            only: MPI_INTEGER, MPI_INTEGER8
      use mpisetup,        only: LAST
      use particle_utils,  only: count_all_particles
#endif /* NBODY_1FILE */

      implicit none

      integer(HID_T),                           intent(in) :: cgl_g_id    !< cg group identifier
      integer(kind=4), dimension(:),   pointer, intent(in) :: cg_n        !< offset for cg group numbering
      integer(kind=4), dimension(:,:), pointer, intent(in) :: cg_all_n_b  !< all cg sizes
      integer(kind=4), dimension(:,:), pointer, intent(in) :: cg_all_n_o  !< all cg sizes, expanded by external boundaries

      integer(HID_T)                                       :: filespace_id, memspace_id
      integer(kind=4)                                      :: error       !< error perhaps should be of type integer(HID_T)
      integer(kind=4), parameter                           :: rank = 3
      integer(HSIZE_T), dimension(:), allocatable          :: dims
      integer(kind=4)                                      :: i, ncg, n
      integer                                              :: ip
      type(grid_container),            pointer             :: cg
      type(cg_list_element),           pointer             :: cgl
      real, dimension(:,:,:),          pointer             :: data_dbl ! double precision buffer (internal default, single precision buffer is the plotfile output default, overridable by h5_64bit)
      type(cg_output)                                      :: cg_desc
#ifdef NBODY_1FILE
      integer(kind=4)                                      :: n_part
      integer(HID_T)                                       :: id
      integer(HID_T), dimension(LAST+1)                    :: tmp
#endif /* NBODY_1FILE */
      character(len=*), parameter :: wrdc_label = "IO_write_data_v2_cg", wrdc1s_label = "IO_write_data_v2_1cg_(serial)", wrdc1p_label = "IO_write_data_v2_1cg_(parallel)"

      call ppp_main%start(wrdc_label, PPP_IO)

      if (nstep - nstep_start < 1) call enable_all_hdf_var  ! just in case things have changed meanwhile

#ifdef NBODY_1FILE
      call cg_desc%init(cgl_g_id, cg_n, nproc_io, gdf_translate(pack(hdf_vars, hdf_vars_avail)), gdf_translate(pdsets))
#else /* !NBODY_1FILE */
      call cg_desc%init(cgl_g_id, cg_n, nproc_io, gdf_translate(pack(hdf_vars, hdf_vars_avail)))
#endif /* !NBODY_1FILE */

      if (cg_desc%tot_cg_n < 1) call die("[data_hdf5:write_cg_to_output] no cg available!")

      ! all arrays are rank 3 here
      allocate(dims(ndims))
      ! Allocate data with the size of first cg
      allocate( data_dbl(cg_all_n_b(xdim, 1), cg_all_n_b(ydim, 1), cg_all_n_b(zdim, 1)) )

      if (nproc_io == 1) then ! perform serial write
         ! write all cg, one by one
         if (cg_desc%tot_cg_n *(ubound(hdf_vars,1, kind=4) + I_ONE) > tag_ub) call die("[data_hdf5:write_cg_to_output] this MPI implementation has too low MPI_TAG_UB attribute")
         do ncg = 1, cg_desc%tot_cg_n
            call ppp_main%start(wrdc1s_label, PPP_IO + PPP_CG)
            dims(:) = [ cg_all_n_b(xdim, ncg), cg_all_n_b(ydim, ncg), cg_all_n_b(zdim, ncg) ]
            call recycle_data(dims, cg_all_n_b, ncg, data_dbl)

            if (master) then
               if (.not. can_i_write) call die("[data_hdf5:write_cg_to_output] Master can't write")

               ip = lbound(hdf_vars,1) - 1
               do i = lbound(hdf_vars,1, kind=4), ubound(hdf_vars,1, kind=4)
                  if (.not.hdf_vars_avail(i)) cycle
                  ip = ip + 1
                  if (cg_desc%cg_src_p(ncg) == proc) then
                     cg => get_nth_cg(cg_desc%cg_src_n(ncg))
                     call get_data_from_cg(hdf_vars(i), cg, data_dbl)
                  else
                     call MPI_Recv(data_dbl(1,1,1), size(data_dbl, kind=4), MPI_DOUBLE_PRECISION, cg_desc%cg_src_p(ncg), ncg + cg_desc%tot_cg_n*i, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err_mpi)
                  endif
                  if (.not.hdf_vars_avail(i)) cycle
                  if (h5_64bit) then
                     call h5dwrite_f(cg_desc%dset_id(ncg, ip), H5T_NATIVE_DOUBLE, data_dbl, dims, error, xfer_prp = cg_desc%xfer_prp)
                  else
                     call h5dwrite_f(cg_desc%dset_id(ncg, ip), H5T_NATIVE_REAL, real(data_dbl, kind=FP_REAL), dims, error, xfer_prp = cg_desc%xfer_prp)
                  endif
               enddo
            else
               if (can_i_write) call die("[data_hdf5:write_cg_to_output] Slave can write")
               if (cg_desc%cg_src_p(ncg) == proc) then
                  cg => get_nth_cg(cg_desc%cg_src_n(ncg))
                  do i = lbound(hdf_vars,1, kind=4), ubound(hdf_vars,1, kind=4)
                     if (hdf_vars_avail(i)) call get_data_from_cg(hdf_vars(i), cg, data_dbl)
                     if (.not.hdf_vars_avail(i)) cycle
                     call MPI_Send(data_dbl(1,1,1), size(data_dbl, kind=4), MPI_DOUBLE_PRECISION, FIRST, ncg + cg_desc%tot_cg_n*i, MPI_COMM_WORLD, err_mpi)
                  enddo
               endif
            endif
            !Serial write for particles
            call ppp_main%stop(wrdc1s_label, PPP_IO + PPP_CG)
         enddo
#ifdef NBODY_1FILE
         n_part = count_all_particles()
         if (n_part .gt. 0) then
            if (is_multicg) call die("[data_hdf5:write_cg_to_output] multicg is not implemented for NBODY_1FILE")
            do i=lbound(pdsets, dim=1, kind=4), ubound(pdsets, dim=1, kind=4)
               tmp(:) = 0
               id=0
               if (master) then
                  tmp(:) = cg_desc%pdset_id(:, i)
               endif
               if (kind(id) == 4) then
                  call MPI_Scatter(tmp, 1, MPI_INTEGER, id, 1, MPI_INTEGER, FIRST, MPI_COMM_WORLD, err_mpi)
               else if (kind(id) == 8) then
                  call MPI_Scatter(tmp, 1, MPI_INTEGER8, id, 1, MPI_INTEGER8, FIRST, MPI_COMM_WORLD, err_mpi)
               else
                  call die("[data_hdf5:write_cg_to_output] no recognized kind of HID_T")
               endif
               if (master) then
                  call nbody_datafields(id, gdf_translate(pdsets(i)), n_part)
               else
                  call nbody_datafields(id, gdf_translate(pdsets(i)), n_part)
               endif
            enddo
         endif
#endif /* NBODY_1FILE */
      else ! perform parallel write
         ! This piece will be a generalization of the serial case. It should work correctly also for nproc_io == 1 so it should replace the serial code
         if (can_i_write) then
            ! write own
            n = 0
            cgl => leaves%first
            do while (associated(cgl))
               call ppp_main%start(wrdc1p_label, PPP_IO + PPP_CG)
               n = n + I_ONE
               ncg = cg_desc%offsets(proc) + n
               dims(:) = [ cg_all_n_b(xdim, ncg), cg_all_n_b(ydim, ncg), cg_all_n_b(zdim, ncg) ]
               call recycle_data(dims, cg_all_n_b, ncg, data_dbl)
               cg => cgl%cg

               ip = lbound(hdf_vars,1) - 1
               do i = lbound(hdf_vars,1, kind=4), ubound(hdf_vars,1, kind=4)
                  if (.not.hdf_vars_avail(i)) cycle
                  ip = ip + 1
                  call get_data_from_cg(hdf_vars(i), cg, data_dbl)
                  if (.not.hdf_vars_avail(i)) cycle
                  if (h5_64bit) then
                     call h5dwrite_f(cg_desc%dset_id(ncg, ip), H5T_NATIVE_DOUBLE, data_dbl, dims, error, xfer_prp = cg_desc%xfer_prp)
                  else
                     call h5dwrite_f(cg_desc%dset_id(ncg, ip), H5T_NATIVE_REAL, real(data_dbl, kind=FP_REAL), dims, error, xfer_prp = cg_desc%xfer_prp)
                  endif
               enddo

#ifdef NBODY_1FILE
               n_part = count_all_particles()
               if (n_part .gt. 0) then
                  do i=lbound(pdsets, dim=1, kind=4), ubound(pdsets, dim=1, kind=4)
                     call nbody_datafields(cg_desc%pdset_id(ncg, i), gdf_translate(pdsets(i)), n_part)
                  enddo
               endif
#endif /* NBODY_1FILE */

               cgl => cgl%nxt
               call ppp_main%stop(wrdc1p_label, PPP_IO + PPP_CG)
            enddo

            ! Behold the MAGIC in its purest form!
            ! Following block of code does exactly *nothing*, yet it's necessary for collective calls of PHDF5
            !>
            !! \deprecated BEWARE, we assume that at least 1 cg exist on a given proc (or at leas we fake it)
            !! \todo there should be something like H5S_NONE as a contradiction to H5S_ALL, yet I cannot find it...
            !<

            dims(:) = [ cg_all_n_b(xdim, 1), cg_all_n_b(ydim, 1), cg_all_n_b(zdim, 1) ]

            ! completely bogus values, only to make HDF5 happy
            call h5screate_simple_f(rank, dims, filespace_id, error)
            call h5screate_simple_f(rank, dims, memspace_id, error)
            ! empty filespace
            call h5sselect_none_f(filespace_id, error)
            ! empty memoryspace
            call h5sselect_none_f(memspace_id, error)

            call recycle_data(dims, cg_all_n_b, I_ONE, data_dbl)
            do ncg = I_ONE, maxval(cg_n)-n
               ip = lbound(hdf_vars,1) - 1
               do i = lbound(hdf_vars, 1, kind=4), ubound(hdf_vars, 1, kind=4)
                  if (.not.hdf_vars_avail(i)) cycle
                  ip = ip + 1
                  ! It is crashing due to FPE when there are more processes than blocks
                  ! because data_dbl contains too large values for single precision.
                  ! If a process doesn't have a block, data_dbl serves just as
                  ! a placeholder to complete collective HDF5 calls.
                  !
                  ! Yes, something stinks here.
                  !
                  ! On uniform grid a process without a cg means that the user made an error and assigned too many processes for too little task.
                  ! In AMR such situation may occur when in a large simulation a massive derefinement happens.
                  ! Usually it will mean that there is something wrong with refinement criteria but still the user
                  ! deserves to get the files, not a FPE crash.
                  if (h5_64bit .or. n < 1) then
                     call h5dwrite_f(cg_desc%dset_id(1, i), H5T_NATIVE_DOUBLE, data_dbl, dims, error, &
                          &          xfer_prp = cg_desc%xfer_prp, file_space_id = filespace_id, mem_space_id = memspace_id)
                  else
                     call h5dwrite_f(cg_desc%dset_id(1, i), H5T_NATIVE_REAL, real(data_dbl, kind=FP_REAL), dims, error, &
                          &          xfer_prp = cg_desc%xfer_prp, file_space_id = filespace_id, mem_space_id = memspace_id)
                  endif
               enddo
            enddo

            call h5sclose_f(memspace_id, error)
            call h5sclose_f(filespace_id, error)
            ! receive (from whom?)
         else
            call die("[data_hdf5:write_cg_to_output] nproc != nproc_io not implemented yet")
            ! send (where?)
         endif
      endif

      ! clean up
      if (allocated(dims)) deallocate(dims)
      if (associated(data_dbl)) deallocate(data_dbl)
      call cg_desc%clean()

      call ppp_main%stop(wrdc_label, PPP_IO)

      if (.false.) i = size(cg_all_n_o, kind=4) ! suppress compiler warning

   contains
      !>
      !! Try to avoid pointless data reallocation for every cg if shape doesn't change
      !<
      subroutine recycle_data(dims, cg_all_n_b, i, data)

         use constants, only: xdim, ydim, zdim
         use hdf5,      only: HSIZE_T

         implicit none

         integer(HSIZE_T), dimension(:)                          :: dims        !< shape of current cg
         integer(kind=4), dimension(:,:), pointer, intent(in)    :: cg_all_n_b  !< all cg sizes
         integer(kind=4),                          intent(in)    :: i           !< no. of cg
         real, dimension(:,:,:), pointer,          intent(inout) :: data        !< temporary storage array used for I/O

         if (associated(data)) then
            if ( any(dims /= shape(data)) ) then
               deallocate(data)
               allocate(data(cg_all_n_b(xdim, i), cg_all_n_b(ydim, i), cg_all_n_b(zdim, i)))
            endif
         endif
      end subroutine recycle_data

   end subroutine write_cg_to_output

   subroutine get_data_from_cg(hdf_var, cg, tab)

      use common_hdf5,      only: cancel_hdf_var
      use dataio_pub,       only: warn, msg, printinfo
      use dataio_user,      only: user_vars_hdf5
      use grid_cont,        only: grid_container
      use mpisetup,         only: master
      use named_array_list, only: qna

      implicit none

      character(len=*),                intent(in)    :: hdf_var
      type(grid_container),   pointer, intent(inout) :: cg
      real, dimension(:,:,:), pointer, intent(inout) :: tab

      integer :: ierrh

      ierrh = 0

      ! Try some default names first
      call datafields_hdf5(hdf_var, tab, ierrh, cg)

      ! Call user routines for user variables or quantities computed in user routines
      if (associated(user_vars_hdf5) .and. ierrh /= 0) call user_vars_hdf5(hdf_var, tab, ierrh, cg)

      ! Check if a given name was registered in named arrays. This is lowest-priority identification.
      if (ierrh /= 0) then  ! All simple scalar named arrays should be handled here
         if (qna%exists(hdf_var)) then
            tab(:,:,:) = real(cg%q(qna%ind(hdf_var))%span(cg%ijkse), kind(tab))
            ierrh = 0
         endif
      endif

      if (ierrh /= 0) then
         write(msg,'(3a)') "[data_hdf5:get_data_from_cg]: '", trim(hdf_var), "' is not recognized as a name of defined variables/fields, not defined in datafields_hdf5 and not found in user_vars_hdf5."
         if (master) then
            call printinfo("", .true.)
            call warn(msg)
         endif
         call cancel_hdf_var(hdf_var)
      endif

   end subroutine get_data_from_cg

#ifdef NBODY_1FILE
   subroutine create_empty_cg_datasets_in_output(cg_g_id, cg_n_b, cg_n_o, Z_avail, n_part, st_g_id)
#else /* !NBODY_1FILE */
   subroutine create_empty_cg_datasets_in_output(cg_g_id, cg_n_b, cg_n_o, Z_avail)
#endif /* !NBODY_1FILE */

      use common_hdf5, only: create_empty_cg_dataset, hdf_vars, hdf_vars_avail, O_OUT
      use hdf5,        only: HID_T, HSIZE_T
#ifdef NBODY_1FILE
      use cg_particles_io, only: pdsets
#endif /* NBODY_1FILE */

      implicit none

      integer(HID_T),                intent(in) :: cg_g_id
      integer(kind=4), dimension(:), intent(in) :: cg_n_b
      integer(kind=4), dimension(:), intent(in) :: cg_n_o
      logical(kind=4),               intent(in) :: Z_avail

      integer :: i
#ifdef NBODY_1FILE
      integer(kind=8)                           :: n_part
      integer(HID_T),                intent(in) :: st_g_id
#endif /* NBODY_1FILE */

      do i = lbound(hdf_vars,1), ubound(hdf_vars,1)
         if (.not.hdf_vars_avail(i)) cycle
         call create_empty_cg_dataset(cg_g_id, gdf_translate(hdf_vars(i)), int(cg_n_b, kind=HSIZE_T), Z_avail, O_OUT)
      enddo

      if (.false.) i = size(cg_n_o) ! suppress compiler warning

#ifdef NBODY_1FILE
      do i = lbound(pdsets,1), ubound(pdsets,1)
         call create_empty_cg_dataset(st_g_id, gdf_translate(pdsets(i)), (/n_part/), Z_avail, O_OUT)
      enddo
#endif /* NBODY_1FILE */

   end subroutine create_empty_cg_datasets_in_output

   subroutine h5_write_to_single_file_v1(fname)

      use cg_level_finest, only: finest
      use cg_list,         only: cg_list_element
      use common_hdf5,     only: hdf_vars, hdf_vars_avail
      use constants,       only: ndims, LO, FP_REAL
      use dataio_pub,      only: die, h5_64bit
      use domain,          only: is_multicg !, is_uneven
      use grid_cont,       only: grid_container
      use hdf5,            only: HID_T, HSIZE_T, H5FD_MPIO_COLLECTIVE_F, H5P_DATASET_CREATE_F, H5P_DATASET_XFER_F, &
           &                     H5S_SELECT_SET_F, H5T_NATIVE_REAL, H5T_NATIVE_DOUBLE, H5F_ACC_RDWR_F, H5P_FILE_ACCESS_F, &
           &                     h5dwrite_f, h5screate_simple_f, h5pcreate_f, h5dcreate_f, h5sclose_f, h5dget_space_f, h5sselect_hyperslab_f, &
           &                     h5pset_dxpl_mpio_f, h5dclose_f, h5open_f, h5close_f, h5fopen_f, h5fclose_f, h5pclose_f, h5pset_fapl_mpio_f !, h5pset_chunk_f
      use MPIF,            only: MPI_INFO_NULL, MPI_COMM_WORLD

      implicit none

      character(len=*), intent(in)      :: fname
      integer(HID_T)                    :: file_id                 !< File identifier
      integer(HID_T)                    :: plist_id, plist_idf     !< Property list identifier
      integer                           :: i
      integer(kind=4)                   :: error                   !< error perhaps should be of type integer(HID_T)
      type(cg_list_element), pointer    :: cgl
      type(grid_container),  pointer    :: cg
      real, dimension(:,:,:), pointer   :: data                    !< Data to write
      integer(kind=4), parameter        :: rank = ndims            !< Dataset rank = 3
      integer(HID_T)                    :: dset_id                 !< Dataset identifier
      integer(HID_T)                    :: filespace               !< Dataspace identifier in file
      integer(HID_T)                    :: memspace                !< Dataspace identifier in memory
      integer(HSIZE_T), dimension(rank) :: count, offset, stride, block, dimsf, chunk_dims

      ! Sometimes the data(:,:,:) is created in an associated state, sometimes not
      nullify(data)

      call h5open_f(error)
      !
      ! Setup file access property list with parallel I/O access.
      !
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_idf, error)
#ifdef MPIF08
      call h5pset_fapl_mpio_f(plist_idf, MPI_COMM_WORLD%mpi_val, MPI_INFO_NULL%mpi_val, error)
#else /* !MPIF08 */
      call h5pset_fapl_mpio_f(plist_idf, MPI_COMM_WORLD, MPI_INFO_NULL, error)
#endif /* !MPIF08 */
      !
      ! Create the file collectively.
      !
      call h5fopen_f(trim(fname), H5F_ACC_RDWR_F, file_id, error, access_prp = plist_idf)
      call h5pclose_f(plist_idf, error)

      !! \todo check if finest is complete, if not then find finest complete level
      dimsf  = finest%level%l%n_d(:)    ! Dataset dimensions
      !
      ! Create the data space for the  dataset.
      !
      call h5screate_simple_f(rank, dimsf, filespace, error)
      do i = 1, size(hdf_vars)
         if (.not.hdf_vars_avail(i)) cycle

         ! Create chunked dataset.
         call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

         ! Cannot use in multiblock
         ! if (.not. is_uneven) call h5pset_chunk_f(plist_id, rank, chunk_dims, error) !> \todo check how much performance it gives (massively parallel I/O is required)

         call h5dcreate_f(file_id, gdf_translate(hdf_vars(i)), H5T_NATIVE_REAL, filespace, dset_id, error, plist_id)
         call h5sclose_f(filespace, error)

         call h5dget_space_f(dset_id, filespace, error)

         ! Create property list for collective dataset write
         call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
         if (.not. is_multicg) call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

         !! \todo if there are fine levels restrict the data first and write warning that v2 should be used instead
         cgl => finest%level%first
         if (.not. associated(cgl)) call die("[data_hdf5:h5_write_to_single_file_v1] I/O v1 cannot handle empty cg lists.")
         do while (associated(cgl))
            cg => cgl%cg

            if (.not.associated(data)) allocate(data(cg%nxb, cg%nyb, cg%nzb))
            call get_data_from_cg(hdf_vars(i), cg, data)
            if (.not.hdf_vars_avail(i)) cycle

            chunk_dims = cg%n_b(:) ! Chunks dimensions
            call h5screate_simple_f(rank, chunk_dims, memspace, error)

            ! Each process defines dataset in memory and writes it to the hyperslab in the file.
            stride(:) = 1
            count(:) =  1
            block(:) = chunk_dims(:)
            offset(:) = cg%my_se(:, LO)

            ! Select hyperslab in the file.
            call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, stride, block)

            ! Write the dataset collectively.
            if (h5_64bit) then
               call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dimsf, error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
            else
               call h5dwrite_f(dset_id, H5T_NATIVE_REAL, real(data, kind=FP_REAL), dimsf, error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
            endif

            ! Close dataspaces.
            call h5sclose_f(memspace, error)

            if (associated(data)) deallocate(data)

            cgl => cgl%nxt
         enddo

         !call h5pclose_f(plist_id, error) ! does it matter?

         ! Close the dataset.
         call h5dclose_f(dset_id, error)
      enddo

      call h5sclose_f(filespace, error)

      ! Close the property list.
      call h5fclose_f(file_id, error)
      call h5close_f(error)

   end subroutine h5_write_to_single_file_v1

   subroutine h5_write_to_multiple_files(sequential)

      use cg_leaves,   only: leaves
      use cg_list,     only: cg_list_element
      use common_hdf5, only: dump_announcement, dump_announce_time, hdf_vars, hdf_vars_avail
      use constants,   only: cwdlen, dsetnamelen, xdim, ydim, zdim, HDF, I_ONE
      use dataio_pub,  only: last_hdf_time, nhdf
      use grid_cont,   only: grid_container
      use h5lt,        only: h5ltmake_dataset_double_f
      use hdf5,        only: H5F_ACC_TRUNC_F, h5fcreate_f, h5open_f, h5fclose_f, h5close_f, HID_T, h5gcreate_f, h5gclose_f, HSIZE_T

      implicit none

      logical,               intent(in) :: sequential
      type(cg_list_element), pointer    :: cgl
      type(grid_container),  pointer    :: cg
      integer(kind=4), parameter        :: rank = 3
      integer(kind=4)                   :: error, i         !< error perhaps should be of type integer(HID_T)
      integer(HID_T)                    :: file_id, grp_id
      integer(kind=8)                   :: ngc              !< current grid index
      integer(HSIZE_T), dimension(rank) :: dims
      character(len=dsetnamelen)        :: gname
      character(len=cwdlen)             :: fname
      real, dimension(:,:,:), pointer   :: data             !< Data to write

      call dump_announcement(HDF, nhdf, fname, last_hdf_time, sequential)

      call h5open_f(error)
      call h5fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, error)
      cgl => leaves%first
      ngc = 0
      nullify(data)
      do while (associated(cgl))
         cg => cgl%cg

         write(gname,'("grid",i8.8)') ngc-1
         call h5gcreate_f(file_id, gname, grp_id, error)

         ! set attributes here
         call h5ltmake_dataset_double_f(grp_id, "fbnd", int(2,kind=4), shape(cg%fbnd,kind=HSIZE_T), cg%fbnd, error)

         if (.not.associated(data)) allocate(data(cg%n_b(xdim),cg%n_b(ydim),cg%n_b(zdim)))
         dims = cg%n_b(:)
         do i = I_ONE, size(hdf_vars, kind=4)
            call get_data_from_cg(hdf_vars(i), cg, data)
            if (.not.hdf_vars_avail(i)) cycle
            call h5ltmake_dataset_double_f(grp_id, hdf_vars(i), rank, dims, data(:,:,:), error)
         enddo
         if (associated(data)) deallocate(data)
         call h5gclose_f(grp_id, error)
         ngc = ngc + 1
         cgl => cgl%nxt
      enddo
      call h5fclose_f(file_id, error)
      call h5close_f(error)

      call dump_announce_time

   end subroutine h5_write_to_multiple_files

end module data_hdf5
