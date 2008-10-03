#define RNG nb+1:nx-nb, nb+1:ny-nb, nb+1:nz-nb
module dataio_hdf5
   use h5lt
   use hdf5

   character(LEN=10), dimension(2) :: dname = (/"fluid","mag"/)

   contains
   subroutine write_restart_hdf5(filename,chdf)
     use types
     use mpi_setup
     use arrays, only : nx,ny,nz,nu, u,b,x,y,z, nxb,nyb,nzb
     use start, only : t,nxd,nyd,nzd,nb, domain, xmin,xmax, &
         ymin,ymax, zmin,zmax, nstep,dt
     use init_problem, only : problem_name, run_id
!     use dataio, only : nhdf, ntsl, nres, nlog, step_hdf, last_hdf_time
     IMPLICIT NONE
     type(hdf) :: chdf
     integer :: llun,fe
     CHARACTER :: filename*128,lfile*128  ! File name

     integer(HID_T) :: file_id       ! File identifier 
     integer(HID_T) :: dset_id       ! Dataset identifier 
     integer(HID_T) :: plist_id      ! Property list identifier 
     integer(HID_T) :: filespace     ! Dataspace identifier in file 
     integer(HID_T) :: memspace      ! Dataspace identifier in memory

     integer(HSIZE_T),  DIMENSION(4) :: count  
     integer(HSSIZE_T), DIMENSION(4) :: offset 
     integer(HSIZE_T),  DIMENSION(4) :: stride
     integer(HSIZE_T),  DIMENSION(4) :: block
     integer(HSIZE_T),  DIMENSION(1) :: acount  
     integer(HSSIZE_T), DIMENSION(1) :: aoffset 
     integer(HSIZE_T),  DIMENSION(1) :: astride
     integer(HSIZE_T),  DIMENSION(1) :: ablock
     integer(HSIZE_T),  DIMENSION(4) :: dimsf, dimsfi, chunk_dims
     integer(HSIZE_T),  DIMENSION(1) :: adimsf, adimsfi, achunk_dims

     integer, dimension(3) :: dims
     integer :: error, drank = 4, arank = 1

     real, dimension(1) :: rbuf
     integer(SIZE_T) :: bufsize = 1

     llun  = chdf%log_lun
     lfile = chdf%log_file

     dims(1) = nx
     dims(2) = ny
     dims(3) = nz

     CALL h5open_f(error)
     CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
     CALL h5pset_fapl_mpio_f(plist_id, comm3d, info, error)

     CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
     CALL h5pclose_f(plist_id, error)

     !----------------------------------------------------------------------------------
     !  WRITE FLUID VARIABLES
     !
     dimsf = (/nu,nx*pxsize,ny*pysize,nz*pzsize/) ! Dataset dimensions
     dimsfi = dimsf
     chunk_dims = (/nu,nx,ny,nz/)                   ! Chunks dimensions

     ! Create the data space for the  dataset. 
     CALL h5screate_simple_f(drank, dimsf, filespace, error)
     CALL h5screate_simple_f(drank, chunk_dims, memspace, error)
     
     ! Create chunked dataset.
     CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
     CALL h5pset_chunk_f(plist_id, drank, chunk_dims, error)
     CALL h5dcreate_f(file_id, dname(1), H5T_NATIVE_DOUBLE, filespace, &
                      dset_id, error, plist_id)
     CALL h5sclose_f(filespace, error)

     ! Each process defines dataset in memory and writes it to the hyperslab
     ! in the file. 
     stride(:) = 1 
     count(:)  = 1 
     block(:)  = chunk_dims(:)

     offset(1)   = 0
     offset(2:4) = pcoords(1:3)*chunk_dims(2:4)

     ! Select hyperslab in the file.
     CALL h5dget_space_f(dset_id, filespace, error)
     CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, &
                                 stride, block)
     CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
     CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
     CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, u(:,:,:,:), dimsfi, error, &
                     file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

     CALL h5sclose_f(filespace, error)
     CALL h5sclose_f(memspace, error)
     CALL h5dclose_f(dset_id, error)
     !----------------------------------------------------------------------------------

     !----------------------------------------------------------------------------------
     !  WRITE MAG VARIABLES
     !
     dimsf = (/3,nx*pxsize,ny*pysize,nz*pzsize/) ! Dataset dimensions
     dimsfi = dimsf
     chunk_dims = (/3,nx,ny,nz/)                 ! Chunks dimensions

     ! Create the data space for the  dataset. 
     CALL h5screate_simple_f(drank, dimsf, filespace, error)
     CALL h5screate_simple_f(drank, chunk_dims, memspace, error)
     
     ! Create chunked dataset.
     CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
     CALL h5pset_chunk_f(plist_id, drank, chunk_dims, error)
     CALL h5dcreate_f(file_id, dname(2), H5T_NATIVE_DOUBLE, filespace, &
                      dset_id, error, plist_id)
     CALL h5sclose_f(filespace, error)

     ! Each process defines dataset in memory and writes it to the hyperslab
     ! in the file. 
     stride(:) = 1 
     count(:)  = 1 
     block(:)  = chunk_dims(:)

     offset(1)   = 0
     offset(2:4) = pcoords(1:3)*chunk_dims(2:4)

     ! Select hyperslab in the file.
     CALL h5dget_space_f(dset_id, filespace, error)
     CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, &
                                 stride, block)
     CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
     CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
     CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, b(:,:,:,:), dimsfi, error, &
                     file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

     CALL h5sclose_f(filespace, error)
     CALL h5sclose_f(memspace, error)
     CALL h5dclose_f(dset_id, error)
     !----------------------------------------------------------------------------------

     !----------------------------------------------------------------------------------
     !  WRITE X Axis
     !
     adimsf  = (/nx*pxsize/) ! Dataset dimensions
     adimsfi = adimsf
     achunk_dims = (/nx/)    ! Chunks dimensions

     ! Create the data space for the  dataset. 
     CALL h5screate_simple_f(arank, adimsf, filespace, error)
     CALL h5screate_simple_f(arank, achunk_dims, memspace, error)
     
     ! Create chunked dataset.
     CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
     CALL h5pset_chunk_f(plist_id, arank, achunk_dims, error)
     CALL h5dcreate_f(file_id, "X axis", H5T_NATIVE_DOUBLE, filespace, &
                      dset_id, error, plist_id)
     CALL h5sclose_f(filespace, error)

     ! Each process defines dataset in memory and writes it to the hyperslab
     ! in the file. 
     astride(:) = 1 
     acount(:)  = 1 
     ablock(:)  = achunk_dims(:)

     aoffset(1) = pcoords(1)*achunk_dims(1)

     ! Select hyperslab in the file.
     CALL h5dget_space_f(dset_id, filespace, error)
     CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, aoffset, acount, error, &
                                 astride, ablock)
     CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
     CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
     CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x(:), adimsfi, error, &
                     file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

     CALL h5sclose_f(filespace, error)
     CALL h5sclose_f(memspace, error)
     CALL h5dclose_f(dset_id, error)
     !----------------------------------------------------------------------------------

     !----------------------------------------------------------------------------------
     !  WRITE Y Axis
     !
     adimsf  = (/ny*pysize/) ! Dataset dimensions
     adimsfi = adimsf
     achunk_dims = (/ny/)    ! Chunks dimensions

     ! Create the data space for the  dataset. 
     CALL h5screate_simple_f(arank, adimsf, filespace, error)
     CALL h5screate_simple_f(arank, achunk_dims, memspace, error)
     
     ! Create chunked dataset.
     CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
     CALL h5pset_chunk_f(plist_id, arank, achunk_dims, error)
     CALL h5dcreate_f(file_id, "Y axis", H5T_NATIVE_DOUBLE, filespace, &
                      dset_id, error, plist_id)
     CALL h5sclose_f(filespace, error)

     ! Each process defines dataset in memory and writes it to the hyperslab
     ! in the file. 
     astride(:) = 1 
     acount(:)  = 1 
     ablock(:)  = achunk_dims(:)

     aoffset(1) = pcoords(2)*achunk_dims(1)

     ! Select hyperslab in the file.
     CALL h5dget_space_f(dset_id, filespace, error)
     CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, aoffset, acount, error, &
                                 astride, ablock)
     CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
     CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
     CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, y(:), adimsfi, error, &
                     file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

     CALL h5sclose_f(filespace, error)
     CALL h5sclose_f(memspace, error)
     CALL h5dclose_f(dset_id, error)
     !----------------------------------------------------------------------------------

     !----------------------------------------------------------------------------------
     !  WRITE Z Axis
     !
     adimsf  = (/nz*pzsize/) ! Dataset dimensions
     adimsfi = adimsf
     achunk_dims = (/nz/)    ! Chunks dimensions

     ! Create the data space for the  dataset. 
     CALL h5screate_simple_f(arank, adimsf, filespace, error)
     CALL h5screate_simple_f(arank, achunk_dims, memspace, error)
     
     ! Create chunked dataset.
     CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
     CALL h5pset_chunk_f(plist_id, arank, achunk_dims, error)
     CALL h5dcreate_f(file_id, "Z axis", H5T_NATIVE_DOUBLE, filespace, &
                      dset_id, error, plist_id)
     CALL h5sclose_f(filespace, error)

     ! Each process defines dataset in memory and writes it to the hyperslab
     ! in the file. 
     astride(:) = 1 
     acount(:)  = 1 
     ablock(:)  = achunk_dims(:)

     aoffset(1) = pcoords(3)*achunk_dims(1)

     ! Select hyperslab in the file.
     CALL h5dget_space_f(dset_id, filespace, error)
     CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, aoffset, acount, error, &
                                 astride, ablock)
     CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
     CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
     CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, z(:), adimsfi, error, &
                     file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

     CALL h5sclose_f(filespace, error)
     CALL h5sclose_f(memspace, error)
     CALL h5dclose_f(dset_id, error)
     !----------------------------------------------------------------------------------
     CALL h5pclose_f(plist_id, error)
     CALL h5fclose_f(file_id, error)
     if(proc == 0) then
        CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error)

        bufsize = 1
        call h5ltset_attribute_double_f(file_id,"/","time", &
           (/t/),bufsize,error)
        call h5ltset_attribute_double_f(file_id,"/","timestep", &
           (/dt/),bufsize,error)
        call h5ltset_attribute_int_f(file_id,"/","nstep", &
           (/chdf%nstep/),bufsize,error)

        call h5ltset_attribute_int_f(file_id,"/","nres", &
           (/chdf%nres+1/),bufsize,error)
        call h5ltset_attribute_int_f(file_id,"/","nhdf", &
           (/chdf%nhdf/),bufsize,error)
        call h5ltset_attribute_int_f(file_id,"/","ntsl", &
           (/chdf%ntsl/),bufsize,error)
        call h5ltset_attribute_int_f(file_id,"/","nlog", &
           (/chdf%nlog/),bufsize,error)
        call h5ltset_attribute_int_f(file_id,"/","step_res", &
           (/chdf%nstep/),bufsize,error)
        call h5ltset_attribute_int_f(file_id,"/","step_hdf", &
           (/chdf%step_hdf/),bufsize,error)
        call h5ltset_attribute_double_f(file_id,"/","last_hdf_time", &
           (/chdf%last_hdf_time/),bufsize,error)

        call h5ltset_attribute_int_f(file_id,"/","nxd", &
           (/nxd/),bufsize,error)
        call h5ltset_attribute_int_f(file_id,"/","nyd", &
           (/nyd/),bufsize,error)
        call h5ltset_attribute_int_f(file_id,"/","nzd", &
           (/nzd/),bufsize,error)
        call h5ltset_attribute_int_f(file_id,"/","nxb", &
           (/nxb/),bufsize,error)
        call h5ltset_attribute_int_f(file_id,"/","nyb", &
           (/nyb/),bufsize,error)
        call h5ltset_attribute_int_f(file_id,"/","nzb", &
           (/nzb/),bufsize,error)
        call h5ltset_attribute_int_f(file_id,"/","nb", &
           (/nb/),bufsize,error)

        call h5ltset_attribute_double_f(file_id,"/","xmin", &
           (/xmin/),bufsize,error)
        call h5ltset_attribute_double_f(file_id,"/","xmax", &
           (/xmax/),bufsize,error)
        call h5ltset_attribute_double_f(file_id,"/","ymin", &
           (/ymin/),bufsize,error)
        call h5ltset_attribute_double_f(file_id,"/","ymax", &
           (/ymax/),bufsize,error)
        call h5ltset_attribute_double_f(file_id,"/","zmin", &
           (/zmin/),bufsize,error)
        call h5ltset_attribute_double_f(file_id,"/","zmax", &
           (/zmax/),bufsize,error)

        bufsize = 3
        call h5ltset_attribute_int_f(file_id,"/","psize", &
           psize,bufsize,error)

        fe = len(trim(problem_name))-1
        call h5ltset_attribute_string_f(file_id,"/","problem name", &
           problem_name(1:fe),error)
        fe = len(trim(domain))-1
        call h5ltset_attribute_string_f(file_id,"/","domain", &
           trim(domain(1:fe)),error)
        call h5ltset_attribute_string_f(file_id,"/","run id", &
           trim(run_id),error)

        CALL h5fclose_f(file_id, error)

        open(llun, file=lfile, position='append')
        write(llun,*) 'Writing restart file: ',trim(filename)
        write(*,*)    'Writing restart file: ',trim(filename)
        close(llun)
     endif
     CALL h5close_f(error)

   end subroutine write_restart_hdf5

   subroutine read_restart_hdf5(chdf)
     use types
     use mpi_setup
     use arrays, only : nx,ny,nz,nu, u,b,x,y,z, nxb,nyb,nzb
     use start, only : t,nxd,nyd,nzd,nb, domain, xmin,xmax, &
         ymin,ymax, zmin,zmax, nstep,dt, new_id
     use init_problem, only : problem_name, run_id
!    use dataio, only : nhdf, ntsl, nres, nlog, step_hdf, &
!        last_hdf_time,step_res, log_lun, log_file
     IMPLICIT NONE
     type(hdf) :: chdf
     integer :: log_lun
     CHARACTER(LEN=128) :: log_file  ! File name
     CHARACTER(LEN=128) :: filename  ! File name

     integer(HID_T) :: file_id       ! File identifier 
     integer(HID_T) :: dset_id       ! Dataset identifier 
     integer(HID_T) :: plist_id      ! Property list identifier 
     integer(HID_T) :: filespace     ! Dataspace identifier in file 
     integer(HID_T) :: memspace      ! Dataspace identifier in memory

     integer(HSIZE_T),  DIMENSION(4) :: count  
     integer(HSSIZE_T), DIMENSION(4) :: offset 
     integer(HSIZE_T),  DIMENSION(4) :: stride
     integer(HSIZE_T),  DIMENSION(4) :: block
     integer(HSIZE_T),  DIMENSION(1) :: acount  
     integer(HSSIZE_T), DIMENSION(1) :: aoffset 
     integer(HSIZE_T),  DIMENSION(1) :: astride
     integer(HSIZE_T),  DIMENSION(1) :: ablock
     integer(HSIZE_T),  DIMENSION(4) :: dimsf, dimsfi, chunk_dims
     integer(HSIZE_T),  DIMENSION(1) :: adimsf, adimsfi, achunk_dims

     integer, dimension(3) :: dims
     integer :: error, drank = 4, arank = 1
     logical file_exist, log_exist

     real, dimension(1) :: rbuf
     integer, dimension(1) :: ibuf
     integer(SIZE_T) :: bufsize = 1

     log_lun  = chdf%log_lun
     log_file = chdf%log_file

     if(proc==0) then
        write (filename,'(a,a1,a3,a1,i4.4,a4)') &
          trim(problem_name),'_', run_id,'_',chdf%nres,'.res'
        write(*,*) 'Reading restart  file: ', trim(filename)
        inquire(file=log_file , exist = log_exist)
          if(log_exist .eqv. .true.) then
              open(log_lun, file=log_file, position='append')
              write(log_lun,*) 'Reading restart  file: ',trim(filename)
              close(log_lun)
          endif
     endif
     call MPI_BCAST(filename, 128,MPI_CHARACTER, 0, comm, ierr)


     inquire(file = filename, exist = file_exist)
     if(file_exist .eqv. .false.) then
        if(log_exist) then
          open(log_lun, file=log_file, position='append')
          write(log_lun,*) 'Restart  file: ', trim(filename), &
               ' does not exist.  ABORTING !!! '
          close(log_lun)
        endif

        write(*,*)       'Restart  file: ', trim(filename), &
            ' does not exist.  ABORTING !!! '
        call MPI_BARRIER(comm3d,ierr)
        call mpistop
        stop
     endif

     CALL h5open_f(error)
     CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
     CALL h5pset_fapl_mpio_f(plist_id, comm3d, info, error)

     CALL h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
     CALL h5pclose_f(plist_id, error)

     !----------------------------------------------------------------------------------
     !  READ FLUID VARIABLES
     !
     dimsf = (/nu,nx*pxsize,ny*pysize,nz*pzsize/) ! Dataset dimensions
     dimsfi = dimsf

     ! Create chunked dataset.
     CALL h5dopen_f(file_id, dname(1), dset_id, error)

     call h5dget_space_f(dset_id, filespace, error)
     call H5sget_simple_extent_ndims_f (filespace,drank,error)
     call H5dget_create_plist_f (dset_id,plist_id,error)
     call h5pget_chunk_f(plist_id, drank, chunk_dims, error)

     ! Each process defines dataset in memory and writes it to the hyperslab
     ! in the file. 
     stride(:) = 1 
     count(:)  = 1 
     block(:)  = chunk_dims(:)

     offset(1)   = 0
     offset(2:4) = pcoords(1:3)*chunk_dims(2:4)

     ! Select hyperslab in the file.
     CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, &
                                 stride, block)
     CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
     CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
     CALL h5screate_simple_f(drank, chunk_dims, memspace, error)
     CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, u, dimsfi, error, &
                     file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

     CALL h5sclose_f(filespace, error)
     CALL h5sclose_f(memspace, error)
     CALL h5pclose_f(plist_id, error)
     CALL h5dclose_f(dset_id, error)
     !----------------------------------------------------------------------------------

     !----------------------------------------------------------------------------------
     !  READ MAG VARIABLES
     !
     dimsf = (/3,nx*pxsize,ny*pysize,nz*pzsize/) ! Dataset dimensions
     dimsfi = dimsf

     ! Create chunked dataset.
     CALL h5dopen_f(file_id, dname(2), dset_id, error)

     call h5dget_space_f(dset_id, filespace, error)
     call H5Sget_simple_extent_ndims_f (filespace,drank,error)
     call H5Dget_create_plist_f (dset_id,plist_id,error)
     call h5pget_chunk_f(plist_id, drank, chunk_dims, error)


     ! Each process defines dataset in memory and writes it to the hyperslab
     ! in the file. 
     stride(:) = 1 
     count(:)  = 1 
     block(:)  = chunk_dims(:)

     offset(1)   = 0
     offset(2:4) = pcoords(1:3)*chunk_dims(2:4)

     ! Select hyperslab in the file.
     CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, &
                                 stride, block)
     CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
     CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
     CALL h5screate_simple_f(drank, chunk_dims, memspace, error)
     CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, b, dimsfi, error, &
                     file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

     CALL h5sclose_f(filespace, error)
     CALL h5sclose_f(memspace, error)
     CALL h5pclose_f(plist_id, error)
     CALL h5dclose_f(dset_id, error)
     !----------------------------------------------------------------------------------
!     CALL h5pclose_f(plist_id, error)
     CALL h5fclose_f(file_id, error)
     if(proc == 0) then
        CALL h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error)
        bufsize = 1
        call h5ltget_attribute_double_f(file_id,"/","time", &
           rbuf,error)
        t = rbuf(1)
        call h5ltget_attribute_double_f(file_id,"/","timestep", &
           rbuf,error)
        dt = rbuf(1)
        call h5ltget_attribute_int_f(file_id,"/","nstep", &
           ibuf,error)
        chdf%nstep = ibuf(1)

        call h5ltget_attribute_int_f(file_id,"/","nres", &
           ibuf,error)
        chdf%nres = ibuf(1)
        call h5ltget_attribute_int_f(file_id,"/","nhdf", &
           ibuf,error)
        chdf%nhdf = ibuf(1)
        call h5ltget_attribute_int_f(file_id,"/","ntsl", &
           ibuf,error)
        chdf%ntsl = ibuf(1)
        call h5ltget_attribute_int_f(file_id,"/","nlog", &
           ibuf,error)
        chdf%nlog = ibuf(1)
        call h5ltget_attribute_int_f(file_id,"/","step_res", &
           ibuf,error)
        chdf%step_res = ibuf(1)
        call h5ltget_attribute_int_f(file_id,"/","step_hdf", &
           ibuf,error)
        chdf%step_hdf = ibuf(1)
        call h5ltget_attribute_double_f(file_id,"/","last_hdf_time", &
           rbuf,error)
        chdf%last_hdf_time = rbuf(1)

        call h5ltget_attribute_string_f(file_id,"/","problem name", &
           problem_name,error)
        call h5ltget_attribute_string_f(file_id,"/","domain", &
           domain,error)
        call h5ltget_attribute_string_f(file_id,"/","run id", &
           new_id,error)

        CALL h5fclose_f(file_id, error)

        open(log_lun, file=log_file, position='append')
        write(log_lun,*) 'Done reading restart file: ',trim(filename)
        write(*,*)    'Done reading restart file: ',trim(filename)
        close(log_lun)
     endif

     call MPI_BCAST(chdf%nstep, 1, MPI_INTEGER, 0, comm3d, ierr)
     call MPI_BCAST(chdf%nres, 1, MPI_INTEGER, 0, comm3d, ierr)
     call MPI_BCAST(chdf%nhdf, 1, MPI_INTEGER, 0, comm3d, ierr)
     call MPI_BCAST(chdf%ntsl, 1, MPI_INTEGER, 0, comm3d, ierr)
     call MPI_BCAST(chdf%nlog, 1, MPI_INTEGER, 0, comm3d, ierr)
     call MPI_BCAST(chdf%step_res, 1, MPI_INTEGER, 0, comm3d, ierr)
     call MPI_BCAST(chdf%step_hdf, 1, MPI_INTEGER, 0, comm3d, ierr)

     call MPI_BCAST(chdf%last_hdf_time, 1, MPI_DOUBLE_PRECISION, 0, comm3d, ierr)
     call MPI_BCAST(t, 1, MPI_DOUBLE_PRECISION, 0, comm3d, ierr)
     call MPI_BCAST(dt, 1, MPI_DOUBLE_PRECISION, 0, comm3d, ierr)

     CALL MPI_BCAST(problem_name, 32, MPI_CHARACTER, 0, comm3d, ierr)
     CALL MPI_BCAST(domain, 16, MPI_CHARACTER, 0, comm3d, ierr)
     CALL MPI_BCAST(new_id, 3, MPI_CHARACTER, 0, comm3d, ierr)
     CALL h5close_f(error)

   end subroutine read_restart_hdf5
!
! ------------------------------------------------------------------------------------
!
   subroutine write_hdf5(chdf)

     use types
     use mpi_setup, only: pcoords, comm3d, proc, info, psize,ierr
     use arrays, only : idna, imxa, imya, imza, iena, nxb,nyb,nzb, u, nx,ny,nz
     use start, only : t,nxd,nyd,nzd,nb, domain, xmin,xmax, &
         ymin,ymax, zmin,zmax, nstep,dt
     use init_problem, only : problem_name, run_id
     
     IMPLICIT NONE
     type(hdf) :: chdf
     integer(HID_T) :: file_id       ! File identifier 
     integer(HID_T) :: plist_id      ! Property list identifier 
     integer :: llun,fe
     CHARACTER :: lfile*128, dd*4 ! File name
     CHARACTER(LEN=32) :: fname

     REAL*4, ALLOCATABLE :: data (:,:,:)  ! Data to write
     integer :: error, error_n, i   
     integer(SIZE_T) :: bufsize = 1

     ! Initialize HDF5 library and Fortran interfaces.
     !
     llun = chdf%log_lun
     lfile = chdf%log_file
     write(dd,'(i4.4)') chdf%nhdf
     fname = trim(problem_name)//"_"//trim(run_id)//"_"//dd//".h5"

     CALL h5open_f(error) 
     ! 
     ! Setup file access property list with parallel I/O access.
     !
     CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
     CALL h5pset_fapl_mpio_f(plist_id, comm3d, info, error)

     !
     ! Create the file collectively.
     !
     CALL h5fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
     CALL h5pclose_f(plist_id, error)
     ALLOCATE (data(nxb,nyb,nzb))

     data = real(u(idna,RNG))
     call write_arr(data,"dens",file_id)
     data = real(u(imxa,RNG) / u(idna,RNG))
     call write_arr(data,"velx",file_id)
     data = real(u(imya,RNG) / u(idna,RNG))
     call write_arr(data,"vely",file_id)
     data = real(u(imza,RNG) / u(idna,RNG))
     call write_arr(data,"velz",file_id)
     data = real(u(iena,RNG))
     call write_arr(data,"ener",file_id)

     DEALLOCATE(data)
     !
     ! Close the property list.
     !
     CALL h5fclose_f(file_id, error)
     if(proc == 0) then
        CALL h5fopen_f (fname, H5F_ACC_RDWR_F, file_id, error)

        bufsize = 1
        call h5ltset_attribute_double_f(file_id,"/","time", &
           (/t/),bufsize,error)
        call h5ltset_attribute_double_f(file_id,"/","timestep", &
           (/dt/),bufsize,error)
        call h5ltset_attribute_int_f(file_id,"/","nstep", &
           (/nstep/),bufsize,error)

        call h5ltset_attribute_int_f(file_id,"/","nxd", &
           (/nxd/),bufsize,error)
        call h5ltset_attribute_int_f(file_id,"/","nyd", &
           (/nyd/),bufsize,error)
        call h5ltset_attribute_int_f(file_id,"/","nzd", &
           (/nzd/),bufsize,error)
        call h5ltset_attribute_int_f(file_id,"/","nxb", &
           (/nxb/),bufsize,error)
        call h5ltset_attribute_int_f(file_id,"/","nyb", &
           (/nyb/),bufsize,error)
        call h5ltset_attribute_int_f(file_id,"/","nzb", &
           (/nzb/),bufsize,error)
        call h5ltset_attribute_int_f(file_id,"/","nb", &
           (/nb/),bufsize,error)

        call h5ltset_attribute_double_f(file_id,"/","xmin", &
           (/xmin/),bufsize,error)
        call h5ltset_attribute_double_f(file_id,"/","xmax", &
           (/xmax/),bufsize,error)
        call h5ltset_attribute_double_f(file_id,"/","ymin", &
           (/ymin/),bufsize,error)
        call h5ltset_attribute_double_f(file_id,"/","ymax", &
           (/ymax/),bufsize,error)
        call h5ltset_attribute_double_f(file_id,"/","zmin", &
           (/zmin/),bufsize,error)
        call h5ltset_attribute_double_f(file_id,"/","zmax", &
           (/zmax/),bufsize,error)

        bufsize = 3
        call h5ltset_attribute_int_f(file_id,"/","psize", &
           psize,bufsize,error)

        call h5ltset_attribute_string_f(file_id,"/","problem name", &
           trim(problem_name),error)
        call h5ltset_attribute_string_f(file_id,"/","domain", &
           domain(1:11),error)
        call h5ltset_attribute_string_f(file_id,"/","run id", &
           trim(run_id),error)


        CALL h5fclose_f(file_id, error)
        open(llun, file=lfile, position='append')
        write(llun,*) 'Writing output   file: ',trim(fname)
        write(*,*)       'Writing output   file: ',trim(fname)
        close(llun)
     endif
     call MPI_BARRIER(comm3d,ierr)
     CALL h5close_f(error)

     end subroutine write_hdf5

     subroutine write_arr(data,dsetname,file_id)
          use start, only : nxd,nyd,nzd,nb
          use arrays, only : nxb,nyb,nzb
          use mpi_setup, only: pcoords

          implicit none
          real*4, dimension(:,:,:) :: data

     integer :: rank = 3 ! Dataset rank 

     CHARACTER(LEN=4) :: dsetname    ! Dataset name

     integer(HID_T) :: file_id       ! Dataset identifier 
     integer(HID_T) :: dset_id       ! Dataset identifier 
     integer(HID_T) :: plist_id       ! Dataset identifier 
     integer(HID_T) :: filespace     ! Dataspace identifier in file 
     integer(HID_T) :: memspace      ! Dataspace identifier in memory

     integer, parameter :: ndims = 3
     integer(HSIZE_T),  DIMENSION(ndims) :: count  
     integer(HSSIZE_T), DIMENSION(ndims) :: offset 
     integer(HSIZE_T),  DIMENSION(ndims) :: stride
     integer(HSIZE_T),  DIMENSION(ndims) :: block
     integer(HSIZE_T), DIMENSION(ndims) :: dimsf, dimsfi, chunk_dims
     integer :: error

     dimsf = (/nxd,nyd,nzd/) ! Dataset dimensions
     dimsfi = dimsf
     chunk_dims = (/nxb,nyb,nzb/) ! Chunks dimensions
     !
     ! Create the data space for the  dataset. 
     !
     CALL h5screate_simple_f(rank, dimsf, filespace, error)
     CALL h5screate_simple_f(rank, chunk_dims, memspace, error)

     !
     ! Create chunked dataset.
     !
     CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
     CALL h5pset_chunk_f(plist_id, rank, chunk_dims, error)
     CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_real, filespace, &
                      dset_id, error, plist_id)
     CALL h5sclose_f(filespace, error)

     !
     ! Each process defines dataset in memory and writes it to the hyperslab
     ! in the file. 
     !
     stride(:) = 1 
     count(:) =  1 
     block(:) = chunk_dims(:)

     offset(:) = pcoords(:)*chunk_dims(:)
     ! 
     ! Select hyperslab in the file.
     !
     CALL h5dget_space_f(dset_id, filespace, error)
     CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, &
                                 stride, block)
     ! 
     ! Initialize data buffer with trivial data.
     !
     !
     ! Create property list for collective dataset write
     !
     CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
     CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
     
     !
     ! Write the dataset collectively. 
     !
     CALL h5dwrite_f(dset_id, H5T_NATIVE_real, data, dimsfi, error, &
                     file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

     !
     ! Close dataspaces.
     !
     CALL h5sclose_f(filespace, error)
     CALL h5sclose_f(memspace, error)
     !
     ! Close the dataset.
     !
     CALL h5dclose_f(dset_id, error)

       end subroutine write_arr

END module dataio_hdf5
