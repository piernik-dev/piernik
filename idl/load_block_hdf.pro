;===============================================================================
;
;  function load_block_hdf reads scientific data produced by "MHD" code in
;  HDF format.
;
;-------------------------------------------------------------------------------
;
;  arguments:
;
;    file - file to read data from
;    var  - variable to be readed
;
function load_block_hdf, file, var, xcoord = xc, ycoord = yc, zcoord = zc, $
                             time = t, tstep =dt, nstep=nstep, dims = dims, $
                             nxf=nx,nyf=ny,nzf=nz, nxb=nxb, nyb=nyb, nzb=nzb, nb=nb, $
                             xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,zmin=zmin,zmax=zmax, $
                             problem=problem_name, gamma=gamma

    on_error, 2

    sd_id = hdf_sd_start(file, /read)

    pdims    = uintarr(3)
    pdims[*]    = 1

    attr_index = hdf_sd_attrfind(sd_id, 'psize1')
    hdf_sd_attrinfo, sd_id, attr_index, data = temp
    pdims[0] = temp
    attr_index = hdf_sd_attrfind(sd_id, 'psize2')
    hdf_sd_attrinfo, sd_id, attr_index, data = temp
    pdims[1] = temp
    attr_index = hdf_sd_attrfind(sd_id, 'psize3')
    hdf_sd_attrinfo, sd_id, attr_index, data = temp
    pdims[2] = temp

    pcoords    = uintarr(3)
    pcoords[*]    = 0

    attr_index = hdf_sd_attrfind(sd_id, 'pcoords1')
    hdf_sd_attrinfo, sd_id, attr_index, data = temp
    pcoords[0] = temp
    attr_index = hdf_sd_attrfind(sd_id, 'pcoords2')
    hdf_sd_attrinfo, sd_id, attr_index, data = temp
    pcoords[1] = temp
    attr_index = hdf_sd_attrfind(sd_id, 'pcoords3')
    hdf_sd_attrinfo, sd_id, attr_index, data = temp
    pcoords[2] = temp

    dims    = uintarr(3)
    dims[*]    = 1

    attr_index = hdf_sd_attrfind(sd_id, 'dims1')
    hdf_sd_attrinfo, sd_id, attr_index, data = temp
    dims[0] = temp
    attr_index = hdf_sd_attrfind(sd_id, 'dims2')
    hdf_sd_attrinfo, sd_id, attr_index, data = temp
    dims[1] = temp
    attr_index = hdf_sd_attrfind(sd_id, 'dims3')
    hdf_sd_attrinfo, sd_id, attr_index, data = temp
    dims[2] = temp

    attr_index = hdf_sd_attrfind(sd_id, 'problem')
    hdf_sd_attrinfo, sd_id, attr_index, data = problem_name
    attr_index = hdf_sd_attrfind(sd_id, 'nstep')
    hdf_sd_attrinfo, sd_id, attr_index, data = nstep
    attr_index = hdf_sd_attrfind(sd_id, 'time')
    hdf_sd_attrinfo, sd_id, attr_index, data = t
    attr_index = hdf_sd_attrfind(sd_id, 'timestep')
    hdf_sd_attrinfo, sd_id, attr_index, data = dt
;mh    attr_index = hdf_sd_attrfind(sd_id, 'gamma')
;mh    hdf_sd_attrinfo, sd_id, attr_index, data = gamma

    sizes = dims

    data = fltarr(sizes[0], sizes[1], sizes[2])
    xc    = fltarr(sizes[0])
    yc    = fltarr(sizes[1])
    zc    = fltarr(sizes[2])

    attr_index = hdf_sd_attrfind(sd_id, 'time')
    hdf_sd_attrinfo, sd_id, attr_index, data = t

    sds_index = hdf_sd_nametoindex(sd_id, var)
    if (sds_index ge 0) then begin
      sds_id = hdf_sd_select(sd_id, sds_index)
      hdf_sd_getdata, sds_id, data
      hdf_sd_endaccess, sds_id
    endif else begin
      print, 'Variable ' + var + ' is not present in file!'
      stop
    endelse

    sds_index = hdf_sd_nametoindex(sd_id, 'xc')
    sds_id = hdf_sd_select(sd_id, sds_index)
    hdf_sd_getdata, sds_id, xc
    hdf_sd_endaccess, sds_id

    sds_index = hdf_sd_nametoindex(sd_id, 'yc')
    sds_id = hdf_sd_select(sd_id, sds_index)
    hdf_sd_getdata, sds_id, yc
    hdf_sd_endaccess, sds_id

    sds_index = hdf_sd_nametoindex(sd_id, 'zc')
    sds_id = hdf_sd_select(sd_id, sds_index)
    hdf_sd_getdata, sds_id, zc
    hdf_sd_endaccess, sds_id


    hdf_sd_end, sd_id


  return, data
end
