;===============================================================================
;
;  function LOADZEUSHDF reads scientific data produced by ZEUS3D in
;  HDF format.
;
;-------------------------------------------------------------------------------
;
;  arguments:
;
;    file - file to read data from
;    var  - variable to be readed
;
pro load_dims_hdf,  file, pdims=pdims, pcoords=pcoords, dims = dims, $
                    nxd,nyd,nzd, nxb,nyb,nzb, nb, $
                    xmin, xmax, ymin, ymax, zmin, zmax 

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

    nxd=0l
    nyd=0l
    nzd=0l
    nxb=0l
    nyb=0l
    nzb=0l
    nb=0l

    attr_index = hdf_sd_attrfind(sd_id, 'nxd')
    hdf_sd_attrinfo, sd_id, attr_index, data = temp
    nxd =temp(0)
    attr_index = hdf_sd_attrfind(sd_id, 'nyd')
    hdf_sd_attrinfo, sd_id, attr_index, data = temp
    nyd =temp(0)
    attr_index = hdf_sd_attrfind(sd_id, 'nzd')
    hdf_sd_attrinfo, sd_id, attr_index, data = temp
    nzd =temp(0)
    attr_index = hdf_sd_attrfind(sd_id, 'nxb')
    hdf_sd_attrinfo, sd_id, attr_index, data = temp
    nxb =temp(0)
    attr_index = hdf_sd_attrfind(sd_id, 'nyb')
    hdf_sd_attrinfo, sd_id, attr_index, data = temp
    nyb =temp(0)
    attr_index = hdf_sd_attrfind(sd_id, 'nzb')
    hdf_sd_attrinfo, sd_id, attr_index, data = temp
    nzb =temp(0)
    attr_index = hdf_sd_attrfind(sd_id, 'nb')
    hdf_sd_attrinfo, sd_id, attr_index, data = temp
    nb =temp(0)

    attr_index = hdf_sd_attrfind(sd_id, 'xmin')
    hdf_sd_attrinfo, sd_id, attr_index, data = xmin
    attr_index = hdf_sd_attrfind(sd_id, 'xmax')
    hdf_sd_attrinfo, sd_id, attr_index, data = xmax
    attr_index = hdf_sd_attrfind(sd_id, 'ymin')
    hdf_sd_attrinfo, sd_id, attr_index, data = ymin
    attr_index = hdf_sd_attrfind(sd_id, 'ymax')
    hdf_sd_attrinfo, sd_id, attr_index, data = ymax
    attr_index = hdf_sd_attrfind(sd_id, 'zmin')
    hdf_sd_attrinfo, sd_id, attr_index, data = zmin
    attr_index = hdf_sd_attrfind(sd_id, 'zmax')
    hdf_sd_attrinfo, sd_id, attr_index, data = zmax

    hdf_sd_end, sd_id

end
