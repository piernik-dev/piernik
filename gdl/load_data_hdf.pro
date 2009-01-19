;===============================================================================
;
;  function load_data_hdf reads scientific data produced by "MHD" code in
;  HDF format.
;
;-------------------------------------------------------------------------------
;
;  arguments:
;
;    file - file to read data from
;    var  - variable to be readed
;
function load_data_hdf, filepref, step, var, xcoord = x, ycoord = y, zcoord = z, $
                             time = t, tstep =dt, nstep=nstep, dims = dims, $
                             nxa=nxa,nya=nya,nza=nza, nxd=nxd, nyd=nyd, nzd=nzd, nb=nb, $ 
                             xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,zmin=zmin,zmax=zmax, $
                             problem=problem_name, gamma=gamma

;  on_error, 2

  frame   = string(0, format = '(i2.2)')+'_' $
           +string(0, format = '(i2.2)')+'_' $
           +string(0, format = '(i2.2)')+'_' $
           +string(step, format = '(i4.4)')

  filename= filepref+'_'+frame+'.hdf'  

  LOAD_DIMS_HDF, filename, pdims=pdims, pcoords=pcoords, dims=dims, $
                           nxd,nyd,nzd, nxb,nyb,nzb, nb, $
                           xmin, xmax, ymin, ymax, zmin, zmax 
;  nproc = pdims(0)*pdims(1)*pdims(2)
  if(nxb eq dims(0)) then begin
    nx=nxb
    ny=nyb
    nz=nzb
    nb = 0
  endif else begin
    nx = dims(0)
    ny = dims(1)
    nz = dims(2)
  endelse

  framedisp   = '0x_0x_0x_'+string(step, format = '(i4.4)')
  filenamedisp= filepref+'_'+framedisp+'.hdf'  
;;;  print, 'file(s) = ',filenamedisp


  for pc0 = 0, pdims(0)-1 do begin
    for pc1 = 0, pdims(1)-1 do begin
      for pc2 = 0, pdims(2)-1 do begin

        frame   = string(pc0, format = '(i2.2)')+'_' $
                 +string(pc1, format = '(i2.2)')+'_' $
                 +string(pc2, format = '(i2.2)')+'_' $
                 +string(step, format = '(i4.4)')


        filename= filepref+'_'+frame+'.hdf'  

     case var of
       'eimt' : begin
                  ener = LOAD_BLOCK_HDF(filename, 'ener',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  dens = LOAD_BLOCK_HDF(filename, 'dens',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  velx = LOAD_BLOCK_HDF(filename, 'velx',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  vely = LOAD_BLOCK_HDF(filename, 'vely',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  velz = LOAD_BLOCK_HDF(filename, 'velz',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  magx = LOAD_BLOCK_HDF(filename, 'magx',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  magy = LOAD_BLOCK_HDF(filename, 'magy',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  magz = LOAD_BLOCK_HDF(filename, 'magz',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  a_block = ener - 0.5*(velx^2+vely^2+velz^2)*dens - 0.5*(magx^2+magy^2+magz^2)
                end
       'eint' : begin
                  ener = LOAD_BLOCK_HDF(filename, 'ener',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  dens = LOAD_BLOCK_HDF(filename, 'dens',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  velx = LOAD_BLOCK_HDF(filename, 'velx',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  vely = LOAD_BLOCK_HDF(filename, 'vely',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  velz = LOAD_BLOCK_HDF(filename, 'velz',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  a_block = (ener - 0.5*(velx^2+vely^2+velz^2)*dens) 
                end
       'temp' : begin
                  ener = LOAD_BLOCK_HDF(filename, 'ener',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  dens = LOAD_BLOCK_HDF(filename, 'dens',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  velx = LOAD_BLOCK_HDF(filename, 'velx',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  vely = LOAD_BLOCK_HDF(filename, 'vely',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  velz = LOAD_BLOCK_HDF(filename, 'velz',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  a_block = (ener - 0.5*(velx^2+vely^2+velz^2)*dens)/dens 
                end
       'ekin' : begin
                  dens = LOAD_BLOCK_HDF(filename, 'dens',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  velx = LOAD_BLOCK_HDF(filename, 'velx',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  vely = LOAD_BLOCK_HDF(filename, 'vely',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  velz = LOAD_BLOCK_HDF(filename, 'velz',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
  		  a_block = 0.5*dens*(velx^2 + vely^2 + velz^2)
                end
       'ekdi' : begin
                  ener = LOAD_BLOCK_HDF(filename, 'ener',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  dens = LOAD_BLOCK_HDF(filename, 'dens',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  velx = LOAD_BLOCK_HDF(filename, 'velx',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  vely = LOAD_BLOCK_HDF(filename, 'vely',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  velz = LOAD_BLOCK_HDF(filename, 'velz',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  a_block = 0.5*dens*(velx^2 + vely^2 + velz^2)/(ener - 0.5*(velx^2+vely^2+velz^2)*dens )
                end
       'eidk' : begin
                  ener = LOAD_BLOCK_HDF(filename, 'ener',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  dens = LOAD_BLOCK_HDF(filename, 'dens',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  velx = LOAD_BLOCK_HDF(filename, 'velx',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  vely = LOAD_BLOCK_HDF(filename, 'vely',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  velz = LOAD_BLOCK_HDF(filename, 'velz',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  a_block = (ener - 0.5*(velx^2+vely^2+velz^2)*dens )/(0.5*dens*(velx^2 + vely^2 + velz^2))
                end
       'emdi' : begin
                  ener = LOAD_BLOCK_HDF(filename, 'ener',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  dens = LOAD_BLOCK_HDF(filename, 'dens',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  velx = LOAD_BLOCK_HDF(filename, 'velx',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  vely = LOAD_BLOCK_HDF(filename, 'vely',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  velz = LOAD_BLOCK_HDF(filename, 'velz',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  magx = LOAD_BLOCK_HDF(filename, 'magx',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  magy = LOAD_BLOCK_HDF(filename, 'magy',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  magz = LOAD_BLOCK_HDF(filename, 'magz',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
	          a_block = 0.5*(magx^2 + magy^2 + magz^2)/(ener - 0.5*dens*(velx^2 + vely^2 + velz^2) - 0.5*(magx^2 + magy^2 + magz^2))
                end
       'eidm' : begin
                  ener = LOAD_BLOCK_HDF(filename, 'ener',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  dens = LOAD_BLOCK_HDF(filename, 'dens',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  velx = LOAD_BLOCK_HDF(filename, 'velx',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  vely = LOAD_BLOCK_HDF(filename, 'vely',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  velz = LOAD_BLOCK_HDF(filename, 'velz',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  magx = LOAD_BLOCK_HDF(filename, 'magx',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  magy = LOAD_BLOCK_HDF(filename, 'magy',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  magz = LOAD_BLOCK_HDF(filename, 'magz',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
  	          a_block = (ener - 0.5*dens*(velx^2 + vely^2 + velz^2) - 0.5*(magz^2 + magz^2 + magz^2))/(0.5*(magz^2 + magz^2 + magz^2))
                end
       'emdk' : begin
                  dens = LOAD_BLOCK_HDF(filename, 'dens',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  velx = LOAD_BLOCK_HDF(filename, 'velx',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  vely = LOAD_BLOCK_HDF(filename, 'vely',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  velz = LOAD_BLOCK_HDF(filename, 'velz',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  magx = LOAD_BLOCK_HDF(filename, 'magx',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  magy = LOAD_BLOCK_HDF(filename, 'magy',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  magz = LOAD_BLOCK_HDF(filename, 'magz',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
	          a_block = 0.5*(magx^2 + magy^2 + magz^2)/(0.5*dens*(velx^2 + vely^2 + velz^2))
                end
       'ekdm' : begin
                  dens = LOAD_BLOCK_HDF(filename, 'dens',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  velx = LOAD_BLOCK_HDF(filename, 'velx',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  vely = LOAD_BLOCK_HDF(filename, 'vely',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  velz = LOAD_BLOCK_HDF(filename, 'velz',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  magx = LOAD_BLOCK_HDF(filename, 'magx',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  magy = LOAD_BLOCK_HDF(filename, 'magy',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  magz = LOAD_BLOCK_HDF(filename, 'magz',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
             a_block = 0.5*dens*(velx^2 + vely^2 + velz^2)/(0.5*(magx^2 + magy^2 + magz^2))
                end
       'btor' : begin
                  magx = LOAD_BLOCK_HDF(filename, 'magx',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  magy = LOAD_BLOCK_HDF(filename, 'magy',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  a_block = magx
                  for iii=0,nx-1 do begin
                  for jjj=0,ny-1 do begin
                  a_block(iii,jjj,*) = (magy(iii,jjj,*)*x_block(iii)-magx(iii,jjj,*)*y_block(jjj))/sqrt(x_block(iii)^2+y_block(jjj)^2)
                  endfor
                  endfor
                end
       'bout' : begin
                  magx = LOAD_BLOCK_HDF(filename, 'magx',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  magy = LOAD_BLOCK_HDF(filename, 'magy',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  a_block = magx
                  for iii=0,nx-1 do begin
                  for jjj=0,ny-1 do begin
                  a_block(iii,jjj,*) = (magx(iii,jjj,*)*x_block(iii)+magy(iii,jjj,*)*y_block(jjj))/sqrt(x_block(iii)^2+y_block(jjj)^2)
                  endfor
                  endfor
                end
       'csnd' : begin
                  ener = LOAD_BLOCK_HDF(filename, 'ener',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  dens = LOAD_BLOCK_HDF(filename, 'dens',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  velx = LOAD_BLOCK_HDF(filename, 'velx',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  vely = LOAD_BLOCK_HDF(filename, 'vely',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  velz = LOAD_BLOCK_HDF(filename, 'velz',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
	          a_block = sqrt((ener-0.5*dens*(velx^2+vely^2+velz^2))/dens)
                end
       'vela' : begin
                  velx = LOAD_BLOCK_HDF(filename, 'velx',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  vely = LOAD_BLOCK_HDF(filename, 'vely',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  velz = LOAD_BLOCK_HDF(filename, 'velz',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
		  a_block = (velx^2 + vely^2 + velz^2)^0.5
                end
       'velr' : begin
                  velx = LOAD_BLOCK_HDF(filename, 'velx',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  vely = LOAD_BLOCK_HDF(filename, 'vely',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
		  a_block = (velx^2 + vely^2)^0.5
                end
       'vlxa' : begin
                  velx = LOAD_BLOCK_HDF(filename, 'velx',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  vely = LOAD_BLOCK_HDF(filename, 'vely',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  velz = LOAD_BLOCK_HDF(filename, 'velz',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
		  a_block = velx/(velx^2 + vely^2 + velz^2)^0.5
                end
       'vlya' : begin
                  velx = LOAD_BLOCK_HDF(filename, 'velx',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  vely = LOAD_BLOCK_HDF(filename, 'vely',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  velz = LOAD_BLOCK_HDF(filename, 'velz',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
		  a_block = vely/(velx^2 + vely^2 + velz^2)^0.5
                end
       'vlza' : begin
                  velx = LOAD_BLOCK_HDF(filename, 'velx',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  vely = LOAD_BLOCK_HDF(filename, 'vely',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  velz = LOAD_BLOCK_HDF(filename, 'velz',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
		  a_block = velz/(velx^2 + vely^2 + velz^2)^0.5
                end
       'maga' : begin
                  magx = LOAD_BLOCK_HDF(filename, 'magx',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  magy = LOAD_BLOCK_HDF(filename, 'magy',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  magz = LOAD_BLOCK_HDF(filename, 'magz',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
	          a_block = (magx^2 + magy^2 + magz^2)^0.5
                end
       'mgxa' : begin
                  magx = LOAD_BLOCK_HDF(filename, 'magx',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  magy = LOAD_BLOCK_HDF(filename, 'magy',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  magz = LOAD_BLOCK_HDF(filename, 'magz',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
	          a_block = magx/(magx^2 + magy^2 + magz^2)^0.5
                end
       'mgya' : begin
                  magx = LOAD_BLOCK_HDF(filename, 'magx',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  magy = LOAD_BLOCK_HDF(filename, 'magy',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  magz = LOAD_BLOCK_HDF(filename, 'magz',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
	          a_block = magy/(magx^2 + magy^2 + magz^2)^0.5
                end
       'mgza' : begin
                  magx = LOAD_BLOCK_HDF(filename, 'magx',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  magy = LOAD_BLOCK_HDF(filename, 'magy',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  magz = LOAD_BLOCK_HDF(filename, 'magz',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
	          a_block = magz/(magx^2 + magy^2 + magz^2)^0.5
                end
       'momx' : begin
                  dens = LOAD_BLOCK_HDF(filename, 'dens',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  velx = LOAD_BLOCK_HDF(filename, 'velx',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
		  a_block = dens*velx
                end
       'momy' : begin
                  dens = LOAD_BLOCK_HDF(filename, 'dens',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  vely = LOAD_BLOCK_HDF(filename, 'vely',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
		  a_block = dens*vely
                end
       'momz' : begin
                  dens = LOAD_BLOCK_HDF(filename, 'dens',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  velz = LOAD_BLOCK_HDF(filename, 'velz',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
		  a_block = dens*velz
                end
       'moma' : begin
                  dens = LOAD_BLOCK_HDF(filename, 'dens',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  velx = LOAD_BLOCK_HDF(filename, 'velx',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  vely = LOAD_BLOCK_HDF(filename, 'vely',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  velz = LOAD_BLOCK_HDF(filename, 'velz',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
		  a_block = dens*(velx^2 + vely^2 + velz^2)^0.5
                end
       'momr' : begin
                  dens = LOAD_BLOCK_HDF(filename, 'dens',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  velx = LOAD_BLOCK_HDF(filename, 'velx',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  vely = LOAD_BLOCK_HDF(filename, 'vely',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
		  a_block = dens*(velx^2 + vely^2)^0.5
                end
       'ccur' : begin
                  curx = LOAD_BLOCK_HDF(filename, 'curx',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  cury = LOAD_BLOCK_HDF(filename, 'cury',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  curz = LOAD_BLOCK_HDF(filename, 'curz',  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
                  a_block = (curx^2+cury^2+curz^2)^0.5
                end

     else :       a_block = LOAD_BLOCK_HDF(filename, var,  xcoord = x_block, ycoord = y_block, zcoord = z_block, time = t)
     endcase
    
        if(pc2 EQ 0) then begin
          a_column=a_block
          x_column=x_block
          y_column=y_block
          z_column=z_block
          nz_column=nz
        endif else begin  
          a_column = [[[a_column]],[[a_block]]] 
          z_column = [z_column,z_block]
          nz_column = nz_column + nz
        endelse
        a   = a_column
        z   = z_column
        nza = nz_column

      endfor

      if(pc1 EQ 0) then begin
        a_plane=a_column
        x_plane=x_column
        y_plane=y_column
        z_plane=z_column
        ny_plane=ny
      endif else begin  
        a_plane = [[a_plane],[a_column]] 
        y_plane = [y_plane,y_column]
        ny_plane = ny_plane + ny
      endelse
      a   = a_plane
      y   = y_plane
      nya = ny_plane

    endfor

    if(pc0 EQ 0) then begin
      a_volume=a_plane
      x_volume=x_plane
      y_volume=y_plane
      z_volume=z_plane
      nx_volume=nx
    endif else begin  
      a_volume = [a_volume,a_plane] 
      x_volume = [x_volume,x_column]
      nx_volume = nx_volume + nx
    endelse
    a   = a_volume
    x   = x_volume
    nxa = nx_volume

  endfor

  return, a
end
