PRO SHOW

data_dir = '../runs/sedov'
prefix = 'sedov_tst'

png_dir = data_dir+'/frames'

DEVICE,DECOMPOSED=0,RETAIN=2

;===============================================================================
;!PATH = '/home/MYHOME/piernik/idl:' + !PATH 
;===============================================================================
data_files = FINDFILE(data_dir+'/'+prefix+'_*.h5', Count=n_files)
IF(n_files EQ 0) THEN GOTO, SKIP
;===============================================================================

;print, data_files

first_frame  	= 0 ; 1

last_frame   	= n_files-1
freq_frame   	= 1

png_output	= 'n'

display_frames	= 'y'		; 'y' or 'n'

win_vert_size 	= 600		; Vertical size (in pixels) 'xz' and 'yz' slices
            				; The other sizes are scaled according to the
				            ; physical sizes 5of the computational box
n_vect_x     	= 16 		; The numbers of points for vector fields
n_vect_y     	= 16 		; in each direction
n_vect_z     	= 16		;

coord_sys 	= 'xyz'         ; 'xyz', 'zrp' or 'rtp'

log_scal        = 'n'

rebin_factor    = 1

load_data       = 'y'

;===============================================================================
SHOW_INIT, s, slice_array, n_slices, plot_par, disp_par, units, batchmode
;===============================================================================

s.sw	= 'on'  
s.panel_name	= 'a'			; If more slices of the same type 
                                        ; are needed use this index

s.type		= 'xy'			; Chose 'yz', 'xz' or 'xy' plane
s.coord		=  0.0         ; Position at the complementary coordinate
s.vect_disp	= 'veli'			; Vector field to display: 'b' or 'v'
s.vect_scaling	= 'free'		; 'fix' or 'free' 
s.vect_scale	=  1.0


s.scal_disp	= 'deni'	        ; Scalar field to display 'd' or 'e'
s.scal_pert     = ''			; inactive
s.scal_scaling	= 'free'		; 'fix' or 'free'
s.scal_scale	= [-4.0,-1.5]

s.name=s.vect_disp+'_'+s.scal_disp+'_'+s.type+'_'+s.panel_name

slice_array = [slice_array,s]
n_slices = n_slices+1

;------------------------------------------------------------------------------


s.sw	= 'on'  
s.panel_name	= 'b'			; If more slices of the same type 
                                        ; are needed use this index

s.type		= 'xy'			; Chose 'yz', 'xz' or 'xy' plane
s.coord		=  0.0         		; Position at the complementary coordinate
s.vect_disp	= 'veli'		; Vector field to display: 'mag', 'veli', 'veln' or 'veld'
s.vect_scaling	= 'free'		; 'fix' or 'free' 
s.vect_scale	=  1.0


s.scal_disp	= 'enei'	        ; Scalar field to display: any valid variable, eg. 'deni', 'vlxi',etc... 
s.scal_pert     = ''			; inactive
s.scal_scaling	= 'free'		; 'fix' or 'free'
s.scal_scale	= [-4.0,-1.5]

s.name=s.vect_disp+'_'+s.scal_disp+'_'+s.type+'_'+s.panel_name

slice_array = [slice_array,s]
n_slices = n_slices+1

;------------------------------------------------------------------------------




;==============================================================================
n_vectors = [n_vect_x,n_vect_y,n_vect_z]
;==============================================================================


  DIRS, data_dir, png_dir, png_output

  READ_DATA, data_dir, prefix, first_frame, data, /noarrays
  COORDS, data.attr.crdsys, csys
  

  PLOT_SIZES, data, win_vert_size, data.attr.crdsys, psizes, plot_par
 
  win_open = 'y'
  first_call = 'y'
  FOR i_frame = first_frame, last_frame, freq_frame DO BEGIN

 
     IF(load_data NE 'n') THEN BEGIN
        READ_DATA, data_dir, prefix, i_frame, data
     ENDIF

     PLOT_SLICES,  data, slice_array, psizes, csys, plot_par, disp_par, units, n_vectors, $
                  display_frames, i_frame, png_output, png_dir, win_open, first_call

  END

SKIP:

END
