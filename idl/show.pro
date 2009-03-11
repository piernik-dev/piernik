PRO SHOW

data_dir = '../runs/sedov/'
prefix = 'sedov_tst'


png_dir = data_dir+'/frames'

DEVICE,DECOMPOSED=0,RETAIN=2

;===============================================================================
!PATH = '/home/MYHOME/piernik/idl:' + !PATH 
;===============================================================================
data_files = FINDFILE(data_dir+'/'+prefix+'_00_00_00_*.hdf', Count=n_files)
IF(n_files EQ 0) THEN GOTO, SKIP
;===============================================================================


first_frame  	= 0
last_frame   	= n_files-1
freq_frame   	= 1

png_output	= 'y'

display_frames	= 'y'		; 'y' or 'n'

win_vert_size 	= 600		; Vertical size (in pixels) 'xz' and 'yz' slices
            				; The other sizes are scaled according to the
				            ; physical sizes of the computational box
n_vect_x     	= 16 		; The numbers of points for vector fields
n_vect_y     	= 16 		; in each direction
n_vect_z     	= 16		;

coord_sys 	= 'xyz'         ; 'xyz', 'zrp' or 'rtp'

log_scal        = 'n'

;===============================================================================
SHOW_INIT, s, slice_array, n_slices
;===============================================================================

s.sw	= 'on'  
s.panel_name	= 'a'			; If more slices of the same type 
                                        ; are needed use this index

s.type		= 'xy'			; Chose 'yz', 'xz' or 'xy' plane
s.coord		=  0.0         ; Position at the complementary coordinate
s.vect_disp	= 'v'			; Vector field to display: 'b' or 'v'
s.vect_scaling	= 'free'		; 'fix' or 'free' 
s.vect_scale	=  1.0


s.scal_disp	= 'd'			; Scalar field to display 'd' or 'e'
s.scal_pert     = ''			; inactive
s.scal_scaling	= 'free'		; 'fix' or 'free'
s.scal_scale	= [0.0,5.0]

s.name=s.vect_disp+s.scal_disp+'_'+s.type+'_'+s.panel_name

slice_array = [slice_array,s]
n_slices = n_slices+1

;-----------------------------------------------------------------------------

s.sw	= 'on' 
s.panel_name	= 'b'			; If more slices of the same type 
                                        ; are needed use this index

s.type		= 'xy'			; Chose 'yz', 'xz' or 'xy' plane
s.coord		=  0.0                  ; Position at the complementary coordinate
s.vect_disp	= 'b'			; Vector field to display: 'b' or 'v'
s.vect_scaling	= 'free'			; 'fix' or 'free' 
s.vect_scale	=  1.0


s.scal_disp	= 'e'			; Scalar field to display 'd' or 'e'
s.scal_pert     = ''			; inactive
s.scal_scaling	= 'free'		; 'fix' or 'free'
s.scal_scale	= [0.0,1.0e3]

s.name=s.vect_disp+s.scal_disp+'_'+s.type+'_'+s.panel_name


slice_array = [slice_array,s]
n_slices = n_slices+1

;==============================================================================
vars = ['den1','vlx1','vly1','vlz1','ene1','magx','magy','magz']
n_vectors = [n_vect_x,n_vect_y,n_vect_z]
;==============================================================================


  COORDS, coord_sys
  DIRS, data_dir, png_dir, png_output
  PLOT_SIZES, data_files(0), win_vert_size, coord_sys
  WINDOWS_OPEN,  slice_array,  display_frames
 
  FOR i_frame = first_frame, last_frame, freq_frame DO BEGIN

    file = data_files(i_frame)


    READ_DATA, data_dir,prefix,i_frame , vars, time, log_scal

    PRINT, FORMAT='(a6,1x,a,1x,a10,e10.3)', $
                'file =', file, '   time =', time

    file_sep = STR_SEP(file,'.')
    n_sep = N_ELEMENTS(file_sep)
    hdf_num = file_sep(n_sep-1)

    PLOT_SLICES,  slice_array,n_vectors, $
                  display_frames, hdf_num, png_output, png_dir, time, i_frame

  END

SKIP:

END
