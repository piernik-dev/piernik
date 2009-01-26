invpal,0 
;set_plot,'Z' 
for i = 401,624,1 do begin
   velmag,16,16,i,cr
   ce = cr + 511 
   str = STRING(i,FORMAT='(I4.4)')+'.png' 
   im = TVRD(TRUE=1) 
   write_png,str,im[*,cr:ce,cr:ce]
endfor
;set_plot,'X' 
end 
