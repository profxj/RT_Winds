; mk_spec_fits, 'Outputs/fiducial_feII.spec', 'Outputs/fiducial_mgII.spec', 'Outputs/fiducial_grid.fits'
; mk_spec_fits, 'ISM/Output/fe_ism.spec',  'ISM/Output/mg_ism.spec', 'ISM/Output/ism_grid.fits'
; mk_spec_fits,  'Radiation/Output/wind_grid_radiative.asc_fe.spec',  'Radiation/Output/wind_grid_radiative.asc_mg.spec',  'Radiation/Output/radiative_grid.fits', /nochk
pro mk_spec_fits, feii_file, mgii_file, file_out, NOCHK=nochk

for qq=0,1 do begin

   case qq of 
      0: filen = feii_file
      1: filen = mgii_file
      else: stop
   endcase

   nlam = 1
   nx   = 1
   ny   = 1
   
   openr,1,filen 
   readf,1,nlam,nx,ny
   readf,1,lam0,dlam,dx
   
   x_arr = dindgen(nx)*dx
   y_arr = x_arr
   wave  = lam0 + dindgen(nlam)*dlam
   
   data  = 0*dblarr(nx,ny,nlam)
   image = 0*dblarr(nx,ny)
   flux  = 0*dblarr(nlam)
   L_tot = 0.0d
   
   w = 0.0d
   f = 0.0d
   x = 0.0d
   y = 0.0d
   c = 0.0d
   
   for i=0,nx-1 do begin
      for j=0,ny-1 do begin
         for k=0,nlam-1 do begin
            
            readf,1,f,c
            
            data(i,j,k) = f
            image(i,j) += f
            flux(k) += f
            L_tot += f
            
         endfor
      endfor
   endfor

   close,1

  ;; renormalize
   nrm = 1.d
   if (nlam gt 1) then begin
      L_tot *= (wave(1) - wave(0))
      nrm *= (wave(1) - wave(0))*nlam
   endif
   if (nx gt 1) then begin
      L_tot *= (x_arr(1) - x_arr(0))^2.0
      nrm *= (x_arr(1) - x_arr(0))^2.0
   endif
   print,'L_tot = ',L_tot
   print,'nrm = ', nrm
   
   x_splot,wave,flux, /bloc
;print,'integrated spectrum, press return to continue'
;junk =''
;read,junk

   loadct,13
   if not keyword_set(NOCHK) then tvimage,bytscl(image)
   
  ;; Write to FITS
   mwrfits, data*nrm, file_out, create=(qq EQ 0)
   mwrfits, wave, file_out
endfor

spawn, 'gzip -f '+file_out

print, 'mk_spec_fits:  All done'

end
  
