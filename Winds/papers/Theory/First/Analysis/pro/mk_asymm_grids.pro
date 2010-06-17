; mk_spec_fits, 'Outputs/fiducial_feII.spec', 'Outputs/fiducial_mgII.spec', 'Outputs/fiducial_grid.fits'
pro mk_asymm_grids

 mg_fils = findfile('Asymmetric/angle_mg*',count=nmg_fil)
 fe_fils = findfile('Asymmetric/angle_fe*',count=nfe_fil)
 file_out = 'Asymmetric/mgII_grids.fits'
 file2_out = 'Asymmetric/feII_grids.fits'

  for ss=0,1 do begin
      case ss of 
         0: begin
            nfil = nmg_fil
            outfile = file_out
            files = mg_fils
         end
         1: begin
            nfil = nfe_fil
            outfile = file2_out
            files = fe_fils
         end
      endcase

      for qq=0,nfil-1 do begin
         
         filen = files[qq]
         print, 'Reading ', filen
         
         nlam = 1
         nx   = 1
         ny   = 1
         
         close, 1
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
         ;; Write to FITS
         if qq EQ 0 then mwrfits, wave, outfile, /create
         mwrfits, data*nrm, outfile
      endfor
      
      spawn, 'gzip -f '+outfile
   endfor

  print, 'mk_asymm_fits:  All done'

end
  
