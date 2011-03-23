; mk_spec_fits, 'Outputs/fiducial_feII.spec', 'Outputs/fiducial_mgII.spec', 'Outputs/fiducial_grid.fits'
pro view_output, spec_fil, SHOW_WAVE=show_wave 

  if not keyword_set(SPEC_FIL) then spec_fil = 'spec.dat'
  if not keyword_set(SHOW_WAVE) then show_wave = 1215.6701

  nlam = 1
  nx   = 1
  ny   = 1
  
  ;; Read
  close, /all
  openr,1,spec_fil 
  readf,1,nlam,nx,ny
  readf,1,lam0,dlam,dx
  close, /all
  rdfloat, spec_fil, f, c, skiplin=2
  
  x_arr = dindgen(nx)*dx
  y_arr = x_arr
  wave  = lam0 + dindgen(nlam)*dlam
  
  data  = 0*dblarr(nx,ny,nlam)
  image = 0*dblarr(nx,ny)
  flux  = 0*dblarr(nlam)
  L_tot = 0.0d
  
  ;; Fill
  w = 0.0d
  x = 0.0d
  y = 0.0d
  
  cnt = 0L
  for i=0,nx-1 do begin
     for j=0,ny-1 do begin
        for k=0,nlam-1 do begin
           
           data(i,j,k) = f[cnt]
           image(i,j) += f[cnt]
           flux(k) += f[cnt]
           L_tot += f[cnt]
           cnt += 1
           
        endfor
     endfor
  endfor
         
         
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

  ;; Show 1D 
;  mg_fil = '..//Output/spec_MgII_lbg_covering.dat'
;  readcol, mg_fil, wv, fx, noscatt_fx, /silen
;  nrm2 = median(fx[where(wv GT 2815)])
  wv = findgen(100)
  fx = fltarr(100)

  spec = total(total(data,1),1)
  pix = where(wave GT 1215.6701)
  nrm = median(spec[pix])
  sz_data = size(data,/dimen)
  x_splot, wave, spec/nrm, xtwo=wv, ytwo=fx, /blo

  ;; Show IFU
  mn = min(abs(wave-show_wave),img_pix)
  ifu_data = data[*,*,img_pix]
  xatv, ifu_data, /blo
  stop

end
  
