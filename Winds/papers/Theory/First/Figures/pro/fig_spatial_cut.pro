pro fig_spatial_cut, wrest, dv, grid_file, psfile, XRNG=xrng, YRNG=yrng, $
                     LBL=lbl, YMRG=ymrg

  if not keyword_set(CSZ) then csz = 1.7
  if not keyword_set(CSZ2) then csz2 = 1.1
  if not keyword_set(lSZ) then lsz = 1.5
  if not keyword_set(XNCOLORS) then xncolors=200L
  if not keyword_set(XRNG) then xrng=[-20, 20]  ; kpc
  if not keyword_set(YRNG) then yrng=[1e-3, 1e3]

  ;; 
  rsum = [ [-10., 10], [-5., 5], [-2, 2], [-1, 1]]
  szr = size(rsum, /dimens)

  ;; Read Data
  if wrest GT 2700. then idx = [2,3] else idx = [0,1]
  raw_data = xmrdfits(grid_file, idx[0], /silent)
  raw_wave = xmrdfits(grid_file, idx[1], /silent)
  dwv = abs(raw_wave[1]-raw_wave[0])
  raw_data = float(raw_data)
  sz_raw = size(raw_data,/dimen)


  xval = xrng[0] + (xrng[1]-xrng[0])*findgen(sz_raw[0])/float(sz_raw[0]-1)

  img_wave = wrest + dv/3e5 * wrest
  ncut = n_elements(dv)

  if keyword_set(PSFILE) then x_psopen, psfile, /maxs
  clr = getcolor(/load)
  clrs = x_setclrs(/dark)
  xmrg = [9,1]
  if not keyword_set(YMRG) then ymrg = [5.0,1]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, ytitle='Relative Flux', $
        xtitle='Distance (kpc)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, /ylog

  ;; Continuum
  spec = total(total(raw_data,1),1)
  conti = spec[0]

  for qq=0,ncut-1 do begin

      ;; Image
      mn = min(abs(raw_wave-img_wave[qq]),img_pix)
      data = raw_data[*,*,img_pix]

      ;; Smash spatially
      spec_spatial = total(data,1)
      oplot, xval, spec_spatial, psym=10, color=clrs[qq]

      xyouts, xrng[0]+0.1*(xrng[1]-xrng[0]), yrng[1]/(2.0^(qq+2)), $
              'v='+strtrim(round(dv[qq]),2)+' km s!u-1!N', $
              color=clrs[qq], charsiz=lsz

      ;; Do some calculations
      tot_flux = total(spec_spatial)
      tot_flux_conti = total(spec_spatial) - conti
;      if tot_flux_conti GT 0. then stop
      for ii=0L,szr[1]-1 do begin
         pix = where(xval GE rsum[0,ii] and xval LE rsum[1,ii])
         fx = total(spec_spatial[pix]) / tot_flux
         if tot_flux_conti GT 0. then $
            fx_c = (total(spec_spatial[pix])-conti) / tot_flux_conti  $
         else fx_c = -9.99
         print, dv[qq], rsum[0,ii], rsum[1,ii], fx, fx_c
      endfor
   endfor

  if keyword_set(LBL) then $
     xyouts, xrng[0]+0.7*(xrng[1]-xrng[0]), yrng[1]/4.0, $
             LBL, color=clr.black, charsiz=lsz
  
  if keyword_set( PSFILE ) then begin 
     x_psclose
     !p.multi = [0,1,1]
  endif

  return

end
