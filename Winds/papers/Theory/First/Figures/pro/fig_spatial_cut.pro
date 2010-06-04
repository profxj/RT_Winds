pro fig_spatial_cut, wrest, dv, grid_file, psfile, XRNG=xrng, YRNG=yrng

  if not keyword_set(CSZ) then csz = 1.7
  if not keyword_set(CSZ2) then csz2 = 1.1
  if not keyword_set(lSZ) then lsz = 1.5
  if not keyword_set(XNCOLORS) then xncolors=200L
  if not keyword_set(XRNG) then xrng=[-20, 20]  ; kpc
  if not keyword_set(YRNG) then yrng=[1e-3, 1e3]

  ;; Read in MgII Data
  if wrest GT 2700. then idx = [2,3] else idx = [0,1]
  raw_data = xmrdfits(grid_file, idx[0], /silent)
  raw_wave = xmrdfits(grid_file, idx[1], /silent)
  dwv = abs(raw_wave[1]-raw_wave[0])
  raw_data = float(raw_data)
  sz_raw = size(raw_data,/dimen)

  xval = xrng[0] + (xrng[1]-xrng[0])*findgen(sz_raw[0])/float(sz_raw[0]-1)

  img_wave = wrest + dv/3e5 * wrest
  ncut = n_elements(dv)

  x_psopen, psfile, /maxs
  clr = getcolor(/load)
  clrs = x_setclrs()
  xmrg = [9,1]
  ymrg = [5.0,1]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, ytitle='Relative Flux', $
        xtitle='r (kpc)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, /ylog

  spec = total(total(raw_data,1),1)
  for qq=0,ncut-1 do begin

      ;; Image
      mn = min(abs(raw_wave-img_wave[qq]),img_pix)
      data = raw_data[*,*,img_pix]

      ;; Smash spatially
      spec_spatial = total(data,1)
      oplot, xval, spec_spatial, psym=10, color=clrs[qq]

      xyouts, xrng[0]+0.1*(xrng[1]-xrng[0]), yrng[1]/(2.0^(qq+1)), $
              'v='+strtrim(round(dv[qq]),2)+' km s!u-1!N', $
              color=clrs[qq], charsiz=lsz
  endfor
  
  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]

  return

end
