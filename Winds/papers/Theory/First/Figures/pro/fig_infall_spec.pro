pro fig_infall_spec, RREAL=rreal

  if not keyword_set( PSFILE ) then psfile = 'fig_infall_spec.ps'
  if not keyword_set(PAD_FRAC) then pad_frac = 0.1
  if not keyword_set(CSZ) then csz = 1.3
  if not keyword_set(CSZ2) then csz2 = 2.3
  if not keyword_set(lSZ) then lsz = 1.5
  if not keyword_set(lSZ2) then lsz2 = 1.3
  if not keyword_set(XNCOLORS) then xncolors=200L
  if not keyword_set(YSTP) then ystp = 0.07

  xlbl = 0.05
  xlbl2 = 0.05
  ylbl = 0.90

  close, /all

  ;;; BEGIN PLOTS
  x_psopen, psfile, /portrait
  !p.multi = [0,1,2]
  clr = getcolor(/load)
  clrs = x_setclrs()

  xmrg = [8,4]
  ymrg = [3.0,3.5]

  ;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; MgII
  !p.multi = [1,1,2]
  ;; Plot MgII
  yrng=[-0.1, 2.8]
  xrng=[2786., 2809.8]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, ytitle='Relative Flux', $
        xtitle='Wavelength (Ang)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=9, psym=1, /nodata

  ;; Fiducial
  Mg_fil = '../Analysis/Fiducial/Output/spec_MgII_fiducial.dat'
  readcol, Mg_fil, wv, fx, /silen
  nrm = median(fx[where(wv GT 2815)])
  oplot, wv, fx/nrm, color=clr.black, psym=10, thick=3

  ;; Radiation
  Mg_fil = '../Analysis/Radiation/Output/spec_rad_MgII.dat'
  readcol, Mg_fil, wv, fx, noscatt_fx, /silen
  nrm = median(fx[where(wv GT 2815)])
  oplot, wv, fx/nrm, color=clr.red, psym=10, thick=3

  nrm = median(noscatt_fx[where(wv GT 2815)])
  oplot, wv, noscatt_fx/nrm, color=clr.black, psym=10, thick=3, linesty=1

  xrng2 = (xrng/2796.352 - 1)*3e5
  axis, xaxis=1, charsiz=csz, xsty=1, xrang=xrng2, xtitl='Velocity (km/s) Relative to MgII 2796'

  oplot, replicate(2796.352,2), yrng, color=clr.orange, linesty=2, thick=2
  oplot, replicate(2803.531,2), yrng, color=clr.orange, linesty=2, thick=2
  oplot, [-9e9,9e9], [0.,0.], color=clr.green, linesty=2, thick=2
  xyouts, xrng[0]+xlbl*(xrng[1]-xrng[0]), yrng[1]*ylbl, $
          'MgII', color=clr.black, charsiz=lsz
;  oplot, xrng, [1., 1.], color=clr.red, linestyle=1, thick=1

  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]
  close, /all

  return

end
