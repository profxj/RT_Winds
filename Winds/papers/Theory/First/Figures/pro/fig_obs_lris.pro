pro fig_obs_lris, RREAL=rreal

  if not keyword_set( PSFILE ) then psfile = 'fig_obs_lris.ps'
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
  x_psopen, psfile, /maxs
  !p.multi = [0,2,1]
  clr = getcolor(/load)
  clrs = x_setclrs()

  xmrg = [8,4]
  ymrg = [3.0,3.5]

  ;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; MgII

  ;; Plot MgII
  yrng=[-0.1, 1.8]
  xrng=[2786., 2809.8]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, ytitle='Relative Flux', $
        xtitle='Wavelength (Ang)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=9, psym=1, /nodata

  ;; ISM
  Mg_fil = '../Analysis/ISM/Output/spec_ISM_MgII.dat'
  readcol, Mg_fil, wv, fx, noscatt_fx, /silen
  nrm = median(fx[where(wv GT 2815)])
  fx = fx/nrm
  dwv = wv[1]-wv[0]
  npix = n_elements(wv)

  ;; Smooth
  fwhm_pix = 2.3 / dwv  ;; Corresponds to 250 km/s FWHM
  nsmooth = fwhm_pix/(2.*sqrt(2*alog(2)))
  kernel = gauss_kernel(nsmooth)
  smth = convol(fx, kernel,/edge_wrap)

  ;; Rebin
  nlow = npix/10
  low_fx = congrid(smth, nlow)
  low_wv = congrid(wv, nlow)

  ;; Noise
  low_fx = x_addnoise(low_fx, 7., seed=-2211L) 

  oplot, low_wv, low_fx, color=clr.black, psym=10, thick=3

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
