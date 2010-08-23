pro fig_obs_lris, FE=fe

  if not keyword_set( PSFILE ) then psfile = 'fig_obs_lris.ps'
  if keyword_set(FE) then psfile = 'fig_obs_lris_fe.ps'
  if not keyword_set(PAD_FRAC) then pad_frac = 0.1
  if not keyword_set(CSZ) then csz = 1.6
  if not keyword_set(CSZ2) then csz2 = 2.3
  if not keyword_set(lSZ) then lsz = 1.7
  if not keyword_set(lSZ2) then lsz2 = 1.3
  if not keyword_set(XNCOLORS) then xncolors=200L
  if not keyword_set(YSTP) then ystp = 0.07
  if not keyword_set(FE) then wrest= 2796.352 else wrest = 2600.173

  xlbl = 0.05
  xlbl2 = 0.05
  ylbl = 0.90

  close, /all

  ;;; BEGIN PLOTS
  x_psopen, psfile, /maxs
  !p.multi = [0,2,1]
  clr = getcolor(/load)
  clrs = x_setclrs()

  xmrg = [8,1]
  ymrg = [4.0,3.5]

  ;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; MgII

  ;; Plot MgII
  if not keyword_set(FE) then begin
     Mg_fil = '../Analysis/ISM/Output/spec_ISM_MgII.dat' 
     xtit='Velocity (km/s) Relative to MgII 2796'
     yrng=[-0.1, 2.8]
     xrng=[2785., 2811.] 
  endif else begin
     stop
     Mg_fil = '../Analysis/Fiducial/Output/spec_FeII_fiducial.dat'
     xtit='Velocity (km/s) Relative to FeII 2600'
     yrng=[-0.1, 1.8]
     xrng = [2580., 2605]
  endelse
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, ytitle='Relative Flux', $
        xtitle='Wavelength (Ang)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=9, psym=1, /nodata

  ;; ISM

  readcol, Mg_fil, wv, fx, noscatt_fx, /silen
  if not keyword_set(FE) then nrm = median(fx[where(wv GT 2815)]) $
  else nrm = median(fx[where(wv GT 2634)])
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
  print, 'dlambda = ', low_wv[1]-low_wv[0]

  ;; Noise
  seed = -2211L
  low_fx = x_addnoise(low_fx, 7., seed=seed)

  ;; Plot
  oplot, wv, fx, color=clr.darkgray, linesty=1, thick=3
  oplot, low_wv, low_fx, color=clr.black, psym=10, thick=3

  xrng2 = (xrng/wrest - 1)*3e5
  axis, xaxis=1, charsiz=csz, xsty=1, xrang=xrng2, xtitl=xtit

  oplot, replicate(2796.352,2), yrng, color=clr.orange, linesty=2, thick=2
  oplot, replicate(2803.531,2), yrng, color=clr.orange, linesty=2, thick=2
  oplot, [-9e9,9e9], [0.,0.], color=clr.green, linesty=2, thick=2
  xyouts, xrng[0]+xlbl*(xrng[1]-xrng[0]), yrng[1]*ylbl, $
          '(a)', color=clr.black, charsiz=lsz
;  oplot, xrng, [1., 1.], color=clr.red, linestyle=1, thick=1
  
  ;; Stack spectrum
  nstack = 100L
  all_spec = fltarr(npix, nstack)
  shft = round(10*randomn(seed, nstack))

  for ii=0L,nstack-1 do all_spec[*,ii] = shift(fx,shft[ii])
  all_spec[*] = convol(all_spec[*], kernel,/edge_wrap)
  all_smooth = reform(rebin(all_spec[*], nlow*nstack), nlow, nstack)

  all_smooth = reform(rebin(all_spec[*], nlow*nstack), nlow, nstack)
  all_smooth = x_addnoise(all_smooth, 2., seed=seed)
  avg_spec = total(all_smooth,2) / nstack

  ;; Plot
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, ytitle='Relative Flux', $
        xtitle='Wavelength (Ang)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=9, psym=1, /nodata

  oplot, wv, fx, color=clr.darkgray, linesty=1, thick=3
  oplot, low_wv, avg_spec, color=clr.black, psym=10, thick=3

  ;; Label
  xrng2 = (xrng/wrest - 1)*3e5
  axis, xaxis=1, charsiz=csz, xsty=1, xrang=xrng2, xtitl=xtit

  if not keyword_set(FE) then begin
     oplot, replicate(2796.352,2), yrng, color=clr.orange, linesty=2, thick=2
     oplot, replicate(2803.531,2), yrng, color=clr.orange, linesty=2, thick=2
  endif else begin
     oplot, replicate(2586.650,2), yrng, color=clr.orange, linesty=2, thick=2
     oplot, replicate(2600.173,2), yrng, color=clr.orange, linesty=2, thick=2
  endelse
  oplot, [-9e9,9e9], [0.,0.], color=clr.green, linesty=2, thick=2

  xyouts, xrng[0]+xlbl*(xrng[1]-xrng[0]), yrng[1]*ylbl, $
          '(b)', color=clr.black, charsiz=lsz
  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]
  close, /all

  return

end
