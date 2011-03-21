pro test_voigts, STRCT=strct

  if not keyword_set( PSFILE ) then psfile = 'test_voigts.ps'
  if not keyword_set(CSZ) then csz = 2.0
  if not keyword_set(lSZ) then lsz = 2.0
  if not keyword_set(XNCOLORS) then xncolors=200L

  c = x_constants()

  ;; Voigt a
  v_doppler = 20*1e5 ; cm/s
  Dnu   = v_doppler * 2.466e15 / 2.9979e10 ;   
  voigt_a   = 6.265e8 / (4 * 3.14159 * Dnu) ;   // gamma/(4 pi Dnu)
  print, 'Voigt a = ', voigt_a

  ;; Plot
  x_psopen, psfile, /maxs
  clr = getcolor(/load)

  ;; Profiles
  xmrg = [9,1]
  ymrg = [4.0,1]
  yrng=[1e-10, 1]
  xrng=[-100, 100]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, $
        ytitle='Voigt', $
        xtitle='x', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, /ylog

  ;; x
  npt = 1000L
  xval = xrng[0] + (xrng[1]-xrng[0])*dindgen(npt)/(npt-1)

  ;; IDL
  voigt_IDL = voigt(voigt_a, xval)
  oplot, xval, voigt_IDL, color=clr.black


  ;; Label
  xlbl = 12.
;  xyouts, xlbl, 1e2, 'v!dr!N', color=clr.blue, charsiz=lsz

  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]

;  x_splot, vel, tau, /bloc
;  stop

  return

end
