pro fig_noemiss, GRID_FILE, PSFILE

  if not keyword_set( PSFILE ) then psfile = 'fig_noemiss.ps'
  if not keyword_set( GRID_FILE ) then $
    grid_file = '../Analysis/Outputs/fiducial_grid.fits'
  if not keyword_set(CSZ) then csz = 1.3
  if not keyword_set(CSZ2) then csz2 = 2.3
  if not keyword_set(lSZ) then lsz = 1.5
  if not keyword_set(lSZ2) then lsz2 = 1.3
  if not keyword_set(XNCOLORS) then xncolors=200L
  if not keyword_set(YRNG) then yrng=[0., 1.5]

  if not keyword_set(MGII_REST) then mgii_rest = 2796.35d

  ;;; BEGIN PLOTS
  x_psopen, psfile, /portrait
  !p.multi = [0,1,2]
  clr = getcolor(/load)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; MgII first
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  mgII_data = xmrdfits(grid_file, 2, /silent)
  mgII_wave = xmrdfits(grid_file, 3, /silent)
  dwv_mgII = abs(mgii_wave[1]-mgii_wave[0])
  mgII_data = float(mgII_data)
  sz_mgii = size(mgII_data,/dimen)

  spec_mgii = total(total(mgii_data,1),1)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Plot MgII
  xrng=[2786., 2810]
  xmrg = [10, 1]
  ymrg = [5, 1]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmarg=xmrg, ymarg=ymrg, $
        ytitle='Normalized Flux', $
        xtitle='Wavelength (Ang)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata

  oplot, mgii_wave, spec_mgii, color=clr.black, psym=10, thick=3

  oplot, replicate(2796.352,2), yrng, color=clr.gray, linesty=2
  oplot, replicate(2803.531,2), yrng, color=clr.gray, linesty=2
  xlbl = 0.5
  ylbl = 0.85
  xyouts, xrng[0]+xlbl*(xrng[1]-xrng[0]), yrng[1]*ylbl, $
          'MgII', color=clr.black, charsiz=lsz
  oplot, xrng, [1., 1.], color=clr.red, linestyle=1, thick=1

  ;; No emission
  fig_nvtau_vs_r, strct=strct
  wrest = strct.wrest
  wave = wrest - wrest * strct.vel / 3e5
  oplot, wave, exp(-1.*(strct.tau < 10)), color=clr.red, psym=10

  ;; 2803
  wrest = 2803.531d
  wave = wrest - wrest * strct.vel / 3e5
  tau = strct.tau / 2
  oplot, wave, exp(-1.*(tau < 10)), color=clr.red, psym=10

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Read in FeII Data
  feII_data = xmrdfits(grid_file, 0, /silent)
  feII_wave = xmrdfits(grid_file, 1, /silent)
  dwv_feII = abs(feii_wave[1]-feii_wave[0])
  feII_data = float(feII_data)
  sz_feii = size(feII_data,/dimen)

  spec_feII = total(total(feii_data,1),1)

  ;; Plot FeII
  yrng=[0., 1.8]
  xrng=[2580, 2605]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        ytitle='Normalized Flux', $
        xtitle='Wavelength (Ang)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata

  oplot, feii_wave, spec_feII, color=clr.black, psym=10, thick=3

  oplot, replicate(2586.650,2), yrng, color=clr.gray, linesty=2
  oplot, replicate(2600.173,2), yrng, color=clr.gray, linesty=2
  xyouts, xrng[0]+xlbl*(xrng[1]-xrng[0]), yrng[1]*ylbl, $
          'FeII', color=clr.black, charsiz=lsz
  oplot, xrng, [1., 1.], color=clr.red, linestyle=1, thick=1

  ;; No emission
  wrest = 2586.650
  getfnam, wrest, f, nam
  wave = wrest - wrest * strct.vel / 3e5
  tau = strct.tau * wrest / strct.wrest * f / strct.fval
  oplot, wave, exp(-1.*(tau < 10)), color=clr.red, psym=10

  wrest = 2600.1729d
  getfnam, wrest, f, nam
  wave = wrest - wrest * strct.vel / 3e5
  tau = strct.tau * wrest / strct.wrest * f / strct.fval
  oplot, wave, exp(-1.*(tau < 10)), color=clr.red, psym=10

  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]

  return

end
