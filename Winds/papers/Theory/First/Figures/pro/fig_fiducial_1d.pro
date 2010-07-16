pro fig_fiducial_1d, RREAL=rreal

  if not keyword_set( PSFILE ) then psfile = 'fig_fiducial_1d.ps'
  if not keyword_set( GRID_FILE ) then $
    grid_file = '../Analysis/Outputs/fiducial_grid.fits'
  if not keyword_set(PAD_FRAC) then pad_frac = 0.1
  if not keyword_set(CSZ) then csz = 1.3
  if not keyword_set(CSZ2) then csz2 = 2.3
  if not keyword_set(lSZ) then lsz = 1.5
  if not keyword_set(lSZ2) then lsz2 = 1.3
  if not keyword_set(XNCOLORS) then xncolors=200L

  if not keyword_set(MGII_REST) then mgii_rest = 2796.35d

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Read in MgII Data
  mgII_data = xmrdfits(grid_file, 2, /silent)
  mgII_wave = xmrdfits(grid_file, 3, /silent)
  dwv_mgII = abs(mgii_wave[1]-mgii_wave[0])
  mgII_data = float(mgII_data)
  sz_mgii = size(mgII_data,/dimen)

  spec_mgii = total(total(mgii_data,1),1)

  ;;; BEGIN PLOTS
  x_psopen, psfile, /maxs
  clr = getcolor(/load)

  ;; Plot MgII
  yrng=[0., 2.5]
  xrng=[2786., 2812]
  pos=[0.08, 0.6, 0.48, 0.95]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        ytitle='Normalized Flux', $
        xtitle='Wavelength (Ang)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, $
        pos=pos

  oplot, mgii_wave, spec_mgii, color=clr.black, psym=10, thick=3

  oplot, replicate(2796.352,2), yrng, color=clr.gray, linesty=2
  oplot, replicate(2803.531,2), yrng, color=clr.gray, linesty=2
  xlbl = 0.5
  ylbl = 0.85
  xyouts, xrng[0]+xlbl*(xrng[1]-xrng[0]), yrng[1]*ylbl, $
          'MgII', color=clr.black, charsiz=lsz
  oplot, xrng, [1., 1.], color=clr.red, linestyle=1, thick=1

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
  xrng=[2580, 2635]
  pos[1]=0.17
  pos[3]=0.52
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        ytitle='Normalized Flux', $
        xtitle='Wavelength (Ang)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, $
        pos=pos, /noeras

  oplot, feii_wave, spec_feII, color=clr.black, psym=10, thick=3

  oplot, replicate(2586.650,2), yrng, color=clr.gray, linesty=2
  oplot, replicate(2600.173,2), yrng, color=clr.gray, linesty=2
  xyouts, xrng[0]+xlbl*(xrng[1]-xrng[0]), yrng[1]*ylbl, $
          'FeII', color=clr.black, charsiz=lsz
  oplot, xrng, [1., 1.], color=clr.red, linestyle=1, thick=1


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Velocity plot
  trans = [2796.352, 2803.531,2586.650, 2600.1729, 2612.6542d]
  ntrans = n_elements(trans)
  !p.multi = [ntrans+1,2,ntrans+1,0,1]

  xmrg = [8,3]
  ymrg = [0,0]
  xrng=[-1100, 500]
  xtit = ''

  for qq=0L,ntrans-1 do begin

     ;; Load up
     if trans[qq] LT 2700. then begin
        wave = feii_wave
        spec = spec_feii
        if trans[qq] GT 2605 then yrng=[0.91, 1.24] else yrng=[0., 1.8]
     endif else begin
        wave = mgii_wave
        spec = spec_mgii
        yrng=[0., 2.4]
     endelse

     if qq EQ (ntrans-1) then delvarx, xspaces $ 
     else xspaces = replicate(' ',30) 
     if qq EQ (ntrans-1) then xtit = 'Relative Velocity (km s!u-1!N)'
     if qq EQ (ntrans-1) then ytki=0.1

     plot, [0], [0], color=clr.black, background=clr.white, charsize=csz2,$
           xmargin=xmrg, ymargin=ymrg, $
           xtitle=xtit, yrange=yrng, thick=4, $
           xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, $
           xtickn=xspaces, ytickinter=ytki

     vel = (wave-trans[qq])/trans[qq] * 3e5  ;; km/s
     oplot, vel, spec, color=clr.black, psym=10, thick=3
     ;; Lines and labels
     oplot, [0., 0.], yrng, color=clr.blue, linestyle=2, thick=1
;     oplot, xrng, [0., 0.], color=clr.green, linestyle=1, thick=1
     oplot, xrng, [1., 1.], color=clr.red, linestyle=1, thick=1

     getfnam, trans[qq], fv, nam 
     nam = strtrim(nam,2)
     xyouts, 0.70*(xrng[1]-xrng[0])+xrng[0], $
             yrng[0]+ (yrng[1]-yrng[0])*0.10, $
             strtrim(nam,2), color=clr.black, charsize=LSZ2
  endfor
  xyouts, 0.53, 0.55, 'Normalized Flux', align=0.5, orient=90, $
          color=clr.black, charsize=LSZ2, /NORM
  
  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]

  return

end
