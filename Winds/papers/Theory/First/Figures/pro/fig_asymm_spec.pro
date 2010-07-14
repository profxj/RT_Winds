pro fig_asymm_spec, RREAL=rreal

  if not keyword_set( PSFILE ) then psfile = 'fig_asymm_spec.ps'
  if not keyword_set(PAD_FRAC) then pad_frac = 0.1
  if not keyword_set(CSZ) then csz = 1.3
  if not keyword_set(CSZ2) then csz2 = 2.3
  if not keyword_set(lSZ) then lsz = 1.5
  if not keyword_set(lSZ2) then lsz2 = 1.3
  if not keyword_set(XNCOLORS) then xncolors=200L

  if not keyword_set(FeII_grid) then FeII_grid =  '../Analysis/Asymmetric/feII_grids.fits'
  if not keyword_set(MgII_grid) then MgII_grid = '../Analysis/Asymmetric/mgII_grids.fits'
  if not keyword_set( GRID_FILE ) then $
    grid_file = '../Analysis/Outputs/fiducial_grid.fits'

  if not keyword_set(ALL_ANG) then all_ang = [0, 120, 150, 180, 30, 60, 90]
  if not keyword_set(ANGLES) then angles = [0, 180, 90, 120]
  if not keyword_set(FE_ANGLES) then fe_angles = [0, 180, 90]

  xlbl = 0.3
  xlbl2 = 0.06
  ylbl = 0.90

  close, /all

  ;;; BEGIN PLOTS
  x_psopen, psfile, /portrait
  !p.multi = [0,1,2]
  clr = getcolor(/load)
  clrs = x_setclrs(/dark)
  tmp = clrs[1]
  clrs[1] = clrs[2]
  clrs[2] = tmp

  xmrg = [8,3]
  ymrg = [3.0,3.5]


  xrng=[2580, 2618]
  xcut = 2605.
  off = 1.
  
  for ss=0,1 do begin
     ;; Plot FeII
     case ss of 
        0: begin
           yrng=[0., 1.8] 
           ysty = 9
           wvmnx = [xrng[0], xcut-off]
        end
        1: begin
              yrng=[0.95,1.25]
              ysty = 5
              wvmnx = [xcut+off,xrng[1]]
           end
        else: stop
     endcase
     
     plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
           xmargin=xmrg, ymargin=ymrg, ytitle='Relative Flux', $
           xtitle='Wavelength (Ang)', yrange=yrng, thick=4, $
           xrange=xrng, ystyle=ysty, xstyle=9, psym=1, /nodata, /NOERASE
     if ss EQ 1 then axis, yaxis=1, charsiz=csz, ysty=1, xrang=yrng, ytickint=0.1
     
     ;; Fiducial
     feII_data = xmrdfits(grid_file, 0, /silent)
     feII_wave = xmrdfits(grid_file, 1, /silent)
     dwv_feII = abs(feii_wave[1]-feii_wave[0])
     feII_data = float(feII_data)
     sz_feii = size(feII_data,/dimen)
     
     spec_feii = total(total(feii_data,1),1)
     pix = where(feii_wave GT wvmnx[0] and feii_wave LT wvmnx[1])
     oplot, feii_wave[pix], spec_feii[pix], color=clr.black, psym=10, thick=3
     ystp = 0.07
     if ss EQ 0 then $
        xyouts, xrng[0]+xlbl2*(xrng[1]-xrng[0]), yrng[1]-ystp*(yrng[1]-yrng[0]), $
                'Fiducial', color=clr.black, charsiz=lsz
     
     nFe = n_elements(fe_angles)
     feII_wave = xmrdfits(FeII_grid, 0, /silent)
     dwv_feII = abs(feii_wave[1]-feii_wave[0])
     pix = where(feii_wave GT wvmnx[0] and feii_wave LT wvmnx[1])
     for kk=0L,nFe-1 do begin
        jj = where(fe_angles[kk] EQ all_ang)
        feII_data = xmrdfits(FeII_grid, jj[0]+1, /silent)
        feII_data = float(feII_data)
        sz_feii = size(feII_data,/dimen)
        
        spec_feii = total(total(feii_data,1),1)
        ;;
        oplot, feII_wave[pix], spec_feii[pix], color=clrs[kk+1], psym=10, thick=3
        ;;
        if ss EQ 0 then $
           xyouts, xrng[0]+xlbl2*(xrng[1]-xrng[0]), yrng[1]-(kk+2)*ystp*(yrng[1]-yrng[0]), $
                   'Angle='+strtrim(fe_angles[kk],2), color=clrs[kk+1], charsiz=lsz
     endfor
     
  endfor
     
  xrng2 = (xrng/2600.173 - 1)*3e5
  axis, xaxis=1, charsiz=csz, xsty=1, xrang=xrng2, xtitl='Velocity (km/s) Relative to FeII 2600'

  oplot, replicate(xcut,2), yrng, color=clr.orange, linesty=2
  oplot, replicate(2586.650,2), yrng, color=clr.orange, linesty=2, thick=2
  oplot, replicate(2600.173,2), yrng, color=clr.orange, linesty=2, thick=2
  oplot, replicate(2612.6542,2), yrng, color=clr.orange, linesty=2, thick=2
  xyouts, xrng[0]+xlbl*(xrng[1]-xrng[0]), yrng[0]+(yrng[1]-yrng[0])*ylbl, $
          'FeII', color=clr.black, charsiz=lsz
;  oplot, xrng, [1., 1.], color=clr.red, linestyle=1, thick=1

  x_curvefill, [xcut-off,xcut+off], [0., 0.], [10., 10], color=clr.tan

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Plot MgII
  !p.multi = [1,1,2]

  yrng=[0., 2.5]
  xrng=[2786., 2809.8]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, ytitle='Relative Flux', $
        xtitle='Wavelength (Ang)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=9, psym=1, /nodata
  
  ;; Fiducial
  mgII_data = xmrdfits(grid_file, 2, /silent)
  mgII_wave = xmrdfits(grid_file, 3, /silent)
  dwv_mgII = abs(mgii_wave[1]-mgii_wave[0])
  mgII_data = float(mgII_data)
  sz_mgii = size(mgII_data,/dimen)
  
  spec_mgii = total(total(mgii_data,1),1)
  oplot, mgii_wave, spec_mgii, color=clr.black, psym=10, thick=3
  ystp = 0.07
  xyouts, xrng[0]+xlbl2*(xrng[1]-xrng[0]), yrng[1]-ystp*(yrng[1]-yrng[0]), $
          'Fiducial', color=clr.black, charsiz=lsz
  
  
  nMg = n_elements(angles)
  mgII_wave = xmrdfits(MgII_grid, 0, /silent)
  dwv_mgII = abs(mgii_wave[1]-mgii_wave[0])
  for kk=0L,nMg-1 do begin
     jj = where(angles[kk] EQ all_ang)
     mgII_data = xmrdfits(MgII_grid, jj[0]+1, /silent)
     mgII_data = float(mgII_data)
     sz_mgii = size(mgII_data,/dimen)
     
     spec_mgii = total(total(mgii_data,1),1)
     ;;
     oplot, mgII_wave, spec_mgii, color=clrs[kk+1], psym=10, thick=3
     ;;
     xyouts, xrng[0]+xlbl2*(xrng[1]-xrng[0]), yrng[1]-(kk+2)*ystp*(yrng[1]-yrng[0]), $
             'Angle='+strtrim(angles[kk],2), color=clrs[kk+1], charsiz=lsz
  endfor
  
  xrng2 = (xrng/2796.352 - 1)*3e5
  axis, xaxis=1, charsiz=csz, xsty=1, xrang=xrng2, xtitl='Velocity (km/s) Relative to MgII 2796'

  oplot, replicate(2796.352,2), yrng, color=clr.orange, linesty=2, thick=2
  oplot, replicate(2803.531,2), yrng, color=clr.orange, linesty=2, thick=2
  xyouts, xrng[0]+xlbl*(xrng[1]-xrng[0]), yrng[1]*ylbl, $
          'MgII', color=clr.black, charsiz=lsz
  oplot, xrng, [1., 1.], color=clr.red, linestyle=1, thick=1

  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]
  close, /all

  return

end
