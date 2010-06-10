pro fig_asymm_spec, RREAL=rreal

  if not keyword_set( PSFILE ) then psfile = 'fig_asymm_spec.ps'
  if not keyword_set(PAD_FRAC) then pad_frac = 0.1
  if not keyword_set(CSZ) then csz = 1.3
  if not keyword_set(CSZ2) then csz2 = 2.3
  if not keyword_set(lSZ) then lsz = 1.5
  if not keyword_set(lSZ2) then lsz2 = 1.3
  if not keyword_set(XNCOLORS) then xncolors=200L

;  if not keyword_set(FeII_grid) then FeII_grid =  '../Analysis/Asymmetric/FeII_grids.fits'
  if not keyword_set(MgII_grid) then MgII_grid = '../Analysis/Asymmetric/mgII_grids.fits'
  if not keyword_set( GRID_FILE ) then $
    grid_file = '../Analysis/Outputs/fiducial_grid.fits'

  xlbl = 0.5
  xlbl2 = 0.06
  ylbl = 0.85

  close, /all

  ;;; BEGIN PLOTS
  x_psopen, psfile, /portrait
  !p.multi = [0,1,2]
  clr = getcolor(/load)
  clrs = x_setclrs()

  xmrg = [9,1]
  ymrg = [4.5,1]

  if keyword_set(FeII_grid) then begin
     nFe = 6

     ;; Plot FeII
     yrng=[0., 1.8]
     xrng=[2580, 2635]
     plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
           xmargin=xmrg, ymargin=ymrg, ytitle='Relative Flux', $
           xtitle='Wavelength (Ang)', yrange=yrng, thick=4, $
           xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata
     
     Fe_fil = ''
     for kk=0L,nFe-1 do begin
        readf, 1, Fe_fil
        readcol, Fe_fil, wv, fx, /silen
        ;;
        oplot, wv, fx, color=clrs[kk], psym=10, thick=3
     endfor
     
     oplot, replicate(2586.650,2), yrng, color=clr.gray, linesty=2
     oplot, replicate(2600.173,2), yrng, color=clr.gray, linesty=2
     xyouts, xrng[0]+xlbl*(xrng[1]-xrng[0]), yrng[1]*ylbl, $
             'FeII', color=clr.black, charsiz=lsz
     oplot, xrng, [1., 1.], color=clr.red, linestyle=1, thick=1
  endif

  if keyword_set(MgII_grid) then begin

     ;; Plot MgII
     yrng=[0., 2.5]
     xrng=[2786., 2810]
     plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
           xmargin=xmrg, ymargin=ymrg, ytitle='Relative Flux', $
           xtitle='Wavelength (Ang)', yrange=yrng, thick=4, $
           xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata

     ;; Fiducial
     mgII_data = xmrdfits(grid_file, 2, /silent)
     mgII_wave = xmrdfits(grid_file, 3, /silent)
     dwv_mgII = abs(mgii_wave[1]-mgii_wave[0])
     mgII_data = float(mgII_data)
     sz_mgii = size(mgII_data,/dimen)
     
     spec_mgii = total(total(mgii_data,1),1)
     oplot, mgii_wave, spec_mgii, color=clr.black, psym=10, thick=3
     ystp = 0.065
     xyouts, xrng[0]+xlbl2*(xrng[1]-xrng[0]), yrng[1]-ystp*(yrng[1]-yrng[0]), $
             'Fiducial', color=clr.black, charsiz=lsz


     angles = [0, 120, 150, 180, 30, 60, 90]
     nMg = n_elements(angles)
     mgII_wave = xmrdfits(MgII_grid, 0, /silent)
     dwv_mgII = abs(mgii_wave[1]-mgii_wave[0])
     for kk=0L,nMg-1 do begin
        mgII_data = xmrdfits(MgII_grid, kk+1, /silent)
        mgII_data = float(mgII_data)
        sz_mgii = size(mgII_data,/dimen)

        spec_mgii = total(total(mgii_data,1),1)
        ;;
        oplot, mgII_wave, spec_mgii, color=clrs[kk+1], psym=10, thick=3
        ;;
        xyouts, xrng[0]+xlbl2*(xrng[1]-xrng[0]), yrng[1]-(kk+2)*ystp*(yrng[1]-yrng[0]), $
                'Angle='+strtrim(angles[kk],2), color=clrs[kk+1], charsiz=lsz
     endfor

     oplot, replicate(2796.352,2), yrng, color=clr.gray, linesty=2
     oplot, replicate(2803.531,2), yrng, color=clr.gray, linesty=2
     xyouts, xrng[0]+xlbl*(xrng[1]-xrng[0]), yrng[1]*ylbl, $
             'MgII', color=clr.black, charsiz=lsz
     oplot, xrng, [1., 1.], color=clr.red, linestyle=1, thick=1
  endif

  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]
  close, /all

  return

end
