pro fig_doppler_spec, RREAL=rreal

  if not keyword_set( PSFILE ) then psfile = 'fig_doppler_spec.ps'
  if not keyword_set(PAD_FRAC) then pad_frac = 0.1
  if not keyword_set(CSZ) then csz = 1.3
  if not keyword_set(CSZ2) then csz2 = 2.3
  if not keyword_set(lSZ) then lsz = 1.5
  if not keyword_set(lSZ2) then lsz2 = 1.3
  if not keyword_set(XNCOLORS) then xncolors=200L
  xlbl = 0.5
  ylbl = 0.85

  close, /all
  openr, 1,'Input/fig_doppler_spec.inp'

  ;;; BEGIN PLOTS
  x_psopen, psfile, /portrait
  !p.multi = [0,1,2]
  clr = getcolor(/load)
  clrs = x_setclrs()

  xmrg = [9,1]
  ymrg = [4.5,1]

  nFe = 0
  readf, 1, nFe

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

  nMg = 0
  readf, 1, nMg

  ;; Plot MgII
  yrng=[0., 2.5]
  xrng=[2786., 2810]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, ytitle='Relative Flux', $
        xtitle='Wavelength (Ang)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata

  Mg_fil = ''
  for kk=0L,nMg-1 do begin
     readf, 1, Mg_fil
     ;; Parse
     i1pos = strpos(Mg_fil, 'b')
     i2pos = strpos(Mg_fil, '.dat')
     if i1pos LT 0 then bval = 15. else $
        bval = float(strmid(Mg_fil, i1pos+1, i2pos-i1pos))
     readcol, Mg_fil, wv, fx, /silen
     nrm = median(fx[where(wv GT 2815)])
     ;; Plot
     oplot, wv, fx/nrm, color=clrs[kk], psym=10, thick=3

     ;; Label
     xyouts, xrng[0] + 0.1*(xrng[1]-xrng[0]), $
             yrng[1] - (0.1 + kk*0.08)*(yrng[1]-yrng[0]), $
             'b = '+string(round(bval),format='(i3)')+' km s!u-1!N', $
             color=clrs[kk], charsi=lsz

  endfor

  oplot, replicate(2796.352,2), yrng, color=clr.gray, linesty=2
  oplot, replicate(2803.531,2), yrng, color=clr.gray, linesty=2
  xyouts, xrng[0]+xlbl*(xrng[1]-xrng[0]), yrng[1]*ylbl, $
          'MgII', color=clr.black, charsiz=lsz
  oplot, xrng, [1., 1.], color=clr.red, linestyle=1, thick=1

  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]
  close, /all

  return

end
