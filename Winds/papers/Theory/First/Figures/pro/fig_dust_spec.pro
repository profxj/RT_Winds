pro fig_dust_spec, RREAL=rreal

  if not keyword_set( PSFILE ) then psfile = 'fig_dust_spec.ps'
  if not keyword_set(PAD_FRAC) then pad_frac = 0.1
  if not keyword_set(CSZ) then csz = 1.3
  if not keyword_set(CSZ2) then csz2 = 2.3
  if not keyword_set(lSZ) then lsz = 1.5
  if not keyword_set(lSZ2) then lsz2 = 1.3
  if not keyword_set(XNCOLORS) then xncolors=200L
  if not keyword_set(YSTP) then ystp = 0.07

  xlbl = 0.3
  xlbl2 = 0.05
  ylbl = 0.90

  close, /all
  openr, 1,'Input/fig_dust_spec.inp'

  ;;; BEGIN PLOTS
  x_psopen, psfile, /portrait
  !p.multi = [0,1,2]
  clr = getcolor(/load)
  clrs = [clr.black, clr.blue, clr.red, clr.darkgreen] ;x_setclrs()

  xmrg = [8,4]
  ymrg = [3.0,3.5]

  nFe = 0
  readf, 1, nFe
  sv_fil = strarr(nFe)

  ;; Plot FeII
  yrng=[0., 1.8]
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
              yrng=[0.95,1.3]
              ysty = 5
              wvmnx = [xcut+off,xrng[1]]
           end
        else: stop
     endcase

     plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
           xmargin=xmrg, ymargin=ymrg, ytitle='Relative Flux', $
           xtitle='Wavelength (Ang)', yrange=yrng, thick=4, $
           xrange=xrng, ystyle=ysty, xstyle=9, psym=1, /nodata, /noerase
     if ss EQ 1 then axis, yaxis=1, charsiz=csz, ysty=1, xrang=yrng, ytickint=0.1
  
     Fe_fil = ''
     for kk=0L,nFe-1 do begin
        if ss EQ 0 then begin
           readf, 1, Fe_fil
           sv_fil[kk] = Fe_fil
        endif else Fe_fil = sv_fil[kk]
        ;; Parse
        i1pos = strpos(Fe_fil, 'tau')
        i2pos = strpos(Fe_fil, '.dat')
        if i1pos LT 0 then tau = 0. else $
           tau = float(strmid(Fe_fil, i1pos+3, i2pos-i1pos+2))
;        print, 'Reading', fe_fil
        readcol, Fe_fil, wv, fx, /silen
        nrm = median(fx[where(wv GT 2605)])
        pix = where(wv GT wvmnx[0] and wv LT wvmnx[1])
        ;; Plot
        oplot, wv[pix], fx[pix]/nrm, color=clrs[kk], psym=10, thick=3
        
        ;; Label
        if ss EQ 0 then  xyouts, xrng[0] + xlbl2*(xrng[1]-xrng[0]), $
                                 yrng[1] - (0.1 + kk*ystp)*(yrng[1]-yrng[0]), $
                                 '!9t!X!ddust!N = '+string(tau,format='(f4.1)'), $
                                 color=clrs[kk], charsi=lsz
        ;; Write
        if ss EQ 0 and (kk GT 0) then begin
           outfil = 'Output/dust_FeII_tau'+strtrim(string(tau,format='(f4.1)'),2)+'.dat'
           writecol, outfil, wv, fx/nrm
        endif
     endfor
  endfor
  xrng2 = (xrng/2600.173 - 1)*3e5
  axis, xaxis=1, charsiz=csz, xsty=1, xrang=xrng2, xtitl='Velocity (km/s) Relative to FeII 2600'

  oplot, replicate(2586.650,2), yrng, color=clr.orange, linesty=2, thick=2
  oplot, replicate(2600.173,2), yrng, color=clr.orange, linesty=2, thick=2
  oplot, replicate(2612.6542,2), yrng, color=clr.orange, linesty=2, thick=2
  xyouts, xrng[0]+xlbl*(xrng[1]-xrng[0]), yrng[0]+(yrng[1]-yrng[0])*ylbl, $
          'FeII', color=clr.black, charsiz=lsz
;     oplot, xrng, [1., 1.], color=clr.red, linestyle=1, thick=1
  x_curvefill, [xcut-off,xcut+off], [0., 0.], [10., 10], color=clr.tan
  
  ;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; MgII
  nMg = 0
  readf, 1, nMg

  !p.multi = [1,1,2]
  ;; Plot MgII
  yrng=[0., 2.7]
  xrng=[2786., 2809.8]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, ytitle='Relative Flux', $
        xtitle='Wavelength (Ang)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=9, psym=1, /nodata

  Mg_fil = ''
  for kk=0L,nMg-1 do begin
     readf, 1, Mg_fil
     ;; Parse
     i1pos = strpos(Mg_fil, 'tau')
     i2pos = strpos(Mg_fil, '.dat')
     if i1pos LT 0 then tau = 0. else $
        tau = float(strmid(Mg_fil, i1pos+3, i2pos-i1pos+2))
     readcol, Mg_fil, wv, fx, /silen
     nrm = median(fx[where(wv GT 2815)])
     ;; Plot
     oplot, wv, fx/nrm, color=clrs[kk], psym=10, thick=3

     ;; Label
     xyouts, xrng[0] + xlbl2*(xrng[1]-xrng[0]), $
             yrng[1] - (0.1 + kk*ystp)*(yrng[1]-yrng[0]), $
             '!9t!X!ddust!N = '+string(tau,format='(f4.1)'), $
             color=clrs[kk], charsi=lsz

     ;; Write
     if kk GT 0 then begin
        outfil = 'Output/dust_MgII_tau'+strtrim(string(tau,format='(f4.1)'),2)+'.dat'
        writecol, outfil, wv, fx/nrm
     endif

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
