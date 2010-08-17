pro fig_plaw_spec, RREAL=rreal

  if not keyword_set( PSFILE ) then psfile = 'fig_plaw_spec.ps'
  if not keyword_set(PAD_FRAC) then pad_frac = 0.1
  if not keyword_set(CSZ) then csz = 1.8
  if not keyword_set(CSZ2) then csz2 = 2.3
  if not keyword_set(lSZ) then lsz = 1.5
  if not keyword_set(lSZ2) then lsz2 = 1.3
  if not keyword_set(XNCOLORS) then xncolors=200L
  if not keyword_set(YSTP) then ystp = 0.07
  lbl = ['A','D','G','B','E','H','C','F','I']
;  lbl = ['A','B','C','D','E','F','G','H','I','J','K']
  root = '../Analysis/PowLaws/Output/spec_plaw_'

  xlbl = 0.2
  xlbl2 = 0.05
  ylbl = 0.90

  close, /all

  ;;; BEGIN PLOTS
  x_psopen, psfile, /maxs
  !p.multi = [0,3,2,0,1]
  clr = getcolor(/load)
  clrs = [clr.black, clr.blue, clr.red, clr.darkgreen]; x_setclrs()

  xmrg = [8,3]
  ymrg = [4.5,0.5]

  for qq=0,2 do begin

     ;; Plot FeII
     yrng=[0., 1.8]
     xrng=[2577, 2618]
     xcut = 2605.
     off = 1.
     !p.multi = [6-qq*2,3,2,0,1]

     for kk=0,2 do begin
        cnt = qq*3 + kk
        Fe_fil = root+'FeII_'+lbl[cnt]+'.dat'
        readcol, Fe_fil, wv, fx, /silen
        nrm = median(fx[where(wv LT 2577)])

        for ss=0,1 do begin
           case ss of 
              0: begin
                 yrng=[0., 2.1] 
                 ysty = 9
                 wvmnx = [xrng[0], xcut-off]
              end
              1: begin
                 yrng=[0.95,1.4]
                 ysty = 5
                 wvmnx = [xcut+off,xrng[1]]
              end
              else: stop
           endcase
           
           plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
                 xmargin=xmrg, ymargin=ymrg, ytitle='Relative Flux', $
                 xtitle='Wavelength (Ang)', yrange=yrng, thick=4, $
                 xrange=xrng, ystyle=ysty, xstyle=1, psym=1, /nodata, /noerase
           if ss EQ 1 then axis, yaxis=1, charsiz=csz, ysty=1, xrang=yrng, ytickint=0.1
           
           pix = where(wv GT wvmnx[0] and wv LT wvmnx[1])
           ;; Plot
           oplot, wv[pix], fx[pix]/nrm, color=clrs[qq+1], psym=10, thick=3, linest=kk+1
           
           ;; Label
           if ss EQ 0 then  xyouts, xrng[0] + xlbl2*(xrng[1]-xrng[0]), $
                                    yrng[1] - (0.1 + kk*ystp)*(yrng[1]-yrng[0]), $
                                    lbl[cnt], color=clrs[qq+1], charsi=lsz
        endfor
     endfor
     oplot, replicate(2586.650,2), yrng, color=clr.gray, linesty=2, thick=1
     oplot, replicate(2600.173,2), yrng, color=clr.gray, linesty=2, thick=1
     oplot, replicate(2612.6542,2), yrng, color=clr.gray, linesty=2, thick=1
     xyouts, xrng[0]+xlbl*(xrng[1]-xrng[0]), yrng[0]+(yrng[1]-yrng[0])*ylbl, $
             'FeII', color=clr.black, charsiz=lsz
     x_curvefill, [xcut-off,xcut+off], [0., 0.], [10., 10], color=clr.tan
     
     ;;;;;;;;;;;;;;;;;;;;;;;;
     ;; MgII

     !p.multi = [5-qq*2,3,2,0,1]
     ;; Plot MgII
     yrng=[0., 4.0]
;     if qq GE 2 then yrng = [0., 2.5]
     xrng=[2785., 2809.8]
     plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
           xmargin=xmrg, ymargin=ymrg, ytitle='Relative Flux', $
           xtitle='Wavelength (Ang)', yrange=yrng, thick=4, $
           xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata

     for kk=0L,2 do begin
        cnt = qq*3 + kk
        Mg_fil = root+'MgII_'+lbl[cnt]+'.dat'
        readcol, Mg_fil, wv, fx, /silen
        nrm = median(fx[where(wv GT 2815)])

        ;; Plot
        oplot, wv, fx/nrm, color=clrs[qq+1], psym=10, thick=3, linest=kk+1

        ;; Label
        xyouts, xrng[0] + xlbl2*(xrng[1]-xrng[0]), $
                yrng[1] - (0.1 + kk*ystp)*(yrng[1]-yrng[0]), $
                lbl[cnt], $
                color=clrs[qq+1], charsi=lsz

     endfor
     
     oplot, replicate(2796.352,2), yrng, color=clr.gray, linesty=2, thick=1
     oplot, replicate(2803.531,2), yrng, color=clr.gray, linesty=2, thick=1
     xyouts, xrng[0]+xlbl*(xrng[1]-xrng[0]), yrng[1]*ylbl, $
             'MgII', color=clr.black, charsiz=lsz
     oplot, xrng, [1., 1.], color=clr.red, linestyle=1, thick=1
  endfor
     
  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]
  close, /all

  return

end
