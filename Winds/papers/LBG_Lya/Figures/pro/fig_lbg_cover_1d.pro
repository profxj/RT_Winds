pro fig_lbg_cover_1d, wv, fx, YRNG=yrng

  
  if not keyword_set( PSFILE ) then psfile = 'fig_lbg_cover_1d.ps'
  if not keyword_set( CSZ ) then csz = 1.7
  if not keyword_set( LSZ ) then lsz = 1.6

  wrest = 1215.6701d

  ;; IGM
  cd, '../Analysis/IGM/pro/', curr=curr
  resolve_routine, 'calc_tau_igm'
  cd, curr

  ;; Plot
  x_psopen, psfile, /portrait
  clr = getcolor(/load)
  !p.multi=[0,2,5,0,1]

  xmrg = [8,2]
  ymrg = [0.0,0.0]
;  xrng=[1205., 1225] ; Ang
  xrng=[-1500., 1500]
;  yrng=[0., 2.7]
  xlbl = 300.

  for qq=0,1 do begin

     ;; Input
     case qq of
        0: file= '../Analysis/LBG/Covering/1D/Output/spec_Lya_lbg_covering_continuum.dat'
        1: file= '../Analysis/LBG/Covering/1D/Output/spec_Lya_lbg_covering_emiss_100A_50kms.dat'
        else: stop
     endcase
     if qq EQ 1 then !p.multi=[5,2,5,0,1]

     readcol, file, wv, fx, fx_noscat, format='D,F,F'
     conti = median(fx)
     vel = (wv-wrest) * 3e5 / wrest
     npix = n_elements(fx)

     for ss=0,3 do begin
        xtit = ''
        ylog = 0
        xspcs = replicate(' ', 30)
        yrng=[0.0, 2.0]
        yconti = conti
        yval = fx
        ylbl = 0.2
        case ss of
           0: begin
              ylog=1
              yrng=[0.1, 300.]
;              if qq EQ 0 then yrng = [0.1, 10] else yrng=[0.1, 300.]
              if qq EQ 0 then yval = replicate(1., npix) else begin
                 ;; 100A line
                 sigma_L = 50. / 3e5 * 1215.6701 ;; Ang
                 A0 = 100. / sqrt(2*!pi) / sigma_L 
                 yval = 1. + A0 * exp(-1*(wv-wrest)^2/2/sigma_L^2)
                 yconti = 1.
              endelse
              lbl = 'Intrinsic'
              ylbl = 0.2
           end
           1: begin
              lbl = 'Wind'
           end
           2: begin
              calc_tau_IGM, vel, tau_IGM=tau_IGM 
              boost = 2.
              yval = fx * exp(-1.*tau_IGM*boost)
              lbl = 'Wind+IGM'
           end
           3: begin
              xtit = 'Relative Velocity (km/s)' 
              xspcs = ''
           end
           else: 
        endcase
        
        plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
              xmargin=xmrg, ymargin=ymrg, xtitle=xtit, ytitle='Relative Flux', $ 
              yrange=yrng, thick=5, xrange=xrng, ystyle=1, xstyle=1, psym=1, $
              /nodata, ylog=ylog, xtickn=xspcs

        oplot, replicate(0,2), yrng, color=clr.gray, lines=2

        ;; Continuum model
        oplot, vel, yval/yconti, psym=10, color=clr.black, thick=3

        ;; Labels
        xyouts, xlbl, ylbl, lbl, color=clr.black, charsi=lsz, align=0.
     endfor
  endfor

  ;; Data
  readcol, '../Data/lbg_stack_2003.dat', wv, fx, format='D,F'
  conti = 2.04d-30
  vel = (wv-wrest) * 3e5 / wrest
;  oplot, vel, fx/conti, psym=10, color=clr.blue, thick=3, linest=1

  ;; EW
;  print, 'EW = ', total(1-fx)*(wv[1]-wv[0]), ' Ang'



  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  print, 'chk_spec:  All done!'
;  stop
       
  return
end
      
      
