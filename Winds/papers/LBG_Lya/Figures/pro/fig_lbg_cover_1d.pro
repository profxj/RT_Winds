pro fig_lbg_cover_1d, wv, fx, YRNG=yrng

  
  if not keyword_set( PSFILE ) then psfile = 'fig_lbg_cover_1d.ps'
  if not keyword_set( CSZ ) then csz = 1.7
  if not keyword_set( LSZ ) then lsz = 1.6

  wrest = 1215.6701d

  ;; Plot
  x_psopen, psfile, /maxs
  clr = getcolor(/load)
  !p.multi=[0,1,1]

  xmrg = [8,2]
  ymrg = [4.0,0.5]
;  xrng=[1205., 1225] ; Ang
  xrng=[-5000., 2200]
;  yrng=[0., 2.7]
  yrng=[0.1, 100.]

  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, xtitle='Relative Velocity (km/s)', $;Rest Wavelength (Ang)', $
        ytitle='Relative Flux', $ 
        yrange=yrng, thick=5, xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, /ylog

  ;; Labels
  lins = [0.]
  for nn=0L,n_elements(lins)-1 do begin
      oplot, replicate(lins[nn],2), yrng, color=clr.gray, lines=2
;      loff = lins[nn] - 200/3e5*lins[nn]
;      oplot, replicate(loff,2), yrng, color=clr.blue, lines=2
  endfor

  ;; Continuum model
  file_conti = '../Analysis/LBG/Covering/1D/Output/spec_Lya_lbg_covering_continuum.dat'
  readcol, file_conti, wv, fx, fx_noscat, format='D,F,F'
  conti = median(fx)
  vel = (wv-wrest) * 3e5 / wrest
  oplot, vel, fx/conti, psym=10, color=clr.black, thick=3

  ;; Continuum model
  file_EW100 = '../Analysis/LBG/Covering/1D/Output/spec_Lya_lbg_covering_emiss_100A_50kms.dat'
  readcol, file_EW100, wv, fx, fx_noscat, format='D,F,F'
  conti = median(fx)
  vel = (wv-wrest) * 3e5 / wrest
  oplot, vel, fx/conti, psym=10, color=clr.red, thick=3


  ;; Data
  readcol, '../Data/lbg_stack_2003.dat', wv, fx, format='D,F'
  conti = 2.04d-30
  vel = (wv-wrest) * 3e5 / wrest
  oplot, vel, fx/conti, psym=10, color=clr.blue, thick=3, linest=1

  ;; EW
;  print, 'EW = ', total(1-fx)*(wv[1]-wv[0]), ' Ang'



  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  print, 'chk_spec:  All done!'
;  stop
       
  return
end
      
      
