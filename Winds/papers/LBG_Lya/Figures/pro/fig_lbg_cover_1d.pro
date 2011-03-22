pro fig_lbg_cover_1d, wv, fx, YRNG=yrng

  
  if not keyword_set( PSFILE ) then psfile = 'fig_lbg_cover_1d.ps'
  if not keyword_set( CSZ ) then csz = 1.7
  if not keyword_set( LSZ ) then lsz = 1.6


  ;; Plot
  x_psopen, psfile, /maxs
  clr = getcolor(/load)
  !p.multi=[0,1,1]

  xmrg = [8,2]
  ymrg = [4.0,0.5]
  xrng=[1210., 1225] ; Ang
  yrng=[0., 1.7]

  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, xtitle='Rest Wavelength (Ang)', $
        ytitle='Relative Flux', $ 
        yrange=yrng, thick=5, xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata

  ;; Continuum model
  file_conti = '../Analysis/LBG/Covering/1D/Output/spec_Lya_lbg_covering_continuum.dat'
  readcol, file_conti, wv, fx, fx_noscat, format='D,F,F'
  conti = median(fx)
  oplot, wv, fx/conti, psym=10, color=clr.black, thick=3

  ;; Continuum model
  file_EW100 = '../Analysis/LBG/Covering/1D/Output/spec_Lya_lbg_covering_emiss_100A_50kms.dat'
  readcol, file_EW100, wv, fx, fx_noscat, format='D,F,F'
  conti = median(fx)
  oplot, wv, fx/conti, psym=10, color=clr.red, thick=3


  ;; EW
;  print, 'EW = ', total(1-fx)*(wv[1]-wv[0]), ' Ang'

  ;; Labels
  lins = [1215.6701]
  for nn=0L,n_elements(lins)-1 do begin
      oplot, replicate(lins[nn],2), yrng, color=clr.gray, lines=1
;      loff = lins[nn] - 200/3e5*lins[nn]
;      oplot, replicate(loff,2), yrng, color=clr.blue, lines=2
  endfor


  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  print, 'chk_spec:  All done!'
;  stop
       
  return
end
      
      
