pro chk_spec

  
  if not keyword_set( PSFILE ) then psfile = 'chk_spec.ps'
  if not keyword_set( CSZ ) then csz = 1.7
  if not keyword_set( LSZ ) then lsz = 1.6


  ;; Read
  file = 'spec.dat'
  print, 'Reading ', file
  readcol, file, wv, fx, fx_noscat
  conti_no = median(fx_noscat[0:50])

  ;; Plot
  x_psopen, psfile, /maxs
  clr = getcolor(/load)
  !p.multi=[0,1,1]

  xmrg = [8,2]
  ymrg = [4.0,0.5]
  xrng=[min(wv,max=mwv),mwv]
  yrng=[0., max(fx)*1.05]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, xtitle='Rest Wavelength (Ang)', $
        ytitle='Relative Flux', $ 
        yrange=yrng, thick=5, xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata

  ;; New spec 
  oplot, wv, fx, psym=10, color=clr.black, thick=3
  oplot, wv, fx_noscat/conti_no, psym=10, color=clr.red, thick=3, linest=1

  ;; EW
  print, 'EW = ', total(1-fx)*(wv[1]-wv[0]), ' Ang'

  ;; Labels
  if max(wv) LT 2700 then lins = [2586.600, 2600.1729, 2612.6542, 2626.4511] $
  else lins = [2796.352, 2803.531]
  for nn=0L,n_elements(lins)-1 do begin
      oplot, replicate(lins[nn],2), yrng, color=clr.gray, lines=1
;      loff = lins[nn] - 200/3e5*lins[nn]
;      oplot, replicate(loff,2), yrng, color=clr.blue, lines=2
  endfor


  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  print, 'chk_spec:  All done!'
       
  return
end
      
      
