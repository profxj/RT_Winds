pro fig_test_static_sphere, wv, fx, YRNG=yrng, XRNG=xrng

  
  if not keyword_set( PSFILE ) then psfile = 'fig_test_static_sphere.ps'
  if not keyword_set( CSZ ) then csz = 1.7
  if not keyword_set( LSZ ) then lsz = 1.6

  c = x_constants()

  ;; Line parameters
  v_doppler    =   40631. ;  // cm/s :: T=10K

  ;; Read
  file = 'Output/spec_test_static_sphere.dat'
  print, 'Reading ', file
  readcol, file, wv, fx, fx_noscat, format='D,F,F'
  conti_no = median(fx_noscat)

  vel = (wv-1215.6701) / 1215.6701 * c.c  ;; cm/s
  xval = -1 * vel/v_doppler

  ;; Plot
  x_psopen, psfile, /maxs
  clr = getcolor(/load)
  !p.multi=[0,1,1]

  xmrg = [8,2]
  ymrg = [4.0,0.5]
  xrng=[-20, 20.]

  if not keyword_set( YRNG ) then yrng=[0., max(fx)*1.05]

  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, xtitle='x', $
        ytitle='Relative Flux', $ 
        yrange=yrng, thick=5, xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata

  ;; New spec 
  oplot, wv, fx, psym=10, color=clr.black, thick=3
  oplot, wv, fx_noscat/conti_no, psym=10, color=clr.red, thick=3, linest=1

  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  print, 'chk_spec:  All done!'
;  stop
       
  return
end
      
      
