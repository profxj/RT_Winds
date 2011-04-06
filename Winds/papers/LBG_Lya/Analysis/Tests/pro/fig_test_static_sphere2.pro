pro fig_test_static_sphere2, wv, fx, YRNG=yrng, XRNG=xrng
  ;;
  
  if not keyword_set( PSFILE ) then psfile = 'fig_test_static_sphere2.ps'
  if not keyword_set( CSZ ) then csz = 1.7
  if not keyword_set( LSZ ) then lsz = 1.6

  c = x_constants()

  ;; Line parameters
  v_doppler    =   18.24*1e5 ; // cm/s (T=20,000K)
  Dnu   = v_doppler * 2.466d15 / 2.9979e10 ;   
  voigt_a   = 6.265d8 / (4 * 3.14159 * Dnu) ;   // gamma/(4 pi Dnu)
  tau_0 = 3.31d-14 * (12.85/18.24) * 2d20

  ;; Read
  file = 'Output/spec_test_static_sphere2.dat' ;; Should be tau_0=10^5
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
  xrng=[-30, 30.]

  if not keyword_set( YRNG ) then yrng=[0., max(fx)*1.05]

  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, xtitle='x', $
        ytitle='Relative Flux', $ 
        yrange=yrng, thick=5, xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata

  ;; New spec 
  oplot, xval, fx, psym=10, color=clr.black, thick=3
  mx_fx = max(fx)

  ;; Neufeld (with modification for a sphere)
  J_neuf = sqrt(6)/(24*sqrt(!pi)) * (xval^2 / (voigt_a*tau_0)) / $
           cosh(sqrt(!pi^3/54) * (abs(xval^3)/(voigt_a*tau_0)))
  mxJ = max(J_neuf)
  J_neuf = J_neuf / mxJ * mx_fx

  oplot, xval, J_neuf, color=clr.red

  print, 'atau = ', voigt_a*tau_0

  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  print, 'chk_spec:  All done!'
;  stop
       
  return
end
      
      
