pro comp_tkrs, FIDDLE=fiddle

  
  if not keyword_set( PSFILE ) then psfile = 'comp_tkrs.ps'
  if not keyword_set( CSZ ) then csz = 1.7
  if not keyword_set( LSZ ) then lsz = 1.6


  ;; Read
  file = 'spec.dat'
  print, 'Reading ', file
  readcol, file, wv, fx
  npix = n_elements(wv)
  tkrs_fx = x_readspec('TKRS4389_b_mask_XF_031609.fits',wav=tkrs_wv,inf=2)
  zgal = 0.6942

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

  ;; TKRS
  rwave = tkrs_wv/(1+zgal)
  wvc = where(rwave GT 2780 and rwave LT 2785)
  nrm = median(tkrs_fx[wvc])
  oplot, tkrs_wv/(1+zgal), tkrs_fx/nrm, psym=10, color=clr.black, thick=3

  ;; Wind
  fwhm_pix = 2.3 / abs(wv[1]-wv[0])  ;; 250 km/s FWHM
  nsmooth = fwhm_pix/(2.*sqrt(2*alog(2)))
  kernel = gauss_kernel(nsmooth)
  smth = convol(fx, kernel,/edge_wrap)
  nlow = npix/10
  low_fx = congrid(smth, nlow)
  low_wv = congrid(wv, nlow)
  oplot, low_wv, low_fx, psym=10, color=clr.red, thick=3

  ;; Labels
  lins = [2796.352, 2803.531]
  for nn=0L,n_elements(lins)-1 do begin
      oplot, replicate(lins[nn],2), yrng, color=clr.gray, lines=1
      loff = lins[nn] - 200/3e5*lins[nn]
      oplot, replicate(loff,2), yrng, color=clr.blue, lines=2
  endfor


  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  if keyword_set(FIDDLE) then $
    x_specplot, fx, fltarr(npix), wav=wv, infl=4, /bloc

  print, 'chk_spec:  All done!'
       
  return
end
      
      
