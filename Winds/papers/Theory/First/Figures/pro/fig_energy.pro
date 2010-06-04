;; Plot g(z) versus z
; fig_gz
pro fig_slitwidth, GRIDFIL

  ;; Get structure if necessary
;  if not keyword_set( PSFILE )  then psfile = 'fig_slitwidth.ps'
  if not keyword_set( GRIDFIL ) then $
    gridfil = '/u/xavier/RadTransfer/Winds/Grid/spec_cube.fits'
  if not keyword_set( LSZ ) then lsz = 1.4

  ;; Initialize
  compile_opt strictarr

  ;; Input
  data = xmrdfits(gridfil, 0, /silent)
  wave = xmrdfits(gridfil, 1, /silent)
  sz = size(data, /dimen)

  slit = (1+findgen(20))/20.
  nslit = n_elements(slit)

  if keyword_set( PSFILE ) then x_psopen, psfile, /maxs
  !p.multi = [0,1,1]
  clr = getcolor(/load)
  xclr = x_setclrs()

  ;; zem
  ymnx = [-0.3,4.0]
  xmrg = [8,1]
  xrng = [0, 1. ]
  ymrg = [4., 0.5]

  ;; Plot
  plot, [0.], [0.], color=clr.black, background=clr.white, charsize=1.8,$
        xmargin=xmrg, ymargin=ymrg, xtitle='Slit Width (Fraction)', $
        ytitle='EW (Ang)',  $
        /nodata, xrange=xrng, ystyle=1, yrange=ymnx, xstyle=1



  sv_EW = fltarr(nslit)
  for qq=0L,nslit-1 do begin

      pix = sz[0]/2 + slit[qq]*sz[0]/2*[-1,1]
      pix = round(pix)
      pix[0] = pix[0] > 0L
      pix[1] = pix[1] < (sz[0]-1)
      dat_slit=data[pix[0]:pix[1],*,*]

      ;; Create spectrum
      spec = total(total(dat_slit,1),1)
      x_splot, wave, spec, /bloc

      EW = median(spec[where(wave GT 2808)])*sz[2] - total(spec)
      EW = EW * abs(wave[1]-wave[0])
      sv_EW[qq] = EW  ;; Ang

  endfor

  oplot, xrng, [0., 0.], color=clr.gray, lines=1
  ;; Plot
  oplot, slit, sv_EW, color=clr.black, psym=10


  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]

  print, 'fig_slitwidth: All done!'
  return
end
      
      
