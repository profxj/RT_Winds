;; Plot g(z) versus z
; fig_gz
pro fig_1dspec, GRIDFIL, TKRS=tkrs

  ;; Get structure if necessary
  if not keyword_set( PSFILE )  then psfile = 'fig_1dspec.ps'
  if not keyword_set( GRIDFIL ) then gridfil = 'spec_cube.fits'
  if not keyword_set( LSZ ) then lsz = 1.4

  ;; Initialize
  compile_opt strictarr

  ;; Input
  data = xmrdfits(gridfil, 0, /silent)
  wave = xmrdfits(gridfil, 1, /silent)
  sz = size(data, /dimen)

  slit = [1.0, 0.1, 0.25, 0.5]
  nslit = n_elements(slit)

  if keyword_set( PSFILE ) then x_psopen, psfile, /portrait
  !p.multi = [0,1,nslit]
  clr = getcolor(/load)
  xclr = x_setclrs()


  ;; zem
  ymnx = [0.,2.5]
  xmrg = [8,1]
  xrng = [min(wave,max=mxw), mxw]
  ymrg = [4., 0.5]

  for qq=0L,nslit-1 do begin

      pix = sz[0]/2 + slit[qq]*sz[0]/2*[-1,1]
      pix = round(pix)
      pix[0] = pix[0] > 0L
      pix[1] = pix[1] < (sz[0]-1)
      dat_slit=data[pix[0]:pix[1],*,*]

      ;; Create spectrum
      spec = total(total(dat_slit,1),1)

      ;; Plot
      plot, [0.], [0.], color=clr.black, background=clr.white, charsize=1.8,$
            xmargin=xmrg, ymargin=ymrg, xtitle='Wavelength (Ang)', ytitle='Flux',  $
            /nodata, xrange=xrng, ystyle=1, yrange=ymnx, xstyle=1

      oplot, wave, spec, color=xclr[qq], psym=10

      if qq EQ 0 then sv_spec = spec $
      else oplot, wave, sv_spec, color=clr.black,psym=10, linesty=1

      ;; Label
      oplot, replicate(2796.352,2), ymnx, color=clr.gray, linest=2
      oplot, replicate(2803.531,2), ymnx, color=clr.gray, linest=2
      xyouts, xrng[0]+1., ymnx[1]-0.3, 'Slit='+string(slit[qq],format='(f4.2)'),$
              color=clr.black, charsiz=lsz

      ;; EW
      EW = median(spec[where(wave GT 2808)])*sz[2] - total(spec)
      EW = EW * abs(wave[1]-wave[0])
      print, 'Slit: ', strtrim(qq,2), ', Total EW = ', EW, ' Ang'

      ;; Low-res
      if qq EQ 0 then begin
          fwhm_pix = 2.3 / abs(wave[1]-wave[0])
          nsmooth = fwhm_pix/(2.*sqrt(2*alog(2)))
          kernel = gauss_kernel(nsmooth)
          smth = convol(sv_spec, kernel,/edge_wrap)
          nlow = n_elements(sv_spec)/10
          low_spec = congrid(smth, nlow)
          low_wave = congrid(wave, nlow)
          oplot, low_wave, low_spec, color=clr.orange, psym=10
          ;; 
          if keyword_set(TKRS) then begin
              ztkrs = 0.694
              tkrs_fx = x_readspec('../Simple/TKRS4389_b_mask_XF_031609.fits', $
                                   wav=tkrs_wv,inf=2)
              tkrs_wv = tkrs_wv/(1+ztkrs)
              conti_tkrs = 0.694  ;; Odd!
              scale = median(spec[where(wave GT 2808)])
              oplot, tkrs_wv, tkrs_fx*scale/conti_tkrs, color=clr.cyan, psym=10
          endif
      endif

  endfor


  ;; Close Ps
  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]

  print, 'fig_1dspec: All done!'
  return
end
      
      
