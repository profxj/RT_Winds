pro fig_2dspec, RREAL=rreal

  if not keyword_set( PSFILE ) then psfile = 'fig_2dspec.ps'
  if not keyword_set( GRIDFIL ) then gridfil = 'spec_cube.fits'
  if not keyword_set(PAD_FRAC) then pad_frac = 0.1
  if not keyword_set(CSZ) then csz = 1.3

  ;; Read in Data
  raw_data = xmrdfits(gridfil, 0, /silent)
  raw_wave = xmrdfits(gridfil, 1, /silent)
  dwv = abs(raw_wave[1]-raw_wave[0])
  raw_data = float(raw_data)
  sz_raw = size(raw_data,/dimen)

  ;; Pad
  if not keyword_set(NOPAD) then begin
      pad_raw = round(pad_frac*sz_raw[2]) > 10
      sv_data = fltarr(sz_raw[0], sz_raw[1], 2*pad_raw+sz_raw[2])
      sz = size(sv_data, /dimen)
      sv_data[*,*,pad_raw:pad_raw+sz_raw[2]-1] = raw_data
      for ii=0L,pad_raw-1 do $
             sv_data[*,*,ii] = raw_data[*,*,0]
      for ii=pad_raw+sz_raw[2],sz[2]-1 do $
             sv_data[*,*,ii] = raw_data[*,*,sz_raw[2]-1]
      wave = [raw_wave[0]-reverse(findgen(pad_raw)+1)*dwv, $
              raw_wave, $
              max(raw_wave)+(findgen(pad_raw)+1)*dwv]
  endif

  ;; Rebinning (spatial and spectral)
  rx = sz[2]/10
  ry = sz[1]/5

  ;; SLIT
  slit = 1.0
  pix = sz[0]/2 + slit*sz[0]/2*[-1,1]
  pix = round(pix)
  pix[0] = pix[0] > 0L
  pix[1] = pix[1] < (sz[0]-1)

  angstrom = '!6!sA!r!u!9 %!6!n'


  set_plot, 'x'
  !p.font = 0
  device, decompose=0
  ps_open, filename=PSFILE, /color, bpp=8, /MAXS;, XSIZE=xsize, YSIZE=ysize
  !p.thick = 6
  !x.thick = 6
  !y.thick = 6
  !p.charthick = 3

  device, /times, isolatin=1
  devicefactor=2540.  ;; cm to inch

  imsize1 = 1.6
  imsize2 = 2.3
  xpos1 = 0.7
  ypos1 = 0.5

  stretch_lo = -0.01
  stretch_hi = 0.02

  for qq=0,2 do begin
  
      ;; Image
      data = sv_data
      if qq EQ 1 then data[33:66,33:66,*] = 0.

      ;;
      if qq EQ 2 then begin
          if keyword_set(RREAL) then spec2d = xmrdfits(RREAL,/silen) else begin
              fwhm_pix = 5
              nx = 4*round(fwhm_pix)
              ny = nx
              sig_pix = fwhm_pix / (2*sqrt(2*alog(2)))
              A = [0., 1., sig_pix, sig_pix, nx/2., ny/2.]
              x = FINDGEN(nx) # REPLICATE(1.0, ny) 
              Y = REPLICATE(1.0, nx) # FINDGEN(ny) 
              ;; Create an ellipse: 
              U = ((X-A[4])/A[2])^2 + ((Y-A[5])/A[3])^2 
              ;; Kernel
              kernel2d = A[0] + A[1] * EXP(-U/2) 
              kernel2d = kernel2d / total(kernel2d[*])
              
              ;; Convolve
              smooth_data = fltarr(sz)
              for ii=0L,sz[2]-1 do $
                     smooth_data[*,*,ii] = convol(sv_data[*,*,ii], kernel2d, /edge_wrap)
              ;;xatv, smooth_data, /blo
              dat_slit=smooth_data[pix[0]:pix[1],*,*]
              spec2d = transpose(total(dat_slit,1))
              sz2 = size(spec2d,/dimen)
              ;; Smooth in spectral
              fwhm_pix = 3. / abs(wave[1]-wave[0])
              nsmooth = fwhm_pix/(2.*sqrt(2*alog(2)))
              kernel = gauss_kernel(nsmooth)
              for ii=0L,sz2[1]-1 do $
                     spec2d[*,ii] = convol(spec2d[*,ii],kernel, /edge_wrap)
              ;; Rebin the 'pixel' scale
              mwrfits, spec2d, 'real_tmp.fits', /create
          endelse
          spec2d = congrid(spec2d, rx, ry)
;          wave = congrid(wave, rx)
          stretch_hi=0.03
      endif else begin
          dat_slit=data[pix[0]:pix[1],*,*]
          spec2d = transpose(total(dat_slit,1))
      endelse

      ;; Log stretch
      disp_flux2d = -1*alog10((spec2d>stretch_lo<stretch_hi)+1.5)
      
      ;; Plot
      loadct,0
      yoffset = (imsize1+0.8)*qq
      if qq EQ 2 then imsize=imsize2 else imsize=imsize1
      tvscl, disp_flux2d, xpos1, ypos1+yoffset, ysize=imsize, /inches

      ;; Axes
      dims = size(disp_flux2d,/dim)
      xlabel = findgen(dims[0]+1)
      ylabel = findgen(dims[1]+1)
      thisPosition = devicefactor*[xpos1, $
                                   ypos1+yoffset, $
                                   xpos1+(imsize*dims[0]/dims[1]), $
                                   ypos1+yoffset+imsize]
      plot, xlabel, ylabel, /nodata, /device, /noerase, position=thisPosition, $
            xrange=[min(xlabel),max(xlabel)], yrange=[min(ylabel),max(ylabel)], $
            xstyle=5, ystyle=5
      ymnx = [-10,10]
      plot, wave, wave, /device, /noerase, xrange=[min(wave),max(wave)], $
            yrange=ymnx, xtitle='Wavelength (Ang)', charsiz=csz, $
            ytitle='kpc', /nodata, xsty=1, $
            position=thisPosition ;, ytickname=['-2','-1','0','+1','+2']
      clr = getcolor(/load)
      oplot, [2796.352, 2796.352], ymnx, color=clr.blue, linesty=1
      oplot, [2803.531, 2803.531], ymnx, color=clr.red, linesty=1
  endfor
  
  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]

end
