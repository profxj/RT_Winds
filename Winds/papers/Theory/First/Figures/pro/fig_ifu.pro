pro fig_ifu, RREAL=rreal

  if not keyword_set( PSFILE ) then psfile = 'fig_ifu.ps'
  if not keyword_set( GRIDFIL ) then $
    gridfil = '/u/xavier/RadTransfer/Winds/Grid/spec_cube.fits'
  if not keyword_set(PAD_FRAC) then pad_frac = 0.1
  if not keyword_set(CSZ) then csz = 1.1
  if not keyword_set(CSZ2) then csz2 = 1.1
  if not keyword_set(lSZ) then lsz = 1.5
  if not keyword_set(XNCOLORS) then xncolors=200L

  ;; Read in Data
  raw_data = xmrdfits(gridfil, 0, /silent)
  raw_wave = xmrdfits(gridfil, 1, /silent)
  dwv = abs(raw_wave[1]-raw_wave[0])
  raw_data = float(raw_data)
  sz_raw = size(raw_data,/dimen)

  set_plot, 'x'
  !p.font = 0
  device, decompose=0
  ps_open, filename=PSFILE, /color, bpp=8, /portrait;, XSIZE=xsize, YSIZE=ysize
  !p.thick = 6
  !x.thick = 6
  !y.thick = 6
  !p.charthick = 3

  device, /times, isolatin=1

  img_wave = [2793.7, 2795.2, 2796.35, 2797.5] ;; Ang
  dv = (img_wave-2796.35)/2796.35 * 3e5  ;; km/s
  ncut = n_elements(dv)

  ;; Spectrum first

  clr = getcolor(/load)
  xmrg = [8,7]
  ymrg = [4.0,1]
  yrng=[0., 2.5]
  xrng=[-600., 1200]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, ytitle='Flux', $
        xtitle='v (km s!u-1!N)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, $
        pos=[0.15, 0.43, 0.85, 0.57]


  spec = total(total(raw_data,1),1)
  vel = (raw_wave-2796.35)/2796.35 * 3e5  ;; km/s
  oplot, vel, spec, color=clr.black, psym=10, thick=3

  for jj=0L,ncut-1 do $
      oplot, replicate(dv[jj],2), yrng, color=clr.gray, lines=1

  ;; IFU shots next

  devicefactor=2540.  ;; cm to inch
  imsize = 2.6  ;; inch
  x0 = 1.2
  y0 = 0.7
;  stretch_lo = 0.01

  ;; Normalize
  irange1 = [0., 15.]
  irange2 = [0., 75.]
  tmp = raw_data
  tmp[33:66,33:66,*] = 0.
  mxd = max(tmp)
  raw_data = irange2[1] * raw_data / mxd

  ;; Color bars
  ctload, 0, ncolors=xncolors;, /rever
  coyote_colorbar, pos=[0.48, 0.66, 0.52, 0.96], /verti, range=irange1, $
                   ncolor=xncolors, /invert

  coyote_colorbar, pos=[0.48, 0.03, 0.52, 0.33], /verti, range=irange2, $
                   ncolor=xncolors, /invert

  xpos1 = [0.5, 4.7, 0.5, 4.7]
  ypos1 = [6., 6., 0.3, 0.3]
  for qq=0,3 do begin

      if qq LE 1 then irange = irange1 else irange = irange2

;      stretch_hi = 0.02 
  
      ;; Image
      mn = min(abs(raw_wave-img_wave[qq]),img_pix)
      data = raw_data[*,*,img_pix]
      ;; Remove inner region
      if qq GT 1 then data[33:66,33:66] = 0.
      print, 'max = ', max(data)

      ;; Log stretch
;      disp_flux = -1*alog10((data>stretch_lo<stretch_hi)+1.5)
      scaled = bytscl(data, min=irange[0], max=irange[1], $
                      top=(xncolors - 1))
      
      ;; Plot
      ctload, 0, ncolors=xncolor, /rever
      tv, scaled, xpos1[qq], ypos1[qq], ysize=imsize, /inches

      ;; Axes
      dims = size(scaled,/dim)
      xlabel = findgen(dims[0]+1)
      ylabel = findgen(dims[1]+1)
      thisPosition = devicefactor*[xpos1[qq], $
                                   ypos1[qq], $
                                   xpos1[qq]+(imsize*dims[0]/dims[1]), $
                                   ypos1[qq]+imsize]
      ctload, 0, ncolors=xncolor
      plot, xlabel, ylabel, /nodata, /device, /noerase, $
            position=thisPosition, $
            xrange=[min(xlabel),max(xlabel)], $
            yrange=[min(ylabel),max(ylabel)], $
            xstyle=5, ystyle=5
      ymnx = [-10,10]
      plot, [0], [0], /device, /noerase, xrange=ymnx, $
            yrange=ymnx, xtitle='kpc', charsiz=csz2, $
            ytitle='kpc', /nodata, xsty=1, $
            position=thisPosition ;, ytickname=['-2','-1','0','+1','+2']
      clr = getcolor(/load)
      if round(dv[qq]) LT 0. then pclr = clr.blue else pclr=clr.red
      if abs(round(dv[qq])) LT 10. then pclr = clr.black
      xyouts, -9., 8., 'v='+strtrim(round(dv[qq]),2),color=pclr, $
              charsiz=lsz
  endfor
  
  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]

  return

end
