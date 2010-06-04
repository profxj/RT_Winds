pro fig_1dspec, RREAL=rreal

  if not keyword_set( PSFILE ) then psfile = 'fig_1dspec.ps'
  if not keyword_set( GRIDFIL ) then $
    gridfil = '../Analysis/Outputs/fiducial_mgII.fits'
  if not keyword_set(PAD_FRAC) then pad_frac = 0.1
  if not keyword_set(CSZ) then csz = 1.1
  if not keyword_set(CSZ2) then csz2 = 2.0
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
  ps_open, filename=PSFILE, /color, bpp=8, /maxs;, XSIZE=xsize, YSIZE=ysize
  !p.thick = 6
  !x.thick = 6
  !y.thick = 6
  !p.charthick = 3

  device, /times, isolatin=1

  img_wave = [2793.7, 2795.2, 2796.35, 2797.5] ;; Ang
  dv = (img_wave-2796.35)/2796.35 * 3e5  ;; km/s
  ncut = n_elements(dv)

  ;; MgII Spectrum 
  clr = getcolor(/load)
  xmrg = [8,7]
  ymrg = [4.0,1]
  yrng=[0., 2.5]
  xrng=[2789., 2810]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, ytitle='Relative Flux', $
        xtitle='Wavelength (Ang)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, $
        pos=[0.08, 0.2, 0.48, 0.5]

  spec = total(total(raw_data,1),1)
  vel = (raw_wave-2796.35)/2796.35 * 3e5  ;; km/s
  oplot, raw_wave, spec, color=clr.black, psym=10, thick=3

  oplot, replicate(2796.352,2), yrng, color=clr.gray, linesty=2
  oplot, replicate(2803.531,2), yrng, color=clr.gray, linesty=2
  xyouts, 2792., 2.2, 'MgII', color=clr.black, charsiz=lsz

  ;; Velocity plot
  trans = [2796.352, 2803.531, 2796.352, 2803.531, 2796.352]
  ntrans = n_elements(trans)
  !p.multi = [ntrans+1,2,ntrans+1,0,1]

  xmrg = [8,3]
  ymrg = [0,0]
  yrng=[0., 2.4]
  xrng=[-600, 500]
  xtit = ''

  for qq=0L,ntrans-1 do begin
     if qq EQ (ntrans-1) then delvarx, xspaces $ 
     else xspaces = replicate(' ',30) 
     if qq EQ (ntrans-1) then xtit = 'Relative Velocity (km s!u-1!N)'


     plot, [0], [0], color=clr.black, background=clr.white, charsize=csz2,$
           xmargin=xmrg, ymargin=ymrg, ytitle='Relative Flux', $
           xtitle=xtit, yrange=yrng, thick=4, $
           xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, $
           xtickn=xspaces

     vel = (raw_wave-trans[qq])/trans[qq] * 3e5  ;; km/s
     oplot, vel, spec, color=clr.black, psym=10, thick=3
     ;; Lines and labels
     oplot, [0., 0.], yrng, color=clr.blue, linestyle=2, thick=1
;     oplot, xrng, [0., 0.], color=clr.green, linestyle=1, thick=1
     oplot, xrng, [1., 1.], color=clr.red, linestyle=1, thick=1

     getfnam, trans[qq], fv, nam 
     nam = strtrim(nam,2)
     xyouts, 0.65*(xrng[1]-xrng[0])+xrng[0], $
             yrng[0]+ (yrng[1]-yrng[0])*0.10, $
             strtrim(nam,2), color=clr.black, charsize=LSZ
  endfor
  
  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]

  return

end
