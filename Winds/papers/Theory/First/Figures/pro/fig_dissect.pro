pro fig_dissect, RREAL=rreal

  if not keyword_set( PSFILE ) then psfile = 'fig_dissect.ps'
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
  ps_open, filename=PSFILE, /color, bpp=8, /maxs;, XSIZE=xsize, YSIZE=ysize
  !p.thick = 6
  !x.thick = 6
  !y.thick = 6
  !p.charthick = 3

  device, /times, isolatin=1

  img_wave = [2793.7, 2795.2, 2796.35, 2797.5] ;; Ang
  dv = (img_wave-2796.35)/2796.35 * 3e5  ;; km/s
  ncut = n_elements(dv)

  ;; IFU shots next

  devicefactor=2540.  ;; cm to inch
  imsize = 3.5  ;; inch
  x0 = 1.2
  y0 = 0.7
;  stretch_lo = 0.01

  ;; Normalize
  tmp = raw_data
  tmp[33:66,33:66,*] = 0.
  mxd = max(tmp)
  raw_data = 100. * raw_data / mxd


  ;; Wind model
  ngrid = 100L
  box_size = 20.           ; Total box size (kpc)
  dl = box_size / ngrid    ; Cell size (kpc) 
  v_wind  = 200.           ; Wind speed (km/s)
  dv_wind = 200.           ; Velocity width of wind (km/s)
  rw_inner = 7.            ; Inner boundary of the wind (kpc)
  rw_outer = 10.           ; Outer boundary of the wind (kpc)
  rg_outer =  3.           ; Outer boundary of galaxy (kpc)
  doppler = 20.            ; Doppler parameter (km/s)
  dens_inner = 1e-2        ; Density of the wind at rw_inner
  rexpon = -1.             ; Exponent of the density law

  rval = dl* (findgen(ngrid) - ngrid/2 + 0.5)  ;; Dimensionless

  ;; Velocity
  speed = v_wind - dv_wind/2. + $
               dv_wind*(rval-rw_inner)/(rw_outer-rw_inner) 
  zero = where(rval LT rw_inner or rval GT rw_outer)
  speed[zero] = 0.

  nalph = 360L
  alpha = findgen(nalph)*!pi/180. ; radians
  v1 = -100. * cos(alpha)
  v2 = -300. * cos(alpha)
  cfac = 1. / cos(alpha)
  zro = where(abs(alpha*180./!pi -90.) LT 1.,nz) 
  if nz GT 0 then cfac[zro] = 1e10

  xpos1 = 1.
  ypos1 = 2.
  for qq=0,ncut-1 do begin

      ctload, 0, ncolors=xncolor
      plot, [0], [0], /nodata, color=255

      ;; Image
      mn = min(abs(raw_wave-img_wave[qq]),img_pix)
      data = raw_data[*,*,img_pix]
      ;; Remove inner region
      if dv[qq] GT -10 then data[33:66,33:66] = 0.
      irange = [0., max(data)]

      ;; Color bars
      ctload, 0, ncolors=xncolors ;, /rever
      coyote_colorbar, pos=[0.48, 0.20, 0.52, 0.76], /verti, range=irange, $
                       ncolor=xncolors, /invert

      ;; Scale
      scaled = bytscl(data, min=irange[0], max=irange[1], $
                      top=(xncolors - 1))
      
      ;; Plot
      ctload, 0, ncolors=xncolor, /rever
      tv, scaled, xpos1, ypos1, ysize=imsize, /inches

      ;; Axes
      dims = size(scaled,/dim)
      xlabel = findgen(dims[0]+1)
      ylabel = findgen(dims[1]+1)
      thisPosition = devicefactor*[xpos1, $
                                   ypos1, $
                                   xpos1+(imsize*dims[0]/dims[1]), $
                                   ypos1+imsize]
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

      ;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Diagnoistic plot
      plot, [0], [0], /noerase, /device, xrange=ymnx, $
            yrange=ymnx, xtitle='kpc', charsiz=csz2, $
            ytitle='kpc', /nodata, xsty=1, $
            position=devicefactor*[6.8, 2., 10.3, 5.5] 

      ;; Wind boundaries
      x_oplotcirc,  7., color=clr.gray
      x_oplotcirc, 10., color=clr.gray

      ;; alpha range
      for kk=0L,nalph-1 do begin
          if (dv[qq]-v1[kk])*(dv[qq]-v2[kk]) LT 0 then begin
              ;; Find the radius
              vr = dv[qq] * cfac[kk]
              mn = min(abs(vr-(-1*speed)), imn)
              x = -1 * rval[imn] * cos(alpha[kk])
              y = rval[imn] * sin(alpha[kk])
              oplot, [x], [y], psym=1, color=clr.red
          endif
      endfor
      ;; 
      
  endfor
  
  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]

  return

end
