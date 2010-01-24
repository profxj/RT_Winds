pro fig_snaps, RREAL=rreal

  if not keyword_set( PSFILE ) then psfile = 'fig_snaps.ps'
  if not keyword_set( GRIDFIL ) then gridfil = 'spec_cube.fits'
  if not keyword_set(PAD_FRAC) then pad_frac = 0.1
  if not keyword_set(CSZ) then csz = 1.3
  if not keyword_set(lSZ) then lsz = 1.5

  ;; Read in Data
  raw_data = xmrdfits(gridfil, 0, /silent)
  raw_wave = xmrdfits(gridfil, 1, /silent)
  dwv = abs(raw_wave[1]-raw_wave[0])
  raw_data = float(raw_data)
  sz_raw = size(raw_data,/dimen)

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

  imsize = 3.0
  x0 = 1.2
  y0 = 0.7
  stretch_lo = -0.01

  img_wave = [2793.7, 2795.2, 2796.35, 2797.5] ;; Ang
  dv = (img_wave-2796.35)/2796.35 * 3e5  ;; km/s

  xpos1 = x0
  for qq=0,3 do begin

      if qq GT 1 then stretch_hi = 0.02 else stretch_hi = 0.02
  
      ;; Position
      if qq EQ 2 then xpos1 = xpos1+(x0+imsize)
      if (qq MOD 2) EQ 1 then ypos1 = 2*y0+imsize else ypos1=y0

      ;; Image
      mn = min(abs(raw_wave-img_wave[qq]),img_pix)
      data = raw_data[*,*,img_pix]

      ;; Remove inner region
      if qq GT 1 then data[33:66,33:66] = 0.

      ;; Log stretch
      disp_flux = -1*alog10((data>stretch_lo<stretch_hi)+1.5)
      
      ;; Plot
      loadct,0
      tvscl, disp_flux, xpos1, ypos1, ysize=imsize, /inches

      ;; Axes
      dims = size(disp_flux,/dim)
      xlabel = findgen(dims[0]+1)
      ylabel = findgen(dims[1]+1)
      thisPosition = devicefactor*[xpos1, $
                                   ypos1, $
                                   xpos1+(imsize*dims[0]/dims[1]), $
                                   ypos1+imsize]
      plot, xlabel, ylabel, /nodata, /device, /noerase, position=thisPosition, $
            xrange=[min(xlabel),max(xlabel)], yrange=[min(ylabel),max(ylabel)], $
            xstyle=5, ystyle=5
      ymnx = [-10,10]
      plot, [0], [0], /device, /noerase, xrange=ymnx, $
            yrange=ymnx, xtitle='kpc', charsiz=csz, $
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

end
