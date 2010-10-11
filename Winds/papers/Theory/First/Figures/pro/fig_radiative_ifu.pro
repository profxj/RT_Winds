pro fig_radiative_ifu, wrest, dv, grid_file, psfile, YMNX=ymnx, XRNG=xrng, YRNG=yrng

  if not keyword_set(PAD_FRAC) then pad_frac = 0.1
  if not keyword_set(CSZ) then csz = 1.1
  if not keyword_set(CSZ2) then csz2 = 1.0
  if not keyword_set(CSZ3) then csz3 = 0.9
  if not keyword_set(lSZ) then lsz = 1.3
  if not keyword_set(XNCOLORS) then xncolors=200L
  if not keyword_set(YMNX) then ymnx = [-50,0]  ; kpc (size of box)
  if not keyword_set(XRNG) then xrng=[-600., 1200]
  if not keyword_set(YRNG) then yrng=[0., 2.5]
  if not keyword_set(PSFILE) then psfile = 'fig_radiative_ifu.ps'

  ;; 'Inputs'
  all_wrest = [2796.35, 2796.35, 2600.173, 2612.654]
  dv = [-100., 100., 0, 0]
  ncut = n_elements(dv)
  grid_file = '../Analysis/Radiation/Output/radiative_grid.fits'
;  ymnx=[-20,20] 

  if ncut GT 4 then stop

  x_psopen, psfile, /portrait
  clr = getcolor(/load)
  xmrg = [8,7]
  ymrg = [4.0,1]

  ;; IFU shots 
  devicefactor=2540.  ;; cm to inch
  imsize = 2.6  ;; inch
  x0 = 1.2
  y0 = 0.7
;  stretch_lo = 0.01

  ;; Normalize
  irange1 = [0.01, 1000.]
  irange2 = [0.01, 1000.]
;  stop

  ;; Color bars
  ctload, 0, ncolors=xncolors;, /rever
  coyote_colorbar, pos=[0.50, 0.66, 0.54, 0.96], /verti, range=alog10(irange1), $
                   ncolor=xncolors, /invert, FORMAT='(f4.1)', $
                   title='Log Surface Brightness',charsiz=csz3 

  if ncut GT 2 then $
     coyote_colorbar, pos=[0.50, 0.285, 0.54, 0.585], /verti, range=alog10(irange2), $
                      ncolor=xncolors, /invert, FORMAT='(f4.1)',$
                   title='Log Surface Brightness',charsiz=csz3 

  xpos1 = [0.5, 4.45, 0.5, 4.45]
  ypos1 = [6., 6., 2.7, 2.7]
  for qq=0,ncut-1 do begin
     wrest = all_wrest[qq]
     img_wave = wrest + dv[qq]/3e5 * wrest
     
     ;; Read in data
     if wrest GT 2700. then idx = [2,3] else idx = [0,1]
     raw_data = xmrdfits(grid_file, idx[0], /silent)
     raw_wave = xmrdfits(grid_file, idx[1], /silent)
     dwv = abs(raw_wave[1]-raw_wave[0])
     raw_data = float(raw_data)
     sz_raw = size(raw_data,/dimen)

     off = 1
     z0 = sz_raw[0]/2 - off 
     z1 = sz_raw[0]/2 + off


     tmp = raw_data
     tmp[z0:z1,z0:z1,*] = 0.
     mxd = max(tmp)
     raw_data = irange2[1] * raw_data / mxd


      if qq LE 1 then irange = irange1 else irange = irange2

;      stretch_hi = 0.02 
  
      ;; Image
      mn = min(abs(raw_wave-img_wave),img_pix)
      print, 'qq: ', qq, img_pix
;      stop
      data = raw_data[*,*,img_pix]
      ;; Remove inner region
;      if qq GT 1 then data[z0:z1,z0:z1] = 0.
      print, 'max = ', max(data)

;      disp_flux = -1*alog10((data>stretch_lo<stretch_hi)+1.5)
      ;; Linear
;      scaled = bytscl(data, min=irange[0], max=irange[1], $
;                      top=(xncolors - 1))
      ;; Log stretch
      scaled = bytscl(alog10(data), min=alog10(irange[0]), max=alog10(irange[1]), $
                      /nan, $
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
      if (qq MOD 2) EQ 1 then begin
         ytit=''
         yspaces = replicate(' ', 30)
      endif  else begin
         ytit = 'kpc'
         yspaces = ''
      endelse
      plot, [0], [0], /device, /noerase, xrange=ymnx, $
            yrange=ymnx, xtitle='kpc', charsiz=csz2, ytickn=yspaces, $
            ytitle=ytit, /nodata, xsty=1, ysty=1, $
            position=thisPosition ;, ytickname=['-2','-1','0','+1','+2']
      if (qq MOD 2) EQ 1 then axis, yaxis=1, charsiz=csz2, ysty=1, xrang=yrng, ytit='kpc'
      clr = getcolor(/load)
      if round(dv[qq]) LT 0. then pclr = clr.blue else pclr=clr.red
      if abs(round(dv[qq])) LT 10. then pclr = clr.black
      case qq of 
         0: wlbl = 'MgII '
         1: wlbl = 'MgII '
         2: wlbl = 'FeII '
         3: wlbl = 'FeII* '
         else: stop
      endcase
      xyouts, ymnx[0]+0.05*(ymnx[1]-ymnx[0]), ymnx[0]+0.90*(ymnx[1]-ymnx[0]),  $
              wlbl+strtrim(fix(wrest),2)+'; v='+strtrim(round(dv[qq]),2)+' km/s', $
              color=clr.black, charsiz=lsz
  endfor
  
  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]

  return

end
