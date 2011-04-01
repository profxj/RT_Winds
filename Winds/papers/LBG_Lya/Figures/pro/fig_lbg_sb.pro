pro fig_lbg_sb, SHOW_INT=show_int

  if not keyword_set(PSFILE) then psfile = 'fig_lbg_sb.ps'
  if not keyword_set(GRID_FILE) then $
     grid_file =  '../Analysis/LBG/Covering/3D/Output/spec_Lya3D_lbg_covering.fits'
  if not keyword_set(CSZ) then csz = 1.7
  if not keyword_set(lSZ) then lsz = 1.5
  if not keyword_set(Z) then z=1.0

  ;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Plot
  x_psopen, psfile, /maxs
  clr = getcolor(/load)
;  xmrg = [9,7]
  xmrg = [9,2]
  ymrg = [3.0,3.5]

  xrng=[0., 100] ; kpc
  yrng=[1e-7, 0.2]

  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, xtitle='Impact Parameter (kpc)', $
        ytitle='Normalized !9m!X (flux per kpc!u2!N)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, /noerase, /ylog

  idx = [0,1] 
  wrest = 1215.6701;; 
  vmnx = [-100, 400.]

  ;; Read in Data
  raw_data = xmrdfits(grid_file, idx[0], /silent)
  raw_wave = xmrdfits(grid_file, idx[1], /silent)
  dwv = abs(raw_wave[1]-raw_wave[0])
  raw_data = float(raw_data)
  sz_raw = size(raw_data,/dimen)

  source_conti = total(raw_data[*,*,0]) ;; Total flux in the source
  
  dr = 100. / sz_raw[0]    ; kpc
  dvec = dr*findgen(sz_raw[0])
  xcell = dvec # replicate(1.,sz_raw[0])
  ycell = replicate(1.,sz_raw[0]) # dvec
  r_cell = sqrt(xcell^2 + ycell^2)
  
;  off = 1
;  z0 = sz_raw[0]/2 - off 
;  z1 = sz_raw[0]/2 + off
  
  ;; Spectrum first
  getfnam, wrest, f, nam
  
  spec = total(total(raw_data,1),1)
  vel = (raw_wave-wrest)/wrest * 3e5 ;; km/s
  gdpix = where(vel GT vmnx[0] and vel LT vmnx[1],ngdp)
  
  source_img = raw_data[*,*,0]*float(ngdp) ;; Total flux in the source
  
  ;; Sum up in these channels
  sum_raw = total(raw_data[*,*,gdpix], 3) 
  wind_img = sum_raw - source_img
  stop

  ;; Calculate the averages
  neval = sz_raw[0]
  avg_r = fltarr(neval-1)
  avg_mu = fltarr(neval-1)
  avg_f = fltarr(neval-1) ;; In a 'shell' that is 1kpc thick
  for qq=0L,neval-2 do begin
     ;; 
     if qq EQ 0 then r0=0. else r0 = r_cell[qq, qq]
;     r1 = r_cell[qq+(sz_raw[0]/2),sz_raw[0]/2]
     r1 = r_cell[qq+1,qq+1]
     ;; Good cells
     gd_cell = where(r_cell GT r0 and r_cell LE (r1*1.01), ngdc)
     if ngdc EQ 0 then stop
     ;; Average mu
     avg_mu[qq] = mean(wind_img[gd_cell]) * dwv ;; Trust me, this works!
     avg_r[qq] = mean(r_cell[gd_cell])
     ;; Average f
     Area = !pi* ( (avg_r[qq]+0.5)^2 - ( (avg_r[qq]-0.5) > 0)^2)
     avg_f[qq] = avg_mu[qq] * Area
  endfor

  avg_mu_kpc = avg_mu / dr^2
  avg_f = avg_f / dr^2
  norm = total(avg_mu_kpc*2*!pi*avg_r)*dr
  print, 'Total EW of line emission = ', norm
  kpc_arcsec = cosm_dist(z, /init, /angular)*1e3
  
  ;; Plot
  oplot,  avg_r, avg_mu_kpc, color=clr.black, psym=10, linesty=lsty
  if keyword_set(SHOW_INT) then $
     oplot,  avg_r, avg_f*1e-2, color=clr.red, linesty=lsty
  
  ;; Label
  oplot, [0.5,0.54]*xrng[1], replicate(yrng[1]/4, 2), color=clr.black, $
         linesty=lsty
  xyouts, 0.55*xrng[1], yrng[1]/4/1.1, 'Wind Model (0 EW; no IGM)', color=clr.black, charsiz=lsz
  oplot, [0.5,0.54]*xrng[1], replicate(yrng[1]/8, 2), color=clr.red, $
         linesty=lsty
  xyouts, 0.55*xrng[1], yrng[1]/8/1.1, 'Observed profile (Steidel+11)', $
          color=clr.red, charsiz=lsz
     
  ;; Data
  xval = findgen(100)
  oplot, xval, 1e-2 * exp(-xval/22.), color=clr.red

  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]

  return

end
