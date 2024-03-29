pro fig_obs_slit, NO_RADIA=no_radia

  if not keyword_set( PSFILE ) then psfile = 'fig_obs_slit.ps'
  if not keyword_set(CSZ) then csz = 2.0
  if not keyword_set(lSZ) then lsz = 1.5
  if not keyword_set(pSZ) then psz = 1.2
  if not keyword_set(XNCOLORS) then xncolors=200L

  if not keyword_set(TAU_RAD) then tau_rad = 0.2

  c = x_constants()

  npt = 100
  slit_ratio = 1.5 * findgen(npt)/npt

  ;; Plot
  x_psopen, psfile, /maxs
  clr = getcolor(/load)

  xmrg = [8,1]
  ymrg = [4.0,1]
  yrng=[0., 2]
  xrng=[0., max(slit_ratio)]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, $
        ytitle='|W!de!N| / W!da!N    [Ang]', $
        xtitle='Slit Width / [2 * r!d!9t!X=0.2!N]', $
        yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Radiative
  if not keyword_set(NO_RADIA) then begin
     fiducial_rtau = 40. ;; kpc  [By hand!]
     
     grid_file='../Analysis/Radiation/Output/radiative_grid.fits'
     idx=[2,3]
     raw_data = xmrdfits(grid_file, idx[0], /silent)
     sz = size(raw_data,/dimen)
     ngrid = sz[0]
     raw_wave = xmrdfits(grid_file, idx[1], /silent)
     raw_data = float(raw_data)
     spec = total(total(raw_data,1),1)
     dwv = abs(raw_wave[1]-raw_wave[0])
     dl = 100. / ngrid      ;; kpc
;     yval = dl*findgen(ngrid) ;; kpc 
     yval = findgen(ngrid) - ngrid/2 + 0.5 ;; Dimensionless
     yval = yval*dl
     
     for qq=0L,2 do begin
        case qq of 
           0: begin
              wrest = 2796.35
              vmnx = [-423., -112.4, 348]
              lsty = 0
           end
           1: begin
              wrest = 2803.531
              vmnx = [-388.6, -142.4, 350]
              lsty = 1
           end
           2: begin
              wrest = 2600.1729
              vmnx = [-1037, -66, 714, 2491, 3363]
              lsty = 2
              idx=[0,1]
              raw_data = xmrdfits(grid_file, idx[0], /silent)
              raw_wave = xmrdfits(grid_file, idx[1], /silent)
              raw_data = float(raw_data)
              spec = total(total(raw_data,1),1)
              dwv = abs(raw_wave[1]-raw_wave[0])
           end
           else: stop
        endcase
        
        vel = (raw_wave-wrest)/wrest * c.c / 1e5 ; km/s
;     if qq EQ 2 then stop
        abs_pix = where(vel GT vmnx[0] and vel LT vmnx[1]) ;; Determined by hand!
        W_abs = total(dwv*(1.-spec[abs_pix])) ;; Ang
        
        em_pix = where(vel GT vmnx[1] and vel LT vmnx[2], n_em)
        if qq GT 1 then begin
           em_pix = [em_pix, where(vel GT vmnx[3] and vel LT vmnx[4], n_em2)]
           n_em = n_em + n_em2
        endif
        W_em  = total(dwv*(1.-spec[em_pix]))
        
        print, wrest, W_abs, W_em
        
        mgii_em = total(total(raw_data[*,*,em_pix],3),1)
        w_ratio = fltarr(npt)
        for jj=0L,npt-1 do begin
           ;; Good cells
           gd = where(abs(yval) LE slit_ratio[jj]*fiducial_rtau, ngd) 
           if ngd GT 0 then begin
              w_ratio[jj] = (dwv*(total(mgii_em[gd])-n_em)) / W_abs
;        print, slit_ratio[jj], w_ratio[jj]*w_abs, w_ratio[jj], ngd
           endif
        endfor
        
        oplot, slit_ratio, w_ratio, color=clr.blue, psym=10, linesty=lsty
     endfor
  endif

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; LBG
  if not keyword_set(fcmax) then fcmax = 0.6 ; Covering fraction
  if not keyword_set(r_min) then r_min = 1.0 ; kpc
  if not keyword_set(r_max) then r_max = 100. ; kpc
  if not keyword_set(gamma) then gamma = 0.5 ; Covering fraction
  if not keyword_set(NPTS1) then npts1 = 10000L                  ; Log steps
  if not keyword_set(NPTS2) then npts2 = 10000L                   ; Log steps

  rcut = 1.2
  rval_lo = r_min * 10.^(alog10(rcut) * findgen(npts1) / npts1) ; kpc
  rval_hi = (r_min*rcut) * 10.^(alog10(r_max/r_min/rcut) * findgen(npts2) / (npts2-1)) ; kpc
  rval = [rval_lo, rval_hi]

  fc_lbg = fcmax * (rval/r_min)^(-1*gamma)
  I_lbg = 1 - fc_lbg
  tau_r = -1*alog(I_lbg)

  mn = min(abs(tau_r - tau_rad), imn)
  fiducial_rtau = rval[imn] ;; kpc
  print, 'LBG: rtau = ', fiducial_rtau

  grid_file='../Analysis/LBG/Covering/3D/Models/lbg_s10_grid_mgii.fits'
  idx=[0,1] ;; MgII
  raw_data = xmrdfits(grid_file, idx[0], /silent)
  sz = size(raw_data,/dimen)
  ngrid = sz[0]
  raw_wave = xmrdfits(grid_file, idx[1], /silent)
  raw_data = float(raw_data)
  spec = total(total(raw_data,1),1)
  dwv = abs(raw_wave[1]-raw_wave[0])
  dl = 100. / ngrid ;; kpc
  yval = dl*findgen(ngrid) ;; kpc 


;  for qq=0L,2 do begin
  for qq=0L,1 do begin
     case qq of 
        0: begin
           wrest = 2796.35
           vmnx = [-750., -65., 180]
           lsty = 0
        end
        1: begin
           wrest = 2803.531
           vmnx = [-570., -70., 680]
           lsty = 1
        end
        2: begin
           wrest = 2600.1729
           vmnx = [-1037, -66, 714, 2491, 3363]
           lsty = 2
           idx=[0,1]
           raw_data = xmrdfits(grid_file, idx[0], /silent)
           raw_wave = xmrdfits(grid_file, idx[1], /silent)
           raw_data = float(raw_data)
           spec = total(total(raw_data,1),1)
           dwv = abs(raw_wave[1]-raw_wave[0])
           stop
        end
        else: stop
     endcase
        
     vel = (raw_wave-wrest)/wrest * c.c / 1e5 ; km/s
;     if qq EQ 2 then stop
     abs_pix = where(vel GT vmnx[0] and vel LT vmnx[1]) ;; Determined by hand!
     W_abs = total(dwv*(1.-spec[abs_pix])) ;; Ang

     em_pix = where(vel GT vmnx[1] and vel LT vmnx[2], n_em)
     if qq GT 1 then begin
        em_pix = [em_pix, where(vel GT vmnx[3] and vel LT vmnx[4], n_em2)]
        n_em = n_em + n_em2
     endif
     W_em  = total(dwv*(1.-spec[em_pix]))

     print, wrest, W_abs, W_em

     mgii_em = total(total(raw_data[*,*,em_pix],3),1)
     w_ratio = fltarr(npt)
     for jj=0L,npt-1 do begin
        ;; Good cells
        gd = where(abs(yval) LE slit_ratio[jj]*fiducial_rtau, ngd) 
        if ngd GT 0 then begin
           w_ratio[jj] = (dwv*(total(mgii_em[gd])-n_em)) / W_abs
;        print, slit_ratio[jj], w_ratio[jj]*w_abs, w_ratio[jj], ngd
        endif
     endfor

     oplot, slit_ratio, w_ratio, color=clr.red, psym=10, linesty=lsty
  endfor


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Fiducial
  fig_nvtau_vs_r, strct=fiducial_strct 
  mn = min(abs(fiducial_strct.tau - tau_rad), imn)
  fiducial_rtau = fiducial_strct.rval[imn] ;; kpc

  grid_file='../Analysis/Outputs/fiducial_grid.fits'
  idx=[2,3]
  raw_data = xmrdfits(grid_file, idx[0], /silent)
  sz = size(raw_data,/dimen)
  ngrid = sz[0]
  raw_wave = xmrdfits(grid_file, idx[1], /silent)
  raw_data = float(raw_data)
  spec = total(total(raw_data,1),1)
  dwv = abs(raw_wave[1]-raw_wave[0])
  dl = 40. / ngrid ;; kpc
  yval = findgen(ngrid) - ngrid/2 + 0.5  ;; Dimensionless
  yval = yval*dl


  for qq=0L,2 do begin
     case qq of 
        0: begin
           wrest = 2796.35
           vmnx = [-1060., -50.4, 295]
           lsty = 0
        end
        1: begin
           wrest = 2803.531
           vmnx = [-470., -43.6, 905]
           lsty = 1
        end
        2: begin
           wrest = 2600.1729
           vmnx = [-1037, -66, 714, 2491, 3363]
           lsty = 2
           idx=[0,1]
           raw_data = xmrdfits(grid_file, idx[0], /silent)
           raw_wave = xmrdfits(grid_file, idx[1], /silent)
           raw_data = float(raw_data)
           spec = total(total(raw_data,1),1)
           dwv = abs(raw_wave[1]-raw_wave[0])
        end
        else: stop
     endcase
        
     vel = (raw_wave-wrest)/wrest * c.c / 1e5 ; km/s
;     if qq EQ 2 then stop
     abs_pix = where(vel GT vmnx[0] and vel LT vmnx[1]) ;; Determined by hand!
     W_abs = total(dwv*(1.-spec[abs_pix])) ;; Ang

     em_pix = where(vel GT vmnx[1] and vel LT vmnx[2], n_em)
     if qq GT 1 then begin
        em_pix = [em_pix, where(vel GT vmnx[3] and vel LT vmnx[4], n_em2)]
        n_em = n_em + n_em2
     endif
     W_em  = total(dwv*(1.-spec[em_pix]))

     print, wrest, W_abs, W_em

     mgii_em = total(total(raw_data[*,*,em_pix],3),1)
     w_ratio = fltarr(npt)
     for jj=0L,npt-1 do begin
        ;; Good cells
        gd = where(abs(yval) LE slit_ratio[jj]*fiducial_rtau, ngd) 
        if ngd GT 0 then begin
           w_ratio[jj] = (dwv*(total(mgii_em[gd])-n_em)) / W_abs
;        print, slit_ratio[jj], w_ratio[jj]*w_abs, w_ratio[jj], ngd
        endif
     endfor

     oplot, slit_ratio, w_ratio, color=clr.black, psym=10, linesty=lsty
  endfor


  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]

;  x_splot, vel, tau, /bloc
;  stop

  return

end
