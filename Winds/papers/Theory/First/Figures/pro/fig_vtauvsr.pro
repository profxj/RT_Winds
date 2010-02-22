pro fig_vtauvsr

  if not keyword_set( PSFILE ) then psfile = 'fig_vtauvsr.ps'
  if not keyword_set(CSZ) then csz = 1.7
  if not keyword_set(lSZ) then lsz = 1.5

  c = x_constants()

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Input values
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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

  ;; tau_0 for 2796 (density)
  eval = where(rval GE rw_inner)
  density = fltarr(ngrid)  ;; Hydrogen (cm^-3)
  density[eval] = dens_inner * (rval[eval]/rw_inner)^rexpon

  fval = 0.6123  ;; MgII 2796
  wave = 2796.352 * 1e-8 ;; cm
  NMgII = density * dl * c.kpc * 10.^(7.53-12)
  tau0 = (sqrt(!pi) * c.e^2 / c.me / c.c) * wave * fval * NMgII / (doppler*1e5)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Plot
  x_psopen, psfile, /maxs
  clr = getcolor(/load)
  !p.multi=[0,1,1]
  xmrg = [8,7]
  ymrg = [4.0,1]
  yrng=[0., 400]
  xrng=[5., 10.]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, xtitle='r (kpc)', $
        ytitle='v!dr!N (km s!u-1!N)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=9, xstyle=1, psym=1, /nodata

  ;; Velocity
  oplot, rval, speed, color=clr.black, psym=10

  ;;;;;;;;;;
  ;; tau0
  yrng2 = [0., 30.]
  plot, [0],  [0],  psym=psym, xrange=xrng, $
    yrange=yrng2, color=clr.black, $
    xmarg=xmrg, ymarg=ymrg, charsiz=csz, $
    background=clr.white, ystyle=5, xstyle=5, $
    /nodata, /noerase;, /ylog
  axis, yaxis=1, xrange=yrng, charsize=csz, ystyle=1;, titl='!9t!X!d0!N'

  oplot, rval, tau0, color=clr.blue, psym=10
  xyouts, 10.4, yrng2[1]/2, '!9t!X!d0!N', color=clr.black,charsi=csz,$
          orient=90.

  if keyword_set( PSFILE ) then x_psclose

  return
end
