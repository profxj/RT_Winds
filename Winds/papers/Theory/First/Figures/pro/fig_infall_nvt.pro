pro fig_infall_nvt, RREAL=rreal, STRCT=strct, SHOW=show

  if not keyword_set( PSFILE ) then psfile = 'fig_infall_nvt.ps'
  if not keyword_set(CSZ) then csz = 2.0
  if not keyword_set(lSZ) then lsz = 2.0
  if not keyword_set(XNCOLORS) then xncolors=200L

  if not keyword_set(dv) then dv = 1.
  if not keyword_set(NPTS) then npts = 1000L  ; Log steps
  if not keyword_set(v_0) then v_0 = 250.  ; Normalization (km/s)
  if not keyword_set(n_0) then n_0 = 0.03  ; Normalization cm^-3
  if not keyword_set(b_val) then b_val = 15. ; km/s
  if not keyword_set(DUST) then dust = 0.1  ; Depletion
  if not keyword_set(METAL) then metal = -0.3  ; [M/H]

  xrng=[1, 100.]
  c = x_constants()
  wrest = 2796.352d
  DVDR_MIN = 0.1

  ;; Radius
  R0 = 1. ; kpc
  Rg = 4.  ; kpc
  rval = xrng[0] * exp( alog(xrng[1]/xrng[0]) * findgen(npts)/float(npts-1) )
  dr = rval - shift(rval,1)  ; kpc
  dr[0] = dr[1] 

  ;; Velocity
  v_r = v_0 * sqrt( Rg * (1./R0 - 1./rval) + alog(R0/rval) )
  print, 'sigma = ', v_0 / 2, 'km/s'
  print, 'R_g = ', Rg, ' kpc'

  ;; dv/dr
  dvdr = (v_0/2.) *abs(Rg/rval^2 - 1./rval) /  $
         sqrt( Rg * (1./R0 - 1./rval) + alog(R0/rval) )

  ;; Density
  n_r = n_0 * v_0 / rval^2 / v_r
  n_r[0] = n_r[1]
  n_r[npts-1] = n_r[npts-2]
  print, 'dM/dt = ', n_0*v_0 * (c.kpc^2 * 1e5) * c.mp / c.msun * c.yr, 'Msun/yr'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; INFALL
  r_inf = 50. ; kpc
  n_inf = 1e-3 ; cm^-3
  v_inf = 50. ; km/s
  inf = where(rval GE r_inf)
  ;; Density
  n_r[inf] = n_inf * (rval[inf]/r_inf)^(-2)
  ;; Velocity
  v_r[inf] = -1. * v_inf * rval[inf] / r_inf

  ;; Sobolev
  ;; v turns over at rval = Rg, so this doesnt make much sense
  getfnam, wrest, fval, nam
  kappa_l = !pi*c.e^2/c.me/c.c * fval
  n_Mg = n_r * 10.^(7.53-12.+METAL) * DUST
  tau_r = n_Mg*kappa_l / ((dvdr>DVDR_MIN)*1e5/c.kpc) * (wrest*1e-8) 
  ;if keyword_set(SHOW) then x_splot, rval, tau_r, /blo

  ;; Optical depth
  mgii = x_setline(wrest)
  lines = replicate(mgii, npts)
  lines.b = b_val
  lines.N = alog10(dr * c.kpc * n_r * 10.^(7.53-12.+METAL) * DUST)
  lines.zabs = v_r/3e5

  npix = 2000L
  wav = 10.^(alog10(2795.) + dindgen(npix)*1.449d-6)
  vel = (wav-wrest)/wrest * 3e5
  dvel = median(vel-shift(vel,1))  ; Should be 1 km/s

  fx = x_voigt(wav, lines, /nosmooth, TAU=tau)
  if keyword_set(SHOW) then x_splot, vel, tau, /blo

  if arg_present(STRCT) then begin
     strct = { $
             wrest: wrest, $
             wave: wav, $
             fval: mgii.f, $
             vel: vel, $
             tau: tau $
             }
     return
  endif

  ;; Plot
  x_psopen, psfile, /maxs
  clr = getcolor(/load)


  ;; MgII Spectrum 
  xmrg = [9,1]
  ymrg = [4.0,1]
  yrng=[-15., 50.]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, $
        ytitle='n!dH!u!N [x10!u2!N, cm!u-3!N];   v!dr!N [x10!u-1!N km s!u-1!N];  ' + $
        '!9t!X!d2796!N', $
        xtitle='Radius (kpc)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, /xlog

  ;; Density
  oplot, rval, n_r*100, color=clr.red, linest=2

  ;; Velocity
  oplot, rval, v_r/10., color=clr.blue, linesty=1

  ;; Tau
  r_tau = fltarr(npix)
  mnv = min(v_r, max=mxv)
  gd = where(vel GT mnv and vel LT mxv, ngd, complement=bad)
  for ii=0L,ngd-1 do begin
     mn = min(abs(vel[gd[ii]]-v_r), imn)
     r_tau[gd[ii]] = rval[imn]
  endfor
  oplot, r_tau[gd], tau[gd], color=clr.black

  ;; Label
  xlbl = 11.
  xyouts, xlbl, 6., 'v!dr!N (x10!u-1!N)', color=clr.blue, charsiz=lsz
  xyouts, xlbl, 0.6, '!9t!X!d2796!N', color=clr.black, charsiz=lsz
  xyouts, xlbl, 0.05, 'n!dH!N (x10!u2!N)', color=clr.red, charsiz=lsz

  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]

;  x_splot, vel, tau, /bloc
;  stop

  return

end
