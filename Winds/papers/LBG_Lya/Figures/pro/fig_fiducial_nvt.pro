pro fig_fiducial_nvt, STRCT=strct

  if not keyword_set( PSFILE ) then psfile = 'fig_fiducial_nvt.ps'
  if not keyword_set(CSZ) then csz = 2.0
  if not keyword_set(lSZ) then lsz = 2.0
  if not keyword_set(XNCOLORS) then xncolors=200L

  if not keyword_set(dv) then dv = 1.
  if not keyword_set(NPTS) then npts = 1000L  ; Log steps
  if not keyword_set(v_0) then v_0 = 50.  ; km/s
  if not keyword_set(n_0) then n_0 = 0.1  ; cm^-3
  if not keyword_set(b_val) then b_val = 20. ; km/s
  if not keyword_set(DUST) then dust = 0.1  ; Depletion
  if not keyword_set(METAL) then metal = 0.0  ; [M/H]

  c = x_constants()
  ;; Radius
  r0 = 1.  ; kpc
  r1 = 100. ; kpc
  rval = r0 * exp( alog(r1/r0) * findgen(npts)/float(npts-1) )
  dr = rval - shift(rval,1)  ; kpc
  dr[0] = dr[1] 

  ;; Velocity
;  v_r = v_0 * sqrt(rval / r0)
  v_r = v_0 * rval / r0

  ;; Density
  n_r = n_0 * (rval/r0)^(-2)

  ;; Optical depth
  wrest = 1215.6701d
  lya = x_setline(wrest)
  lines = replicate(lya, npts)
  lines.b = b_val
  lines.N = alog10(dr * c.kpc * n_r ) ;* 10.^(7.53-12.+METAL) * DUST)
  lines.zabs = v_r/3e5

  npix = 2000L
  wav = 10.^(alog10(1215.6701) + dindgen(npix)*1.449d-6)
  vel = (wav-wrest)/wrest * 3e5
  dvel = median(vel-shift(vel,1))  ; Should be 1 km/s

  fx = x_voigt(wav, lines, /nosmooth, TAU=tau)

  ;; Sobolev
  dvdr = v_0/r0
  getfnam, wrest, fval, nam
  kappa_l = !pi*c.e^2/c.me/c.c * fval
  n_Mg = n_r ; * 10.^(7.53-12.+METAL) * DUST
  tau_r = n_Mg*kappa_l / (dvdr*1e5/c.kpc) * (wrest*1e-8) 

  if arg_present(STRCT) then begin
     strct = { $
             rval: rval, $
             wrest: wrest, $
             wave: wav, $
             fval: lya.f, $
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
  yrng=[1, 1e7]
  xrng=[r0, r1]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, $
        ytitle='n!dH!u!N [x10!u6!N, cm!u-3!N];   v!dr!N;  ' + $
        '!9t!X!S!uLy!9a!X!R!N!dS!N', $
        xtitle='Radius (kpc)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, /ylog, /xlog

  ;; Density
  oplot, rval, n_r*1e6, color=clr.red, linest=2

  ;; Velocity
  oplot, rval, v_r, color=clr.blue, linesty=1

  ;; Tau
  r_tau = fltarr(npix)
  mnv = min(v_r, max=mxv)
  gd = where(vel GT mnv and vel LT mxv, ngd, complement=bad)
  for ii=0L,ngd-1 do begin
     mn = min(abs(vel[gd[ii]]-v_r), imn)
     r_tau[gd[ii]] = rval[imn]
  endfor
  oplot, rval, tau_r, color=clr.black
  oplot, r_tau[gd], tau[gd], color=clr.black, linesty=1

  ;; Label
  xlbl = 12.
  xyouts, xlbl, 1e2, 'v!dr!N', color=clr.blue, charsiz=lsz
  xyouts, xlbl, 1e6, '!9t!X!S!uLy!9a!X!R!N!dS!N', color=clr.black, charsiz=lsz
  xyouts, xlbl, 1e4, 'n!dH!N (x10!u6!N)', color=clr.red, charsiz=lsz

  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]

;  x_splot, vel, tau, /bloc
;  stop

  return

end
