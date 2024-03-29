pro fig_ism_diagn, RREAL=rreal, STRCT=strct

  if not keyword_set( PSFILE ) then psfile = 'fig_ism_diagn.ps'
  if not keyword_set(CSZ) then csz = 2.0
  if not keyword_set(lSZ) then lsz = 2.0
  if not keyword_set(XNCOLORS) then xncolors=200L

  if not keyword_set(dv) then dv = 1.
  if not keyword_set(NPTS) then npts = 500L                  ; Log steps
  if not keyword_set(v_0) then v_0 = 50.  ; km/s
  if not keyword_set(n_0) then n_0 = 0.1  ; cm^-3
  if not keyword_set(b_val) then b_val = 15. ; km/s
  if not keyword_set(DUST) then dust = 0.1  ; Depletion
  if not keyword_set(METAL) then metal = -0.3  ; [M/H]

  c = x_constants()
  r0 = 1.  ; kpc
  r1 = 20. ; kpc
  n_ISM    =  1.           ; Density of the ISM (cm^-3)
  r_ISM    =  0.5          ; Inner radius of ISM
  b_ISM    =  40.          ; Doppler within the ISM

  ;; Radius
  rval = r0 * exp( alog(r1/r0) * findgen(npts)/float(npts-1) )
  rval = [r0*findgen(50)/50, rval] > 1e-3
  ISM_pix = where(rval GE r_ISM and rval LT r0)
  npts = npts + 50
  dr = rval - shift(rval,1)  ; kpc
  dr[0] = dr[1] 


  ;; Velocity
;  v_r = v_0 * sqrt(rval / r0)
  v_r = v_0 * rval * (rval GE r0) / r0

  ;; Density
  n_r = n_0 * (rval/r0)^(-2)
  n_r[ISM_pix] = n_ISM
  n_r[0:ISM_pix[0]-1] = 0.

  ;; Optical depth
  wrest = 2796.352d
  mgii = x_setline(wrest)
  lines = replicate(mgii, npts)
  lines.b = b_val
  lines[ISM_pix].b = b_ISM
  lines.N = alog10(dr * c.kpc * n_r * 10.^(7.53-12.+METAL) * DUST)
  lines.zabs = v_r/3e5
  gdlin = where(lines.N GT 0.)
  lines = lines[gdlin]

;  printcol, lines[0:50].N, lines.zabs, lines.b

  npix = 2000L
  wav = 10.^(alog10(2795.) + dindgen(npix)*1.449d-6)
  vel = (wav-wrest)/wrest * 3e5
  dvel = median(vel-shift(vel,1))  ; Should be 1 km/s

  fx = x_voigt(wav, lines, /nosmooth, TAU=tau)

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
  yrng=[0.01, 250.]
  xrng=[r_ISM, r1]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, $
        ytitle='n!dH!u!N [x10!u2!N, cm!u-3!N];   v!dr!N [x10!u-2!N km s!u-1!N];  ' + $
        '!9t!X!d2796!N', $
        xtitle='Radius (kpc)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, /ylog, /xlog

  ;; Density
  oplot, rval, n_r*100, color=clr.red, linesty=2

  ;; Velocity
  oplot, rval, v_r/100., color=clr.blue, linesty=1
  oplot, [1., 1.], [min(v_r/100), yrng[0]], color=clr.blue, lines=1

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
  xlbl = 10.
  xyouts, xlbl, 12., 'v!dr!N (x10!u-2!N)', color=clr.blue, charsiz=lsz
  xyouts, xlbl, 0.7, '!9t!X!d2796!N', color=clr.black, charsiz=lsz
  xyouts, xlbl, 0.02, 'n!dH!N (x10!u2!N)', color=clr.red, charsiz=lsz

  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]

  return

end
