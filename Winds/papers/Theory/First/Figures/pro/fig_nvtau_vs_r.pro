pro fig_nvtau_vs_r, RREAL=rreal

  if not keyword_set( PSFILE ) then psfile = 'fig_nvtau_vs_r.ps'
  if not keyword_set(CSZ) then csz = 1.9
  if not keyword_set(lSZ) then lsz = 1.5
  if not keyword_set(XNCOLORS) then xncolors=200L

  if not keyword_set(dv) then dv = 1.
  if not keyword_set(NPTS) then npts = 100L                  ; Log steps
  if not keyword_set(v_0) then v_0 = 70.  ; km/s
  if not keyword_set(n_0) then n_0 = 0.1  ; cm^-3
  if not keyword_set(b_val) then b_val = 15. ; km/s
  if not keyword_set(DUST) then dust = 0.1  ; Depletion
  if not keyword_set(METAL) then metal = 0.  ; [M/H]

  c = x_constants()
  ;; Radius
  r0 = 1.  ; kpc
  r1 = 10.
  rval = r0 * exp( alog(r1/r0) * findgen(npts)/float(npts-1) )
  dr = rval - shift(rval,1)  ; kpc
  dr[0] = dr[1] 

  ;; Velocity
;  v_r = v_0 * sqrt(rval / r0)
  v_r = v_0 * rval / r0

  ;; Density
  n_r = n_0 * (rval/r0)^(-2)

  ;; Optical depth
  wrest = 2796.352d
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
  

  ;; Plot
  x_psopen, psfile, /maxs
  clr = getcolor(/load)


  ;; MgII Spectrum 
  xmrg = [8,1]
  ymrg = [4.0,1]
  yrng=[0.1, 1000.]
  xrng=[1., 10]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, $
        ytitle='n!dH!u!N [x100, cm!u-3!N];   v [km s!u-1!N];  !9t!X!d2796!N', $
        xtitle='Radius (kpc)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, /ylog

  ;; Density
  oplot, rval, n_r*100, color=clr.black

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
  oplot, r_tau[gd], tau[gd], color=clr.red, psym=10, linesty=2

  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]

  return

end
