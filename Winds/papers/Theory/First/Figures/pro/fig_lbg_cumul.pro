pro fig_lbg_cumul, RREAL=rreal, STRCT=strct

  if not keyword_set( PSFILE ) then psfile = 'fig_lbg_cumul.ps'
  if not keyword_set(CSZ) then csz = 1.3
  if not keyword_set(lSZ) then lsz = 2.0
  if not keyword_set(XNCOLORS) then xncolors=200L

  if not keyword_set(dv) then dv = 1.
  if not keyword_set(NPTS1) then npts1 = 10000L                  ; Log steps
  if not keyword_set(NPTS2) then npts2 = 10000L                   ; Log steps

  ;; LBG stuff
  if not keyword_set(r_min) then r_min = 1.0 ; kpc
  if not keyword_set(r_max) then r_max = 100. ; kpc
  if not keyword_set(gamma) then gamma = 0.5 ; Covering fraction
  if not keyword_set(fcmax) then fcmax = 0.6 ; Covering fraction
  if not keyword_set(alpha) then alpha = 1.3 ; Covering fraction
  if not keyword_set(vmax) then vmax = 800. ; Covering fraction
  if not keyword_set(WREST) then wrest = 2796.352d


  ;; Begin
  c = x_constants()
  rcut = 1.2
  rval_lo = r_min * 10.^(alog10(rcut) * findgen(npts1) / npts1) ; kpc
  rval_hi = (r_min*rcut) * 10.^(alog10(r_max/r_min/rcut) * findgen(npts2) / (npts2-1)) ; kpc
  rval = [rval_lo, rval_hi]
  npts = npts1+npts2
  dr = rval - shift(rval,1)  ; kpc
  dr[0] = dr[1] 

  ;;;;;;;;;;;;;;;;
  ;; LBG
  A_lbg = vmax^2 * (1-alpha)
  v_lbg = -1. * sqrt(A_lbg/(1-alpha)) * sqrt(r_min^(1-alpha) - rval^(1-alpha)) ; km/s
  fc_lbg = fcmax * (rval/r_min)^(-1*gamma)
  I_lbg = 1 - fc_lbg
;  I_v = 1. - fcmax * (1 - (1-alpha)*v_lbg^2/A_lbg)^(-1*gamma/(1-alpha))

  ;;;;;;;;;;;;;;;;
  ;; Model (i)  [Sobolev]

  ;; Key stuff
  tau_r = -1*alog(I_lbg)
  dvdr = sqrt(A_lbg/(1-alpha)) * 0.5 / sqrt(1-rval^(1-alpha)) * (alpha-1) * rval^(-1*alpha)
  dvdr[0] = 2*dvdr[1] ; Kludge

  ;; Density
  getfnam, wrest, fval, nam
  kappa_l = !pi*c.e^2/c.me/c.c * fval
  n_Mg = tau_r * (dvdr*1e5/c.kpc) / (wrest*1e-8) / kappa_l 
  METAL = -0.3
  DUST = 0.1
  b_val = 5.
  nr_i =  n_Mg /  (10.^(7.53-12.+METAL) * DUST)  ;; Hydrogen

  ;; Column
  N_HI = total(nr_i * dr) * c.kpc
  print, 'N_HI = ', N_HI

  ;; Mass
  mass = 4*!pi*total(nr_i * rval^2 * dr, /cumul) * (c.kpc^3) * c.mp / c.msun

  ;; Energy
  energy = 0.5 * 4*!pi*total(nr_i * v_lbg^2 * rval^2 * dr, /cumul) * (c.kpc^3) * c.mp * 1d10

  ;; Momentum
  mom =  4*!pi*total(nr_i * abs(v_lbg) * rval^2 * dr, /cumul) * c.mp * 1e5

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Plot
  x_psopen, psfile, /maxs
  !p.multi = [0,2,2]
  clr = getcolor(/load)

  xmrg = [8,7]
  ymrg = [3.7,0.5]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Density and Velocity
  xrng=[1e-4, 100]
  yrng=[1e-6, 1.]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, $
        ytitle='n!dH!N (cm!u-3!N)', $
        xtitle='[r-1] (kpc)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=9, xstyle=1, psym=1, /nodata, /ylog, /xlog, /noeras, $
        xtickformat='x_logticks'

  ;; Density
  oplot, rval-1, nr_i, color=clr.black

  ;; Velocity
  yrng=[0., 800]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, $
        xtitle='[r-1] (kpc)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=4, xstyle=4, psym=1, /nodata, /xlog, /noeras
  axis, yaxis=1, charsiz=csz, ysty=1, yrang=yrng, ytit='|v!dr!N| (km s!u-1!N)'

  oplot, rval-1, abs(v_lbg), color=clr.black, linesty=1

  ;;;;;;;;;;;;;;;;;;;;;;;
  ;; Integrated mass
  xrng=[1e0, 100]
  yrng=[2e5, 5e8]
  !p.multi = [3,2,2]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, $
        ytitle='Cumulative Mass (M!dSun!N)', $
        xtitle='Radius (kpc)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, /xlog, /ylog

  ;; Mass
  oplot, rval, mass, color=clr.black

  ;;;;;;;;;;;;;;;;;;;;;;;
  ;; Integrated Energy
  xrng=[1e0, 100]
  yrng=[1d51, 2d57]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, $
        ytitle='Cumulative Energy (ergs)', $
        xtitle='Radius (kpc)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, /xlog, /ylog

  ;; Mass
  oplot, rval, energy, color=clr.black

  ;;;;;;;;;;;;;;;;;;;;;;;
  ;; Momentum
  xrng=[1e0, 100]
  yrng=[1d-19, 2d-15]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, $
        ytitle='Cumulative Momentum (g cm s!u-1!N)', $
        xtitle='Radius (kpc)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, /xlog, /ylog

  ;; Mass
  oplot, rval, mom, color=clr.black


  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]

  return
end
