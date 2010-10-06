pro fig_lbg_clump_eandp, RREAL=rreal, STRCT=strct

  if not keyword_set( PSFILE ) then psfile = 'fig_lbg_clump_eandp.ps'
  if not keyword_set(CSZ) then csz = 1.9
  if not keyword_set(lSZ) then lsz = 2.0

  if not keyword_set(dv) then dv = 1.
  if not keyword_set(NPTS1) then npts1 = 10000L                  ; Log steps
  if not keyword_set(NPTS2) then npts2 = 10000L                   ; Log steps

  ;; LBG stuff
  if not keyword_set(r_min) then r_min = 1.0 ; kpc
  if not keyword_set(r_max) then r_max = 100. ; kpc
  if not keyword_set(gamma) then gamma = 0.5 ; Covering fraction
  if not keyword_set(fcmax) then fcmax = 0.6 ; 
  if not keyword_set(alpha) then alpha = 1.3 ;
  if not keyword_set(vmax) then vmax = 800. ; 
  if not keyword_set(WREST) then wrest = 2796.352d


  ;; Begin
  c = x_constants()
  rval = r_min + findgen(100)
  dr = 1. ; kpc

  ;;;;;;;;;;;;;;;;
  ;; LBG
  A_lbg = vmax^2 * (1-alpha)
  v_lbg = -1. * sqrt(A_lbg/(1-alpha)) * sqrt(r_min^(1-alpha) - rval^(1-alpha)) ; km/s
  fc_lbg = fcmax * (rval/r_min)^(-1*gamma)
  I_lbg = 1 - fc_lbg
  I_v = 1. - fcmax * (1 - (1-alpha)*v_lbg^2/A_lbg)^(-1*gamma/(1-alpha))

  ;; Mass
  rel_mass = rval^(1.5)

  ;; Mass
  rel_energy = rval^(1.5) * v_lbg^2
  rel_mo = rval^(1.5) * v_lbg

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Plot
  x_psopen, psfile, /maxs
  !p.multi = [0,2,1]
  clr = getcolor(/load)

  xmrg = [9,2]
  ymrg = [3.7,0.5]

  ;;;;;;;;;;;;;;;;;;;;;;;
  ;; Integrated Energy
  xrng=[1e0, 100]
  yrng=[1., 1e5]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, $
        ytitle='Cumulative Energy (Relative to Mass at 2kpc)', $
        xtitle='Radius (kpc)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, /xlog, /ylog

  ;; Energy
  oplot, rval[1:*], total(rel_energy[1:*],/cumul)/rel_energy[1], color=clr.black

  ;;;;;;;;;;;;;;;;;;;;;;;
  ;; Integrated Momentum
  xrng=[1e0, 100]
  yrng=[1., 1e5]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, $
        ytitle='Cumulative Momentum (Relative to Mass at 2kpc)', $
        xtitle='Radius (kpc)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, /xlog, /ylog

  ;; Momentum
  oplot, rval[1:*], total(rel_mo[1:*],/cumul)/rel_mo[1], color=clr.black

  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]

  spawn, 'ps2pdf '+psfile

  return
end
