pro fig_lbg_clump_deltar, RREAL=rreal, STRCT=strct

  if not keyword_set( PSFILE ) then psfile = 'fig_lbg_clump_deltar.ps'
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
  if not keyword_set(Dv) then Dv = 10. ; km/s


  ;; Begin
  c = x_constants()
  rval = r_min + findgen(10000L)/100
  dr = 1. ; kpc

  ;;;;;;;;;;;;;;;;
  ;; LBG
  A_lbg = vmax^2 * (1-alpha)
  v_lbg = -1. * sqrt(A_lbg/(1-alpha)) * sqrt(r_min^(1-alpha) - rval^(1-alpha)) ; km/s
  fc_lbg = fcmax * (rval/r_min)^(-1*gamma)
  I_lbg = 1 - fc_lbg
  I_v = 1. - fcmax * (1 - (1-alpha)*v_lbg^2/A_lbg)^(-1*gamma/(1-alpha))

  dvdr = sqrt(A_lbg/(1-alpha)) * 0.5 / sqrt(1-rval^(1-alpha)) * (alpha-1) * rval^(-1*alpha)
  dvdr[0] = 2*dvdr[1] ; Kludge

  ;; Delta r
  delta_r = Dv / dvdr

  ;; Sum up to 100kpc from 2kpc
  xval = fltarr(10000L)
  yval = fltarr(10000L)
  xval[0] = 2. ; kpc
  cnt = 0
  while(xval[cnt] LT 100.) do begin
     mn = min(abs(rval-xval[cnt]),imn)
     yval[cnt] = delta_r[imn]
     cnt = cnt+1
     xval[cnt] = xval[cnt-1] + delta_r[imn]
  endwhile

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Plot
  x_psopen, psfile, /maxs
  !p.multi = [0,1,1]
  clr = getcolor(/load)

  xmrg = [10,2]
  ymrg = [3.7,0.5]

  ;;;;;;;;;;;;;;;;;;;;;;;
  ;; Integrated mass
  xrng=[1e0, 100]
  yrng=[1e-3, 10]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, $
        ytitle='!9D!X r (kpc)', $
        xtitle='Radius (kpc)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, /xlog, /ylog

  ;; Delta r
  oplot, rval, delta_r, color=clr.black

;  oplot, xval[0:cnt-1], yval[0:cnt-1], color=clr.red, psym=10

  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]

  return
end
