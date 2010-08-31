pro fig_lbg_aclump, RREAL=rreal, STRCT=strct

  if not keyword_set( PSFILE ) then psfile = 'fig_lbg_aclump.ps'
  if not keyword_set(CSZ) then csz = 1.5
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
  rval = r_min + findgen(10000L)/100

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
  delta_r = Dv / dvdr ; kpc


  ;; Sum up to 100kpc from 2kpc
  xval = fltarr(10000L)
  drval = fltarr(10000L)
  mval = fltarr(10000L)
  nval = fltarr(10000L)
  fval = fltarr(10000L)  ;; Flux of clumps
  xval[0] = 2. ; kpc
  cnt = 0
  sigma = 0.01 ; kpc^2
  mc = 250. ; Msun
  while(xval[cnt] LT 100.) do begin
     ;; Dr
     mn = min(abs(rval-xval[cnt]),imn)
     drval[cnt] = delta_r[imn]
     ;; Density
     nval[cnt] = fc_lbg[imn]/sigma / drval[cnt]  ; kpc^-3
     ;; Flux
     fval[cnt] = nval[cnt] * abs(v_lbg[imn]*1e5/c.kpc) * xval[cnt]^2 *c.yr ; yr^-1
     ;; Mass (in the Shell)
     mval[cnt] = 4*!pi*xval[cnt]^2 * fc_lbg[imn] / sigma * mc
     ;; Increment
     cnt = cnt+1
     xval[cnt] = xval[cnt-1] + delta_r[imn]
  endwhile
  xval = xval[0:cnt-1]
  mval = mval[0:cnt-1]
  nval = nval[0:cnt-1]
  fval = fval[0:cnt-1]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Plot
  x_psopen, psfile, /maxs
  !p.multi = [0,2,2]
  clr = getcolor(/load)

  ymrg = [3.7,0.5]

  ;;;;;;;;;;;;;;;;;;;;;;;
  ;; Delta r
  xmrg = [10,2]
  xrng=[1e0, 100]
  yrng=[1e-3, 10]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, $
        ytitle='!9D!X r (kpc)', $
        xtitle='Radius (kpc)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, /xlog, /ylog

  ;; Delta r
  oplot, rval, delta_r, color=clr.black

  ;;;;;;;;;;;;;;;;;;;;;;;
  ;; Density of clumps
  xmrg = [8,2]
  xrng=[1e0, 100]
  yrng=[1, 1e5]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, $
        ytitle='Clump Density (kpc!u-3!N)', $
        xtitle='Radius (kpc)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, /xlog, /ylog

  ;; Density
  oplot, xval, nval, color=clr.black, psym=10

  ;;;;;;;;;;;;;;;;;;;;;;;
  ;; Mass (cumulative)
  xrng=[1e0, 100]
  yrng=[1e5, 1e11]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, $
        ytitle='Cumulative Mass (Msun)', $
        xtitle='Radius (kpc)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, /xlog, /ylog

  ;; Density
  oplot, xval, total(mval,/cumul), color=clr.black

  ;;;;;;;;;;;;;;;;;;;;;;;
  ;; Flux of clumps
  xrng=[1e0, 100]
  yrng=[1e-3, 0.02]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, $
        ytitle='Clump Flux (nvr!u2!N; clumps/yr)', $
        xtitle='Radius (kpc)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, /xlog

  ;; Density
  oplot, xval, fval, color=clr.black


  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]

  return
end
