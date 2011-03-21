pro emit_lya_line, STRCT=strct

  if not keyword_set( PSFILE ) then psfile = 'emit_lya_line.ps'
  if not keyword_set(CSZ) then csz = 2.0
  if not keyword_set(lSZ) then lsz = 2.0

  c = x_constants()

  ;; Gaussian
  wrest = 1215.6701d
  sigma = 20. ; km/s
  EW = 100. ; Ang
  sigma_A = sigma * wrest / 3e5 ; Ang

  nrm = 1.
  A0 = EW / sigma_A / sqrt(!pi*2)

  ;; Evaluate
  lmin = 1205. ; Ang
  lmax = 1225. ; Ang
  npt = 10000L
  lval = lmin + (lmax-lmin)*dindgen(npt)/(npt-1)
  dl = lval[1]-lval[0]
  feval = nrm + A0*exp(-1*(lval-wrest)^2 / (2 * sigma_A^2))
  cumul = total(feval,/cumul)*dl 
  mxc = max(cumul)
  cumul = cumul / mxc

  ;; Plot
  x_psopen, psfile, /maxs
  clr = getcolor(/load)

  ;; Profiles
  xmrg = [9,1]
  ymrg = [4.0,1]
  yrng=[1., 1e5]
  xrng=[lmin, lmax]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, $
        ytitle='f_lam', $
        xtitle='!9l!X (Ang)', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, /ylog

  ;; x
  nran = 100000L
  seed = -2244L
  ranv = randomu(seed, nran)
  ranl = dblarr(nran)
  
  ;; SLOW
  for qq=0L,nran-1 do begin
     mn = min(abs(ranv[qq]-cumul),imn)
     ranl[qq] = lval[imn]
  endfor

  ;; Plot
  plothist, ranl, bin=0.1, /overpl, color=clr.black

  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]

  stop
  return

end
