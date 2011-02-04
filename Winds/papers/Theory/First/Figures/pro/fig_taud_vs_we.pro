pro fig_taud_vs_we, RREAL=rreal

  if not keyword_set( PSFILE ) then psfile = 'fig_taud_vs_we.ps'
  if not keyword_set(PAD_FRAC) then pad_frac = 0.1
  if not keyword_set(CSZ) then csz = 1.8
  if not keyword_set(lSZ) then lsz = 1.5
  if not keyword_set(lSZ2) then lsz2 = 1.3
  if not keyword_set(XNCOLORS) then xncolors=200L
  if not keyword_set(YSTP) then ystp = 0.07
  if not keyword_set(VMNX) then vmnx = [-40., 692]
  if not keyword_set(WREST) then wrest = 2803.531

  xlbl = 0.3
  xlbl2 = 0.05
  ylbl = 0.90

  c = x_constants()

  ;; Files
  files = findfile('../Analysis/Dust/Output/spec_MgII_tau*.dat', count=nfil)
  taud = fltarr(nfil+1)
  W_e = fltarr(nfil+1)

  ;; Analyze fiducial
  Mg_fil = '../Analysis/Fiducial/Output/spec_MgII_fiducial.dat'
  readcol, Mg_fil, wv, fx, noscatt_fx, /silen
  nrm = median(fx[where(wv GT 2815)])
  fx = fx/nrm
  dwv = wv[1] - wv[0]
  vel = (wv-wrest)/wv*c.c/1e5
  mn = min(abs(vel-vmnx[0]),e0)
  mn = min(abs(vel-vmnx[1]),e1)
  W_e[0] = dwv*total( (1.-fx[e0:e1]) )

  ;; Analyze dust models
  for qq=0L,nfil-1 do begin
     ;; Set taud
     i1pos = strpos(files[qq], 'tau')
     i2pos = strpos(files[qq], '.dat')
     taud[qq+1] = float(strmid(files[qq], i1pos+3, i2pos-i1pos+2))
     ;; Analyze
     readcol, files[qq], wv, fx, /silen
     nrm = median(fx[where(wv GT 2815)])
     fx = fx/nrm
     dwv = wv[1] - wv[0]
     vel = (wv-wrest)/wv*c.c/1e5
     mn = min(abs(vel-vmnx[0]),e0)
     mn = min(abs(vel-vmnx[1]),e1)
     W_e[qq+1] = dwv*total( (1.-fx[e0:e1]) )
  endfor
  srt = sort(taud)
  taud = taud[srt]
  W_e = W_e[srt]

  ;;; BEGIN PLOTS
  x_psopen, psfile, /maxs
  !p.multi = [0,1,1]
  clr = getcolor(/load)

  xmrg = [8,2]
  ymrg = [4.0,1]

  ;; Plot FeII
  yrng=[0., -2.5]
  xrng=[0., 10.]

  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, ytitle='W!de!N (Ang)', $
        xtitle='!9t!X!dd!N', yrange=yrng, thick=4, $
        xrange=xrng, ystyle=ysty, xstyle=1, psym=1, /nodata, /noerase

  oplot, taud, W_e, color=clr.black, thick=5
  oplot, taud, min(W_e)/(1+taud), color=clr.red, thick=3, linest=1
  oplot, taud, min(W_e)*exp(-1.*taud), color=clr.green, thick=3, linest=2

  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]

  return

end
