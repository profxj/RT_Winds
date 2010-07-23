pro fig_obs_ew

  if not keyword_set( PSFILE ) then psfile = 'fig_obs_ew.ps'
  if not keyword_set(CSZ) then csz = 2.0
  if not keyword_set(lSZ) then lsz = 1.5
  if not keyword_set(pSZ) then psz = 1.2
  if not keyword_set(XNCOLORS) then xncolors=200L

  if not keyword_set(dv) then dv = 1.
  if not keyword_set(NPTS) then npts = 1000L  ; Log steps
  if not keyword_set(v_0) then v_0 = 250.  ; Normalization (km/s)
  if not keyword_set(v_1) then v_1 = 1.  ; Scaling in velocity law
  if not keyword_set(n_0) then n_0 = 0.1  ; Normalization cm^-3
  if not keyword_set(b_val) then b_val = 15. ; km/s
  if not keyword_set(DUST) then dust = 0.1  ; Depletion
  if not keyword_set(METAL) then metal = -0.3  ; [M/H]

  c = x_constants()

  ;; Read in the structures
  restore, '../Tables/fiducial_strct.idl'
  strct = all_strct
  restore, '../Tables/other_strct.idl'
  strct = [strct,all_strct[1:*]]  ;; The '1' avoids duplicating the fiducial model
  nmodel = n_elements(strct)

  ;; Plot
  x_psopen, psfile, /maxs
  clr = getcolor(/load)


  ;; MgII Spectrum 
  xmrg = [8,1]
  ymrg = [4.0,1]
  yrng=[-3., 3.]
  xrng=[0.1, 20.]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, $
        ytitle='W!da!N - W!di!N    [Ang]', $
        xtitle='W!di!N    [Ang]', $
        yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, /xlog

;  seed = -1244L
;  ranv = 0.2

  ;; FeII 2600
  feiib= where( fix(strct.anly_strct.wrest) EQ 2600 AND $)
         strct.anly_strct.W_int GT 0.3, ngd)
  Wi_2600 = (strct.anly_strct.W_int)[feiib]
  Wa_2600 = (strct.anly_strct.W_abs)[feiib]
  oplot, Wi_2600, Wa_2600-Wi_2600, color=clr.red, psym=2, symsiz=psz

  ;; FeII 2586
  feiia= where( fix(strct.anly_strct.wrest) EQ 2586 AND $)
         strct.anly_strct.W_int GT 0.3, ngd)
  Wi_2586 = (strct.anly_strct.W_int)[feiia]
  Wa_2586 = (strct.anly_strct.W_abs)[feiia]
  oplot, Wi_2586, Wa_2586-Wi_2586, color=clr.blue, psym=4, symsiz=psz

  ;; MgII 2796
  mgii = where( round(strct.anly_strct.wrest) EQ 2796 AND $)
         strct.anly_strct.W_int GT 0.3, ngd)
  Wi_2796 = (strct.anly_strct.W_int)[mgii]
  Wa_2796 = (strct.anly_strct.W_abs)[mgii]
  oplot, Wi_2796, Wa_2796-Wi_2796, color=clr.black, psym=1, symsiz=psz


  ;; Label
  xlbl = 0.2
  ystp = 0.15
  ooff = 2.09
  toff = 0.05
  xyouts, xlbl, yrng[1]-ooff-(ystp*1), 'MgII 2796', color=clr.black, charsiz=lsz
  oplot, [xlbl-10], [yrng[1]-toff-(ystp*1)], psym=1, color=clr.black

  xyouts, xlbl, yrng[1]-ooff-(ystp*2), 'FeII 2600', color=clr.red, charsiz=lsz
  oplot, [xlbl-10], [yrng[1]-toff-(ystp*2)], psym=2, color=clr.red

  xyouts, xlbl, yrng[1]-ooff-(ystp*3), 'FeII 2586', color=clr.blue, charsiz=lsz
  oplot, [xlbl-10], [yrng[1]-toff-(ystp*3)], psym=4, color=clr.blue

  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]

;  x_splot, vel, tau, /bloc
;  stop

  return

end
