pro fig_obs_peaktau

  if not keyword_set( PSFILE ) then psfile = 'fig_obs_peaktau.ps'
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
  yrng=[0., 3.9]
  xrng=[-900, 100]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, $
        ytitle='Observed !9t!X!dpk!N',$ 
        xtitle='v!d!9t!X!N  [km s!u-1!N]', $
        yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata

  seed = -1244L
  ranv = 0.2

  if keyword_set(FEII) then begin
     ;; FeII 2600
     feiib= where( fix(strct.anly_strct.wrest) EQ 2600 AND $)
            strct.anly_strct.W_abs GT 0.3, ngd)
     vt_2600 = (strct.anly_strct.tau_vel)[feiib]
     tau_2600 = (strct.anly_strct.tau_peak)[feiib]
     tau_2600 = tau_2600 + ranv*(randomu(seed,ngd)-0.5) 
     oplot, vt_2600, tau_2600, color=clr.red, psym=2, symsiz=psz
     
     ;; FeII 2586
     feiia= where( fix(strct.anly_strct.wrest) EQ 2586 AND $)
            strct.anly_strct.W_abs GT 0.3, ngd)
     vt_2586 = (strct.anly_strct.tau_vel)[feiia]
     tau_2586 = (strct.anly_strct.tau_peak)[feiia]
     tau_2586 = tau_2586 + ranv*(randomu(seed,ngd)-0.5) 
     oplot, vt_2586, tau_2586, color=clr.blue, psym=4, symsiz=psz
  endif

  ;; MgII 2796
  mgii = where( round(strct.anly_strct.wrest) EQ 2796 AND $)
         strct.anly_strct.W_abs GT 0.3, ngd)
  vt_2796 = (strct.anly_strct.tau_vel)[mgii]
  tau_2796 = (strct.anly_strct.tau_peak)[mgii]
  valu = where(tau_2796 LT 2.5,nv, complement=lim)
  tau_2796[valu] = tau_2796[valu] + ranv*(randomu(seed,nv)-0.5) 
  oplot, vt_2796[valu], tau_2796[valu], color=clr.black, psym=1, symsiz=psz

  plotsym, 2, 2.0, thick=4
  oplot, vt_2796[lim], tau_2796[lim], color=clr.black, psym=8, symsiz=psz

  ;; Label
  xlbl = -800
  ystp = 0.15
  ooff = 0.09
  toff = 0.05
  if keyword_set(FEII) then begin
     xyouts, xlbl, yrng[1]-ooff-(ystp*1), 'MgII 2796', color=clr.black, charsiz=lsz
     oplot, [xlbl-10], [yrng[1]-toff-(ystp*1)], psym=1, color=clr.black

     xyouts, xlbl, yrng[1]-ooff-(ystp*2), 'FeII 2600', color=clr.red, charsiz=lsz
     oplot, [xlbl-10], [yrng[1]-toff-(ystp*2)], psym=2, color=clr.red
     
     xyouts, xlbl, yrng[1]-ooff-(ystp*3), 'FeII 2586', color=clr.blue, charsiz=lsz
     oplot, [xlbl-10], [yrng[1]-toff-(ystp*3)], psym=4, color=clr.blue
  endif

  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]

;  x_splot, vel, tau, /bloc
;  stop

  return

end
