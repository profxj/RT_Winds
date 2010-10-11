pro fig_obs_edelv

  if not keyword_set( PSFILE ) then psfile = 'fig_obs_edelv.ps'
  if not keyword_set(CSZ) then csz = 2.0
  if not keyword_set(lSZ) then lsz = 2.0
  if not keyword_set(XNCOLORS) then xncolors=200L

  if not keyword_set(dv) then dv = 1.

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
  yrng=[0, 1000]
  xrng=[-900, 0]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, $
        xtitle='v!d!9t!X!N (FeII 2586) [km s!u-1!N]', $
        ytitle='!9D!Xv!de!N (FeII* 2612) [km s!u-1!N]', $
        yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata

  ;; Grab out MgII 2796
  feii = where( fix(strct.anly_strct.wrest) EQ 2586 AND $)
         strct.anly_strct.W_abs GT 0.2, ngd)
  vt_2586 = (strct.anly_strct.vel_tau)[feii]
  sz = size(strct.anly_strct,/dimen)
  s_idx = feii / sz[0]

  ;; Grab out Dv
  feiis = where( fix(strct[s_idx].anly_strct.wrest) EQ 2612 )
  dv_2612 = (strct[s_idx].anly_strct.flux_dv)[feiis]

  ;; Density
  oplot, vt_2586, dv_2612, color=clr.black, psym=4, symsiz=1.4

  ;; Label
  xlbl = 11.
;  xyouts, xlbl, 6., 'v!dr!N (x10!u-1!N)', color=clr.blue, charsiz=lsz

  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]

;  x_splot, vel, tau, /bloc
;  stop

  return

end
