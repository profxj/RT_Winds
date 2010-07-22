pro fig_obs_vf

  if not keyword_set( PSFILE ) then psfile = 'fig_obs_vf.ps'
  if not keyword_set(CSZ) then csz = 2.0
  if not keyword_set(lSZ) then lsz = 2.0
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
  yrng=[-200, 200]
  xrng=[-100, 300]
  plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
        xmargin=xmrg, ymargin=ymrg, $
        ytitle='v!df!N (FeII* 2612) [km s!u-1!N]', $
        xtitle='v!df!N (MgII 2796) [km s!u-1!N]', $
        yrange=yrng, thick=4, $
        xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata

  ;; Grab out MgII 2796
  mgii = where( round(strct.anly_strct.wrest) EQ 2796 )
  vf_2796 = (strct.anly_strct.vel_flux)[mgii]

  feii = where( fix(strct.anly_strct.wrest) EQ 2612 )
  vf_2612 = (strct.anly_strct.vel_flux)[feii]

  ;; Add a random piece for plotting convenience
  seed = -1244L
  ranv = 10.
  vf_2796 = vf_2796 + ranv*(randomu(seed,nmodel)-0.5) 
  vf_2612 = vf_2612 + ranv*(randomu(seed,nmodel)-0.5) 

  ;; Density
  oplot, vf_2796, vf_2612, color=clr.black, psym=1, symsiz=1.4

  ;; Guides
  xypt = -1000 + findgen(2001)
  oplot, xypt, xypt, color=clr.gray, lines=2

  oplot, xrng, [0., 0], color=clr.gray, lines=1

  ;; Label
  xlbl = 11.
;  xyouts, xlbl, 6., 'v!dr!N (x10!u-1!N)', color=clr.blue, charsiz=lsz

  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]

;  x_splot, vel, tau, /bloc
;  stop

  return

end
