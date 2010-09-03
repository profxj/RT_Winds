pro mkall_fig

  if keyword_set(ALL) then  begin
     ifu = 1
     setup = 1
     oned = 1
  endif

;;;; FIDUCIAL ;;;;
  ;; Setup
  fig_nvtau_vs_r

  ;; 1D spec
  fig_fiducial_1d
  fig_noemiss

  ;; Spatial
  ;; MgII
  fig_ifu, 2796.35, [-250., -100., 0, 200], $
           '../Analysis/Outputs/fiducial_grid.fits', 'fig_fiducial_ifu_mgii.ps', $
           ymnx=[-20,20], xrng=[-600, 500], yrng=[0., 2.5]

  fig_feii_ifu
;  fig_ifu, 2600.173, [-100., 0.], $
;           '../Analysis/Outputs/fiducial_grid.fits', 'fig_fiducial_ifu_feii2600.ps', $
;           ymnx=[-20,20], xrng=[-600, 500], yrng=[0., 1.8]

;  fig_ifu, 2612.654, [-100., 0.], $
;           '../Analysis/Outputs/fiducial_grid.fits', 'fig_fiducial_ifu_feii2612.ps', $
;           ymnx=[-20,20], xrng=[-600, 500], yrng=[0.9, 1.3]

;  fig_spatial_cut, 2796.35, [-200., -100., 0, 200], $
;           '../Analysis/Outputs/fiducial_grid.fits', 'fig_scut_fiducial_mgii.ps', $
;           xrng=[-20,20], yrng=[1e-4, 1e0]
  fig_fiducial_cut

;;;;;;;;;;;;;;;;
;; Complexity

  ;; Asymmetry
  fig_asymm_spec
  fig_bicon_spec

  ;; Dust
  fig_doppler_spec

  ;; Dust
  fig_dust_spec
  fig_taud_vs_we

  ;; ISM
  fig_ism_diagn
  fig_ism_spec
  fig_ism_spec, /dust
  fig_ism_ifu

  ;; Normalize
  fig_norm_spec

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Other Laws

  ;; Varying the density and velocity laws
  fig_vary_profiles
  fig_plaw_spec

  ;; Radiation driven
  fig_radiation_nvt
  fig_rad_spec

  ;; LBG
  fig_lbg_sobolev
  fig_lbg_spec
  fig_lbg_cumul
  fig_lbg_ifu

  ;; Clump
  fig_lbg_clump_mass
  fig_lbg_clump_eandp

;;;;;;;;;;;;;;;;
;; Observations
  fig_obs_vf
  fig_obs_vtau
  fig_obs_peaktau
  fig_obs_ew
  fig_obs_edelv
  fig_obs_slit
  fig_obs_sb
  fig_obs_lris
  fig_obs_lris, /fe

;  fig_slitwidth


  return
end
