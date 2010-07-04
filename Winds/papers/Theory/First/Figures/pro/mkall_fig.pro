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
  fig_ifu, 2796.35, [-200., -100., 0, 200], $
           '../Analysis/Outputs/fiducial_grid.fits', 'fig_fiducial_ifu_mgii.ps', $
           ymnx=[-20,20], xrng=[-600, 500], yrng=[0., 2.5]

  fig_ifu, 2600.173, [-150., 0.], $
           '../Analysis/Outputs/fiducial_grid.fits', 'fig_fiducial_ifu_feii2600.ps', $
           ymnx=[-20,20], xrng=[-600, 500], yrng=[0., 1.8]

  fig_ifu, 2612.654, [-150., 0.], $
           '../Analysis/Outputs/fiducial_grid.fits', 'fig_fiducial_ifu_feii2612.ps', $
           ymnx=[-20,20], xrng=[-600, 500], yrng=[0.9, 1.2]

  fig_spatial_cut, 2796.35, [-200., -100., 0, 200], $
           '../Analysis/Outputs/fiducial_grid.fits', 'fig_scut_fiducial_mgii.ps', $
           xrng=[-20,20], yrng=[1e-4, 1e0]

;;;;;;;;;;;;;;;;
;; Complexity

  ;; Asymmetry
  fig_asymm_spec

  ;; Dust
  fig_doppler_spec

  ;; Dust
  fig_dust_spec

  ;; ISM
  fig_ism_diagn
  fig_ism_spec

  ;; LBG
  fig_lbg_sobolev
  fig_lbg_profiles
  fig_lbg_cumul

  ;; Varying the density and velocity laws
  fig_vary_profiles
  fig_plaw_spec

;;;;;;;;;;;;;;;;
;; Observations

;  fig_slitwidth


  return
end
