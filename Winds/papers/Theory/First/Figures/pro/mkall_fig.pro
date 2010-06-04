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

  ;; IFU
  ;; MgII
  fig_ifu, 2796.35, [-200., 0., 125, 300], $
           '../Analysis/Outputs/fiducial_grid.fits', 'fig_fiducial_ifu_mgii.ps', $
           ymnx=[-20,20], xrng=[-600, 500], yrng=[0., 2.5]

  fig_ifu, 2600.173, [-150., 0.], $
           '../Analysis/Outputs/fiducial_grid.fits', 'fig_fiducial_ifu_feii2600.ps', $
           ymnx=[-20,20], xrng=[-600, 500], yrng=[0., 1.8]

  fig_ifu, 2612.654, [-150., 0.], $
           '../Analysis/Outputs/fiducial_grid.fits', 'fig_fiducial_ifu_feii2612.ps', $
           ymnx=[-20,20], xrng=[-600, 500], yrng=[0.9, 1.2]


;;;;;;;;;;;;;;;;

  fig_slitwidth


  return
end
