pro fig_fiducial_cut

  psfile = 'fig_fiducial_cut.ps'
  x_psopen, psfile, /portrait

  !p.multi = [0,1,2]
  ymrg = [4., 0.5]

  ;; MgII
  fig_spatial_cut, 2796.35, [-1500., -200., -100., 0, 200], $
                   '../Analysis/Outputs/fiducial_grid.fits', $
                   xrng=[-20,20], yrng=[1e-4, 1e0], LBL='(a) MgII 2796', $ 
                   ymrg=ymrg

  ;; FeII
  fig_spatial_cut, 2600.173, [-100., 0], $
                   '../Analysis/Outputs/fiducial_grid.fits', $
                   xrng=[-20,20], yrng=[1e-4, 1e0], LBL='(b) FeII 2600', $
                   ymrg=ymrg

  x_psclose
  !p.multi = [0,1,1]

  return
end
