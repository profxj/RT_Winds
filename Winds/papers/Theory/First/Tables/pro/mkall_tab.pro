pro mkall_tab

  ;; Fiducial
  mktab_fiducial_ew

  ;; Big Tables
  mktab_measure_models, TITLE='Line Diagnostics for the Fiducial Model and Variants', $
                        LBL='tab:line_diag'
  mktab_measure_models, 'tab_meas_plaws.tex', INFIL='Input/tab_meas_plaws.inp', $
                        TITLE='Line Diagnostics for the Fiducial Model and Variants', $
                        LBL='tab:plaw_diag'

  return
end
