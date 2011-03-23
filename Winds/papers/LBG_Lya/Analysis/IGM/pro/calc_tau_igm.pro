pro calc_tau_igm, tau=tau

  c = x_constants()

  ;; Cosmology
  z = 3.
  h = 0.72
  hubb = cosm_hubble(z, /w05map)
  nH_bar = 0.04 * (c.rhoc * h^2) * (1+z)^3 / c.mp ;; Should correct for He

  ;;
  nval = 100L

  ;; Model the Halo
  rvir = 200.  ; pMpc
  r = rvir + 9*rvir*findgen(nval)/(nval-1)

  ;; Model density
  delta = 20. / (r/rvir) - 1
  nH_r = (1+delta)*nH_bar ; cm^-3

  ;; Velocity  (Following Fig 3 of Santos)
  ;; This includes Hubble flow!!!
  v_1rvir = -150. ; km/s
  v_r = v_1rvir * (1-alog10(r/rvir))   
  stop

  return
end
