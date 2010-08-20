model_file = "Models/param_lbg_s10.mod"
line_file = "mg2.line"

-- number of photons to use
n_photons = 2e7

-- parameters of output spectrum/image
l_start = 2770  -- start wavelength
l_stop  = 2830  -- stop wavelength
l_delta = 0.25  -- wavelength bin size
l_nx    = 100   -- spatial pixels per dimension
l_xmax  = 100.0  -- spatial extent of the image (kpc)

-- observer viewing angle
O_theta = 0.0
O_phi   = 0.0
