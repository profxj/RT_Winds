model_file = "Models/param_lbg_s10.mod"
line_file = "fe_uv1.lines"

-- number of photons to use
n_photons = 1e8

-- parameters of output spectrum/image
l_start = 2575  -- start wavelength
l_stop  = 2635  -- stop wavelength
l_delta = 0.25  -- wavelength bin size
l_nx    = 100   -- spatial pixels per dimension
l_xmax  = 100.0  -- spatial extent of the image (kpc)

-- observer viewing angle
O_theta = 0.0
O_phi   = 0.0
