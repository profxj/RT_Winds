model_file = "Models/param_lbg_s10_si2.mod"
line_file = "si2_1526.lines"

-- number of photons to use
n_photons = 1e8

-- parameters of output spectrum/image
l_start = 1515  -- start wavelength
l_stop  = 1540  -- stop wavelength
l_delta = 0.1   -- wavelength bin size
l_nx    = 100   -- spatial pixels per dimension
l_xmax  = 100.0  -- spatial extent of the image (kpc)

-- observer viewing angle
O_theta = 0.0
O_phi   = 0.0
