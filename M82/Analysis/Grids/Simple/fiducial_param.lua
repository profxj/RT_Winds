-- model_file = "const_vr_grid.asc"
-- model_file = "wind_grid_asymmetric.asc"
model_file = "wind_grid_fiducial.asc"
line_file = "mg2.line"

-- number of photons to use
n_photons = 1e7

--  factor to rescale opacity by (metallicity and density)
opac_factor = 0.06

-- parameters of output spectrum/image
l_start = 2770  -- start wavelength
l_stop  = 2830  -- stop wavelength
l_delta = 0.25  -- wavelength bin size
l_nx    = 300   -- spatial pixels per dimension
l_xmax  = 40.0  -- spatial extent of the image

-- observer viewing angle
O_theta = 0.0
O_phi   = 0.0

-- uv background (ignore for now)
n_uv =  0
L_uv = 1.25e39
iter_uv = 1
-- T_bb = 40000
