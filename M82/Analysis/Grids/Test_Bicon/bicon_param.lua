-- model_file = "const_vr_grid.asc"
-- model_file = "wind_grid_asymmetric.asc"
-- model_file = "wind_grid_fiducial.asc"
model_file = "../tst_biconical.asc"
line_file = "mg2.line"

-- number of photons to use
n_photons = 1e6

--  factor to rescale opacity by (metallicity and density)
opac_factor = 0.06

-- parameters of output spectrum/image
l_start = 2770  -- start wavelength
l_stop  = 2830  -- stop wavelength
l_delta = 0.25  -- wavelength bin size
l_nx    = 100   -- spatial pixels per dimension
l_xmax  = 40.0  -- spatial extent of the image

-- observer viewing angle
O_phi   = 0.0 --  One the sky rotation (I think)
O_theta = 0.0  -- Views biconical down the barrel (90=sideways)

-- uv background (ignore for now)
n_uv =  0
L_uv = 1.25e39
iter_uv = 1
-- T_bb = 40000
