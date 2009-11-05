pro grid_wind

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Input values
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ngrid = 30L
  box_size = 20.           ; Total box size (kpc)
  dl = box_size / ngrid    ; Cell size (kpc) 
  v_wind  = 200.           ; Wind speed (km/s)
  dv_wind = 200.           ; Velocity width of wind (km/s)
  rw_inner = 7.            ; Inner boundary of the wind (kpc)
  rw_outer = 10.           ; Outer boundary of the wind (kpc)
  rg_outer =  3.           ; Outer boundary of galaxy (kpc)
  doppler = 20.            ; Doppler parameter (km/s)
  density = 1e-2           ; Density within the wind

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Construct the Grid
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  xval = findgen(ngrid) - ngrid/2 + 0.5  ;; Dimensionless
  cell_center = fltarr(ngrid,ngrid,ngrid,3)
  xx = xval # replicate(1., ngrid) 
  yy = replicate(1., ngrid) # xval
  for k=0L,ngrid-1 do cell_center[*,*,k,0] = xx
  for i=0L,ngrid-1 do cell_center[i,*,*,1] = xx
  for i=0L,ngrid-1 do cell_center[i,*,*,2] = yy

  ;;;;;;;
  ;; Corners
  x_corners = fltarr(ngrid,ngrid,ngrid,8)
  y_corners = fltarr(ngrid,ngrid,ngrid,8)
  z_corners = fltarr(ngrid,ngrid,ngrid,8)

  ;; x
  xc = 0.5*[-1,1,-1,1,-1,1,-1,1]
  tmp = (cell_center[*,*,*,0])[*] 
  for qq=0,7 do x_corners[*,*,*,qq] = tmp + xc[qq]
  ;; y
  yc = 0.5*[1,1,-1,-1,1,1,-1,-1]
  tmp = (cell_center[*,*,*,1])[*] 
  for qq=0,7 do y_corners[*,*,*,qq] = tmp + yc[qq]
  ;; z
  zc = 0.5*[-1,-1,-1,-1,1,1,1,1]
  tmp = (cell_center[*,*,*,2])[*] 
  for qq=0,7 do z_corners[*,*,*,qq] = tmp + zc[qq]

  ;;;;;;;;;;;;;;;;;;;;;
  ;; Density
  ;;;;;;;;;;;;;;;;;;;;;
  r_corners = dl * sqrt( x_corners^2 + y_corners^2 + z_corners^2)
  delvarx, x_corners, y_corners, z_corners
  msk_corners = intarr(ngrid,ngrid,ngrid,8)
  gd_corn = where(r_corners GE rw_inner and r_corners LE rw_outer)
  msk_corners[gd_corn] = 1B

  frac_cell = total(msk_corners,4)
  rho_grid = density * float(frac_cell) / 8.

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Velocity
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  r_cell = dl * reform( sqrt( cell_center[*,*,*,0]^2 + $
                              cell_center[*,*,*,1]^2 + $
                              cell_center[*,*,*,2]^2), $
                        ngrid,ngrid,ngrid)
  speed_cell = v_wind - dv_wind/2. + dv_wind*(r_cell-rw_inner)/(rw_outer-rw_inner) 
  vx_grid = fltarr(ngrid,ngrid,ngrid)
  vy_grid = fltarr(ngrid,ngrid,ngrid)
  vz_grid = fltarr(ngrid,ngrid,ngrid)
  vx_grid[*] = (speed_cell)[*] * (cell_center[*,*,*,0])[*] / r_cell 
  vy_grid[*] = (speed_cell)[*] * (cell_center[*,*,*,1])[*] / r_cell 
  vz_grid[*] = (speed_cell)[*] * (cell_center[*,*,*,2])[*] / r_cell 

  gd_cell = where(rho_grid GT 1e-3, complement=bad_cell, ngd)
  vx_grid[bad_cell] = 0.
  vy_grid[bad_cell] = 0.
  vz_grid[bad_cell] = 0.
  dopp = fltarr(ngrid,ngrid,ngrid)
  dopp[gd_cell] = doppler

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Emissivity
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Output
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  outfil = 'wind_grid.asc'
  close, /all
  openw, outfil, 11
  for ii=0L,ngrid-1 do begin
      for jj=0L,ngrid-1 do begin
          writecol, 'bah', rho_grid[ii,jj,*], $
                    vx_grid[ii,jj,*], $
                    vy_grid[ii,jj,*], $
                    vz_grid[ii,jj,*], $
                    emiss[ii,jj,*], $
                    dopp[ii,jj,*], $
                    FMT='(f6.4,1x,f6.1,1x,f6.1,1x,f6.1,1x,f4.2)'
  
  close, /all
  return
end
