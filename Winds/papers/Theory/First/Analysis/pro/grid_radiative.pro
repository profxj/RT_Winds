;; Makes the Grid for the fiducial wind model
;;  v, n follow the Murray et al 2005 model 
pro grid_radiative

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Input values
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  if not keyword_set(DUST) then dust = 0.1  ; Depletion
;  if not keyword_set(METAL) then metal = -0.3  ; [M/H]

  ngrid = 300L
  v_wind  = 250.           ; Peak wind speed (km/s) 
  rw_inner = 1.            ; Inner boundary of the wind (kpc)
  rw_outer = 50.           ; Outer boundary of the wind (kpc)
  box_size = rw_outer*2           ; Total box size (kpc)
  dl = box_size / ngrid    ; Cell size (kpc) 
  rg_outer = 0.2           ; Outer boundary of galaxy (kpc)
  doppler =  15.            ; Doppler parameter (km/s)
  n0      =  0.03           ; Density at inner radius (Hydrogen; cm^-3)
  R0      =  1.00           ; Launch radius (kpc)
  Rg      =  4.00           ; Other radial parameter (kpc)

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
  ;; Radius (kpc)
  ;;;;;;;;;;;;;;;;;;;;;
  r_cell = dl * reform( sqrt( cell_center[*,*,*,0]^2 + $
                              cell_center[*,*,*,1]^2 + $
                              cell_center[*,*,*,2]^2), $
                        ngrid,ngrid,ngrid)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Velocity (needs to be done first)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  speed_cell = v_wind * sqrt( (Rg * (1./R0 - 1./r_cell) + alog(R0/r_cell)) > 0. )

  vx_grid = fltarr(ngrid,ngrid,ngrid)
  vy_grid = fltarr(ngrid,ngrid,ngrid)
  vz_grid = fltarr(ngrid,ngrid,ngrid)
  vx_grid[*] = (speed_cell)[*] * (dl*cell_center[*,*,*,0])[*] / r_cell 
  vy_grid[*] = (speed_cell)[*] * (dl*cell_center[*,*,*,1])[*] / r_cell 
  vz_grid[*] = (speed_cell)[*] * (dl*cell_center[*,*,*,2])[*] / r_cell 

  print, 'Velocity: ', max(vx_grid, min=mn), mn

  ;;;;;;;;;;;;;;;;;;;;;
  ;; Density
  ;;;;;;;;;;;;;;;;;;;;;
  r_corners = dl * sqrt( x_corners^2 + y_corners^2 + z_corners^2)
  delvarx, x_corners, y_corners, z_corners
  msk_corners = intarr(ngrid,ngrid,ngrid,8)
  gd_corn = where(r_corners GE rw_inner and r_corners LE rw_outer)
  msk_corners[gd_corn] = 1B

;  density = n0 * ((r_cell>1e-5)/rw_inner)^(-2)
  density = n0 * v_wind / r_cell^2 / (speed_cell > 1e-3)

  frac_cell = total(msk_corners,4)
  rho_grid = density * float(frac_cell) / 8.
  delvarx, msk_corners

  ;;;;;;;;;;;;;;;;;
  ;; Doppler and zeroing out
  ;;;;;;;;;;;;;;;;;

  gd_cell = where(rho_grid GT 0., $
                  complement=bad_cell, ngd)
  vx_grid[bad_cell] = 0.
  vy_grid[bad_cell] = 0.
  vz_grid[bad_cell] = 0.
  dopp = fltarr(ngrid,ngrid,ngrid)
  dopp[gd_cell] = doppler

  gd_cell = where(speed_cell GT 0., $
                  complement=bad_cell, ngd)
  rho_grid[bad_cell] = 0.


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Emissivity
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  msk_corners = intarr(ngrid,ngrid,ngrid,8)
  gd_corn = where(r_corners LT rg_outer)
  msk_corners[gd_corn] = 1B

  frac_cell = total(msk_corners,4)
  emiss = float(frac_cell) / 8.

  stop

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Output
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  outfil = 'wind_grid_radiative.asc'
  close, /all
  openw, 11, outfil
  printf, 11, ngrid, dl
  for ii=0L,ngrid-1 do begin
      for jj=0L,ngrid-1 do begin
          writecol, 'bah', rho_grid[ii,jj,*], $
                    vx_grid[ii,jj,*], $
                    vy_grid[ii,jj,*], $
                    vz_grid[ii,jj,*], $
                    emiss[ii,jj,*], $
                    dopp[ii,jj,*], $
                    FMT='(e12.4,1x,f7.1,1x,f7.1,1x,f7.1,1x,f4.2,1x,f4.1)', $
                    FILNUM=11
      endfor
  endfor
  spawn, 'gzip -f '+outfil
  
  close, /all
  return
end
