pro plot_wind,filen


range = [0.5*10^37d,4*10^40d]

nlam = 1
nx   = 1
ny   = 1

openr,1,filen 
readf,1,nlam,nx,ny
readf,1,lam0,dlam,dx


x_arr = dindgen(nx)*dx
y_arr = x_arr
wave  = lam0 + dindgen(nlam)*dlam

data  = 0*dblarr(nx,ny,nlam)
image = 0*dblarr(nx,ny)
flux  = 0*dblarr(nlam)
L_tot = 0.0d


w = 0.0d
f = 0.0d
x = 0.0d
y = 0.0d
c = 0.0d

for i=0,nx-1 do begin
 for j=0,ny-1 do begin
  for k=0,nlam-1 do begin

   readf,1,f,c

   data(i,j,k) = f
   image(i,j) += f
   flux(k) += f
   L_tot += f

  endfor
 endfor
endfor

close,1

;; renormalize
if (nlam gt 1) then L_tot *= (wave(1) - wave(0))
if (nx gt 1) then L_tot *= (x_arr(1) - x_arr(0))^2.0
print,'L_tot = ',L_tot


plot,wave,flux
print,'integrated spectrum, press return to continue'
junk =''
read,junk



loadct,13
tvimage,bytscl(image)



stop

;; Write to FITS
mwrfits, data, 'spec_cube.fits', /create
mwrfits, wave, 'spec_cube.fits' 

end
  
