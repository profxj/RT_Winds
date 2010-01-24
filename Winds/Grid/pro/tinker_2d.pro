;; READ
if not keyword_set( GRIDFIL ) then gridfil = 'spec_cube.fits'
raw_data = xmrdfits(gridfil, 0, /silent)
raw_wave = xmrdfits(gridfil, 1, /silent)
dwv = abs(raw_wave[1]-raw_wave[0])
raw_data = float(raw_data)
sz_raw = size(raw_data,/dimen)

;; PAD with continuum
pad_frac = 0.2
pad_raw = round(pad_frac*sz_raw[2]) > 10
sv_data = fltarr(sz_raw[0], sz_raw[1], 2*pad_raw+sz_raw[2])
sz = size(sv_data, /dimen)
sv_data[*,*,pad_raw:pad_raw+sz_raw[2]-1] = raw_data
for ii=0L,pad_raw-1 do $
       sv_data[*,*,ii] = raw_data[*,*,0]
for ii=pad_raw+sz_raw[2],sz[2]-1 do $
       sv_data[*,*,ii] = raw_data[*,*,sz_raw[2]-1]
wave = [raw_wave[0]-reverse(findgen(pad_raw)+1)*dwv, $
        raw_wave, $
        max(raw_wave)+(findgen(pad_raw)+1)*dwv]
;;
rx = sz[2]/10
ry = sz[1]/4

;; SLIT
slit = 1.0
;slit = 0.25
pix = sz[0]/2 + slit*sz[0]/2*[-1,1]
pix = round(pix)
pix[0] = pix[0] > 0L
pix[1] = pix[1] < (sz[0]-1)

;;;;;;;;;;;;;;;;
;; IDEAL
;IDEAL = 1
if keyword_set(IDEAL) then begin
    print, 'Ideal: All data'
    data = sv_data
    dat_slit=data[pix[0]:pix[1],*,*]
    spec2d = transpose(total(dat_slit,1))
    xatv, spec2d, /block, max=max(spec2d)*1.1, min=0., /inver
    
    print, 'Ideal: No central'
    data = sv_data
    data[33:66,33:66,*] = 0.
    dat_slit=data[pix[0]:pix[1],*,*]
    spec2d = transpose(total(dat_slit,1))
    xatv, spec2d, /block, max=max(spec2d)*1.1, min=0., /inver
endif
stop

;;;;;;;;;;;;;;;;
;; 2D Spectral
if keyword_set(SPECTRAL) then begin
    print, 'Spectral: All data'
    data = sv_data
    dat_slit=data[pix[0]:pix[1],*,*]
    spec2d = transpose(total(dat_slit,1))
    print, 'Tot1: ', total(spec2d[*])
    sz2 = size(spec2d,/dimen)
    ;; Smooth in spectral
    fwhm_pix = 3. / abs(wave[1]-wave[0])
    nsmooth = fwhm_pix/(2.*sqrt(2*alog(2)))
    kernel = gauss_kernel(nsmooth)
    for ii=0L,sz2[1]-1 do begin
        tmp = convol(spec2d[*,ii],kernel, /edge_wrap)
        spec2d[*,ii] = tmp
    endfor
    print, 'Tot2: ', total(spec2d[*])
    xatv, spec2d, /block, max=max(spec2d)*1.1, min=0., /inver
endif

;;;;;;;;;;;;;;;;
;; 2D Spatial

;; Create Gaussian kernel
fwhm_pix = 5
nx = 4*round(fwhm_pix)
ny = nx
sig_pix = fwhm_pix / (2*sqrt(2*alog(2)))
A = [0., 1., sig_pix, sig_pix, nx/2., ny/2.]
x = FINDGEN(nx) # REPLICATE(1.0, ny) 
Y = REPLICATE(1.0, nx) # FINDGEN(ny) 
;; Create an ellipse: 
U = ((X-A[4])/A[2])^2 + ((Y-A[5])/A[3])^2 
;; Kernel
kernel2d = A[0] + A[1] * EXP(-U/2) 
kernel2d = kernel2d / total(kernel2d[*])

if keyword_set(SPATIAL) or keyword_set(FULL) then begin
    print, 'Spatial: All data'
    
    ;; Convolve
    smooth_data = fltarr(sz)
    for qq=0L,sz[2]-1 do $
           smooth_data[*,*,qq] = convol(sv_data[*,*,qq], kernel2d, /edge_wrap)
    ;;xatv, smooth_data, /blo
    dat_slit=smooth_data[pix[0]:pix[1],*,*]
    spec2d = transpose(total(dat_slit,1))
    r_spec2d = congrid(spec2d, rx, ry)
    if keyword_set(SPATIAL) then $
      xatv, r_spec2d, /block, max=max(r_spec2d)*1.1, min=0.
endif
    
;;;;;;;;;;;;;;;;
;; FULL
if keyword_set(FULL) then begin
    ;; Use spec2d from SPATIAL
    sz2 = size(spec2d,/dimen)
    ;; Smooth in spectral
    fwhm_pix = 3. / abs(wave[1]-wave[0])
    nsmooth = fwhm_pix/(2.*sqrt(2*alog(2)))
    kernel = gauss_kernel(nsmooth)
    for ii=0L,sz2[1]-1 do $
           spec2d[*,ii] = convol(spec2d[*,ii],kernel, /edge_wrap)
    r_spec2d = congrid(spec2d, rx, ry)
    xatv, r_spec2d, /block, max=max(r_spec2d)*1.1, min=0.
    spec_1d = total(r_spec2d,2)
    r_wave = congrid(wave, rx)
    x_splot, r_wave, spec_1d, /blo
endif

end
