;; calc_ew_values, 'Outputs/fiducial_grid.fits', strct
pro calc_all_values, mgII_strct, feII_strct, strct, TRANS=trans

  if not keyword_set(MINI) then MINI = 0.05
  if not keyword_set(TRANS) then $
     trans = [2796.352,2803.531, 2586.650d, 2600.1729, 2612.6542d, $
              2626.4511, 2632.1081]
  ntrans = n_elements(trans)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Read in MgII Data
  spec_mgII = mgII_strct.spec
  mgII_wave = mgII_strct.wave
  dwv_mgII = abs(mgii_wave[1]-mgii_wave[0])

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Read in FeII Data
  spec_feII = feII_strct.spec
  feII_wave = feII_strct.wave
  dwv_feII = abs(feii_wave[1]-feii_wave[0])

;  x_specplot, spec_feII, fltarr(n_elements(spec_feII)), wave=feII_wave, inf=4, /blo
;  stop
  ;;;;;;;;;
  ;;  
  tmp = { $
        wrest: 0.d, $
        flg_trans: 0, $      ;; 0=Resonance, 1=Fine-structure
        fval: 0., $
        j: 0., $             ;; Angular momentum
        vmnx_abs: fltarr(2), $   ;; End points of W integration
        vmnx_em: fltarr(2), $   ;; End points of W integration
        W_abs: 0., $
        W_em: 0., $
        tau_peak: 0., $   ;; Maximum optical depth of the line
        tau_vel: 0., $   ;; Velocity of maximum optical depth
        vel_tau: 0., $   ;; Tau-weighted velocity centroid for absorption
        flux_peak: 0., $   ;; Peak flux of the line
        flux_vel: 0., $   ;; Velocity of maximum flux
        vel_flux: 0. $   ;; Flux-weighted velocity centroid
        }
  strct = replicate(tmp, ntrans)

  for qq=0L,ntrans-1 do begin
     ;; Load up
     if trans[qq] LT 2700. then begin
        wave = feii_wave
        spec = spec_feii
        dwv = dwv_feii
     endif else begin
        wave = mgii_wave
        spec = spec_mgii
        dwv = dwv_mgii
     endelse

     vel = (wave-trans[qq])/trans[qq] * 3e5

     ;; Read info
     getfnam, trans[qq], fv, nam 
     abslin = x_setline(trans[qq], /overrid)

     strct[qq].wrest = trans[qq]
     strct[qq].fval = fv
     strct[qq].j = abslin.j

     ;; Fine?
     if strpos(nam,'*') GT 0 then strct[qq].flg_trans = 1

     mn = min(abs(wave-trans[qq]),imn)

     ;; Absorption
     if strct[qq].flg_trans EQ 0 then begin
        ;; Find closest pixel to the rest wavelength

        ;; Find the endpoints
        if spec[imn] LT 1. then begin
           ;; Absorption already?
           a = where(spec[imn+1:*] GT 0.95, na)
           if na EQ 0 then stop
           a1 = a[0] + imn
           ;; Returns to the continuum
           a = where(spec[0:imn] GT 0.95, na)
           a0 = a[na-1]
        endif else begin
           ;; Last pixel of absorption
           a = where(spec[0:imn] LT 0.95, na)
           if na EQ 0 then a1 = imn-1 else a1 = a[na-1]
           ;; First pixel of absorption
           a = where(spec[0:a1] GT 0.95, na)
           if na EQ 0 then a0 = a1-1 else a0 = a[na-1]
        endelse
        ;; Sum it up
        strct[qq].W_abs = dwv*total( (1.-spec[a0:a1]) )
        strct[qq].vmnx_abs = [vel[a0],vel[a1]]

        ;; Peak optical depth
        strct[qq].tau_peak = -1*alog(min(spec[a0:a1],imnt) > MINI)
        strct[qq].tau_vel = vel[a0+imnt]
        tau_vals = -1*alog(spec[a0:a1] > MINI)
        strct[qq].vel_tau = total(vel[a0:a1]*tau_vals) / total(tau_vals) ;; km/s
     endif

     ;;;;;;;;;;;;;;;;
     ;; Emission
     if spec[imn] LT 1.02 then begin
        e = where(spec[imn:*] GT 1.02, na)
        e0 = e[0]+imn
        e = where(spec[e0:*] LT 1.02, na)
        e1 = e[0]+e0
     endif else begin
        e = where(spec[0:imn] LT 1.02, na)
        e0 = e[na-1]
        e = where(spec[imn:*] LT 1.02, na)
        e1 = e[0]+imn
     endelse
     ;; Sum it up
     strct[qq].W_em = dwv*total( (1.-spec[e0:e1]) )
     strct[qq].vmnx_em = [vel[e0],vel[e1]]

     ;; Peak flux
     strct[qq].flux_peak = max(spec[e0:e1],imx)
     strct[qq].flux_vel = vel[e0+imx]
     strct[qq].vel_flux= total(vel[e0:e1]*spec[e0:e1]) / total(spec[e0:e1]) ;; km/s

  endfor
  

end
