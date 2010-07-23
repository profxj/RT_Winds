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
  intr_mgII = mgII_strct.noscatt_spec
  mgII_wave = mgII_strct.wave
  dwv_mgII = abs(mgii_wave[1]-mgii_wave[0])

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Read in FeII Data
  spec_feII = feII_strct.spec
  intr_feII = feII_strct.noscatt_spec
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
        W_abs: 0., $       ;; Absorption EW (resonance lines only)
        W_int: 0., $   ;; (intrinsic) EW in the absence of scattering/re-emission
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
        intr_spec = intr_feii
        dwv = dwv_feii
     endif else begin
        wave = mgii_wave
        spec = spec_mgii
        intr_spec = intr_mgii
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
           ;; This failed in Models D and G of the Power-Laws
           ;; Kludges
           case round(trans[qq]) of 
              2796: begin 
                 a = where(spec[0:imn] LT 0.95, na)
                 a0 = a[0]
              end
              2803: begin 
                 ;; Find peak 2796 emission
                 mn = min(abs(wave-2796.352),imn2)
                 mx = max(spec[imn2:imn-20],imx2)
                 ;; Grab min from there
                 a = where(spec[imn2+imx2:imn] LT 0.95, na)
                 a0 = imn2+imx2+a[0]
              end
              else: begin
                 a = where(spec[0:imn] GT 0.95, na)
                 a0 = a[na-1]
              end
           endcase
        endif else begin
           ;; Last pixel of absorption
           a = where(spec[0:imn] LT 0.95, na)
           if na EQ 0 then a1 = imn-1 else a1 = a[na-1]
           ;; Find peak absorption
           kdg_pix = 50
           mn = min(spec[a1-kdg_pix:a1],imn2)
           isrch = a1-kdg_pix+imn2
           ;; First pixel of absorption
           a = where(spec[0:isrch] GT 0.95, na)  ;; The '5' is a kludge
           if na EQ 0 then a0 = a1-1 else a0 = a[na-1]
        endelse
        ;; Sum it up
        strct[qq].W_abs = dwv*total( (1.-spec[a0:a1]) )
        strct[qq].vmnx_abs = [vel[a0],vel[a1]]

        ;; Peak optical depth
        strct[qq].tau_peak = -1*alog(min(spec[a0:a1],imnt) > MINI)  > 0.
        if strct[qq].tau_peak GT 0. then begin 
           strct[qq].tau_vel = vel[a0+imnt]
           tau_vals = -1*alog(spec[a0:a1] > MINI)
           strct[qq].vel_tau = total(vel[a0:a1]*tau_vals) / total(tau_vals) ;; km/s
        endif 

        ;; Intrinsic EW
        if fix(trans[qq]) NE 2803 then pix = where(vel GT -1000 and vel LT 200) $
        else pix = where(vel GT -789 and vel LT 200) 
        strct[qq].W_int = dwv*total( (1.-intr_spec[pix]) )
     endif

     ;;;;;;;;;;;;;;;;
     ;; Emission
     if spec[imn] LT 1.02 then begin
        if strct[qq].flg_trans EQ 0 then begin
           ;; The following assumes the line is to the red 
           e = where(spec[imn:*] GT 1.02, na)
           e0 = e[0]+imn
           e = where(spec[e0:*] LT 1.02, na)
           e1 = e[0]+e0
        endif else begin
           ;; Find the peak (hopefully there is one that is positive)
           mx = max(spec[imn-20+lindgen(41)], imx)
           imx = imn-20+imx
           e = where(spec[0:imx] LT 1.02, na)
           e0 = e[na-1]
           e = where(spec[imx:*] LT 1.02, na)
           e1 = e[0]+imn
        endelse
     endif else begin
        e = where(spec[0:imn] LT 1.02, na)
        e0 = e[na-1]
        e = where(spec[imn:*] LT 1.02, na)
        e1 = e[0]+imn
        ;; Hard code 2626 and 2632 to deal with blending
        case round(trans[qq]) of 
           2626: begin
              wcen = (2626.4511+2632.1081d)/2
              mn = min(abs(wave-wcen), imn3)
              e1 = e1 < imn3
           end
           2632: begin
              wcen = (2626.4511+2632.1081d)/2
              mn = min(abs(wave-wcen), imn3)
              e0 = e0 > imn3
           end
           else:
        endcase
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
