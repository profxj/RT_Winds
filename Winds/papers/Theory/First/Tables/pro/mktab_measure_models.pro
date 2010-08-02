;; mktab_summ
pro mktab_measure_models, outfil, all_strct=all_strct, INFIL=infil, TITLE=title, $
                          LBL=lbl, SAV_FIL=sav_fil

  if not keyword_set( INFIL ) then infil = 'Input/tab_meas_models.inp'
  if not keyword_set( OUTFIL ) then outfil = 'tab_meas_models.tex'
  if not keyword_set( SAV_FIL ) then sav_fil = 'measures.idl'
  if not keyword_set(ALL_ANG) then all_ang = [0, 120, 150, 180, 30, 60, 90]
  if not keyword_set(TRANS) then $
     trans = [2796.352,2803.531, 2586.650d, 2600.1729, 2612.6542d, $
              2626.4511, 2632.1081]
  ntran = n_elements(trans)

  ;; Initialize
  compile_opt strictarr

  ;; Analyze
  close, /all
  openr, 1, infil
  readf, 1, nmodels
  model_nm = ''
  data_fil = ''
  flg_anly = 0
  if not keyword_set(ALL_STRCT) then begin
     for qq=0L,nmodels-1 do begin
        ;; Read name
        readf, 1, model_nm
        ;; Read in data
        readf, 1, flg_anly
        case flg_anly of
           1: begin  ;; IFU data cube (standard mode)
              stop  ;; Dont use this mode!!  JXP, 7/20/2010
              readf, 1, data_fil
              mgII_data = xmrdfits(data_fil, 2, /silent)
              mgII_wave = xmrdfits(data_fil, 3, /silent)
              mgII_data = float(mgII_data)
              spec_mgii = total(total(mgii_data,1),1)
              mgII_strct = { $
                           spec: spec_mgii, $
                           wave: mgII_wave $
                           }
              feII_data = xmrdfits(data_fil, 0, /silent)
              feII_wave = xmrdfits(data_fil, 1, /silent)
              feII_data = float(feII_data)
              spec_feii = total(total(feii_data,1),1)
              feII_strct = { $
                           spec: spec_feii, $
                           wave: feII_wave $
                           }
           end
           2: begin ;; outflow.cc output
              ;; MgII
              readf, 1, data_fil
              readcol, data_fil, wv, fx, noscatt_fx, /silen
              nrm = median(fx[where(wv GT 2815)])
              fx = fx/nrm
              nrm2 = median(noscatt_fx[where(wv GT 2815)])
              noscatt_fx = noscatt_fx/nrm2
              mgII_strct = { $
                           fil: data_fil, $
                           spec: fx, $
                           noscatt_spec: noscatt_fx, $
                           wave: wv $
                           }
              ;; FeII
              readf, 1, data_fil
              readcol, data_fil, wv, fx, noscatt_fx, /silen
              nrm = median(fx[where(wv LT 2577)])
              fx = fx/nrm
              nrm2 = median(noscatt_fx[where(wv LT 2577)])
              noscatt_fx = noscatt_fx/nrm2
              feII_strct = { $
                           fil: data_fil, $
                           spec: fx, $
                           noscatt_spec: noscatt_fx, $
                           wave: wv $
                           }
           end
           else: begin ; Asymmetry data cube
              ;; MgII
              readf, 1, data_fil
              jj = where(flg_anly EQ all_ang)
              mgII_wave = xmrdfits(data_fil, 0, /silent)
              mgII_data = xmrdfits(data_fil, jj[0]+1, /silent)
              mgII_data = float(mgII_data)
              spec_mgii = total(total(mgii_data,1),1)
              mgII_strct = { $
                           spec: spec_mgii, $
                           noscatt_spec: replicate(2., n_elements(spec_mgii)), $
                           wave: mgII_wave $
                           }
              ;; FeII
              readf, 1, data_fil
              jj = where(flg_anly EQ all_ang)
              feII_wave = xmrdfits(data_fil, 0, /silent)
              feII_data = xmrdfits(data_fil, jj[0]+1, /silent)
              feII_data = float(feII_data)
              spec_feii = total(total(feii_data,1),1)
              feII_strct = { $
                           spec: spec_feii, $
                           noscatt_spec: replicate(2., n_elements(spec_feii)), $
                           wave: feII_wave $
                           }
           end
        endcase

;        if strmid(model_nm,0,4) EQ 'Radi' then stop
        ;; Analyze
        cd, '../Analysis/pro/', curr=curr
        RESOLVE_ROUTINE, 'calc_all_values'
        calc_all_values, mgII_strct, feII_strct, strct, TRANS=trans
        cd, curr
        if qq EQ 0 then begin
           tmp = { $
                 name: model_nm, $
                 data_fil: data_fil, $
                 anly_strct: strct $
                 }
           all_strct = replicate(tmp, nmodels)
        endif else begin
           all_strct[qq].name = model_nm
           all_strct[qq].data_fil = data_fil
           all_strct[qq].anly_strct = strct
        endelse
     endfor
  endif

  ;; 
  c = x_constants()

  ;; Table header
  close, /all
  openw, 91, outfil

  printf, 91,  ' '
;      printf, 1,  '\clearpage'
  printf, 91,  ' '
  printf, 91,  '\begin{deluxetable}{ccrccccccccccc}'
  printf, 91,  '\rotate'
  printf, 91,  '\tablewidth{0pc}'
  printf, 91,  '\tablecaption{'+TITLE
  printf, 91,  '\label{'+LBL+'}}'
  printf, 91,  '\tabletypesize{\footnotesize}'
  printf, 91,  '\tablehead{\colhead{Transition} & \colhead{Model} & ' + $
          '\colhead{$v_{\rm int}^a$} & ' + $
          '\colhead{$W_{\rm i}$} & ' + $
          '\colhead{$W_{\rm a}$} & \colhead{$\tau_{\rm pk}$} & \colhead{$v_\tau$} '
  printf, 91,  '& \colhead{$v_{\bar \tau}$}'
  printf, 91,  '& \colhead{$v_{\rm int}^b$} & ' + $
    '\colhead{$W_{\rm e}$} & \colhead{$f_{\rm pk}$} & \colhead{$v_f$} '
  printf, 91,  '& \colhead{$v_{\bar f}$}'
  printf, 91,  '\\'
  printf, 91,  '&& (\kms) & (\AA) & (\AA) && (\kms) & (\kms) & (\kms) & (\AA) & & (\kms) & (\kms)}'
  printf, 91, '\startdata'

  ;; LOOP
  for kk=0L,ntran-1 do begin
     ;; Name
     getfnam, trans[kk], f, nam
     lin = nam
     ;; PRINT
     lin = lin+'\\'
     printf, 91, lin
     
     for qq=0,nmodels-1 do begin
        lin = '&'
        lin = lin+all_strct[qq].name
        lin = lin+'&'

        ;; Absorption?
        if all_strct[qq].anly_strct[kk].flg_trans EQ 0 then begin 
           ;; v interval
           lin = lin+'[$'
           lin = lin+strtrim(round(all_strct[qq].anly_strct[kk].vmnx_abs[0]),2)
           lin = lin+','
           lin = lin+strtrim(round(all_strct[qq].anly_strct[kk].vmnx_abs[1]),2)
           lin = lin+'$]'
           lin = lin+'&'
           ;; Intrinsic EW
           if all_strct[qq].anly_strct[kk].W_int GT 0. then $
              lin = lin+string(all_strct[qq].anly_strct[kk].W_int,format='(f5.2)') $
           else lin = lin + '$\dots$'
           lin = lin+'&'
           ;; EW
           lin = lin+string(all_strct[qq].anly_strct[kk].W_abs,format='(f5.2)')
           lin = lin+'&'
           ;; tau
           lin = lin+string(all_strct[qq].anly_strct[kk].tau_peak,format='(f4.2)')
           lin = lin+'&$'
           lin = lin+string(round(all_strct[qq].anly_strct[kk].tau_vel),format='(i5)')
           lin = lin+'$&$'
           lin = lin+string(round(all_strct[qq].anly_strct[kk].vel_tau),format='(i5)')
           lin = lin+'$&'
        endif else lin = lin+'&&&&&&'

        ;; Emission
        ;; v interval
        lin = lin+'[$'
        lin = lin+strtrim(round(all_strct[qq].anly_strct[kk].vmnx_em[0]),2)
        lin = lin+','
        lin = lin+strtrim(round(all_strct[qq].anly_strct[kk].vmnx_em[1]),2)
        lin = lin+'$]'
        lin = lin+'&'
        ;; EW
        lin = lin+'$'
        lin = lin+string(all_strct[qq].anly_strct[kk].W_em,format='(f5.2)')
        lin = lin+'$&'
        ;; Flux
        lin = lin+string(all_strct[qq].anly_strct[kk].flux_peak,format='(f5.2)')
        lin = lin+'&$'
        lin = lin+string(round(all_strct[qq].anly_strct[kk].flux_vel),format='(i5)')
        lin = lin+'$&$'
        lin = lin+string(round(all_strct[qq].anly_strct[kk].vel_flux),format='(i5)')
        lin = lin+'$&'

        ;; PRINT
        lin = lin+'\\'
        printf, 91, lin
     endfor
  endfor

  ;; Write the structure
  save, all_strct, filenam=SAV_FIL

  ;; End of Table
  printf, 91, '\enddata'
  printf, 91, '\tablecomments{{L}isted are the equivalent widths (intrinsic, absorption, and emission), the peak optical depth for the absorption'
  printf, 91, '$\tau_{\rm pk} \equiv -\ln(I_{\rm min})$, the velocity where the optical depth peaks $v_\tau$, the optical depth-weighted velocity centroid '
  printf, 91, '$v_{\bar \tau} \equiv \int dv \, v \ln[I(v)] / \int dv \ln[I(v)]$, the peak flux $f_{\rm pk}$ in emission, the velocity where the flux peaks '
  printf, 91, '$v_f$, and the flux-weighted velocity centroid of the emission line $v_{\bar f}$.}'
  printf, 91, '\end{deluxetable}' 

  close, /all

  ;; All done
  print, 'mktab_summ: All done'
  return

  return
end
  
