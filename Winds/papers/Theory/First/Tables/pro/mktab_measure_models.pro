;; mktab_summ
pro mktab_measure_models, outfil, all_strct=all_strct, INFIL=infil, TITLE=title, $
                          LBL=lbl, SAV_FIL=sav_fil

  if not keyword_set( INFIL ) then infil = 'Input/tab_meas_models.inp'
  if not keyword_set( OUTFIL ) then outfil = 'tab_meas_models.tex'
  if not keyword_set( SAV_FIL ) then sav_fil = 'measures.idl'
  if not keyword_set(ALL_ANG) then all_ang = [0, 120, 150, 180, 30, 60, 90]
  if not keyword_set(ALL_THETA) then all_theta = [10, 45, 60]
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
        flg_int = 0
        ;; Read name
        readf, 1, model_nm
        print, qq, model_nm
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
           else: begin ; Asymmetry data cube (includes biconical)
              ;; MgII
              readf, 1, data_fil
              if flg_anly GE 0 then jj = where(flg_anly EQ all_ang) $
              else jj = where(abs(flg_anly) EQ all_theta)
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
              feII_wave = xmrdfits(data_fil, 0, /silent)
              feII_data = xmrdfits(data_fil, jj[0]+1, /silent)
              feII_data = float(feII_data)
              spec_feii = total(total(feii_data,1),1)
              feII_strct = { $
                           spec: spec_feii, $
                           noscatt_spec: replicate(2., n_elements(spec_feii)), $
                           wave: feII_wave $
                           }
              ;; Intrinsic
              case abs(flg_anly) of 
                 0: flg_int = 0 
                 180: flg_int = 1
                 else: flg_int = 1  ;; Biconical
              endcase
           end
        endcase

;        if strmid(model_nm,0,4) EQ 'Radi' then stop
        ;; Analyze
        cd, '../Analysis/pro/', curr=curr
        RESOLVE_ROUTINE, 'calc_all_values'
        calc_all_values, mgII_strct, feII_strct, strct, TRANS=trans;, DEBUG=(qq EQ 2)
        cd, curr
        if qq EQ 0 then begin
           tmp = { $
                 name: model_nm, $
                 data_fil: data_fil, $
                 anly_strct: strct $
                 }
           all_strct = replicate(tmp, nmodels)
        endif else begin
           ;; Asymmetry kludge
           if flg_int EQ 1 then strct.W_int = all_strct[0].anly_strct.W_int
           ;; Save
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
  printf, 91,  '& \colhead{$v_{\rm int}^e$} & ' + $
    '\colhead{$W_{\rm e}$} & \colhead{$f_{\rm pk}$} & \colhead{$v_f$} '
  printf, 91,  '& \colhead{$v_{\bar f}$} & \colhead{$\Delta v_{\rm e}$}'
  printf, 91,  '\\'
  printf, 91,  '&& (\kms) & (\AA) & (\AA) && (\kms) & (\kms) & (\kms) & (\AA) & & (\kms) ' + $
          '& (\kms) & (\kms)}'
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
           if all_strct[qq].anly_strct[kk].W_int LE 0. then begin
              for ii=0L,5 do begin
                 lin = lin + '$\dots$'
                 lin = lin+'&'
              endfor
           endif else begin
              ;; v interval
              lin = lin+'[$'
              lin = lin+strtrim(round(all_strct[qq].anly_strct[kk].vmnx_abs[0]),2)
              lin = lin+','
              lin = lin+strtrim(round(all_strct[qq].anly_strct[kk].vmnx_abs[1]),2)
              lin = lin+'$]'
              lin = lin+'&'
              ;; Intrinsic EW
              lin = lin+string(all_strct[qq].anly_strct[kk].W_int,format='(f5.2)') 
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
           endelse
        endif else lin = lin+'&&&&&&'

        ;; Emission
        if all_strct[qq].anly_strct[kk].W_em GE 0. then begin
           for ii=0L,5 do begin
              lin = lin + '$\dots$'
              if ii NE 5 then lin = lin+'&'
           endfor
        endif else begin
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
           lin = lin+string(round(all_strct[qq].anly_strct[kk].flux_dv),format='(i4)')
        endelse

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
  printf, 91, '$v_f$, the flux-weighted velocity centroid of the emission line $v_{\bar f}$ (occasionally affected by blends with neighboring emission lines), and the $90\%$ width $\Delta v_{\rm e}$.' 
  printf, 91, 'The $v^a_{\rm int}$ and $v^e_{\rm int}$ columns give the velocity range used to calculate the absorption and emission characteristics, respecitvely.  These were defined by the velcoities where the profile crossed 0.95 in the normalized flux.}'
  printf, 91, '\end{deluxetable}' 

  close, /all

  ;; All done
  print, 'mktab_summ: All done'
  return

  return
end
  
