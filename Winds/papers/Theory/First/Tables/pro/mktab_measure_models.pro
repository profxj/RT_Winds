;; mktab_summ
pro mktab_measure_models, outfil, all_strct=all_strct, INFIL=infil

  if not keyword_set( INFIL ) then infil = 'Input/tab_meas_models.inp'
  if not keyword_set( OUTFIL ) then outfil = 'tab_meas_models.tex'
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
           1: begin
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
           2: begin
              ;; MgII
              readf, 1, data_fil
              readcol, data_fil, wv, fx, /silen
              nrm = median(fx[where(wv GT 2815)])
              fx = fx/nrm
              mgII_strct = { $
                           spec: fx, $
                           wave: wv $
                           }
              ;; FeII
              readf, 1, data_fil
              readcol, data_fil, wv, fx, /silen
              nrm = median(fx[where(wv GT 2605)])
              fx = fx/nrm
              feII_strct = { $
                           spec: fx, $
                           wave: wv $
                           }
           end
           else: begin ; Asymmetry
              ;; MgII
              readf, 1, data_fil
              jj = where(flg_anly EQ all_ang)
              mgII_wave = xmrdfits(data_fil, 0, /silent)
              mgII_data = xmrdfits(data_fil, jj[0]+1, /silent)
              mgII_data = float(mgII_data)
              spec_mgii = total(total(mgii_data,1),1)
              mgII_strct = { $
                           spec: spec_mgii, $
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
                           wave: feII_wave $
                           }
           end
        endcase

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
  printf, 91,  '\begin{deluxetable}{ccrcccccccccc}'
  printf, 91,  '\rotate'
  printf, 91,  '\tablewidth{0pc}'
  printf, 91,  '\tablecaption{Line Diagnostics \label{tab:line_diag}}'
  printf, 91,  '\tabletypesize{\footnotesize}'
  printf, 91,  '\tablehead{\colhead{Transition} & \colhead{Model} & \colhead{$v_{\rm int}^a$} & ' + $
    '\colhead{$W_{\lambda}$} & \colhead{$\tau_{\rm pk}$} & \colhead{$\tau_{\rm vel}$} '
  printf, 91,  '& \colhead{$<v(\tau)>$}'
  printf, 91,  '& \colhead{$v_{\rm int}^b$} & ' + $
    '\colhead{$W_{\lambda}$} & \colhead{$f_{\rm pk}$} & \colhead{$f_{\rm vel}$} '
  printf, 91,  '& \colhead{$<v(f)>$}'
  printf, 91,  '\\'
  printf, 91,  '&& (\kms) & (\AA) }'
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
        endif else lin = lin+'&&&&&'

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

  ;; End of Table
  printf, 91, '\enddata'
;  printf, 91, '\tablecomments{Unless specified otherwise, all quantities refer to the \sna=2 threshold.  The cosmology assumed has $\Omega_\Lambda = 0.7, ' + $
;          '\Omega_m = 0.3$, and $H_0 = 72 \mkms \rm Mpc^{-1}$.}'
  printf, 91, '\tablenotetext{a}{Velocity interval over which the equivalent width is calculated.}'
  printf, 91, '\end{deluxetable}'

  close, /all

  ;; All done
  print, 'mktab_summ: All done'
  return

  return
end
  
