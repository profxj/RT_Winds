;; mktab_summ
pro mktab_fiducial_ew, outfil 

  if not keyword_set( OUTFIL ) then outfil = 'tab_fiducial_ew.tex'

  ;; Initialize
  compile_opt strictarr

  cd, '../Analysis/pro/', curr=curr
  RESOLVE_ROUTINE, 'calc_ew_values'
  calc_ew_values, '../Outputs/fiducial_grid.fits', strct
  cd, curr
  nlin = n_elements(strct)

  ;; 
  c = x_constants()

  ;; Table header
  close, /all
  openw, 91, outfil

  printf, 91,  ' '
;      printf, 1,  '\clearpage'
  printf, 91,  ' '
  printf, 91,  '\begin{deluxetable}{ccr}'
;  printf, 91,  '\rotate'
  printf, 91,  '\tablewidth{0pc}'
  printf, 91,  '\tablecaption{Equivalent Widths for the ' + $
          'Fiducial Model\label{tab:fiducial_EW}}'
  printf, 91,  '\tabletypesize{\footnotesize}'
  printf, 91,  '\tablehead{\colhead{Transition} & \colhead{$v_{\rm int}^a$} & ' + $
    '\colhead{$W_{\lambda}$} \\ '
  printf, 91,  '& (\kms) & (\AA) }'
  printf, 91, '\startdata'

  ;; LOOP
  for kk=0L,nlin-1 do begin
     ;; Name
     getfnam, strct[kk].wrest, f, nam
     lin = nam
     lin = lin+'&'

     ;; Absorption?
     if strct[kk].flg_trans EQ 0 then begin 
        ;; v interval
        lin = lin+'[$'
        lin = lin+strtrim(round(strct[kk].vmnx_abs[0]),2)
        lin = lin+','
        lin = lin+strtrim(round(strct[kk].vmnx_abs[1]),2)
        lin = lin+'$]'
        lin = lin+'&'

        ;; EW
        lin = lin+string(strct[kk].W_abs,format='(f5.2)')

        ;; PRINT
        lin = lin+'\\'
        printf, 91, lin

        ;; Start new line
        lin = '&'
     endif

     ;; Emission
     ;; v interval
     lin = lin+'[$'
     lin = lin+strtrim(round(strct[kk].vmnx_em[0]),2)
     lin = lin+','
     lin = lin+strtrim(round(strct[kk].vmnx_em[1]),2)
     lin = lin+'$]'
     lin = lin+'&'
     
     ;; EW
     lin = lin+'$'
     lin = lin+string(strct[kk].W_em,format='(f5.2)')
     lin = lin+'$'

     ;; PRINT
     lin = lin+'\\'
     printf, 91, lin
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
  
