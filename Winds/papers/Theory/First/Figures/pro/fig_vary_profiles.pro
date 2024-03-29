pro fig_vary_profiles, RREAL=rreal, STRCT=strct

  if not keyword_set( PSFILE ) then psfile = 'fig_vary_profiles.ps'
  if not keyword_set(CSZ) then csz = 1.4
  if not keyword_set(lSZ) then lsz = 1.4
  if not keyword_set(XNCOLORS) then xncolors=200L

  if not keyword_set(dv) then dv = 1.
  if not keyword_set(NPTS) then npts = 500L                  ; Log steps
  if not keyword_set(b_val) then b_val = 15. ; km/s
  if not keyword_set(DUST) then dust = 0.1  ; Depletion
  if not keyword_set(METAL) then metal = -0.3  ; [M/H]

  if not keyword_set(VLAWS) then vlaws = [-2, -1, 0.5]
  if not keyword_set(VNORM) then vnorm = [2., 50., 100.]
  n_vlaw = n_elements(vlaws)
  if not keyword_set(NLAWS) then nlaws = [-3, 0., 2]
  n_nlaw = n_elements(nlaws)
  nmodel = n_vlaw*n_nlaw + 1
  if not keyword_set(N_0) then n_0 = [0.4, 0.5, 0.3, 0.01, 0.01, 0.02, 0.01, 1e-3, 0.0001, 0.1]
  lbl = ['A','B','C','D','E','F','G','H','I','J','K']
  if n_elements(N_0) NE nmodel then stop
  ;; Blue-dotted, blue-dashed, blue-dot-dash
  ;; Red-dotted, red-dashed, red-dot-dash
  ;; Green-dotted, green-dashed, green-dot-dash
  ;; Black

  c = x_constants()

  ;;;;;;;;;;;;;;;;;;
  ;; Radius
  r0 = 1.  ; kpc
  r1 = 20. ; kpc
  rval = r0 * exp( alog(r1/r0) * findgen(npts)/float(npts-1) )
  dr = rval - shift(rval,1)  ; kpc
  dr[0] = dr[1] 

  ;;;;;;;;;;;;;;;;;;;;
  ;; Spectrum
  npix = 2000L
  wrest = 2796.352d
  wav = 10.^(alog10(2795.) + dindgen(npix)*1.449d-6)
  vel = (wav-wrest)/wrest * 3e5
  dvel = median(vel-shift(vel,1))  ; Should be 1 km/s

  ;; Plot
  close, /all
  x_psopen, psfile, /portrait
  !p.multi = [0,1,2]
  clr = getcolor(/load)
  xclrs = x_setclrs(/dark)
  xmrg = [11,2]
  ymrg = [4.0,1]

  openw, 11, 'plaws_parms.dat'

  for qq=0L,nmodel-1 do begin
     if qq EQ (nmodel-1) then begin
        ;; FIDUCIAL ;;;;;;;;;;;;;;;;
        ;; Velocity
        v_0 = 50. ; km/s
        v_r = v_0 * rval / r0

        ;; Density
        n_r = n_0[qq] * (rval/r0)^(-2)
        idx_v = -1
        idx_n = -1
     endif else begin
        ;; Velocity
        idx_v = qq MOD 3
        if vlaws[idx_v] LT 0 then v_r = vnorm[idx_v] * (rval / r1)^vlaws[idx_v] $
        else v_r = vnorm[idx_v] * (rval / r0)^vlaws[idx_v]
        ;; Density
        idx_n = qq / 3
        n_r = n_0[qq] * (rval/r0)^nlaws[idx_n]
     endelse
     if (qq NE (nmodel-1)) then begin
        print, 'Density: ', idx_n, nlaws[idx_n], n_0[qq]
        print, 'Velocity: ', idx_v, vlaws[idx_v], vnorm[idx_v]
     endif
        
     ;; Optical depth
     mgii = x_setline(wrest)
     lines = replicate(mgii, npts)
     lines.b = b_val
     lines.N = alog10(dr * c.kpc * n_r * 10.^(7.53-12.+METAL) * DUST)
     lines.zabs = v_r/3e5

     fx = x_voigt(wav, lines, /nosmooth, TAU=tau)
;     if qq EQ (nmodel-1) then stop

     ;; Table
     if qq NE (nmodel-1) then $
        printf, 11, lbl[qq], ' & ', '$r^{'+string(nlaws[idx_n],format='(i2)')+'}$& '+$
               string(n_0[qq],format='(f6.4)')+' & '+$
               '$r^{'+string(vlaws[idx_v],format='(f4.1)')+'}$& '+$
               string(vnorm[idx_v],format='(f5.1)')+' & '+$
               string(alog10(max(tau)),format='(f3.1)')+' \\ '

     ;; Density and Velocity
     !p.multi = [2,1,2]
     yrng=[1e-3, 1e3]
     xrng=[r0, r1]
     plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
           xmargin=xmrg, ymargin=ymrg, ytickformat='x_logticks', $
           ytitle='n!dH!u!N [x10!u2!N, cm!u-3!N];   v!dr!N [x10!u-2!N km s!u-1!N]', $
           xtitle='Radius (kpc)', yrange=yrng, thick=4, $
           xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, /ylog, /noerase

     if qq EQ 0 then xyouts, xrng[0]+0.2*(xrng[1]-xrng[0]), $
                             yrng[1]/4., $
                             '(a)', charsi=lsz, color=clr.black

     ;; Density
     if idx_v EQ 0 then oplot, rval, n_r*100., color=clr.gray, linesty=idx_n+1
     if idx_v EQ -1 then oplot, rval, n_r*100, color=clr.black, thick=2

     ;; Velocity
     if idx_n EQ 0 then oplot, rval, v_r/100, color=xclrs[idx_v+1]
     if idx_n EQ -1 then oplot, rval, v_r/100., color=clr.black, thick=2, linesty=4

     ;; tau
     !p.multi = [1,1,2]
     yrng=[0.01, 3e3]
     xrng=[0., 1000.]
     plot, [0], [0], color=clr.black, background=clr.white, charsize=csz,$
           ytitle='!9t!X!d2796!N', ytickformat='x_logticks', $
           xtitle='v!N (km s!u-1!N)', yrange=yrng, thick=4, $
           xmargin=xmrg, ymargin=ymrg, $
           xrange=xrng, ystyle=1, xstyle=1, psym=1, /nodata, /noerase, /ylog

     if qq EQ 0 then xyouts, xrng[0]+0.2*(xrng[1]-xrng[0]), $
                             yrng[1]/4., $
                             '(b)', charsi=lsz, color=clr.black

     if idx_n EQ -1 then thk = 2 else thk = 4
     ;; Tau
     oplot, vel, tau, color=xclrs[idx_v+1], linesty=idx_n+1, thick=thk

     ;; Label
     ylbl = replicate(yrng[1]/(2^(idx_n+1)),2)
     oplot, 400.+idx_v*200+[30., 110], 1.1*[ylbl,ylbl],$
            color=xclrs[idx_v+1], linesty=idx_n+1
     if qq LT (nmodel-1) then $
         xyouts, 400 + idx_v*200 + 130., ylbl, lbl[qq], color=xclrs[idx_v+1], charsiz=lsz,$
                 align=0.
  endfor
     
  if keyword_set( PSFILE ) then x_psclose
  !p.multi = [0,1,1]
  close, /all

  return

end
