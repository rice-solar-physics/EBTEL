   pro ebtel2, ttime, heat, length, t, n, p, v, ta, na, pa, c11, dem_tr, dem_cor,  $
       logtdem, f_ratio, rad_ratio, cond, rad_cor,$
       classical=classical, dynamic=dynamic, dem_old=dem_old, $
       flux_nt=flux_nt, energy_nt=energy_nt, rtv=rtv
     ;
     ; NAME:  Enthalpy-Based Thermal Evolution of Loops (EBTEL)
     ;
     ; PURPOSE:
     ;   Compute the evolution of spatially-averaged (along the field) loop quantities
     ;   using simplified equations.  The instantaneous differential emission measure of
     ;   the transition region is also computed. This version incorporates all the modifications
     ;   from Cargill et al (2012) and is written in modular form for clarity.
     ;   DEM parts unchanged except for TR pressure correction (see
     ;   Paper 2). Original list of JK modifications removed from preamble.
     ;
     ; INPUTS:
     ;   ttime  = time array (s)
     ;   heat   = heating rate array (erg cm^-3 s^-1)   (direct heating only)
     ;              (the first element, heat(0), determines the initial static equilibrium)
     ;   length = loop half length (top of chromosphere to apex) (cm)
     ;
     ; OPTIONAL KEYWORD INPUTS:
     ;   classical = set to use the UNsaturated classical heat flux
     ;   dynamic   = set to use dynamical r1 and r2 (NOT recommended, especially when T > 10 MK). Now redundant.
     ;   dem_old   = set to use old technique of computing DEM(T) in the trans. reg.
     ;               (weighted average of demev, demcon, and demeq)
     ;   flux_nt   = energy flux array for nonthermal electrons impinging the chromosphere
     ;               (input as a positive quantity; erg cm^-2 s^-1)
     ;   energy_nt = mean energy of the nonthermal electrons in keV (default is 50 keV)
     ;   rtv       = set to use Rosner, Tucker, Vaiana radiative loss function (Winebarger form)
     ;
     ; OUTPUTS:
     ;   t (t_a) = temperature array corresponding to time (avg. over coronal section of loop / apex)
     ;   n (n_a) = electron number density array (cm^-3) (coronal avg. / apex)
     ;   p (p_a) = pressure array (dyn cm^-2) (coronal avg. /apex)
     ;   v = velocity array (cm s^-1) (r4 * velocity at base of corona)
     ;   c11 = C1 (or r3 in this code)
     ;   dem_tr = differential emission measure of transition region, dem(time,T), both legs
     ;             (dem = n^2 * ds/dT  cm^-5 K^-1)
     ;             (Note:  dem_tr is not reliable when a nonthermal electron flux is used.)
     ;   dem_cor = differential emission measure of corona, dem(time,T), both legs
     ;             (Int{dem_cor+dem_tr dT} gives the total emission measure of a
     ;             loop strand having a cross sectional area of 1 cm^2)
     ;   logtdem = logT array corresponding to dem_tr and dem_cor
     ;   f_ratio = ratio of heat flux to equilibrium heat flux
     ;             (ratio of heat flux to tr. reg. radiative loss rate)
     ;   rad_ratio = ratio of tr. reg. radiative loss rate from dem_tr and from r3*(coronal rate)
     ;   cond = conductive loss from corona
     ;   rad_cor =  coronal radiative loss
     ;
     ; CORRESPONDENCE WITH VARIABLES IN ASTROPHYSICAL JOURNAL ARTICLES:
     ;          (Klimchuk et al., 2008, ApJ, 682, 1351; Cargill et al., ApJ 2012)
     ;   r1 = c_3
     ;   r2 = c_2
     ;   r3 = c_1
     ;   f, ff = F_0
     ;   f_eq, ff_eq = - R_tr
     ;   dem_eq = DEM_se
     ;
     ; USAGE:
     ;   To include the transition region DEM:
     ;      IDL> ebtel2, ttime, heat, t, n, p, v, dem_tr, dem_cor, logtdem
     ;   To exclude the transition region DEM (faster):
     ;      IDL> ebtel2, ttime, heat, t, n, p, v
     ;   To include a nonthermal electron energy flux:
     ;      IDL> ebtel2, ttime, heat, t, n, p, v, flux_nt=flux_nt
     ;   To compute rad_ratio:
     ;      IDL> ebtel2, ttime, heat, t, n, p, v, dem_tr, dem_cor, logtdem, f_ratio, rad_ratio
     ;      (Takes 25% more computing time.)
     ;
     ; INTENSITIES:
     ;   For observations in which temperature response function, G(T), has units of
     ;   DN s^-1 pix^-1 cm^5 and the loop diameter, d, is larger than the pixel dimension, l_pix:
     ;
     ;      I_cor_perp = d/(2L)* Int{G(T)*dem_cor(T)*dT}
     ;      I_tr_perp = d/l_pix * Int{G(T)*dem_tr(T)*dT}
     ;      I_tr_parallel = Int{G(T)*dem_tr(T)*dT} ,
     ;
     ;   for lines-of-sight perpendicular and parallel to the loop axis.  I_tr_perp assumes that
     ;   the transition region is thinner than l_pix.
     ;
     ; MISCELLANEOUS COMMENTS:
     ;   A 1 sec time is generally adequate, but in applications where exceptionally strong conductive
     ;      cooling is expected (e.g., intense short duration heating events, especially in short loops),
     ;      a shorter time step may be necessary. If there is a question, users should compare runs with
     ;      different time steps and verify that there are no significant differences in the results.
     ;   Runs much more quickly if the transition region DEM is not computed.
     ;   Speed can be increased by increasing the minimum DEM temperature from 10^4 to, say, 10^5 K
     ;      or by decreasing the maximum DEM temperature from 10^8.5 to, say, 10^7.5
     ;      (search on 450 and 451).
     ;   The equilibrium base heat flux coefficient of 2/7 is appropriate for uniform heating;
     ;      a coefficient of 4/7 is more appropriate for apex heating.
     ;   To have equal amounts of thermal and nonthermal heating:  flux_nt = heat*length.
     ;   It is desirable to have a low-level background heating during the cooling phase so that the
     ;      coronal temperature does not drop below values at which the corona DEM is invalid.
     ;   r1 = c_3 = 0.7 gives more accurate coronal evolution than the original 0.5, especially in
     ;      the late phase of cooling.  However, it produces excess DEM at the very hottest temperatures
     ;      during impulsive events, and the transition region DEM is somewhat elevated.  We have
     ;      therefore introduced r1_tr = 0.5, which provides a more accurate transition region DEM at
     ;      the same time that r1 = 0.7 provides a more accurate radiative cooling.
     ;   v = (c_2/c_3)*(t_tr/t)*v_0 = (r2/r1)*(t_tr/t)*(v/r4) at temperature t_tr in the transition
     ;      region, where t is the average coronal temperature.
     ;

     ; HISTORY:
     ; May 2012. PC version. Modular form.
     ; See original ebtel.pro for many additional comments.
     ; 2013 Jan 15, JAK, Fixed a bug in the compution of the radiation function at 10^4 K;
     ;      important for computing radiation losses based on dem_tr;
     ;      ge vs. gt in computing rad;  lt vs. le in computing rad_dem
     ; ---------------------------------------------------------------------
     common params, k_b, mp, kappa_0

     ntot = n_elements(ttime)

     ; Physical constants Can comment out Hydrad lines if needed.
     k_b = 1.38e-16
     mp = 1.67e-24
     n_he_n_p = 0.075;   He/p abundance.
     z_avg = (1 + 2*n_he_n_p)/(1 + n_he_n_p); Include Helium
     z_avg = 1.; For Hydrad comparison.
     kb_fact = 0.5*(1.+1./z_avg)
     k_b = k_b*kb_fact; Modify equation of state for non-e-p plasma
     m_fact = (1 + n_he_n_p*4.)/(2 + 3.*n_he_n_p); Include Helium
     m_fact = (1 + n_he_n_p*4.)/2.; For Hydrad comparison
     mp = mp*m_fact*(1 + z_avg)/z_avg; Average ion mass

     kappa_0 = 1.e-6
     kappa_0 = 8.12e-7; Hydrad value

     ; Set ksX in loss function (T-breaks).

     radloss, rad, 1.e6, 0, rtv=rtv

     ; Calculate initial Cs (rs in this code).

     calc_c3, r1
     r1=0.6

     ; Ratio of temperature at top of transition region to apex temperature,
     ;      used only for calculation of transition region DEM, dem_tr
     ;      See miscellaneous comments above.
     ; r1_tr = 0.5
     r1_tr = r1

     ; Ratio of average to apex temperature (c_2 in ApJ paper)

     calc_c2, r2
     r2 = 0.9

     ; Ratio of conductive to radiative losses in equilibrium (c_1 in ApJ paper)
     ; Initial value of C1. Fix value now then later iteration on initial state sets it properly

     r3 = 2

     ; Ratio of average to base velocity
     r4 = 1.0

     ; Invalid heat array
     if (heat(0) eq 0.) then begin
       print, '* * * * No initial loop heating * * * *'
       goto, jump99
     endif

     q = heat
     t = fltarr(ntot)
     n = fltarr(ntot)
     ta = fltarr(ntot)
     na = fltarr(ntot)
     p = fltarr(ntot)
     pa = fltarr(ntot)
     v = fltarr(ntot)
     c11 = fltarr(ntot)
     cond = fltarr(ntot)
     rad_cor = fltarr(ntot)
     dz = fltarr(ntot)

     ; Set up thermal conduction parameters
     c1 = -2./7.*kappa_0  ;

     ;  conduction coefficient
     c_sat = -1.5*k_b^1.5/9.1e-28^0.5
     sat_limit = 1./6.

     ; Set up nonthermal electrons
     if not keyword_set(flux_nt) then flux_nt = ttime*0.
     flux_nt = -flux_nt

     if not keyword_set(energy_nt) then energy_nt = 50.   ; 50 keV
     energy_nt = 1.602e-9*energy_nt
     j_nt = flux_nt/energy_nt

     dlogt_cor = -alog10(r2)
     dj = fix(100*dlogt_cor)
     nj = 2*dj + 1

     ; Set up DEM in transition region
     if n_params() gt 12 then begin

       logtdem = findgen(451)/100. + 4.
       tdem = 10.^logtdem
       root_tdem = tdem^0.5
       fourth_tdem = tdem^0.25
       dem_tr = fltarr(ntot,451)
       dem_cor = fltarr(ntot,451)
       rad_dem = fltarr(451)
       demev = fltarr(ntot,451)
       demcon = demev
       demeq = demev
       rad_ratio = fltarr(ntot)

       f_ratio = fltarr(ntot)
       root_c2 = (kappa_0/(20.*k_b))^0.5/k_b    ; root_c2 to avoid overflow in dem_ev
       c3 = -5.*k_b
       c4 = (kappa_0/14.)^0.5/k_b

       ;    Radiation in transition region

       for i = 0, 450 do begin
         radloss, rad,tdem(i),1, rtv=rtv
         rad_dem(i)=rad
;         if tdem(i) le 1.e4 then begin
         if tdem(i) lt 1.e4 then begin
           rad_dem(i) = 1.
         endif
       endfor
       root_rad_dem = rad_dem^0.5

     endif

     ; ---------------------------
     ; Initial static equilibrium
     ; 2 methods. (a) Use EBTEL eqm. (b) Use scalings laws. (a) recommended.
     ; ---------------------------

     ; Set up trial values for C1 = 2

     tt_old = r2*(3.5*r3/(1. + r3)*length*length*q(0)/kappa_0)^(2./7.)
     radloss, rad, tt_old, 1, rtv=rtv
     nn = (q(0)/((1. + r3)*rad))^0.5
     nn_old = nn

     ; Iterate on TT and r3

     for i=1,100 do begin
       calc_c1, tt_old, nn, length, rad, r3
       tt_new = r2*(3.5*r3/(1. + r3)*length*length*q(0)/kappa_0)^(2./7.)
       radloss, rad, tt_new, 1, rtv=rtv
       nn = (q(0)/((1. + r3)*rad))^0.5
       err = tt_new - tt_old
       err_n = nn_old - nn
       if abs(err) lt 1e3 then  begin
         i=100
       endif
       print,r3,tt_new,tt_old,err_n,i
       tt_old = tt_new
       nn_old = nn
     endfor

     tt = tt_old
     nn = (q(0)/((1. + r3)*rad))^0.5

     ;   If want to fix out of eqm start, e.g. cooling flare. Section 4.2, Paper 3.
;        tt=1.e7*r2
;        nn = 1e9/r2

     print, ' '
     print, 'Model parameters'
     print, '  r1 = ', r1
     print, '  r2 = ', r2
     print, '  R3 = ', r3
     print, ' '

     t(0) = tt
     n(0) = nn
     p(0) = 2.*k_b*n(0)*t(0)
     v(0) = 0.
     ta(0) = t(0)/r2
     calc_lambda, t(0), sc
     na(0) = n(0)*r2*exp(-2*length*(1.-sin(3.14159/5.))/3.14159/sc);
     pa(0) = 2*k_b*na(0)*ta(0)

     print, 'Initial static equilibrium'
     print, '  L = ', length
     print, '  Q = ', q(0)
     print, '  T, T_a = ', tt,ta(0)
     print, '  n, n_a = ', nn, na(0)

     ; Scaling law values

     lambda_0 = 1.78e-19
     lambda_0 = 1.95e-18                           ;  lambda = lambda_0*T^bb
     bb = -0.5
     bb = -2./3.
     q_0 = heat(0)                                     ;  heating rate (erg cm^-3 s^-1)
     t_0 = r2*(3.5/kappa_0*q_0)^(2./7.)*length^(4./7.) ;  temperature (K) (avg., not apex)
     p_0 = r2^(-7./4.)*(8./7.*kappa_0/lambda_0)^0.5*k_b       $
       *t_0^((11.-2.*bb)/4.)/length                ;  total pressure (dyn cm^-2)
     n_0 = 0.5*p_0/(t_0*k_b)                           ;  electron number density (cm^-3)
     v_0 = 0.                                          ;  velocity

     print, ' '
     print, 'Scaling law values'
     print, '  T = ', t_0
     print, '  P = ', p_0
     print, '  n = ', n_0

     ; ----------------------
     ; Time-dependent heating
     ; ----------------------

     h_tot=0.

     ; Loop over time steps

     for i = 0, ntot-2 do begin
       dt = ttime(i+1) - ttime(i)

       ; Thermal conduction flux at base

       f_cl = c1*(t(i)/r2)^3.5/length

       if keyword_set(classical) then begin
         f = f_cl
       endif else begin
         f_sat = sat_limit*c_sat*n(i)*t(i)^1.5
         f = -f_cl*f_sat/(f_cl*f_cl + f_sat*f_sat)^0.5
       endelse

       radloss, rad, t(i), 1, rtv=rtv

       ; Evaluate C1 - C3

       calc_c1, t(i), n(i), length, rad, r3
       calc_c2, r2
       calc_c3, r1

       r12 = r1/r2
       r12_tr = r1_tr/r2

       c11(i) = r3; Store for output

       ; Equilibrium thermal conduction flux at base (-R_tr in ApJ paper)
       f_eq = -r3*n(i)*n(i)*rad*length

       ;      pv = 0.4*(f_eq - f)
       pv = 0.4*(f_eq - f - flux_nt(i))
       ;      dn = pv*0.5/(r12*k_b*t(i)*length)*dt
       dn = (pv*0.5/(r12*k_b*t(i)*length) + j_nt(i)/length)*dt

       n(i+1) = n(i) + dn

       ;      dp = 2./3.*(q(i) + (1. + 1./r3)*f_eq/length)*dt
       dp = 2./3.*(q(i) + (1. + 1./r3)*f_eq/length    $
         - (1. - 1.5*k_b*t(i)/energy_nt)*flux_nt(i)/length)*dt
       p(i+1) = p(i) + dp

       t(i+1) = p(i+1)/(n(i+1)*2.*k_b)

       v(i+1) = pv/p(i+1)

       h_tot = h_tot+heat(i)

       ; Calculate scale height
       calc_lambda, t(i+1), sc

       ; calculate apex quantities
       ;
       ta(i+1) = t(i+1)/r2;
       na(i+1) = n(i+1)*r2*exp(-2.*length*(1.-sin(3.14159/5.))/3.14159/sc);
       pa(i+1) = 2*k_b*na(i+1)*ta(i+1)

       ; Differential emission measure
       if n_params() gt 12 then begin

         ;   Transition Region
         if (r12_tr*t(i) gt tdem(450)) then begin
           print, 'Transition region temperature outside DEM range'
           return
         endif

         if (f ne f_eq) then        $
           cf = f*f_eq/(f - f_eq)   $
         else                       $
           cf = 1.e10*f

         for j = 0, 450 do begin
           if (tdem(j) lt r12_tr*t(i)) then begin
             ;
             if not keyword_set(dem_old) then begin

               ;    New method
               aaa = kappa_0*tdem(j)^1.5
               bbb = -5.*k_b*n(i)*v(i)
               p2kt2 = (p(i)/(2.*k_b*tdem(j)))^2
               ; Next line has TR fix
               p2kt2 = (p(i)*exp(2*sin(3.14159/5.)*length/3.14159/sc)/(2.*k_b*tdem(j)))^2
               ccc = -p2kt2*rad_dem(j)
               dtds1 = (-bbb + (bbb*bbb - 4.*aaa*ccc)^0.5)/(2.*aaa)
               dtds2 = (-bbb - (bbb*bbb - 4.*aaa*ccc)^0.5)/(2.*aaa)
               dtds = max(dtds1, dtds2)
               dem_tr(i,j) = 2.*p2kt2/dtds  ; factor 2 for both legs

             endif else begin

               ;    Old method
               ;           approximation to tr. reg. dem when evaporation dominates
               dem_ev = (root_c2/n(i))*(root_c2*p(i)*p(i)  $
                 /root_tdem(j))/v(i)
               ;               vabs = abs(v(i))
               ;               vmin = max([vabs, abs(f_eq/(2.5*pv))])
               ;               dem_ev = (root_c2/n(i))*(root_c2*p(i)*p(i)  $
               ;                        /root_tdem(j))/vmin*v(i)/vabs

               ;           approximation to tr. reg. dem when condensation dominates
               dem_con = c3*n(i)*v(i)/rad_dem(j)

               ;           approximation to tr. reg. dem under equilibrium conditions
               dem_eq = c4*p(i)/(root_rad_dem(j)*fourth_tdem(j))

               dem_tr(i,j) = 2.*(f*dem_ev + cf*dem_eq - f_eq*dem_con)   $
                 /(f + cf - f_eq)              ; factor 2 for both legs

               demev(i,j) = 2.*dem_ev
               demcon(i,j) = 2.*dem_con
               demeq(i,j) = 2.*dem_eq
             endelse

             f_ratio(i) = f/f_eq
             if keyword_set(classical) then begin
               f_ratio(i) = f/f_eq
             endif else begin
               f_ratio(i) = f_cl/f_sat
             endelse

           endif

         endfor

         ss = where(dem_tr(i,*) lt 0., nneg)
         if (nneg gt 0) then begin
           print, ' '
           print, '***** Negative DEM;   i = ', i
           print, ' '
           if (i ne 0) then begin
             for j = 0, 450 do dem_tr(i,j) = dem_tr(i-1,j)
           endif else begin
             dem_tr(0,*) = 0.  ; to avoid problems at start of run with saturated heat flux
           endelse
         endif

         ;   Corona   (EM distributed uniformly over temperture interval [t_min, t_max])
         t_max = max([t(i)/r2, 1.1e4])
         t_min = max([t(i)*(2. - 1/r2), 1.e4])
         j_max = fix((alog10(t_max) - 4.0)*100)
         j_min = fix((alog10(t_min) - 4.0)*100)

         em = 2.*n(i)*n(i)*length            ; factor of 2 for both legs

         ;         dem0 = em/(t_max - t_min)
         delta_t = 10.^4*(10.^((j_max+0.5)/100.)    $
           - 10.^((j_min-0.5)/100.))
         dem0 = em/delta_t

         for j = j_min, j_max do   $
           dem_cor(i,j) = dem0

         ;   Transition region radiation losses based on DEM
         if n_params() gt 11 then begin
           rad_loss = 0.

           for j = 0, 450 do begin
             if (tdem(j) lt r12_tr*t(i)) then $
               rad_loss = rad_loss + dem_tr(i,j)*rad_dem(j)*tdem(j)*0.01*2.3  ; 2.3=ln(10)

           endfor

           rad_ratio(i) = -rad_loss/f_eq
         endif

       endif

       cond(i)=f
       rad_cor(i)=f_eq/r3

     endfor

     v = r4*v

     if n_params() gt 12 then begin
       dem_tr(ntot-1,*) = dem_tr(ntot-2,*)
       dem_cor(ntot-1,*) = dem_cor(ntot-2,*)
     endif
     
     jump99:

     return
   end
   