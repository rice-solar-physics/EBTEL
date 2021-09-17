pro calc_c1, temp, den, length, rad, c1
    common params, k_b, mp, kappa_0
    common area, a0_ac, atr_ac, l_c, l_tr, l_star

    ; Calculate scale height
    calc_lambda, temp, sc

    r3_rad_0 = 0.6  ; Radiative phase value, no gravity
    r3_eqm_0 = 2.  ; Value in eqm with no gravity, -2/3 loss power law.
    l_fact_eq  = 5.
    l_fact_rad = 5.

    calc_c2, r2

    f_eq_1 = -den*den*rad*length;

    ; Adjust values for gravity
    ;
    r3_eqm_g =   r3_eqm_0*exp(2*2*sin(3.14159/l_fact_eq)*l_c/3.14159/sc)
    r3_radn_g = r3_rad_0*exp(2*2*sin(3.14159/l_fact_rad)*l_c/3.14159/sc)
    ;
    ; Adjust for loss function
    ;
    r3_eqm  =  r3_eqm_g*1.95e-18/temp^(2./3.)/rad
    r3_radn = r3_radn_g*1.95e-18/temp^(2./3.)/rad
    ; These lines turn off loss function fix
;   r3_eqm = r3_eqm_g
;   r3_radn = r3_radn_g
    ;
    ;Calculate over/under density
    ;    
    n_eq_2 = (a0_ac/atr_ac)*kappa_0/3.5/r3_eqm/rad*(l_star/length/l_c^2)*(temp/r2)^(7./2.)/(1.-l_tr/l_c/r3_eqm)
        
    noneq2=den^2./n_eq_2
    
    r3_cond = 6.
    ;     r3_cond = 2.  ; uncomment this to turn off the conductive phase fix
    if (noneq2 lt 1.) then begin
        r3=(2.*r3_eqm+r3_cond*(1./noneq2-1.))/(1+1./noneq2)
    endif else begin
        ; Radiative loops: transition from eqm.
        r3=(2.*r3_eqm+r3_radn*(noneq2-1.))/(1+noneq2)
    endelse

    ; Final value
    c1=r3

    ; Short cut to use constant C1
;   c1=2. ; EBTEL-2 value
;   c1=4. ; Paper 1 value
    return
end