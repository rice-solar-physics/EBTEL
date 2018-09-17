pro calc_lambda, temp, sc
    common params, k_b, mp, kappa_0
    sc = (2.*k_b*temp/mp)/2.74e4
;   sc = 1e15 ;Include to kill gravity
    return
end