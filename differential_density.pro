     pro dif_den, logtdem, dem_cor, dem_tr, t, n, p, length, dd_cor, dd_tr, helium=helium
	 
;  Procedure to compute the differential density distribution from EBTEL output

;  Input:
;    t = temperature array corresponding to time (avg. over coronal section)
;    n = electron number density array (cm^-3) (coronal avg.)
;    p = pressure array (dyn cm^-2) (coronal avg.)
;    length = loop half length (top of chromosphere to apex) (cm)
;    dem_tr = differential emission measure of transition region, DEM(time,T), both legs
;              (DEM = n^2 * ds/dT  cm^-5 K^-1)
;    dem_cor = differential emission measure of corona, DEM(time,T), both legs
;              (Int{dem_cor+dem_tr dT} gives the total emission measure of a
;              loop strand having a cross sectional area of 1 cm^2)
;    logtdem = logT array corresponding to dem_tr, dem_cor, dd_tr, and dd_cor
;  Output:
;    dd_tr = differential density distribution of transition region, DD(time,T), both legs
;              (DD = n * ds/dT  cm^-2 K^-1)
;    dd_cor = differential density distribution of corona, DD(time,T), both legs
;              (Int{dd_cor+dd_tr dT} gives the column density (Int{n ds}) of a loop strand)
;  Optional keyword input:
;    helium = set to use a plasma that includes Helium


;  History:
;    18-Jun-2020, writted, J. Klimchuk

     k_b = 1.38e-16
	 mp = 1.67e-24
	
; If Helium is included:
     if keyword_set(helium) then begin
	    print, 'Helium included'
	    n_he_n_p = 0.075;   He/p abundance.
        z_avg = (1 + 2*n_he_n_p)/(1 + n_he_n_p); Include Helium
        kb_fact = 0.5*(1.+1./z_avg)
        k_b = k_b*kb_fact; Modify equation of state for non-e-p plasma
        m_fact = (1 + n_he_n_p*4.)/(2 + 3.*n_he_n_p); Include Helium
        mp = mp*m_fact*(1 + z_avg)/z_avg; Average ion mass
	 endif
	
     t_tr = 10.^logtdem
	 
	 dim = size(dem_cor)
	 ntime = dim(1)
	 ntemp = dim(2)
	 
	 dd_cor = dem_cor*0.
	 dd_tr = dem_tr*0.
	 
	 for i = 0, ntime - 1 do begin
	 
	    dd_cor(i,*) = dem_cor(i,*)/n(i)
		
		sc = (2.*k_b*t(i)/mp)/2.74e4 ;   gravitational scale height in corona
        p2kt = p(i)*exp(2*sin(3.14159/5.)*length/3.14159/sc)/(2.*k_b*t_tr)
;  Original EBTEL (no stratification)
;		p2kt = p(i)/(2.*k_b*t_tr)
		  
		n_tr = p2kt
		dd_tr(i,*) = dem_tr(i,*)/n_tr
	 
	 endfor
	 
	 return
	 end