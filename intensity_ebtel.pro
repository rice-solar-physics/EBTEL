   pro intensity_ebtel, dem, logtdem, g_array, t_array, int, $
        t, n, length, int_avg
;
;  PURPOSE:
;     Compute coronal intensity as a function of time for an EBTEL simulation
;
;  INPUT:
;     dem(i,j) = array of DEM at time i
;     logtdem(j) = corresponding logT array
;     t = temperature array
;     n = electron density array
;     g_array = array of G(T) values corresponding to t_array 
;       (e.g., from  
;          IDL> gofnt, 'fe_12', 190, 200, t_array, g_array, density=1.e9
;          IDL> g_array = trace_t_resp('195ao', t_array)
;     t_array = array of temperatures corresponding to g_array
;     length = loop halflength
;
;  OUTPUT:
;     int = intensity based on dem 
;     int_avg = intensity based on T and n (the instantaneous coronal averages from EBTEL).
;
;  HISTORY:
;     08-may-01, written, J.A. Klimchuk
;     08-oct-08, JAK, changed name of input variable from dem_cor to dem
;

   nsize = size(dem)
   ntime = nsize(1)
   int = fltarr(ntime)
   int_avg = fltarr(ntime)

   temp = 10.^logtdem
   g = 0.*temp
   tmin = min(t_array)
   tmax = max(t_array)
   ss = where((temp ge tmin) and (temp le tmax))
   g(ss) = 10.^dspline(alog10(t_array), alog10(g_array), logtdem(ss))
   dlogt = logtdem(1) - logtdem(0)

   for i = 0, ntime-1 do begin
      dem_i = dem(i,*)
      int(i) = alog(10.)*dlogt*total(g*temp*dem_i)
   endfor

   g_avg = 0.*t
   ss = where((t ge tmin) and (t le tmax))
   g_avg(ss) = 10.^dspline(alog10(t_array), alog10(g_array), alog10(t(ss)))
   int_avg = n*n*g_avg*2.*length

   return
   end


