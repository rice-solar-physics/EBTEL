pro radloss, rad, tt, ij, rtv=rtv

     ; Rad loss routine
     ; Parameters: Losses, temp, flags
     ; ij = 0 sets the Ks. ij=1 skips this bit. Ks in common block
     ;
     common ks, kt0,kt1,kt2,kt3,kt4,kt5,kt6

     ; Set the Ks

     if ij eq 0 then begin
       if not keyword_set(rtv) then begin
         ;    Raymond-Klimchuk loss function
         print, ' Klimchuk losses'
         kt0 = 1.e4
         kt1 = 9.3325e4
         kt2 = 4.67735e5
         kt3 = 1.51356e6
         kt4 = 3.54813e6
         kt5 = 7.94328e6
         kt6 = 4.28048e7

       endif else begin
         ;    RTV loss function
         print, ' RTV losses'
         kt0 = 10.^4.3
         kt1 = 10.^4.6
         kt2 = 10.^4.9
         kt3 = 10.^5.4
         kt4 = 10.^5.75
         kt5 = 10.^6.3
       endelse
     end

     if not keyword_set(rtv) then begin
       ;    Raymond-Klimchuk loss function

         if (tt gt kt6) then rad = 1.96e-27*sqrt(tt) else $
         if (tt gt kt5) then rad = 5.4883e-16/tt else $
         if (tt gt kt4) then rad = 3.4629e-25*(tt^(0.3333333)) else $
         if (tt gt kt3) then rad = 3.5300e-13/(tt^(1.5)) else $
         if (tt gt kt2) then rad = 1.8957e-22 else $
         if (tt gt kt1) then rad = 8.8669e-17/tt else $
;         if (tt gt kt0) then rad = 1.0909e-31*tt*tt $
         if (tt ge kt0) then rad = 1.0909e-31*tt*tt $
       else                rad = 0.0

     endif else begin
       ;    RTV loss function

         if (tt gt kt5) then rad = 10^(-17.73)/tt^(.666) else $
         if (tt gt kt4) then rad = 10.^(-21.94) else $
         if (tt gt kt3) then rad = 10.^(-10.4)/tt^2. else $
         if (tt gt kt2) then rad = 10.^(-21.2) else $
         if (tt gt kt1) then rad = 10^(-31.0)*tt^2. else $
;         if (tt gt kt0) then rad = 10^(-21.85) $
         if (tt ge kt0) then rad = 10^(-21.85) $
       else                rad = 0.0
     endelse
 ;      One part fit if needed
 ;      rad = 1.95e-18/tt^(2./3.)
     return
   end