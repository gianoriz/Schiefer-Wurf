program RungeKutta

  implicit none   !RungeKutta 4ter Ordnung

  !Definition der verwendeten Variablen:
  real :: y, v_y, yneu, v_yneu   !y-Koord. [m/s] & Geschw. in y-Richtung [m]
  real :: x, v_x, xneu, v_xneu   !x-Koord. [m/s] & Geschw. in x-Richtung [m] 
  real :: z, v_z, zneu, v_zneu   !z-Koord. [m/s] & Geschw. in z-Richtung [m]
  real :: K11, K12, K13, K14, K15, K16
  real :: K21, K22, K23, K24, K25, K26
  real :: K31, K32, K33, K34, K35, K36
  real :: K41, K42, K43, K44, K45, K46
  real :: v_0                    !Anfangsgeschwindig. [m/s]
  real :: t, dt                  !Flugzeit [s] & Zeitintervall [s]
  real :: g = 9.81               !Erdbeschleunigung   [m/s²]
  real :: alpha                  !Abwurfwinkel        [rad]    
  integer :: i, m                !Laufvariablen
  real :: Pi = 3.1415926535898

  !Setzte Parameter:
  alpha = 45.0 * (Pi/180.0)
  v_0 = 20.0
  dt = 0.2
  m = 14

  !Startwerte: 
  t = 0.0 
  x = 0.0 
  y = 0.0
  write(*,*) t, y                !Ausgeben des Anfangswertes y=0 bei t=0    

  !Geschw. der Komponenten:
  v_x = v_0 * cos(alpha)         !bleibt waehrend der Flugzeit konstant
  v_y = v_0 * sin(alpha)         !aendert sich waehrend der Flugzeit  

  do i = 0, m   

     K11 = v_x * dt
     K12 = v_y * dt 
     K13 = v_z * dt 
     K14 = 0
     K15 = - g * dt 
     K16 = 0

     K21 = v_x * dt + 0.5 * (K14) * dt
     K22 = v_y * dt + 0.5 * (K15) * dt
     K23 = v_z * dt + 0.5 * (K16) * dt
     K24 = 0
     K25 = - g * dt
     K26 = 0

     K31 = v_x * dt + 0.5 * (K24) * dt
     K32 = v_y * dt + 0.5 * (K25) * dt
     K33 = v_z * dt + 0.5 * (K26) * dt

     K34 = 0
     K35 = - g * dt
     K36 = 0

     K41 = v_x * dt + (K34) * dt
     K42 = v_y * dt + (K35) * dt 
     K43 = v_z * dt + (K36) * dt
     K44 = 0
     K45 = - g * dt
     K46 = 0

     xneu = x + (K11)/6 + (K21)/3 + (K31)/3 + (K41)/6     ! x-Werte
     yneu = y + (K12)/6 + (K22)/3 + (K32)/3 + (K42)/6     ! y-Werte
     zneu = z + (K13)/6 + (K23)/3 + (K33)/3 + (K43)/6         ! z-Werte
     v_xneu = v_x + (K14)/6 + (K24)/3 + (K34)/3 + (K44)/6 ! v_x-Werte
     v_yneu = v_y + (K15)/6 + (K25)/3 + (K35)/3 + (K45)/6 ! v_y-Werte
     v_zneu = v_z + (K16)/6 + (K26)/3 + (K36)/3 + (K46)/6 ! v_z-Werte

     v_y = v_yneu                !Puffer
     y = yneu                    !Puffer

     t = (i+1) * dt              !i+1 gerechnet um Offset-Probleme zu vermeiden
     write(*,*) t, y

  end do

end program RungeKutta
