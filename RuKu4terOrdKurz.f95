program RungeKutta

  !PROGRAMM: RUNGE-KUTTA-4TER-ORDNUNG

  implicit none   

  !DEFINITION DER VARIABLEN:
  real :: y, v_y, yneu, v_yneu  !y-Koord. [m] & Geschw. in y-Richt.[m/s]
  real :: x, v_x, xneu, v_xneu  !x-Koord. [m] & Geschw. in x-Richt.[m/s] 
  real :: z, v_z, zneu, v_zneu  !z-Koord. [m] & Geschw. in z-Richt.[m/s]
  real :: t, dt                 !Flugzeit [s] & Zeitintervall [s]
  real :: v_0                   !Anfangsgeschw.    [m/s]
  real :: g = 9.81              !Erdbeschleunigung [m/s²]
  real :: alpha                 !Abwurfwinkel      [rad]    
  integer :: q, n, i, j         !Laufvariablen
  real :: Pi = 3.1415926535898  !Ludolphsche Zahl 
  real, dimension(4,6) :: K     !K = (K(1,1), ..., K(1,6),..., ..., K(4,6))  Matrix
  real, dimension(6) ::   W     !W = (/W(1)= x, W(2)= y, W(3)= z, W(4)= v_x, W(5)= v_y, W(6)= v_z/) 6-dim Vektor  
  real :: Mj 

  !ZUSAETZLICH SOLLEN DIE ERGEBNISSE ALS TEXTDATEI IM FOLGENDEM VERZEICHNIS GESPEICHERT WERDEN:
  open(7, File='/afs/.geo.uni-koeln.de/usr/rizzo/my/f95/NumerischeMethoden/RungeKutta/RuKu4terOrdn/dt0k222222.txt')


  !SETZE PARAMETER:
  alpha = 45.0 * (Pi/180.0)
  v_0 = 20.0
  dt = 0.2
  n = 14


  !STARTWERTE: 
  t = 0.0 
  W = 0.0
  x = 0.0
  y = 0.0
  z = 0.0                       !kann man 0 setzen, wenn man das Problem 2dimensional betrachtet
  v_x = v_0 *cos(alpha)
  v_y = v_0 *sin(alpha)
  v_z = 0.0                     !kann man 0 setzen, wenn man das Problem 2dimensional betrachtet



  !DEFINITION DES VEKTORS W (ORTS- & GESCHWINDIGKEITSKOMPONENTE)N:  
  W(1) = x
  W(2) = y
  W(3) = 0
  W(4) = v_0 *cos(alpha)        !v_x = v_0 * cos(alpha) bleibt waehrend der Flugzeit konstant
  W(5) = v_0 *sin(alpha)        !v_y = v_0 * sin(alpha) aendert sich waehrend der Flugzeit  
  W(6) = 0

  write(*,*) "------------------------------------------------------------------------------------------------------------"
  write(*,*) "Programm zur Berechnung der Bahn des schiefen Wurfes nach Runge-Kutta-4ter-Ordnung" 
  write(*,*) "Erstellt am 01.09.2010."
  write(*,*) "."
  write(*,*) " Time dt[s]         x [m]           y [m]           z[m]          v_x [m/s]       v_y [m/s]      v_z [m/s]" 
  write(*,*) "------------------------------------------------------------------------------------------------------------"
  write(*,*) t, x, y, z, v_x, v_y, v_z        !Ausgeben der Anfangswerte 

  write(7,*) " Time dt[s]         x [m]           y [m]           z[m]          v_x [m/s]       v_y [m/s]      v_z [m/s]" 
  write(7,*) "------------------------------------------------------------------------------------------------------------"
  write(7,*) t, x, y, z, v_x, v_y, v_z   


  !BERECHNUNG:

  do q = 0, n !Grosse do-Zeitschleife

     !BESTIMMUNG DER MATRIXELEMENTE: 

     K(1,1) = W(4) * dt                      !v_x * dt
     K(1,2) = W(5) * dt                      !v_y * dt 
     K(1,3) = W(6) * dt                      !v_z * dt  
     K(1,4) = 0
     K(1,5) = - g * dt
     K(1,6) = 0
     K(2,1) = W(4) * dt + 0.5 * K(1,4) * dt  !v_x * dt + 0.5 * K(1,4) * dt
     K(2,2) = W(5) * dt + 0.5 * K(1,5) * dt  !v_y * dt + 0.5 * K(1,5) * dt
     K(2,3) = W(6) * dt + 0.5 * K(1,6) * dt  !v_z * dt + 0.5 * K(1,6) * dt 
     K(2,4) = 0
     K(2,5) = - g * dt
     K(2,6) = 0
     K(3,1) = W(4) * dt + 0.5 * K(2,4) * dt  !v_x * dt + 0.5 * K(2,4) * dt 
     K(3,2) = W(5) * dt + 0.5 * K(2,5) * dt  !v_y * dt + 0.5 * K(2,5) * dt 
     K(3,3) = W(6) * dt + 0.5 * K(2,6) * dt  !v_z * dt + 0.5 * K(2,6) * dt
     K(3,4) = 0
     K(3,5) = - g * dt
     K(3,6) = 0
     K(4,1) = W(4) * dt + K(3,4) * dt        !v_x * dt + K(3,4) * dt
     K(4,2) = W(5) * dt + K(3,5) * dt        !v_y * dt + K(3,5) * dt 
     K(4,3) = W(6) * dt + K(3,6) * dt        !v_z * dt + K(3,6) * dt
     K(4,4) = 0
     K(4,5) = - g * dt
     K(4,6) = 0


     do i = 1, 6                             !Schleife fuer die Koordinaten x,y,z,v_x,v_y,v_z 
        do j = 1, 4                          !Schleife fuer die Matrixelemente Spalte



           if ((j == 1) .or. (j == 4)) then  !Die if-Schleife setzt den Nenner Mj auf 6 bzw. 3
              Mj = 6
           else
              Mj = 3
           endif


           W(i) = W(i) + K(j,i)/Mj           !Runge-Kutta-4ter-Ordnung-Algorithmus


           t = (q+1) * dt                    !q+1 gerechnet um Offset-Probleme zu vermeiden


        end do                               !Ende der i-Schleife

     end do                                  !Ende der j-Schleife


     write(*,*) t, W
     write(7,*) t, W



  end do                                     !Ende der grossen do-Zeit-Schleife
  close(7)

end program RungeKutta
