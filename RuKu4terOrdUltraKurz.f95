program RungeKuttaUltraKurz

  !PROGRAMM: RUNGE-KUTTA-4TER-ORDNUNG

  implicit none   

  !DEFINITION DER VARIABLEN:
  double precision y, v_y, yneu, v_yneu  !y-Koord. [m] & Geschw. in y-Richt.[m/s]
  double precision x, v_x, xneu, v_xneu  !x-Koord. [m] & Geschw. in x-Richt.[m/s] 
  double precision z, v_z, zneu, v_zneu  !z-Koord. [m] & Geschw. in z-Richt.[m/s]
  double precision t, dt                 !Flugzeit [s] & Zeitintervall [s]
  double precision v_0                   !Anfangsgeschw.    [m/s]
  real :: g = 9.81                       !Erdbeschleunigung [m/s²]
  double precision alpha                 !Abwurfwinkel      [rad]    
  integer :: q, n, i, j                  !Laufvariablen
  real :: Pi = 3.1415926535898           !Ludolphsche Zahl 
  double precision, dimension(4,6) :: K  !K(j,i) = (K(1,1), ..., K(1,6),..., ..., K(4,6))  Matrix
  double precision, dimension(9) :: W    !W(i) = (/W(1)= x, W(2)= y, W(3)= z, W(4)= v_x, W(5)= v_y, W(6)= v_z/) 6-dim Vektor  
  double precision Mj, Mk                !Hilfsgroessen Mj steht fuer 3 bzw. 6 & Mk steht fuer 2 bzw. 1
  double precision A, B                  !Hilfsgroessen um Matrixelemente die es nicht gibt zu umgehen


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
  z = 0.0                                !kann man 0 setzen, wenn man das Problem 2dimensional betrachtet
  v_x = v_0 *cos(alpha)
  v_y = v_0 *sin(alpha)
  v_z = 0.0                              !kann man 0 setzen, wenn man das Problem 2dimensional betrachtet



  !DEFINITION DES VEKTORS W (ORTS- & GESCHWINDIGKEITSKOMPONENTEN):  
  W(1) = x                               !x
  W(2) = y                               !y
  W(3) = 0.0                             !z
  W(4) = v_0 *cos(alpha)                 !v_x = v_0 * cos(alpha) bleibt waehrend der Flugzeit konstant
  W(5) = v_0 *sin(alpha)                 !v_y = v_0 * sin(alpha) aendert sich waehrend der Flugzeit  
  W(6) = 0.0                             !v_z
  W(7) = 0.0                             !a_x 
  W(8) = (-1) * g                        !a_y
  W(9) = 0.0                             !a_z



  !AUSGABE AUF KONSOLE:
  write(*,*) "------------------------------------------------------------------------------------------------------------"
  write(*,*) "Programm zur Berechnung der Bahn des schiefen Wurfes nach Runge-Kutta-4ter-Ordnung" 
  write(*,*) "Erstellt am 01.09.2010."
  write(*,*) "."
  write(*,*) " Time dt[s]         x [m]           y [m]           z[m]          v_x [m/s]       v_y [m/s]      v_z [m/s]" 
  write(*,*) "------------------------------------------------------------------------------------------------------------"
  write(*,522) t, x, y, z, v_x, v_y, v_z !Ausgeben der Anfangswerte
522 format(F7.3, 10X, F7.3, 10X, F7.3, 10X, F7.3, 10X, F7.3, 10X, F7.3, 10X, F7.3, 10X)



  !AUSGABE ALS TEXTDATEI:
  !write(7,*) " Time dt[s]         x [m]           y [m]           z[m]          v_x [m/s]       v_y [m/s]      v_z [m/s]" 
  !write(7,*) "------------------------------------------------------------------------------------------------------------"
   write(7,*) t, y  


  !BERECHNUNG:

  do q = 0, n                            !Grosse do-Zeitschleife



     do j = 1, 4
        do i = 1, 6 


           if (j == 1) then   !Eingefuehrt um Probleme mit Matrixelementen K(0,4),K(0,5),...,K(0,6) zu vermeiden die es nicht gibt  
              A = 0 
           else 
              A = 1
           end if


           if (i > 3) then    !Eingefuehrt um Probleme mit Matrixelementen K(1,7),K(1,8),...,K(1,10) zu vermeiden die es nicht gibt 
              B = 0
           else 
              B= 1
           end if


           if ((j == 2) .or. (j == 3)) then                    !Beruecksichtigt den Nenner Mk, der 2 bzw. 1 sein soll
              Mk = 2
           else
              Mk = 1
           end if




           K(j,i) = dt * W(i+3) + K(j-1,i+3) * (dt/Mk) * A * B ! BERECHNET DIE EINZELNEN MATRIXELEMENTE!!!




        end do
     end do




     do i = 1, 6                                               !Schleife fuer die Koordinaten x,y,z,v_x,v_y,v_z 
        do j = 1, 4                                            !Schleife fuer die Matrixelemente Spalte



           if ((j == 1) .or. (j == 4)) then                    !if-Schleife setzt den Nenner Mj auf 6 bzw. 3
              Mj = 6
           else
              Mj = 3
           endif


           W(i) = W(i) + K(j,i)/Mj                             !RUNGE-KUTTA-4TER-ORDNUNG-ALGORITHMUS!!! 


           t = (q+1) * dt                                      !q+1 gerechnet um Offset-Probleme zu vermeiden


        end do                                                 !Ende der i-Schleife

     end do                                                    !Ende der j-Schleife




 
     write(*,523) t, W(1), W(2), W(3), W(4), W(5), W(6)        !Ausgeben der Anfangswerte
523  format(F7.3, 10X, F7.3, 10X, F7.3, 10X, F7.3, 10X, F7.3, 10X, F7.3, 10X, F7.3, 10X)
     write(7,*) t, W(2)                                        !Ausgabe als Textdatei

  end do                                                       !Ende der grossen do-Zeit-Schleife
  close(7)

end program RungeKuttaUltraKurz
