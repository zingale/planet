!Compute Hydrostatic Equilibrium of a Hot Jupiter 

program CleanHotJupiter
implicit none

  real :: Lad, Lin, mu, m, k, grav !Misc Variables
  real :: P1, P2, Pc, dP !Pressures
  real :: den1, den2 !Densities
  real :: T1, T2, Td !Temperatures
  real :: R, X !Gas Constant, Model Parameter
  real :: z, dz !Position
  real :: exp1, exp2 !Exponents
  integer :: alpha , beta

  integer :: count = 0 !iteration counter
  real :: denMin = 1e-7 !Minimum density to stop iteration



  !Declaration of constants
  Td = 1500.0
  Pc = 1.0e9
  dP = 2e3
  Lad = 2.0/7.0
  Lin = 1.0/2.0
  grav = -1.0e3
  X = 1.0
  
  !Exponents
  alpha = 1
  beta = 0
  exp1 = 1.0 + alpha
  exp2 = 1.0 / (4.0+beta)

  !Gas Constant stuff
  mu = 2.34
  m = 1.66e-24
  k = 1.38e-16
  R = k/(mu*m)
  
  !Initial Conditions
  P1 = 1.0e9
  z = 0
  
  
  open(unit=18,file="newmodelf90.hse",status="replace", action="write") !Open data file to modify
  
  do !Begin do loop to calculate variables of cells

     P2 = P1 - dP !Pressure of adjacent cell

     T1 = Td * (1.0 + (Lad /(Lin - Lad)) * ((P1/Pc)**exp1))**exp2 !Temp of first cell
     T2 = Td * (1.0 + (Lad /(Lin - Lad)) * ((P2/Pc)**exp1))**exp2 !Temp of adjacent cell

     den1 = P1/(R*T1) !Density of first cell
     den2 = P2/(R*T2) !Density of adjacent cell

     dz = -(P1-P2)*(1.0/grav)*(1.0/2.0)*((1.0/den1)+(1.0/den2))
     
     if(den1 < denMin)then !Check if density is significantly small, if so, exit loop
        exit
     endif 

    write(*,*) z , den1 , T1, P1, X 
   

     z = z +dz !Iterate z 
     P1 = P2 !Shift cells right by one
     count = count + 1 !Increase iteration count

  enddo 

  !Go back to beginning of file and write header
  rewind(unit=18)
  write(18,'(a9,i6,/,a22,/,a9,/,a13,/,a10,/,a3)') "# npts = ",count, & 
 "# num of variables = 4" , "# density", "# temperature", "# pressure","# X"

  close(unit=18) !Close data file
endprogram CleanHotJupiter
