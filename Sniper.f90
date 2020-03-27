!!!!!!!!!!!!!!!!!!!! RUNGE-KUTTA PROGRAM !!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine Runge_Kutta(x1,y1,x,y,dt,b)   
Implicit None
Real*8 :: x1,y1,x,y,dt,f0,f1,f2,f3,g0,g1,g2,g3,F,G,b

f0 =  F(x1,y1,b)    ;     g0 = G(x1,y1,b)    

f1 = F(x1 + 0.5*dt*f0, y1 + 0.5*dt*g0, b) ; g1 = G(x1 + 0.5*dt*f0, y1 + 0.5*dt*g0, b) 

f2 = F(x1 + 0.5*dt*f1, y1 + 0.5*dt*g1, b) ; g2 = G(x1 + 0.5*dt*f1, y1 + 0.5*dt*g1, b) 

f3 = F(x1 + dt*f2 , y1 + dt*g2, b) ;  g3 = G(x1 + dt*f2 , y1 + dt*g2, b) 

y = y1 + (dt*(g0 + 2*g1 + 2*g2 + g3)/6.0) ; x = x1 + (dt*(f0 + 2*f1 + 2*f2 + f3)/6.0) 
End subroutine

Function F(x,y,b) result(val)
Implicit none
Real*8 :: x,y,b,val
val = x*(1 - (x**2) - (y**2)) + y*(x - b)
End function F

Function G(x,y,b) result(val)
Implicit None
Real*8 :: x,y,b,val
val = y*(1 - (x**2) - (y**2)) - x*(x - b)
End function G

!!!!!!!!!!!!!!!!!!! MAIN PROGRAM !!!!!!!!!!!!!!!!!!!!!!!
Program DE
Implicit None
Real*8 :: x,y,t,dt,b,x0,xi,yi,y0
Integer :: i,n
!!!!!!!!!!!!! OPEN FILES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
open(unit=14, file = 'x-curve.dat')
open(unit=12, file = 'y-curve.dat')
open(unit=13, file = 'phase.dat')
!!!!!!!!!!!! INPUT PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!
Print*, "number of iterations"
Read*, n
Print*, "time step"
read*, dt

x0 = -0.1
y0 = -0.1
b = 0.7
t = 0.0
xi = x0
yi = y0
Do i = 1,n,1
  t = t + dt
  Call Runge_Kutta(xi,yi,x,y,dt,b)
  !print*, xi,vi
  xi = x 
  yi = y
  
  !print*, t,x,v
  Write(14,*) t, x
  Write(13,*) y, x
  write(12,*) t, y
End Do

End Program DE





