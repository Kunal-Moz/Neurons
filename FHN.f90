!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! FITZHUGH-NAGUMO NEURON MODEL !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        dv/dt = (v - ((v**3)/3) - w)/epsilon
!        dw/dt = v + a - gamma*y
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! v - EQUATION 
function F(v,w,epsi,ie) result(val)
implicit none
real*8 :: v,w,epsi,val,ie
val = (v - ((v**3)/3.0) - w - Ie)
end function f

! w - EQUATION
function G(v,w,a,gam,epsi) result(val)
implicit none
real*8 :: v,a,val,w,gam,epsi
val = (v + a - gam*w)
end function G

! RUNGE-KUTTA METHOD 
Subroutine Runge_Kutta(x1,y1,x,y,dt,a,epsi,gam,ie)   
Implicit None
Real*8 :: x1,y1,x,y,dt,f0,f1,f2,f3,g0,g1,g2,g3,F,G,a,epsi,gam,ie

f0 =  F(x1,y1,epsi,ie)    ;     g0 = G(x1,y1,a,gam,epsi)    

f1 = F(x1 + 0.5*dt*f0, y1 + 0.5*dt*g0, epsi,ie) ; g1 = G(x1 + 0.5*dt*f0, y1 + 0.5*dt*g0, a,gam,epsi) 

f2 = F(x1 + 0.5*dt*f1, y1 + 0.5*dt*g1, epsi,ie) ; g2 = G(x1 + 0.5*dt*f1, y1 + 0.5*dt*g0,  a,gam,epsi) 

f3 = F(x1 + dt*f2 , y1 + dt*g2, epsi,ie) ;  g3 = G(x1 + dt*f2 , y1 + dt*g0, a,gam,epsi) 

y = y1 + (dt*(g0 + 2*g1 + 2*g2 + g3)/6.0) ; x = x1 + (dt*(f0 + 2*f1 + 2*f2 + f3)/6.0) 
End subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! MAIN PROGRAM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Program Fitzhugh_Nagumo
Implicit None
Real*8 :: v,w,t,dt,a,wi,vi,epsi,gam,ie
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

gam = 0.005
epsi = 1.0
a = -0.97
ie = 10

t = 0.0
vi = 10
wi = -0.7
Do i = 1,n,1
  t = t + dt
  Call Runge_Kutta(vi,wi,v,w,dt,a,epsi,gam,ie)
  !print*, xi,vi
  vi = v 
  wi = w
  
  !print*, t,x,v
  Write(14,*) t, v
  Write(13,*) v, w 
  write(12,*) t, w
End Do

End Program Fitzhugh_Nagumo


