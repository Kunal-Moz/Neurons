!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!! RUNGE-KUTTA METHOD !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine RK(x1,y1,dt,y,x,mu)
Implicit None
Real*8 :: x1,y1,dt,f0,f1,f2,f3,g0,g1,g2,g3,x,y,F,G,mu,t

                 f0 = F(x1,y1,t,mu)                 ;           g0 = G(x1,y1,t,mu)

        f1 = F(x1 + 0.5*dt*f0,y1 + 0.5*dt*g0,t + 0.5*dt,mu)  ;   g1 = G(x1 + 0.5*dt*f0,y1 + 0.5*dt*g0,t + 0.5*dt,mu)

        f2 = F(x1 + 0.5*dt*f1,y1 + 0.5*dt*g1,t + 0.5*dt,mu)  ;   g2 = G(x1 + 0.5*dt*f1,y1 + 0.5*dt*g1,t + 0.5*dt,mu)
 
        f3 = F(x1 + dt*f2,y1 + dt*g2,t + dt,mu)          ;   g3 = G(x1 + dt*f2,y1 + dt*g2,t + dt,mu)
  y = y1 + (dt*(g0 + 2*g1 + 2*g2 + g3)/6.0)    ;       x = x1 + (dt*(f0 + 2*f1 + 2*f2 + f3)/6.0)
End subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! LIST OF FOUR INITIAL CONDITONS !!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine var_list(x0,y0,xvar,yvar,dt,mu)
Implicit none
real*8,dimension(1:4) :: xvar,yvar
real*8 :: x0,y0,dt,x,y,mu,t,xi,yi
integer :: ii
xvar(1) = x0
yvar(1) = y0
t = 0.0
xi = x0
yi = y0
Do ii = 1,3,1
  t = t + dt
  call RK(xi,yi,dt,y,x,mu)
  xvar(ii+1) = x
  yvar(ii+1) = y
  xi = x
  yi = y
End do
ENd Subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! ADAM-BASHFORTH-MOULTON METHOD !!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine ABM(x1,x2,x3,x4,y1,y2,y3,y4,x,y,dt,t,mu)
Implicit None
Real*8 :: F,G,f0,f1,f2,f3,g0,g1,g2,g3,t,dt,px,py,x,y,mu,x1,x2,x3,x4,y1,y2,y3,y4 


f0 = F(x1,y1,t,mu)           ;     g0 = G(x1,y1,t,mu)
f1 = F(x2,y2,t+dt,mu)        ;     g1 = G(x2,y2,t+dt,mu)
f2 = F(x3,y3,t+2.0*dt,mu)    ;     g2 = G(x3,y3,t+2.0*dt,mu)
f3 = F(x4,y4,t+3.0*dt,mu)    ;     g3 = G(x4,y4,t+3.0*dt,mu)


px = x4 + (dt/24.0)*(-9.0*f0 + 37.0*f1 - 59.0*f2 + 55.0*f3) ; py = y4 + (dt/24.0)*(-9.0*g0 + 37.0*g1 - 59.0*g2 + 55.0*g3)

x = x4 + (dt/24.0)*(f1 - 5.0*f2 + 19.0*f3 + 9.0*F(px,py,t+4.0*dt,mu))
y = y4 + (dt/24.0)*(g1 - 5.0*g2 + 19.0*g3 + 9.0*G(px,py,t+4.0*dt,mu))

End Subroutine ABM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! V - EQUATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Function F(x,v,t,mu) result(val)
Implicit None
Real*8 :: x,v,t,val,mu
val = v
End function F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!! X - EQUATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Function G(x,v,t,mu) result(value)
Implicit none
Real*8 :: x,v,t,value,mu
value = mu*(1 - (x**2))*v - x
End function G

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! MAIN PROGRAM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Program Van_der_pol
Implicit none
Real*8,Dimension(1:4) :: xvar,vvar
Real*8,Dimension(:),Allocatable :: X,V
Real*8 :: x0,v0,mu,t,dt,xo,vo,x1,x2,x3,x4,v1,v2,v3,v4
Integer :: ii,n

!!!!!!!!!!!!! OPEN FILES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Open(unit=11, file = 'X_curve.dat')
Open(unit=12, file = 'V_curve.dat')
open(unit=13, file = 'Phase.dat')

!!!!!!!!!!!!!!!!! ENTER PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!
Print*, "Step Size:"
Read*, dt
Print*, "Number of Iterations:"
Read*, n

!!!!!!!!!!!!!!!!!! INITIAL CONDITIONS !!!!!!!!!!!!!!!!!!
Allocate(X(n),V(n))
x0 = 0.1
v0 = 0.1
mu = 5.0
t = 0.0

call var_list(x0,v0,xvar,vvar,dt,mu)


Do ii = 1,4,1
  t = t + dt
  X(ii) = xvar(ii)
  V(ii) = Vvar(ii)
End do
x1 = xvar(1) ; v1 = vvar(1)
x2 = xvar(2) ; v2 = vvar(2)
x3 = xvar(3) ; v3 = vvar(3)
x4 = xvar(4) ; v4 = vvar(4)

Do ii =5,n,1
  t = t + dt
  
  call ABM(x1,x2,x3,x4,v1,v2,v3,v4,xo,vo,dt,t,mu)
  
  x1 = x2 ; v1 = v2
  x2 = x3 ; v2 = v3
  x3 = x4 ; v3 = v4
  x4 = xo ; v4 = vo
  X(ii) = xo   ; V(ii) = vo
End do


Do ii = 1,n,1
  t = t + dt
  write(11,*) t, X(ii)
  write(12,*) t, V(ii)
  write(13,*) X(ii),V(ii)
End do 

End Program van_der_pol
  
  
  

















