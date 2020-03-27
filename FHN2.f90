!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!! RUNGE-KUTTA METHOD !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine RK(x1,y1,dt,y,x,a,gam,epsi,ie)
Implicit None
Real*8 :: x1,y1,dt,f0,f1,f2,f3,g0,g1,g2,g3,x,y,F,G,a,gam,t,epsi,ie

                 f0 = F(x1,y1,t,epsi,ie)                 ;           g0 = G(x1,y1,t,a,gam)

        f1 = F(x1 + 0.5*dt*f0,y1 + 0.5*dt*g0,t + 0.5*dt,epsi,ie)  ;   g1 = G(x1 + 0.5*dt*f0,y1 + 0.5*dt*g0,t + 0.5*dt,a,gam)

        f2 = F(x1 + 0.5*dt*f1,y1 + 0.5*dt*g1,t + 0.5*dt,epsi,ie)  ;   g2 = G(x1 + 0.5*dt*f1,y1 + 0.5*dt*g1,t + 0.5*dt,a,gam)
 
        f3 = F(x1 + dt*f2,y1 + dt*g2,t + dt,epsi,ie)          ;   g3 = G(x1 + dt*f2,y1 + dt*g2,t + dt,a,gam)
  y = y1 + (dt*(g0 + 2*g1 + 2*g2 + g3)/6.0)    ;       x = x1 + (dt*(f0 + 2*f1 + 2*f2 + f3)/6.0)
End subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! LIST OF FOUR INITIAL CONDITONS !!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine var_list(x0,y0,xvar,yvar,dt,a,gam,epsi,ie)
Implicit none
real*8,dimension(1:4) :: xvar,yvar
real*8 :: x0,y0,dt,x,y,t,xi,yi,epsi,ie,a,gam
integer :: ii
xvar(1) = x0
yvar(1) = y0
t = 0.0
xi = x0
yi = y0
Do ii = 1,3,1
  t = t + dt
  call RK(xi,yi,dt,y,x,a,gam,epsi,ie)
  xvar(ii+1) = x
  yvar(ii+1) = y
  xi = x
  yi = y
End do
ENd Subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! ADAM-BASHFORTH-MOULTON METHOD !!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine ABM(x1,x2,x3,x4,y1,y2,y3,y4,x,y,dt,t,a,gam,epsi,ie)
Implicit None
Real*8 :: F,G,f0,f1,f2,f3,g0,g1,g2,g3,t,dt,px,py,x,y,x1,x2,x3,x4,y1,y2,y3,y4,a,gam,epsi,ie 


f0 = F(x1,y1,t,epsi,ie)           ;     g0 = G(x1,y1,t,a,gam)
f1 = F(x2,y2,t+dt,epsi,ie)        ;     g1 = G(x2,y2,t+dt,a,gam)
f2 = F(x3,y3,t+2.0*dt,epsi,ie)    ;     g2 = G(x3,y3,t+2.0*dt,a,gam)
f3 = F(x4,y4,t+3.0*dt,epsi,ie)    ;     g3 = G(x4,y4,t+3.0*dt,a,gam)


px = x4 + (dt/24.0)*(-9.0*f0 + 37.0*f1 - 59.0*f2 + 55.0*f3) ; py = y4 + (dt/24.0)*(-9.0*g0 + 37.0*g1 - 59.0*g2 + 55.0*g3)

x = x4 + (dt/24.0)*(f1 - 5.0*f2 + 19.0*f3 + 9.0*F(px,py,t+4.0*dt,epsi,ie))
y = y4 + (dt/24.0)*(g1 - 5.0*g2 + 19.0*g3 + 9.0*G(px,py,t+4.0*dt,a,gam))

End Subroutine ABM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! V - EQUATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function F(v,w,t,epsi,ie) result(val)
implicit none
real*8 :: v,w,epsi,val,ie,t
val = (v - ((v**3)/3.0) - w - Ie)
end function f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! W - EQUATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function G(v,w,t,a,gam) result(val)
implicit none
real*8 :: v,a,val,w,gam,t
val = (v + a - gam*w)
end function G

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! MAIN PROGRAM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Program Fitzhugh_Nagumo
Implicit none
Real*8,Dimension(1:4) :: vvar,wvar
Real*8,Dimension(:),Allocatable :: v,w
Real*8 :: v0,w0,mu,t,dt,vo,wo,v1,v2,v3,v4,w1,w2,w3,w4,a,gam,epsi,ie
Integer :: ii,n

!!!!!!!!!!!!! OPEN FILES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Open(unit=11, file = 'pot_curve.dat')
Open(unit=12, file = 'var_curve.dat')
open(unit=13, file = 'Phase.dat')

!!!!!!!!!!!!!!!!! ENTER PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!
Print*, "Step Size:"
Read*, dt
Print*, "Number of Iterations:"
Read*, n

!!!!!!!!!!!!!!!!!! INITIAL CONDITIONS !!!!!!!!!!!!!!!!!!
Allocate(V(n),W(n))
v0 = 10
w0 = -0.7
gam = 0.005
epsi = 1.0
a = -0.97
ie = 10
t = 0.0

call var_list(v0,w0,vvar,wvar,dt,a,gam,epsi,ie)


Do ii = 1,4,1
  t = t + dt
  v(ii) = vvar(ii)
  w(ii) = wvar(ii)
End do
v1 = vvar(1) ; w1 = wvar(1)
v2 = vvar(2) ; w2 = wvar(2)
v3 = vvar(3) ; w3 = wvar(3)
v4 = vvar(4) ; w4 = wvar(4)

Do ii =5,n,1
  t = t + dt
  
  call ABM(v1,v2,v3,v4,w1,w2,w3,w4,vo,wo,dt,t,a,gam,epsi,ie)
  
  v1 = v2 ; w1 = w2
  v2 = v3 ; w2 = w3
  v3 = v4 ; w3 = w4
  v4 = vo ; w4 = wo
  V(ii) = vo ; w(ii) = wo
End do


Do ii = 1,n,1
  t = t + dt
  write(12,*) t, W(ii)
  write(11,*) t, V(ii)
  write(13,*) V(ii), w(ii)
End do 

End Program Fitzhugh_Nagumo
  
  
  



