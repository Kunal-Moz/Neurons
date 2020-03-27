!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! HINDMARSH-ROSE NEURON MODEL !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      dx/dt = y + phi(x) - z + I = f(x)    !!!! membrane potential
!      dy/dt = psi(x) - y = g(x)            !!!! transport of sodium
!      dz/dt = r*(s*(x - xr) - z) = h(x)    !!!! transport of pottassium
!      ! phi(x) = -a*x**3 + b*x**2
!      ! psi(x) = c - d*x**2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Function phi(x,a,b) result(val)
implicit none
real*8 :: a,b,x,val
val = -a*(x**3) + b*(x**2)
end function phi

Function psi(x,c,d) result(val)
implicit none
real*8 :: c,d,x,val
val = c - d*(x**2)
end function psi
!!!!!!!!!!!!!!!!!!!!!! X equation !!!!!!!!!!!!!!!!!!!!!!
Function F(x,y,z,I,a,b) result(val)
implicit none
real*8 :: x,y,z,phi,val,I,a,b
val = y + phi(x,a,b) - z + I
End function F
!!!!!!!!!!!!!!!!!!!!! Y equation !!!!!!!!!!!!!!!!!!!!!!
Function G(x,y,c,d) result(val)
implicit none
real*8 :: x,y,psi,val,c,d
val = psi(x,c,d) - y
end function G
!!!!!!!!!!!!!!!!!!!! Z equation !!!!!!!!!!!!!!!!!!!!!!!
Function H(x,z,s,r,xr) result(val)
implicit none
real*8 :: x,z,r,s,xr,val
val = r*(s*(x - xr) - z)
end function H

!!!!!!!!!!!!!!!!!! RUNGE-KUTTA ALGORITHM !!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine Runge_Kutta(x1,y1,z1,x,y,z,dt,r,xr,s,i,a,b,c,d)   
Implicit None
Real*8 :: x1,y1,z1,x,y,z,dt,r,xr,f0,f1,f2,f3,g0,g1,g2,g3,h0,h1,h2,h3,i,phi,psi,f,g,h,s,a,b,c,d

f0 =  F(x1,y1,z1,I,a,b)    ;     g0 = G(x1,y1,c,d)     ;    h0 = H(x1,z1,s,r,xr)

f1 = F(x1 + 0.5*dt*f0, y1 + 0.5*dt*g0, z1 + 0.5*dt*h0, I,a,b) ;
g1 = G(x1 + 0.5*dt*f0, y1 + 0.5*dt*g0, c,d) ;
h1 = H(x1 + 0.5*dt*f0, z1 + 0.5*dt*h0 ,s,r,xr)

f2 = F(x1 + 0.5*dt*f1, y1 + 0.5*dt*g1, z1 + 0.5*dt*h1, I,a,b) ; 
g2 = G(x1 + 0.5*dt*f1, y1 + 0.5*dt*g1, c,d) ;
h2 = H(x1 + 0.5*dt*f1, z1 + 0.5*dt*h1, s,r,xr)

f3 = F(x1 + dt*f2 , y1 + dt*g2 , z1 + dt*h2, I,a,b) ;  
g3 = G(x1 + dt*f2 , y1 + dt*g2,c,d) ;
h3 = H(x1 + dt*f2 , y1 + dt*h2 , s, r, xr)

y = y1 + (dt*(g0 + 2*g1 + 2*g2 + g3)/6.0) ; x = x1 + (dt*(f0 + 2*f1 + 2*f2 + f3)/6.0) ; z = z1 + (dt*(h0 + 2*h1 + 2*h2 + h3)/6.0)
End subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! MAIN PROGRAM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program Hindmarsh_Rose
implicit none
real*8 :: x,y,z,t,dt,psi,phi,I,a,b,c,d,s,xr,r,x0,xi,yi,zi
integer :: n,j

!!!!!!!!!!!!!! OPEN FILES !!!!!!!!!!!!!!!!!!!!!!!!!
open(unit=1, file = 'Membrane_pot.dat')
open(unit=2, file = 'fc.dat')
open(unit=3, file = 'sc.dat')
open(unit=14, file = 'phase-zx.dat')
open(unit=12, file = 'phase-xy.dat')
open(unit=13, file = 'phase-yz.dat')
!!!!!!!!!!!!!! INPUT PARAMETERS !!!!!!!!!!!!!!!!!!!!

Print*, "number of iterations"
Read*, n
Print*, "time step"
read*, dt

I = 5.0
x0 = -0.5
s = 4.0
r = 0.01 ; xr = -1.6
a = 1.0 ; c = 1.0 ; d = 5.0
b = 2.7 ! Square wave Burst
!b = 2.52 ! Plateau-like Burst

xi = x0
yi = 0.0
zi = 0.0
t = 0.0
Do j = 1,n,1
  t = t + dt
  Call Runge_Kutta(xi,yi,zi,x,y,z,dt,r,xr,s,i,a,b,c,d)
  !print*, xi,vi
  xi = x 
  yi = y
  zi = z
  !print*, t,x,v
  Write(1,*) t, x
  Write(12,*) x, y
  Write(13,*) y, z
  write(14,*) z, x
  write(2,*) t, y
  write(3,*) t, z
End Do

End Program Hindmarsh_Rose

