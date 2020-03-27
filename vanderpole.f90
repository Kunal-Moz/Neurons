Program van_der_Pole
Implicit None
Real*8:: RKf,RKg,F,G,xi,vi,x0,v0,x,v,mu,E,dt,t
Integer :: n,i

!!!!!!!!!!!!!!!!!!!!! FILES !!!!!!!!!!!!!!!!!!!!!!!!!!!!
Open(unit=11, file = 'Displacement.dat')
Open(unit=12, file = 'Velocity.dat')
Open(unit=13, file = 'Energy.dat')
open(unit=14, file = 'Phase.dat')

!!!!!!!!!!!!!!!!!!!!! INITIAL CONDITION !!!!!!!!!!!!!!!!!!!!!!!
Print*, "Intial Velocity: "
Read*, v0
Print*, "Intial Displacement: "
Read*, x0
Print*, "Step Size:"
Read*, dt
Print*, "Number of Iterations:"
Read*, n
!!!!!!!!!!!!!!!!!!!!! DIFFERENTIAL EQUATIONs !!!!!!!!!!!!!!!!
!! dv/dt = -k*x/m
!! dx/dt = v

!!!!!!!!!!!!!!!!!!!!! FUNCTION VALUES !!!!!!!!!!!!!!!!!!!!!!!'

Mu = 5
xi = x0
vi = v0
t = 0

Do i = 1,n,1
  t = t + dt
  Call RK(xi,vi,dt,v,x,mu)
  !print*, xi,vi
  xi = x 
  vi = v
  !E = 0.5*(k*(x**2) + m*(v**2))
  !print*, t,x,v
  Write(11,*) t, x
  Write(12,*) t, v
  !Write(13,*) t, E
  write(14,*) x, v
End Do

End Program van_der_Pole  
  

!!!!!!!!!!!!!!!!!!!!! RUNGE KUTTA-4th ORDER !!!!!!!!!!!!!!!!!!!!!!!!

Function F(x,v) result(value)
Implicit None
Real*8 :: x,v,value
value = v   !!dx/dt = v = f(x,v)
End Function F

Function G(x,v,mu) result(val)
Implicit None
Real*8 :: x,v,mu,val
val = mu*(1 - (x**2))*v - x   !!dv/dt = mu*(1-x**2)*v - x = g(x,v)
End Function

Subroutine RK(x1,v1,dt,v,x,mu)
Implicit None
Real*8 :: x1,v1,x2,v2,x3,v3,x4,v4,dt,f0,f1,f2,f3,g0,g1,g2,g3,x,v,F,G,mu

           f0 = F(x1,v1)                 ;           g0 = G(x1,v1,mu)

        f1 = F(x1 + 0.5*dt*f0,v1 + 0.5*dt*g0)  ;   g1 = G(x1 + 0.5*dt*f0,v1 + 0.5*dt*g0,mu)

        f2 = F(x1 + 0.5*dt*f1,v1 + 0.5*dt*g1)  ;   g2 = G(x1 + 0.5*dt*f1,v1 + 0.5*dt*g1,mu)
 
        f3 = F(x1 + dt*f2,v1 + dt*g2)          ;   g3 = G(x1 + dt*f2,v1 + dt*g2,mu)
    v = v1 + (dt*(g0 + 2*g1 + 2*g2 + g3)/6.0)        ;       x = x1 + (dt*(f0 + 2*f1 + 2*f2 + f3)/6.0)
End subroutine

