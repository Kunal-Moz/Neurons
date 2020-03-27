Function f(x,r) result(val)
Implicit none
Real*8 :: x,r,val
val = r*x - x**3
End Function f

Program Runge_Kutta
Implicit None
Real*8 :: f1,f2,f3,f4,x,dt,x0,f 
Integer :: t,n,rn,r
!Real*8,Dimension(:),allocatable :: 
open(unit=1, file = 'bif.dat')

Print*, "Intial conditon:"
Read*, x0

Print*, "Parameters limit : "
Read*, rn

!Print*, "Time iterations:"
!Read*, N

Print*, "Time step:"
Read*, dt

x = x0
Do r = 1,rn,1
Do t = 1,20,1
  f1 = f(x,real(r)*0.01)*dt
  f2 = f(x + 0.5*f1,real(r)*0.01)*dt
  f3 = f(x + 0.5*f2,real(r)*0.01)*dt
  f4 = f(x + f3,real(r)*0.01)*dt
  x = x + (f1 + 2.0*f2 + 2.0*f3 + f4)/6.0
End do
Write(1,*) r,x
End do


End Program Runge_Kutta
