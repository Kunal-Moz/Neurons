Program Rulkov_model
implicit none
Real*8 :: alpha,mu,sigma,Ie,F
Integer :: t,n,i
Real*8,Dimension(:),allocatable :: X,Y

!!!!!!!!!!!!!!!!!!!!!!! Output Files !!!!!!!!!!!!!!!!!!!!!!!!!!!!
open(unit=1, file = 'Mem_voltage.dat')
open(unit=2, file = 'slow_var.dat')
open(unit=3, file = 'phase.dat')
!!!!!!!!!!!!!!!! Input Parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Print*, "Number of Values:"
Read*, n
allocate(x(n),y(n))

alpha = 6.0
mu = 0.005
sigma = -0.1
x(1) = -0.1
y(1) = 0.0

Do t=1,n-1,1
  x(t+1) = F(x(t),(y(t) + 0.0),alpha)  
  y(t+1) = y(t) - mu*(x(t) + 1.0 - sigma) 
End do

!Print*, x,y
Do i=1,n,1
  write(1,*) i, X(i)
  write(2,*) i, y(i)
End do 

Do i=1,n-1,1
  write(3,*) X(i),X(i+1)
End do 

End Program 


!!!!!!!!!!!!!!!!!!!!!! EXTERNAL APPLIED CURRENT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
function Ie(t) result(val)
implicit none
real*8 :: val
integer :: t
  if (t < 25 .or. t > 75 ) then
    val = 0.0  
  else 
    val = 10.0
  end if
end function Ie

function F(x,y,alpha) result(value)
implicit none
real*8 :: value,x,y,alpha
  if (x .le. 0.0) then
    value = ( alpha/(1.0 - x) ) + y
  else if ( x > 0.0 .and. x < (alpha + y) ) then
    value = alpha + y
  else 
    value = -1.0
  end if
end function F



























