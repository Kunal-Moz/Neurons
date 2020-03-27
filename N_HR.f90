!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Code for a chain of HR neurons with a delayed coupling and time-scale mismatch !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! HINDMARSH-ROSE NEURON MODEL !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      dxi/dt = yi + phi(x1) - z1 + I + gs*((V - x1)/(1.0 + exp(-lambda*(xj(t- tau) - K)))) = f(x)    !!!! membrane potential
!      dyi/dt = psi(xi) - yi = g(xi)            !!!! transport of sodium
!      dzi/dt = epsi*(s*(xi - xr) - zi) = h(xi)    !!!! transport of pottassium
!      ! phi(xi) = -a*xi**3 + b*xi**2
!      ! psi(xi) = c - d*xi**2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!! X equation !!!!!!!!!!!!!!!!!!!!!!
Function F(x,y,z,Ie,gs,coupl_sum) result(val)
implicit none
real*8 :: x,y,z,val,Ie,V,coupl_sum,gs
V = 2.0
val = y - 1.0*x**3 + 3.0*x**2 - z + Ie + gs*(V - x)*coupl_sum
End function F

!!!!!!!!!!!!!!!!!!!!! Y equation !!!!!!!!!!!!!!!!!!!!!!
Function G(x,y) result(val)
implicit none
real*8 :: x,y,val
val = 1.0 - 5.0*x**2 - y
end function G

!!!!!!!!!!!!!!!!!!!! Z equation !!!!!!!!!!!!!!!!!!!!!!!
Function H(x,z,epsi) result(val)
implicit none
real*8 :: x,z,epsi,val
val = epsi*(4.0*(x + 1.6) - z)
end function H

!!!!!!!!!!!!!!!!!!!!! Nonlinear Coupling Function !!!!!!!!!!!!!!!!!!!!!!!!
function Cupling(x) result(cop)
implicit none
real*8 :: x,K,lambda,cop
K = -0.25
lambda = 1
cop = 1.0/(1.0 + exp(-lambda*(x - K)))
end function

!!!!!!!!!!!!!!!!!!!!!!! Uniform random number generator !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function random_uniform(seed) result(ran)
implicit none
integer :: seed
real*8 :: ran
seed = mod(7**5*seed,2147483647)
ran = seed/2147483647.
end function random_uniform

!!!!!!!!!!!!!!!!!!!!!!!!!! Initialize connectivity matrix !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine Adj_Matrix_chain(A,N,seed)
implicit None
real*8,Dimension(N,N) :: A
integer :: i,j,n,seed,ir,il
real*8 :: p,random_uniform

Do i=1,n,1
  Do j=1,n,1
    !ir=i+1 ; if (ir > n) ir = 1
    il=i-1 ; if (il < 1) il = n
    if (il .eq. j) then
      a(i,j) = 1.0
    else
      a(i,j) = 0.0
    end if  
  End do
End do  
  
End Subroutine


!!!!!!!!!!!!!!!!!!!! Calculate the derivatives !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine differetial(x,y,z,Ie,epsi,gs,tau,A,N,xdot,ydot,zdot,eta)
implicit none
real*8 :: dt,epsi,gs,tau,F,G,H,cupling,Sum1
integer :: n,j,i
real*8,dimension(1:n) :: X,Y,Z,Xdot,Ydot,Zdot,eta,Ie
real*8,dimension(N,N) :: A
!Sum1 = 0
do i=1,n,1
  Sum1 = 0
  do j=1,n,1
    Sum1 = Sum1 + a(i,j)*cupling(x(j))
  End do

  Xdot(i) = eta(i)*F(X(i),Y(i),Z(i),Ie(i),gs,sum1)
  ydot(i) = eta(i)*G(X(i),Y(i))
  zdot(i) = eta(i)*H(X(i),Z(i),epsi)
End do

End Subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!! Runge-Kutta for networks !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine Rk_network(X,Y,Z,Ie,epsi,gs,eta,tau,A,N,dt,niter)
implicit none
real*8 :: xd,dt,F,G,H,gs,tau,epsi,t
integer :: N,j,i,niter
real*8,dimension(1:n) :: X,Y,Z,eta,Ie
real*8,dimension(N:N) :: A
real*8,Dimension(:,:),allocatable :: Kmx,Kmy,Kmz
real*8,Dimension(:),allocatable :: Xdot,Ydot,Zdot,Xt,Yt,Zt
Allocate(Kmx(n,4),Kmy(n,4),Kmz(n,4))
allocate(Xt(n),Yt(n),Zt(n),xdot(n),ydot(n),zdot(n)) 
t = 0.0

Xt = X ; Yt = Y ; Zt = Z

Do j = 1,niter,1
  t = t + dt
  
  Do i = 1,N,1
    x(i) = xt(i) ;
    y(i) = yt(i) ;
    z(i) = zt(i)
  End do    
  call differetial(xt,yt,zt,Ie,epsi,gs,tau,A,N,xdot,ydot,zdot,eta)
  Do i = 1,N,1
    Kmx(i,1) = xdot(i) ;
    Kmy(i,1) = ydot(i) ;
    Kmz(i,1) = zdot(i)
  End do
  
  Do i = 1,N,1
    xt(i) = xt(i) + 0.5*dt*Kmx(i,1) ;
    yt(i) = yt(i) + 0.5*dt*Kmy(i,1) ;
    zt(i) = zt(i) + 0.5*dt*Kmz(i,1)
  End do 
  call differetial(xt,yt,zt,Ie,epsi,gs,tau,A,N,xdot,ydot,zdot,eta)
  Do i = 1,N,1
    Kmx(i,2) = xdot(i) ;
    Kmy(i,2) = ydot(i) ;
    Kmz(i,2) = zdot(i)
  End do
  
  Do i = 1,N,1
    xt(i) = xt(i) + 0.5*dt*Kmx(i,2) ;
    yt(i) = yt(i) + 0.5*dt*Kmy(i,2) ;
    zt(i) = zt(i) + 0.5*dt*Kmz(i,2)
  End do  
  call differetial(xt,yt,zt,Ie,epsi,gs,tau,A,N,xdot,ydot,zdot,eta)
  Do i = 1,n,1
    Kmx(i,3) = xdot(i) ;
    Kmy(i,3) = ydot(i) ;
    Kmz(i,3) = zdot(i)
  End do
  
  Do i = 1,n,1
    xt(i) = xt(i) + dt*Kmx(i,3) ;
    yt(i) = yt(i) + dt*Kmy(i,3) ;
    zt(i) = zt(i) + dt*Kmz(i,3)
  End do
  call differetial(xt,yt,zt,Ie,epsi,gs,tau,A,N,xdot,ydot,zdot,eta)
  Do i = 1,n,1
    Kmx(i,4) = xdot(i) ;
    Kmy(i,4) = ydot(i) ;
    Kmz(i,4) = zdot(i)
  End do
  
  Do i=1,n,1
    xt(i) = x(i) + (1.0/6.0)*(kmx(i,1) + 2.0*kmx(i,2) + 2.0*kmx(i,3) + kmx(i,4))*dt ;
    yt(i) = y(i) + (1.0/6.0)*(kmy(i,1) + 2.0*kmy(i,2) + 2.0*kmy(i,3) + kmy(i,4))*dt ;
    zt(i) = z(i) + (1.0/6.0)*(kmz(i,1) + 2.0*kmz(i,2) + 2.0*kmz(i,3) + kmz(i,4))*dt
  End do
  
  Do i = 1,n,1
    if(xt(i).ne.xt(i))then
      print*, "Instability"
      goto 101
    endif
  End Do

  if ( j > (Niter - 2**18) ) then
    
  End if 
  !Print*, "one loop"
  
End do
101 continue

End Subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PRGORAM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM Network_HR_Neurons
implicit none
real*8,dimension(:),allocatable :: X,Y,Z,eta
real*8,dimension(:),allocatable :: Ie,gg
real*8,dimension(:,:),allocatable :: A,Xdata
real*8 :: dt,random_uniform,tau,epsi,sumav,sumdv,gs
integer :: N, Niter, seed,i,j,seed2

!!!!! Open Files !!!!!
open(unit = 1, file = "ts1.dat")
open(unit = 2, file = "ts2.dat")
open(unit = 3, file = "ts3.dat")
open(unit = 4, file = "ts4.dat")
open(unit = 10, file = "lfp.dat")
!!!!! No. of Neurons !!!!!
N = 4

!!!!! Parameters !!!!!
tau = 0
epsi = 0.006
Print*, "Coupling Strength: "
Read*, gs

!!!!! No of iterations and Allocation !!!!!
Niter = 300000
dt = 0.01
seed = 34234
seed2 = 85635
allocate(X(N),Y(N),Z(N),eta(n),Ie(n))!Xdelay(N))
allocate(A(n,n))

!!!!! Initiate adjacency matrix !!!!!
call Adj_Matrix_chain(A,n,seed)
Do i=1,n,1
  print*, A(i,1),A(i,2),A(i,3),A(i,4)
End do
!!!!! Initial condition !!!!! 
Do i = 1,n,1
  x(i) = -0.5 + 0.01*random_uniform(seed)
  y(i) = 0.1 + 0.01*random_uniform(seed)
  z(i) = 0.1 + 0.01*random_uniform(seed)
End do

!Do i=1,n,1
!  Do j=1,n,1
!    open(unit=44,file="Coupling.txt")
!    read*, gs(i,j) 
!  End Do
!Do End

Do i=1,n,1
  open(unit=42,file="Slowness.txt")
  read(42,*) eta(i)
  !eta(i) = 1.0
End do
Print*, "Time scale mismatch: ", eta  
 close(42)

Do i=1,N,1 
  open(unit=43,file="Current.txt")
  Read(43,*) Ie(i)
  !Ie(i) = 2.9
End do
Print*, "Current Input to each neuron:", Ie
 close(43)

call Rk_network(X,Y,Z,Ie,epsi,gs,eta,tau,A,N,dt,niter)
 
  
End program Network_HR_Neurons


