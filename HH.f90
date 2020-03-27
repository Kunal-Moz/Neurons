!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! HODGKIN-HUXLEY NEURON MODEL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!            dV/dt = (gl*(v - el) + gk*(n**4)*(v - ek) + gna*(m**3)*h*(v - ena) - i)*(-1/Cm)              !! 
!!                            dn/dt = alphan(v,t)*( 1 - n) - betan(v,t)*n                                  !!
!!                            dm/dt = alpham(v,t)*( 1 - m) - betam(v,t)*m                                  !!
!!                            dh/dt = alphah(v,t)*( 1 - h) - betah(v,t)*h                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Voltage equation !
function F(v,n,m,h,t,gl,el,gk,ek,gna,ena,i,cm) result(fres)
implicit none
Real*8 :: v,n,m,h,t,gl,el,gk,ek,gna,ena,i,cm,fres
fres = (gl*(v - el) + gk*(n**4)*(v - ek) + gna*(m**3)*h*(v - ena) - i)*(-1/Cm)
end function F

! n equation !
function p(v,n,t) result(gres)
implicit none
real*8 :: v,n,t,alphan,betan,gres
gres = alphan(v,t)*(1 - n) - betan(v,t)*n
end function p

! m equation !
function q1(v,m,t) result(hres)
implicit none
real*8 :: v,m,t,alpham,betam,hres
hres = alpham(v,t)*(1 - m) - betam(v,t)*m
end function q1

! h equation !
function q2(v,m,t) result(lres)
implicit none
real*8 :: v,h,t,alphah,betah,lres
lres = alphah(v,t)*(1 - h) - betah(v,t)*h
end function q2



