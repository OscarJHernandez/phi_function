! This module contains all of the operators that are needed to employ Saori's cut-off dependant M1 operators

module Saori_MEC
use deuteronMEC_params
use matrixElements
use NLO_MEC
implicit none

real(8),parameter:: smallMomentum = 0.05d0 ! This parameter determines the small p-limit of the integrals
real(8),allocatable:: I_Lambdafunc(:)
real(8),allocatable:: dI_Lambdafunc(:)
real(8),allocatable:: d2I_Lambdafunc(:)
real(8),allocatable:: Y_Lambdafunc(:)
real(8),allocatable:: dY_Lambdafunc(:)
real(8),parameter::epsabs = 1.0d-8 ! Absolute tolerance for the numerical oscillatory integrals

contains

subroutine init_Saori_currents()
implicit none
integer::i
real(8)::ri

! Here we initalize all arrays that store the functions evaluated at one point
allocate(I_Lambdafunc(NquadPosition))
allocate(dI_Lambdafunc(NquadPosition))
allocate(d2I_Lambdafunc(NquadPosition))
allocate(Y_Lambdafunc(NquadPosition))
allocate(dY_Lambdafunc(NquadPosition))

    do i=1,NquadPosition
        ri = rMax*dr(i)
            I_Lambdafunc(i) = I_Lambda(ri)
            dI_Lambdafunc(i) = dI_Lambda(ri)
            d2I_Lambdafunc(i) = d2I_Lambda(ri)
            Y_Lambdafunc(i) = Y_Lambda(ri)
            dY_Lambdafunc(i) = dY_Lambda(ri)
    end do
end subroutine

! <2||MEC||1>
real(8) function MEC_NLO_Lambda(n1,l1,ss1,jj1,n2,l2,ss2,jj2,t1,mt1,t2,mt2)
implicit none
integer::n1,n2,l1,l2
integer::jj1,jj2,t1,t2,mt1,mt2
integer::tt1,tt2,mmt1,mmt2
integer::ss1,ss2
real(8)::factor1
real(8)::radial,isospin,angular
real(8),external:: CG
real(8)::s

tt1 = 2*t1
tt2 = 2*t2
mmt1 = 2*mt1
mmt2 = 2*mt2 

! 1/[MeV**2]
factor1 = ((gA**(2)))/(Fpi*Fpi)

! The radial Matrix element in the HO basis with cutt-off
angular = MEC_M1_Operator_Lambda(n1,l1,jj1,ss1,n2,l2,jj2,ss2)

! This part seems ok.
isospin =  dsqrt(2.d0)*full_matrixElement_MEC_Isospin(t1,mt1,t2,mt2)

MEC_NLO_Lambda = factor1*angular*isospin
end function

! The Operator for the MEC currents
! With the cutt-off
! This consists of two parts
! <2|MEC_M1 |1>
real(8) function MEC_M1_Operator_Lambda(n1,l1,jj1,ss1,n2,l2,jj2,ss2)
implicit none
integer::n1,n2
integer::l1,jj1,ss1,l2,jj2,ss2
integer::h,hh
real(8)::s
real(8),external:: SIXJ,ANINEJ,THREEJ
real(8)::factor1
real(8)::radial1
real(8)::radial2

! The minus sign in radial one is fine
! Factor 1 has been corrected
factor1 = (2.d0)*dsqrt(2.d0*pi/3.d0)
radial1 = -1.d0*MEC_RadialME1_Lambda(n1,l1,n2,l2)
radial2 = MEC_RadialME2_Lambda(n1,l1,n2,l2)


s=0.d0
do h=0,2
hh=2*h
s = s+dsqrt(dfloat(2*h+1))*THREEJ(2,2,hh,0,0,0)*reduced_ME_Omega(l1,jj1,ss1,l2,jj2,ss2,1,h,1)*radial1
end do

MEC_M1_Operator_Lambda = s*factor1+dsqrt(2.d0)*radial2*reduced_ME_Sigma12_J(l1,jj1,ss1,1,l2,jj2,ss2)

end function

! The radial ME with the regulator 
! the Function f must be in units of MeV
real(8) function MEC_RadialME1_Lambda(n1,l1,n2,l2)
integer::n2,l2,n1,l1
real(8)::s,ri
integer::i
real(8)::f
real(8)::mpiMeV

! 1/fm conversion
mpiMeV = mpi/(hbarc)

s=0.d0

do i =1,NquadPosition
ri = rMax*dr(i)
f = 0.25d0*ri*dY_Lambdafunc(i)*hbarc
f =f -0.5d0*(1.d0/((2.d0*pi)**3))*(d2I_Lambdafunc(i)-(1.d0/ri)*dI_Lambdafunc(i))*(hbarc**3)

! This is the exact limiting form of the function for testing
!f = -1.d0*(mpi/(8.d0*pi))*(1.d0+(1.d0/(mpiMeV*ri)))*dexp(-1.d0*mpiMeV*ri)

!f= f*hbarc
s = s+LagExpArray(n1,l1,i)*LagExpArray(n2,l2,i)*dwR(i)*f
end do

s=rMax*s

MEC_RadialME1_Lambda =s


end function

! The radial ME with the regulator 
! The Function f must be in units of MeV
real(8) function MEC_RadialME2_Lambda(n1,l1,n2,l2)
integer::n2,l2,n1,l1
real(8)::s,ri
integer::i
real(8)::f
real(8)::mpiMeV

! 1/fm conversion
mpiMeV = mpi/(hbarc)



s=0.d0

do i =1,NquadPosition
ri = rMax*dr(i)

f =  0.5d0*(1.d0/((2.d0*pi)**3))*(d2I_Lambdafunc(i)+(1.d0/ri)*dI_Lambdafunc(i))*(hbarc**3)
f = f-0.25*ri*dY_Lambdafunc(i)*hbarc

! This is the exact function in the Infinite cutoff, for testing purposes
!f = (mpi/(8.d0*pi))*dexp(-mpi*ri/hbarc)

!f=f*hbarc
s = s+LagExpArray(n1,l1,i)*LagExpArray(n2,l2,i)*dwR(i)*f
end do

s=rMax*s

MEC_RadialME2_Lambda =s


end function

! This Controls the form of the regulator that we use
! The user can substitute any form
real(8) function regulator(p)
implicit none
real(8)::p
real(8)::pl
pl = p/Lambda
regulator = dexp(-1.d0*pl*pl)
!regulator =1.d0

end function regulator


! This function uses quadpack routine for the Fourier integral over p
real(8) function A11(r)
implicit none
real(8)::r
real(8),parameter::a=0.d0 ! the starting point of the integration
integer,parameter:: integr = 2 ! for Sin(w*x) function
real(8)::abserr
integer::ier
integer::neval
real(8)::result

call qawf(dA11,a,r/(hbarc),integr,epsabs, result, abserr, neval, ier )

A11 =result
!print *, r,result,ier

end function

! The function that we will integrate without Sin(p*r) function
real(8) function dA11(p)
implicit none
real(8)::p

dA11 = p/(mpi*mpi+p*p)
dA11 = dA11*regulator(p)

end function

! This function uses quadpack routine for the Fourier integral over p
real(8) function A12(r)
implicit none
real(8)::r
real(8),parameter::a=0.d0 ! the starting point of the integration
integer,parameter:: integr = 2 ! for Sin(w*x) function
real(8)::abserr
integer::ier
integer::neval
real(8)::result

call qawf(dA12,a,r/(hbarc),integr,epsabs, result, abserr, neval, ier )

A12 =result

!print *, r,result,ier

end function

! The function that we will integrate without Sin(p*r) function
real(8) function dA12(p)
implicit none
real(8)::p

dA12 = (p)/((mpi*mpi+p*p)**2)
dA12 = dA12*regulator(p)

end function

! This function uses quadpack routine for the Fourier integral over p
real(8) function A32(r)
implicit none
real(8)::r
real(8),parameter::a=0.d0 ! the starting point of the integration
integer,parameter:: integr = 2 ! for Sin(w*x) function
real(8)::abserr
integer::ier
integer::neval
real(8)::result

call qawf(dA32,a,r/(hbarc),integr,epsabs, result, abserr, neval, ier )

A32 =result

!print *, r,result,ier

end function

! The function that we will integrate without Sin(p*r) function
real(8) function dA32(p)
implicit none
real(8)::p

dA32 = (p*p*p)/((mpi*mpi+p*p)**2)
dA32 = dA32*regulator(p)

end function

! This function uses quadpack routine for the Fourier integral over p
real(8) function B22(r)
implicit none
real(8)::r
real(8),parameter::a=0.d0 ! the starting point of the integration
integer,parameter:: integr = 1 ! for Cos(w*x) function
real(8)::abserr
integer::ier
integer::neval
real(8)::result

call qawf(dB22,a,r/(hbarc),integr,epsabs, result, abserr, neval, ier )

B22 =result

!print *, r,result,ier

end function

! The function that we will integrate without Sin(p*r) function
real(8) function dB22(p)
implicit none
real(8)::p

dB22 = (p*p)/((mpi*mpi+p*p)**2)
dB22 = dB22*regulator(p)

end function

! \int_{0}^{\infinity} p^k Cos(p*r)*Regulator(p)/((m^2+p^2)^s)
real(8) function Aks(k,s,r)
implicit none
integer::k,s
real(8)::r,pj,s0,dpj
integer::j

s0=0.d0
do j=1,NquadMomentum
pj = dp(j)/(1.d0-dp(j))
dpj = 1.d0/((1.d0-dp(j))**(2.d0))

!pj = Pmax*dp(j)
s0 = s0+((pj**k)*dsin((pj*r/hbarc))/((mpi*mpi+pj*pj)**s))*regulator(pj)*dwP(j)*dpj
!s0 = s0+((pj**k)*dsin((pj*r/hbarc))/((mpi*mpi+pj*pj)**s))*dwP(j)*dpj
end do

Aks =s0
!Bks = s0*Pmax

end function

! \int_{0}^{\infinity} p^k Cos(p*r)*Regulator(p)/((m^2+p^2)^s)
real(8) function Bks(k,s,r)
implicit none
integer::k,s
real(8)::r,pj,s0,dpj
integer::j

s0=0.d0
do j=1,NquadMomentum
pj = dp(j)/(1.d0-dp(j))
dpj = 1.d0/((1.d0-dp(j))**(2.d0))

!pj = Pmax*dp(j)
s0 = s0+((pj**k)*dcos((pj*r/hbarc))/((mpi*mpi+pj*pj)**s))*regulator(pj)*dwP(j)*dpj
!s0 = s0+((pj**k)*dcos((pj*r/hbarc))/((mpi*mpi+pj*pj)**s))*dwP(j)*dpj
end do

Bks =s0
!Bks = s0*Pmax

end function

! Second Derivative of I_Lambda
real(8) function d2I_Lambda(r)
implicit none
real(8)::r

!d2I_Lambda = (dI_Lambda(r-2.d0*dh)-8.d0*dI_Lambda(r-dh)+8.d0*dI_Lambda(r+dh)-dI_Lambda(r+2.d0*dh) )/(12.d0*dh)


! High order derivative method
!d2I_Lambda = ((469.d0/90.d0)*I_Lambda(r)-(223.d0/10.d0)*I_Lambda(r+dh)+(879.d0/20.d0)*I_Lambda(r+2.d0*dh))/dh**2
!d2I_Lambda = d2I_Lambda + (-1.d0*(949.d0/18.d0)*I_Lambda(r+3.d0*dh)+41.d0*I_Lambda(r+4.d0*dh))/dh**2
!d2I_Lambda =d2I_Lambda +(-1.d0*(201.d0/10.d0)*I_Lambda(r+5.d0*dh)+(1019.d0/180.d0)*I_Lambda(r+6.d0*dh))/dh**2
!d2I_Lambda = d2I_Lambda -(7.d0/10.d0)*I_Lambda(r+7.d0*dh)/dh**2
!d2I_Lambda = d2I_Lambda/(dh**2)

!d2I_Lambda = 2.d0*I_Lambda(r)-5.d0*I_Lambda(r+dh)+4.d0*I_Lambda(r+2.d0*dh)-I_Lambda(r+3.d0*dh)
!d2I_Lambda = d2I_Lambda/(dh**2)

! Second derivative with correct units
!d2I_Lambda = (8.d0*pi/r**3)*Aks(1,2,r)-(4.d0*pi/r)*(Aks(3,2,r)/(hbarc**2))-(8.d0*pi/r**2)*(Bks(2,2,r)/hbarc)
!d2I_Lambda = (8.d0*pi/r**3)*A12(r)-(4.d0*pi/r)*(A32(r)/(hbarc**2))-(8.d0*pi/r**2)*(B22(r)/hbarc)

if(r.le.(smallMomentum)) then
d2I_Lambda = (8.d0*pi/r**3)*Aks(1,2,r)-(4.d0*pi/r)*(Aks(3,2,r)/(hbarc**2))-(8.d0*pi/r**2)*(Bks(2,2,r)/hbarc)
else
d2I_Lambda = (8.d0*pi/r**3)*A12(r)-(4.d0*pi/r)*(A32(r)/(hbarc**2))-(8.d0*pi/r**2)*(B22(r)/hbarc)
end if



! This is the analytical second derivative
!d2I_Lambda = mpi*pi*pi*dexp(-mpi*r/(hbarc))*(1.d0/hbarc**3)

end function

! This is the regulated Yukawa function that is used by Sauri
real(8) function Y_Lambda(r)
implicit none
real(8)::r

!Y_Lambda = (1.d0/(2.d0*pi*pi*r))*Aks(1,1,r)

Y_Lambda = (1.d0/(2.d0*pi*pi*r))*A11(r)

!Y_Lambda = (1.d0/(2.d0*pi*pi*r))*(dExp(-mpi*r/(hbarc))*(0.5d0*pi))

! The INF-cut-off limit
!Y_Lambda = dexp(-mpi*r/(hbarc))/(4.d0*pi*r)

end function

! This is the regulated I_Lambda(r) function used by Sauri
real(8) function I_Lambda(r)
implicit none
real(8)::r

if(r.le.smallMomentum) then
I_Lambda = (4.d0*pi/r)*Aks(1,2,r)
else
I_Lambda = (4.d0*pi/r)*A12(r)
end if

!=============================================================
! The infinite limit
!I_Lambda = (4.d0*pi/r)*dExp(-mpi*r/(hbarc))*((pi*r)/(4.d0*hbarc*mpi))
!I_Lambda=((pi**2)/(mpi*hbarc))*Dexp(-(mpi*r)/hbarc)
!=============================================================

!print *,r,I_Lambda ,((pi**2)/(mpi*hbarc))*Dexp(-(mpi*r)/hbarc)

end function

! The partial dervative with respect to r of I_L(r)
real(8) function dI_Lambda(r)
implicit none
real(8)::r

!dI_Lambda = (I_Lambda(r+dh)-I_Lambda(r))/dh
!dI_Lambda = (I_Lambda(r-2.d0*dh)-8.d0*I_Lambda(r-dh)+8.d0*I_Lambda(r+dh)-I_Lambda(r+2.d0*dh) )/(12.d0*dh)

! Higher order finite diff method
dI_Lambda = (-1.d0*I_Lambda(r-3.d0*dh)/60.d0)+I_Lambda(r-2.d0*dh)*(3.d0/20.d0)-(3.d0/4.d0)*I_Lambda(r-dh)
dI_Lambda = dI_Lambda+(3.d0/4.d0)*I_Lambda(r+dh)-(3.d0/20.d0)*I_Lambda(r+2.d0*dh)+(I_Lambda(r+3.d0*dh)/60.d0)
dI_Lambda = dI_Lambda/dh


!dI_Lambda = (-25.d0/12.d0)*I_Lambda(r)+4.d0*I_Lambda(r+dh)-3.d0*I_Lambda(r+2.d0*dh)
!dI_Lambda = dI_Lambda+(4.d0/3.d0)*I_Lambda(r+3.d0*dh)-(1.d0/4.d0)*I_Lambda(r+4.d0*dh)
!dI_Lambda = dI_Lambda/dh

!dI_Lambda = -((4.d0*pi)/r**2)*Aks(1,2,r)+(4.d0*pi/r)*Bks(2,2,r)

!=============================================================
! This is the analytical first derivative
!dI_Lambda = -1.d0*pi*pi*dexp(-(mpi*r)/hbarc)*(1.d0/hbarc**2)
!=============================================================

end function

! First Derivative of I_Lambda
real(8) function dY_Lambda(r)
implicit none
real(8)::r

!dY_Lambda = (Y_Lambda(r+dh)-Y_Lambda(r))/dh
!dY_Lambda = (Y_Lambda(r-2.d0*dh)-8.d0*Y_Lambda(r-dh)+8.d0*Y_Lambda(r+dh)-Y_Lambda(r+2.d0*dh) )/(12.d0*dh)
!dY_Lambda = (1.d0/(2.d0*pi*pi*r))*Bks(2,1,r)-(1.d0/(2.d0*pi*pi*r**2))*Aks(1,1,r)

! higher order derivative
dY_Lambda = (-1.d0*Y_Lambda(r-3.d0*dh)/60.d0)+Y_Lambda(r-2.d0*dh)*(3.d0/20.d0)-(3.d0/4.d0)*Y_Lambda(r-dh)
dY_Lambda = dY_Lambda+(3.d0/4.d0)*Y_Lambda(r+dh)-(3.d0/20.d0)*Y_Lambda(r+2.d0*dh)+(Y_Lambda(r+3.d0*dh)/60.d0)
dY_Lambda = dY_Lambda/dh

!=============================================================
! Exact derivative
!dY_Lambda = (dexp(-mpi*r/(hbarc))/(4.d0*pi*r))*(-1.d0*mpi/hbarc)
!dY_Lambda = dY_Lambda - dexp(-mpi*r/(hbarc))/(4.d0*pi*r**2)
!=============================================================

end function

end module
