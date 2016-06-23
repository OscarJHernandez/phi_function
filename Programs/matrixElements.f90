! This modules will contain all of the pieces needed to generate the matrix elements of an Operator
! < 1 || O || 2 >
module matrixElements
use phiFuncParams
implicit none

real(8),allocatable::daP(:),dbP(:),dp(:),dwP(:),eP(:)
real(8),allocatable::daR(:),dbR(:),dr(:),dwR(:),eR(:)
integer,parameter:: ipolyMomentum=2 !  Legendre polynomial from 0 to 1 
integer,parameter:: ipolyPosition=2 !  Legendre polynomial from 0 to 1 

! Here we have arrays that store calculated points
real(8),allocatable:: LagExpArray(:,:,:) ! The array that stores the associated Laguerre Exp(-2*v*x*x)L(n,l+1/2,x) polynomials
real(8),allocatable:: SphericalBessel(:,:) ! Array that stores the spherical Bessel functions evaluated at quadrature points
real(8),allocatable:: PhiFunction(:,:,:,:) ! Array that stores the Phi functions evaluated at the quadrature points
real(8),allocatable:: RadialSeagullME(:,:,:,:,:) ! Array that stores the radial ME for the seagull term

contains 

! Initialize all of the guassian quadrature points and weights for momentum and position
subroutine gauleg_init_1()
implicit none

    allocate(daP(NquadMomentum),dbP(NquadMomentum),dp(NquadMomentum),dwP(NquadMomentum),eP(NquadMomentum))
    allocate(daR(NquadPosition),dbR(NquadPosition),dr(NquadPosition),dwR(NquadPosition),eR(NquadPosition))
    
    al = 0.d0
    
    call drecur(NquadMomentum,ipolyMomentum,al,dbe,daP,dbP,iderr)
    call dgauss(NquadMomentum,daP,dbP,depsma,dp,dwP,ierr,eP)
    
    call drecur(NquadPosition,ipolyPosition,al,dbe,daR,dbR,iderr)
    call dgauss(NquadPosition,daR,dbR,depsma,dr,dwR,ierr,eR)

end subroutine


! This subroutine initializes all of the Norm(n,l)*Laguerre*Exp(-0.5*v*r*r)*r**(L+1) functions in the basis
subroutine initExpLaguerre()
implicit none
integer::ni,li,i
real(8)::xi,alphai,ri

allocate(LagExpArray(0:nPolyMax,0:lPolyMax,NquadPosition))

do ni=0,nPolyMax
    do li=0,lPolyMax
        do i=1,NquadPosition
            ri = rMax*dr(i)
            xi = 2.d0*v*ri*ri
            alphai = 0.5d0+dfloat(li)
            LagExpArray(ni,li,i) = Norm(ni,li,v)*LagExp(ni,alphai,xi)*(ri**(li+1))
        end do
    end do
end do


end subroutine


! This function tests the orthognonality of the Basis functions, with respect to n1 and n2
real(8) function radialMeTest(n1,l1,n2,l2)
implicit none
integer::n1,n2,l1,l2
real(8)::s
integer::i

s =0.d0

do i=1,NquadPosition
s = s+LagExpArray(n1,l1,i)*LagExpArray(n2,l2,i)*dwR(i)
end do

s = rMax*s

radialMeTest =s

end function



! Calculates the radial matrix elements of the operator r^p 
! Note, we must divide the operator by r/2 
real(8) function radial_matrixElement_rp(n1,l1,n2,l2,k)
implicit none
integer::n1,n2,l1,l2,k
real(8)::ri
real(8)::s
integer::i

s=0.d0

do i =1,NquadPosition
ri = rMax*dr(i)
s = s+LagExpArray(n1,l1,i)*LagExpArray(n2,l2,i)*dwR(i)*((ri/2.d0)**k)
end do

s=rMax*s

radial_matrixElement_rp =s

end function





! This function evaluates the associated Laguerre polynomial with an exponential at a given point.
! This function behaves quite well even for large values of x. The exponential helps keep the polynomial from
! overflow problems.
! LagExp = Exp(-1/2*x)*Lag[n,alpha,x]
real*8 function LagExp(n,alpha,x)
implicit none
integer,intent(in):: n
real*8,intent(in):: alpha
real*8,intent(in):: x
real*8::L,Lm,Lp
integer:: i 

! Initialize the values
Lm = 1.d0*dexp(-0.5d0*x)
L = (1.d0+ alpha -x)*dexp(-0.5d0*x)
Lp = 0.d0


    if(n.eq.0) then
        Lp = Lm
    else if (n.eq.1) then
        Lp = L
    else
    
        do i=1,n-1
        Lp = ((2.d0*dfloat(i)+1.d0+alpha-x)/(dfloat(i)+1.d0)*L - ((dfloat(i)+alpha)/(dfloat(i)+1.d0))*Lm)
        Lm = L
        L = Lp
        end do
    
    end if
    
    LagExp = Lp

end function

! This is the Normalization coefficients for the harmonic Oscillator
! This function computes the normalization for the Harmonic Oscillator wave function
! N(n,l,v) = sqrt( Gamma(n+1)/Gamma(n+l+3/2) * 2 *(2*v)^(l+3/2))
real*8 function Norm(n,l,v)
implicit none
integer:: n,l,i
integer:: Nm
real*8::v
real*8:: frac
real*8:: pi


pi = datan(1.d0)*4.d0


frac = 1.d0/dsqrt(pi)

do i=1,n
    frac =frac*(dfloat(i)/(dfloat(i)-0.5d0))
end do

Nm = n+l+1

do i=n+1,Nm
    frac = frac*(1.d0/(dfloat(i) -0.5d0))
end do

frac = frac*2.d0*((2.d0*v)**(dfloat(l)+1.5d0))
Norm = dsqrt(frac)

return
end function Norm




end module

