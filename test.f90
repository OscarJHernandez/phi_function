program test
use phiFunction
use phiFuncParams
implicit none
integer,parameter::nu = 0
integer,parameter::sigma=0
integer,parameter::L=0
real(8),parameter::dh =1.d0
real(8),parameter::q=500.0d0
real(8)::r
integer::i

!Initialize the gau_leg quadrature stuff
call gauleg_init()
print *, 'r',r
print *, 'q',q
print *, 'Nquad',Nquad1
do i=1,25
r = dh*dfloat(i)
print *, r,phi(q,r,nu,sigma,L), (pi/(4.d0*mpi))*dexp(-mpi*r/hbarc)
print *, r,(phi(q,r,2,0,0)/(hbarc**2)),  (pi/(4.d0*r*hbarc))*dexp(-mpi*r/hbarc)*(2.d0-((mpi*r)/hbarc))
print *, r,(phi(q,r,2,2,0)/(hbarc**2)),  (pi/(4.d0*r*hbarc))*dexp(-mpi*r/hbarc)*(1.d0+((mpi*r)/hbarc))
end do


!print *, phi(q,r,nu,sigma,L), phi_0_LL(q,r,sigma), (pi/(4.d0*mpi))*dexp(-mpi*r/hbarc)
!print *,  phi_0_LL(q,r,sigma), (pi/(4.d0*mpi))*dexp(-mpi*r/hbarc)

!do i=0,50
!p = dh*dfloat(i)
!print *, p,dFunc(q,r,x,nu,sigma,p)
!end do

end program test
