program test
use phiFunction
use phiFuncParams
use qDependent_Seagull_MEC
use qDependent_PionInFlight_MEC

implicit none
integer,parameter::nu = 2
integer,parameter::sigma=2
integer,parameter::L=2
real(8),parameter::q=1.d0
integer::n1,n2,l1,l2
integer,parameter:: n=0
integer,parameter:: intkey=1
real(8)::Time1,Time2


!Initialize the gau_leg quadrature stuff
call gauleg_init()
l1 =5
l2=5
!CALL CPU_TIME(TIME1) 
!print *, sphericalHankelFunction1(n,z)

print *, 'Nquad1',Nquad1
print *, 'Nquad2',Nquad2
do n1=30,50
do n2=30,50
!print *, 'seagull', n1,n2,l1,l2,seagull_radial_me_spherical_bessel(n1,n2,l1,l2,q,v,n,intkey)
!print *, 'seagull', n1,n2,l1,l2,seagull_radial_me_spherical_bessel(n1,n2,l1,l2,q,v,n,intkey-1)
print *, 'inFlight',n1,n2,l1,l2,pionInFlight_radial_me_spherical_bessel(n1,n2,l1,l2,q,v,nu,sigma,L,intkey)
print *, 'inFlight',n1,n2,l1,l2,pionInFlight_radial_me_spherical_bessel(n1,n2,l1,l2,q,v,nu,sigma,L,intkey-1)
end do
end do
!CALL CPU_TIME(TIME2)

print *, 'TIME', (Time2-Time1)


end program test
