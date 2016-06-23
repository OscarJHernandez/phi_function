module phiFuncParams
implicit none
real(8),parameter::hbarc = 197.327053d0
real(8),parameter::mpi=140.d0 ! Mass of pion
real(8),parameter::pi= dAtan(1.d0)*4.d0 ! Pi
integer,parameter:: Nquad1 = 200 ! the number of quadrature points for Legendre quad for integrating over momentum
real*8:: da(Nquad1),db(Nquad1),dx(Nquad1),dw(Nquad1),e(Nquad1)

integer,parameter:: Nquad2 = 300 ! the number of quadrature points for Legendre quad for integrating over momentum
real*8:: da2(Nquad2),db2(Nquad2),dx2(Nquad2),dw2(Nquad2),e2(Nquad2)


integer,parameter:: ipoly1=2 !  Legendre polynomial from 0 to 1 
real*8,parameter:: depsma=1.0d-18
real*8:: al,ierr,iderr,dbe
real*8:: al2,ierr2,iderr2,dbe2

! From input supposedly 
integer,parameter::NquadMomentum=500 ! These are the number of quadrature points to use in p
integer,parameter::NquadPosition=500 ! This r-space quadrature points
real(8),parameter::rMax=23.d0
real(8),parameter::pMax=5000.d0
integer,parameter::nPolyMax = 10
integer,parameter::lPolyMax = 10
real(8),parameter::hw =8.d0
real(8),parameter :: mp=938.27231D0 ! [Mev]
real(8),parameter :: mn=939.56563D0 ! [Mev]
real(8),parameter:: M_N_reduced = (mp*mn)/(mp+mn)
real(8),parameter::v = (M_N_reduced/(hbarc**(2)))*0.5d0*hw



contains

! Initialize all of the guassian quadrature points and weights
subroutine gauleg_init()
implicit none
integer::i

al = 0.d0

call drecur(Nquad1,ipoly1,al,dbe,da,db,iderr)
call dgauss(Nquad1,da,db,depsma,dx,dw,ierr,e)

call drecur(Nquad2,ipoly1,al2,dbe2,da2,db2,iderr2)
call dgauss(Nquad2,da2,db2,depsma,dx2,dw2,ierr2,e2)
  


end subroutine



end module
