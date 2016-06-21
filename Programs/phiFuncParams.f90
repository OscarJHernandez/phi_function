module phiFuncParams
implicit none
real(8),parameter::hbarc = 197.327d0
real(8),parameter::mpi=140.d0 ! Mass of pion
real(8),parameter::pi= dAtan(1.d0)*4.d0 ! Pi
integer,parameter:: Nquad1 = 1000 ! the number of quadrature points for Legendre quad for integrating over momentum
real*8:: da(Nquad1),db(Nquad1),dx(Nquad1),dw(Nquad1),e(Nquad1)
integer,parameter:: ipoly1=2 !  Legendre polynomial from 0 to 1 
real*8,parameter:: depsma=1.0d-18
real*8:: al,ierr,iderr,dbe

contains

! Initialize all of the guassian quadrature points and weights
subroutine gauleg_init()
implicit none

al = 0.d0

call drecur(Nquad1,ipoly1,al,dbe,da,db,iderr)
call dgauss(Nquad1,da,db,depsma,dx,dw,ierr,e)

end subroutine



end module
