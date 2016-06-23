module phiFuncParams
implicit none
real(8),parameter::pMax =500
real(8),parameter:: m=1.d0
integer,parameter:: Nquad1 = 1000 ! the number of quadrature points for Legendre quad for integrating over momentum
integer,parameter:: Nquad2 = 100 ! The second quadrature point scheme (integral over x)
real*8:: da(Nquad1),db(Nquad1),dx(Nquad1),dw(Nquad1),e(Nquad1)
real*8:: da2(Nquad2),db2(Nquad2),dx2(Nquad2),dw2(Nquad2),e2(Nquad2)
integer,parameter:: ipoly1=2 !  Legendre polynomial from 0 to 1 
integer,parameter:: ipoly2=1 !  Legendre polynomial from -1 to 1 
real*8,parameter:: depsma=1.0d-18
real*8:: al,ierr,iderr,dbe


contains

! Initialize all of the guassian quadrature points and weights
subroutine gauleg_init()
implicit none

al = 0.d0

call drecur(Nquad1,ipoly1,al,dbe,da,db,iderr)
call dgauss(Nquad1,da,db,depsma,dx,dw,ierr,e)

call drecur(Nquad2,ipoly2,al,dbe,da2,db2,iderr)
call dgauss(Nquad2,da,db,depsma,dx2,dw2,ierr,e)

end subroutine



end module
