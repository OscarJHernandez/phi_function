! These are the q-dependent Seagull current Matrix elements


module qDependent_Seagull_MEC
use matrixElements
use phiFuncParams
implicit none
real*8,parameter::epsabs = 1d-15 !The Requested accuracy for spherical bessel ME
real*8,parameter::epsrel = 1d-15 !The Requested accuracy
integer,parameter:: key=6 !
real*8,parameter::stepsize= 1.d0 ! the step size used to determine the range of the function
real(8),parameter:: absfunc = 1.d-10 ! This is the |f(x)|. Which determines the range of integration condition

real*8:: abserr ! estimate of modulus of absolute error
integer:: neval ! number of integrand evaluations
integer:: ier ! error message

contains



! This program computes the radial matrix element corresponding to the <n1,l1|jn(qr)| n2 l2 > in the Harmonic Oscillator basis
! This program uses the quadpack routine 'qag' to integrate over a finite range to specified tolerance levels set in
! the module mod_quadpack_params.
!The QAG algorithm is a simple adaptive integration procedure. The integration region is divided into subintervals,
! and on each iteration the subinterval with the largest estimated error is bisected. This reduces the overall error rapidly, 
! as the subintervals become concentrated around local difficulties in the integrand. 
real*8 function seagull_radial_me_spherical_bessel(n1,n2,l1,l2,q,v,n,intkey)
implicit none
integer,intent(in):: n1,n2,l1,l2
integer,intent(in)::intkey
integer,intent(in):: n ! Spherical bessel function order
real*8::a ! the lower limit of the integral
real*8::b ! the upper limit of the integral
real*8:: q,v
real*8:: output! the result

! initiallize parameters for function
a = 0.d0


! Determines the range of the function
b = determineMaxRange(n1,n2,l1,l2,n,q,v)

if(intkey.eq.0) then
call quadpackIntegration(n1,n2,l1,l2,n,q,v,a,b,output)
else
call legendreQuadrature(n1,n2,l1,l2,n,q,v,a,b,output)
end if

seagull_radial_me_spherical_bessel = output

end function seagull_radial_me_spherical_bessel


subroutine quadpackIntegration(n1,n2,l1,l2,n,q,v,a,b,output)
integer::n1,n2,l1,l2,n
real(8)::q,v,a,b
real(8)::output,results     
   
   call qag ( wrapper,a, b, epsabs, epsrel, key, results, abserr, neval, ier )
   
   output = results
       
   contains 
   
   real(8) function wrapper(x)
   implicit none
   real(8)::x
   wrapper = func(x,n1,n2,l1,l2,n,q,v)
   end function
   

end subroutine

subroutine legendreQuadrature(n1,n2,l1,l2,n,q,v,a,b,output)
implicit none
integer::n1,n2,l1,l2,n
real(8)::q,v,a,b
real(8)::output,f,s ,x 
integer::i

  s=0.d0
  
  
    do i=1,Nquad2
        x = b*dx2(i)
        f = func(x,n1,n2,l1,l2,n,q,v)
        s = s +   dw2(i)*f*b
    end do
   
   output =s

end subroutine


! The Radial matrix elements written out in the H.O basis
! <1||3*j_n(1/2*q*r)*(1+1/mr)*(Exp(-mr)/(mr)) ||2>
! q is in units of MeV
real*8 function func(x,ni,nj,li,lj,bessel_n,q,v)
implicit none
integer,intent(in):: ni,nj,li,lj,bessel_n
real*8::v,q
real*8,external:: Bessel
real*8::x,alphai,alphaj,y
real(8):: func2
real(8)::mpiMeV

! Convert from MeV to 1/fm
mpiMeV = (mpi/hbarc)


alphai = dfloat(li)+ 0.5d0
alphaj = dfloat(lj)+0.5d0
y=q*x*0.5d0*(1.d0/hbarc)


    func = 3.d0*(1.d0/(mpiMeV*x))*(1.d0+(1.d0/(mpiMeV*x)))*dexp(-1.d0*mpiMeV*x)
    func = func*BESSEL(bessel_n,y)
    func = func*(LagExp(ni,alphai,2.d0*v*x*x)*Norm(ni,li,v))
    func=func*(LagExp(nj,alphaj,2.d0*v*x*x)*Norm(nj,lj,v))
    func=func*(x**(alphai+alphaj+1.d0))
    
         
end function func

! This function will determine the maximum range of the radial matrix elements before integrating
! over that range.
! uses four points , f1,f2,f3,f4
real(8) function determineMaxRange(ni,nj,li,lj,bessel_n, q, v)
implicit none
integer,intent(in):: ni,nj,li,lj,bessel_n
real*8,intent(in)::v,q
real*8::x1,x2,x3,x4
real*8:: f1,f2,f3,f4
integer::i

! initialize some parameters
f1 = 10.d0
f2 = 10.d0
f3= 10.d0
f4 =10.d0
i=1


do while((abs(f1).gt.absfunc).and.(abs(f2).gt.absfunc).and.(abs(f3).gt.absfunc).and.(abs(f4).gt.absfunc))

x1 = dfloat(i)*stepsize
f1 = abs(func(x1,ni,nj,li,lj,bessel_n,q,v))

x2 = dfloat(i+1)*stepsize
f2 = abs(func(x2,ni,nj,li,lj,bessel_n,q,v))

x3 = dfloat(i+2)*stepsize
f3 = abs(func(x3,ni,nj,li,lj,bessel_n,q,v))

x4 = dfloat(i+3)*stepsize
f4 = abs(func(x4,ni,nj,li,lj,bessel_n,q,v))
!print *,x1,f1!,f2,f3,f4
!print *,x2,f2
!print *,x3,f3
!print *,x4,f4

i=i+4

end do

!print *, f1,f2,f3,f4

determineMaxRange = x4

end function determineMaxRange


end module 
