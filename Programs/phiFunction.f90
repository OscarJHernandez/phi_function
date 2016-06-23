module phiFunction
use phiFuncParams
contains



! This is the spherical Hankel function of the first kind
Complex*16 function sphericalHankelFunction1(l,cz)
implicit real(8) (a,b,d-h,o-z)
implicit complex*16 (c)
complex*16::cz
integer::ll
integer::l
integer::n,ifail,kfn,mode1
      
      
      cl=dcmplx(1.0d0*ll,0.0d0)
      ceta=dcmplx(0.0d0,0.0d0)
      ifail=1
      n=1
      kfn=1
      mode1=1

     call coulcc(cz,ceta,cl,n,cf,col,cfp,colp,csig,mode1,kfn,ifail)
     
sphericalHankelFunction1 = cf+dcmplx(0.d0,1.d0)*col

end function

! This is the Legendre polynomial
real(8) function LegendreP(n,x)
integer::n
real(8)::x
real ( kind = 8 ) pd(0:n+2)
real ( kind = 8 ) pn(0:n+2)

call lpn( n+2, x, pn, pd )
LegendreP = pn(n)

end function


! This function works only for nu+sigma=even,L=even 
! These functions appear to work accurately, for r<= 23.0 fm
! The units of the Phi(q,r) functions are MeV^(nu-1)
real(8) function phi(q,r,nu,sigma,L)
implicit none
integer,intent(in)::nu,sigma,L
real(8)::r
real(8)::q
real(8)::s
integer::i


if((mod(nu+sigma,2)==0).and.(mod(L,2)==0)) then
    
    !! Integrate from x=0, to x=1, like in the notes
    s=0.d0
    do i=1,Nquad1
        s = s +   0.5d0*(1.d0/q)*Gtilda(dx(i),q,sigma,nu,r)*dw(i)*LegendreP(L,dx(i))
    end do
    
else    
print *, 'ERROR: PHI FUNCTION ARGUMENTS, NU+SIGMA AND L MUST BE EVEN'
s = 0.d0
end if

phi = s
end function


!!==============================
!! Here we define Gtilda(x,q,m)
!!==============================
real(8) function Gtilda(x,q,sigma,nu,r)
implicit none
real(8)::x 
real(8)::q ! the momentum transferred
integer:: sigma ! 
integer:: nu
complex(16)::s0
complex(16)::s1
complex(16)::s2 
complex(16)::h1,h0
real*8::r
complex(16):: G
complex*16::z

s0 = dCMPLX(0.d0,1.d0)*dsqrt(0.25d0*q*q+mpi*mpi)
s1 = dCMPLX(0.d0, 1.d0)*dsqrt(0.25d0*(q*q)*(1.d0-x*x)+mpi*mpi)
s2 = dCMPLX(0.5d0*q*x,dsqrt(0.25d0*(q*q)*(1.d0-x*x)+mpi*mpi))

z = (s2*r/hbarc)
h1 = sphericalHankelFunction1(sigma,z)


!z = (s0*r/hbarc)
!h0 = sphericalHankelFunction1(sigma,z)

    G = (s2**(nu+1))*(h1/s1)   ! -(s0**nu)*h0
    G = G*dCMPLX(0.d0,1.d0)*(pi/x)

Gtilda = REALPART(G)

    if(isnan(Gtilda)) then
    print *, ''
    print *, 'sigma',sigma
    print *, 's2',s2
    print *, 'z',z
    print *, 'x',x
    print *, 's1',s1
    print *, 'h1',h1
    stop 'Gtilda is NAN'
    
    end if
end function



end module
