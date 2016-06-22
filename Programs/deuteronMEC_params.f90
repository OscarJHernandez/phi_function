! In this routine, we will define all of the parameters that are needed by the deuteron
! MEC code.

module deuteronMEC_params
implicit none

! Here we define external functions
integer,external::dimensionOfMatrix


! Now we have arrays that store the Hamiltonian
real(8),allocatable:: H0(:,:) ! This Hamiltonian contains just the Ground states of the deuteron
real(8),allocatable:: psi0(:) ! This is the ground state of the system
real(8),allocatable:: eigGS(:) ! The array that will hold the Eigenvalues of the HGS
real(8),allocatable::eigES(:) ! Stores the excited states Eigenvalues
real(8),allocatable:: HE(:,:) ! This Hamiltonian contains the Excited states of the Deuteron
real(8),allocatable::Udagger(:,:) ! The Unitary matrix that diagonalizes the HE matrix
real(8):: E0 ! The ground state of the Deuteron

! Here we define the physical constants for Program
real(8),parameter :: hbarc=197.327053d0 ! MeV*fm
real(8),parameter:: alpha = 1.d0/137.035999139d0
real(8),parameter :: mp=938.27231D0 ! [Mev]
real(8),parameter :: mn=939.56563D0 ! [Mev]
real(8):: gp = 5.585694713d0 ! g factor of the proton
real(8):: gn = -3.82608545d0 ! g factor of the neutron
real(8):: pi
real(8):: M_N_reduced ! [Mev] the reduced mass of the proton and neutron


integer::st
character(50)::word,potential
character*50::input ! The name of the input file
character(50)::string1,string2,string
integer:: Nmax ! The parameter Nmax which is taken from the User
integer::nPolyMax ! The maximum order of the Laguerre Polynomial
integer::lPolyMax ! The maxim alpha of the associated Laguerre Polynomial
real(8):: hw ! Hbar Omega parameter that is given by the user
real(8)::v
integer:: jjMax ! Twice the Maximum needed Jmax for the program
integer::matrixE ! The dimension of ground state matrix for the deuteron
integer::matrixE_oper ! The dimensions of the excited state matrix for the detueron

integer,parameter:: ssMax =2 ! The maximum allowed spin =0,1
integer,parameter:: tMax =1 

! Here we set the deuteron ground state parameters
integer,parameter:: ssMax_0=2  ! The Spin can be either 0->1, but will be restricted to just s=1 by ssfixed
integer,parameter:: jjMax_0=2 ! jj=2, which means j=1 
integer,parameter:: lParity_restrict=1  ! The parity of the angular component of the wave function is even
integer,parameter:: ssfixed=2   ! The spin of the deuteron is fixed to s=1
integer,parameter:: jjfixed=2   ! The total angular momentum j of the deuteron is fixed to j=1
integer,parameter:: tMax_0=0
integer,allocatable:: QGS(:,:) ! This array will store the quantum numbers as a table

integer,parameter::jjfixed_E=-1 ! jj is allowed to vary
integer,parameter::ssfixed_E=-1 ! ss is allowed to vary
integer,parameter:: lParity_restrict_E=0 ! No restiction on the parity of states other than they have to be fermionic
integer,allocatable:: QES(:,:) ! This array will store the quantum numbers as a table

integer:: n ! Used for Lappack subroutine

integer:: constants ! The parameter that defines which constants to use in MEC currents [0=Friar, 1= Sauri]
integer::NLO_current ! Do you want the NLO-current? No=0, NLO(Infty)=1, NLO(Sauri)=2, NLO(Lambda)-Friar =3
integer::DELTA_current ! no=0, yes =1 (Friars Delta Current)
real(8)::cutOff ! The cut-off for the currents
integer::NquadMomentum ! These are the number of quadrature points to use in p
integer::NquadPosition ! This r-space quadrature points
real(8)::rMax
real(8)::pMax


contains


! This subroutine sets all of the physical constants
subroutine setConstants()
pi = datan(1.d0)*4.d0
M_N_reduced = (mp*mn)/(mp+mn)
end subroutine


! This subroutine reads the input file
subroutine readInput(fname)
implicit none
character(50)::fname

    ! Now we open the input from the file: H_input.txt
    open(unit =1, file= fname,iostat=st)
    
    ! Reads from an input file
    read(1,*) word, Nmax
    read(1,*) word,jjMax
    read(1,*) word,input
    read(1,*) word,potential
    read(1,*) word,hw
    read(1,*) word,constants
    read(1,*) word,NLO_current
    read(1,*) word,DELTA_current
    read(1,*) word,cutOff
    read(1,*) word,NquadMomentum
    read(1,*) word,NquadPosition
    read(1,*) word,Rmax
    read(1,*) word,Pmax
    close(1)
    
    ! Now we determine what the HO parameter v is in the correct units
    v= (M_N_reduced/(hbarc**(2)))*0.5d0*hw


end subroutine



! This subroutine initializes the size of the ground state matrix
subroutine initializeGroundState()

    ! Calculate the dimension of the ground state matrix
    matrixE = dimensionOfMatrix(matrixE,Nmax,jjMax_0,ssMax_0,tMax_0,jjfixed,ssfixed,lParity_restrict) 
    
    ! Allocate memory for the table of Quantum numbers 
    allocate(QGS(matrixE,6))
        
    ! Fill the Q_table with quantum numbers with Ground state restrictions
    call generateTable(QGS,matrixE,Nmax,jjMax_0,ssMax_0,tMax_0,jjfixed,ssfixed,lParity_restrict)
    
    ! Determine the Max order of n an L
    nPolyMax = Nmax/2 ! max order of N
    lPolyMax = (jjMax+2)/2 ! Max order of L
    
    ! determine the Max order of l LL=JJ+S
    
    ! Allocate memory for the Ground state Matrix
    allocate(H0(matrixE,matrixE))
    
    ! Allocate memory for the set of Eigenvalues
    allocate(eigGS(matrixE))

end subroutine

! This subroutine initializes the size of the Excited state matrix elements 
subroutine initializeExcitedStates()

        ! Calculate the Dimensions of the Excited State Matrix
        matrixE_oper = dimensionOfMatrix(matrixE_oper,Nmax,jjMax,ssMax,tMax,jjfixed_E,ssfixed_E,lParity_restrict_E) 
    
        print *, 'Length of the New Vector: ',matrixE_oper
        
        ! Now we allocate Memory for the Quantum Numbers of the new Vector
        allocate(QES(matrixE_oper,6))
        
        ! We Generate the table of new Quantum Numbers (Q_table2)
        call generateTable(QES,matrixE_oper,Nmax,jjMax,ssMax,tMax,jjfixed_E,ssfixed_E,lParity_restrict_E)
        
        ! Now we Generate a larger Matrix H
        print *, 'Generating ', matrixE_oper, ' x ', matrixE_oper , 'matrix'
        allocate(HE(matrixE_oper,matrixE_oper))
        allocate(eigES(matrixE_oper))
        
        ! Initialize the Matrix
        HE(:,:)=0.0
    

end subroutine

subroutine readInESMatrixElements()
implicit none
character(50)::input_ES
real(8)::me
integer::i,j,z

write(string,"(I3.3)") Int(hw)
write(string2,"(I3.3)") Int(Nmax)

input_ES = "H_matrix"//"_"//trim(potential)//"_ES_hw"//trim(string)//"_Nmax"//trim(string2)//".txt"
print *, 'Input File to be Read for ES: ', input_ES

open(1,file=input_ES,status='unknown')
read(1,*) Nmax,hw

do z=1,matrixE_oper*(matrixE_oper+1)/2
read(1,*), i,j,me
HE(i,j)= me
end do

close(1)



end subroutine


subroutine readInGSMatrixElements()
! Auxilary variables used for reading in files
implicit none
character(50)::input_GS
real(8)::me
integer::i,j,z


write(string,"(I3.3)") Int(hw)
write(string2,"(I3.3)") Int(Nmax)
   
input_GS = "H_matrix"//"_"//trim(potential)//"_GS_hw"//trim(string)//"_Nmax"//trim(string2)//".txt"
print *, 'Input File to be Read for GS: ', input_GS
open(1,file=input_GS,status='unknown')
    
    ! this line overwrites hw and Nmax
    read(1,*) Nmax,hw
    
    do z=1,matrixE*(matrixE+1)/2
    read(1,*), i,j,me
    H0(i,j)= me
    end do
    
    close(1)

end subroutine


! This subroutine copies the ground state of H, to the array psi0
subroutine extractGS()
implicit none
integer::i

allocate(psi0(matrixE))
psi0(:) = 0.d0

do i=1,matrixE
psi0(i)= H0(i,1)
end do

! Extract the Ground State 
E0 = eigGS(1)


end subroutine


!---------------------------------------------------------!
!Calls the LAPACK diagonalization subroutine DSYEV        !
!input:  a(n,n) = real symmetric matrix to be diagonalized!
!            n  = size of a                               !
!output: a(n,n) = orthonormal eigenvectors of a           !
!        eig(n) = eigenvalues of a in ascending order     !
!---------------------------------------------------------!
SUBROUTINE diagonalize(M,n,W)
INTEGER :: info, lwork, lwmax = 2000000
integer::n
REAL*8 :: M(n,n)
real*8::W(n)
REAL*8, DIMENSION(:), ALLOCATABLE :: work
INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV
!REAL*8, DIMENSION(:), ALLOCATABLE :: W !Vector storing eigenvalues
ALLOCATE(work(lwmax))
ALLOCATE(IPIV(n))
!ALLOCATE(W(n))

! Optimizing lwork
  LWORK = -1
  CALL dsyev( 'Vectors', 'Upper', n, M, n, W, work, lwork, info )
  lwork = min( lwmax, int( work( 1 ) ) )
  CALL dsyev( 'Vectors', 'Upper', n, M, n, W, work, lwork, info )

  IF( info.GT.0 ) THEN
  
  print *, 'Algorithm Failed to Converge'
  
  else if(info.LT.0) THEN
  
  print *, 'the ith argument had an illegal Value', INFO*(-1)
  
  END IF

DEALLOCATE(IPIV,work)
RETURN
END SUBROUTINE diagonalize



! Implementation of matrix Vector multiplication
! Multiplies an A(n by m ) matrix by a m by 1 vector, and produces C(n by 1)
Subroutine MatrixVectorMultiplication(A,n,m,V,C)
implicit none
integer:: n,m,i,j
real*8:: A(n,m)
real*8:: V(m)
real*8:: C(n)
real*8:: s


  do i=1, n
    s=0.0
     do j=1, m
     s = s + A(i,j)*V(j) 
     end do
     
     C(i) =s
  end do 


end subroutine MatrixVectorMultiplication


end module
