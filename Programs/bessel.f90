!--- spherical bessel  functions of first kind j_l(qr) ---
        REAL(8) FUNCTION BESSEL(L,X)
! THE SPH. BESSEL FUNCTION
        IMPLICIT real(8) (A-H,O-Z)
        IF(X.EQ.0.d0.AND.L.EQ.0) THEN
          BESSEL=1.d0
          RETURN
        ENDIF 
        IF(X.EQ.0.d0.AND.L.GT.0) THEN
          BESSEL=0.d0
          RETURN
        ENDIF
        IF(X.GT.L+1.D0) GO TO 2
        EPS=1.D-14
        F=1.d0  
        DO I=1,L
           F=F*X/(2.d0*I+1.d0)
        enddo
        S=1.d0
        T=1.d0
        D=0.5d0*X*X
        DO 4 I=1,40
           T=-T*D/(I*(2.d0*(L+I)+1.d0))
           S=S+T
           IF(ABS(T/S).LT.EPS)THEN
              BESSEL=F*S
              RETURN
           ENDIF
4       enddo
2       continue
        A=SIN(X)/X
        IF(L.EQ.0)THEN  
        BESSEL=A
        RETURN
        ENDIF
        B=(A-COS(X))/X
        IF(L.EQ.1) THEN 
        BESSEL=B
        RETURN
        ENDIF
        DO I=2,L
           C=B/X*(2.d0*I-1.d0)-A
           A=B
           B=C
        enddo
        BESSEL=C  
        RETURN
        END FUNCTION BESSEL
