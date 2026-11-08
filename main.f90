      PROGRAM RSD
      INTEGER N,NPROB,K,J,ITER
      LOGICAL CONV
      DOUBLE PRECISION X(50),XN(50),G(50),RG(50),F,FN,T,GN2,XNORM
      N = 3
      NPROB = 1
      CALL INITPT(N,X,NPROB,1.0D0)
      CALL PROJINIT(N,X)
      CONV = .FALSE.
      ITER = 5000
      DO 100 K = 1,5000
          CALL OBJFCN(N,X,F,NPROB)
          CALL RGRAD(N,X,G,RG,NPROB)
          GN2 = 0.0D0
          DO 10 J = 1,N
              GN2 = GN2 + RG(J)**2
   10     CONTINUE
          IF (DSQRT(GN2).LT.1.0D-6) THEN
              CONV = .TRUE.
              ITER = K - 1
              GO TO 200
          END IF
          DO 15 J = 1,N
              G(J) = -RG(J)
   15     CONTINUE
          T = 1.0D0
   20     CONTINUE
          CALL RETRACT(N,X,G,T,XN)
          CALL OBJFCN(N,XN,FN,NPROB)
          IF (FN.LE.F - 1.0D-4*T*GN2) GO TO 40
          T = 0.5D0*T
          IF (T.GT.1.0D-16) GO TO 20
          ITER = K - 1
          GO TO 200
   40     CONTINUE
          DO 50 J = 1,N
              X(J) = XN(J)
   50     CONTINUE
  100 CONTINUE
  200 CONTINUE
      XNORM = 0.0D0
      DO 60 J = 1,N
          XNORM = XNORM + X(J)**2
   60 CONTINUE
      XNORM = DSQRT(XNORM)
      IF (CONV) THEN
          WRITE(*,*) 'CONVERGIO'
      ELSE
          WRITE(*,*) 'NO CONVERGIO'
      END IF
      WRITE(*,'(A,I6)') ' ITER = ',ITER
      WRITE(*,'(A,ES16.8)') ' F = ',F
      WRITE(*,'(A,ES16.8)') ' |GRAD| = ',DSQRT(GN2)
      WRITE(*,*) 'X ='
      WRITE(*,'(5F16.8)') (X(J),J=1,N)
      WRITE(*,'(A,F16.8)') ' NORMA DE X = ',XNORM
      END
