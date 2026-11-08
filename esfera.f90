      SUBROUTINE PROJINIT(N,X)
C     PROYECTA X A LA ESFERA UNIDAD
      INTEGER N,J
      DOUBLE PRECISION X(N),S
      S = 0.0D0
      DO 10 J = 1,N
          S = S + X(J)**2
   10 CONTINUE
      S = DSQRT(S)
      DO 20 J = 1,N
          X(J) = X(J)/S
   20 CONTINUE
      RETURN
      END

      SUBROUTINE RGRAD(N,X,G,RG,NPROB)
C     GRADIENTE RIEMANNIANO = PROYECCION TANGENTE DEL EUCLIDEO
      INTEGER N,NPROB,J
      DOUBLE PRECISION X(N),G(N),RG(N),S
      CALL GRDFCN(N,X,G,NPROB)
      S = 0.0D0
      DO 10 J = 1,N
          S = S + X(J)*G(J)
   10 CONTINUE
      DO 20 J = 1,N
          RG(J) = G(J) - S*X(J)
   20 CONTINUE
      RETURN
      END

      SUBROUTINE RETRACT(N,X,XI,T,XNEW)
C     XNEW = (X + T*XI)/NORM
      INTEGER N,J
      DOUBLE PRECISION X(N),XI(N),XNEW(N),T,S
      S = 0.0D0
      DO 10 J = 1,N
          XNEW(J) = X(J) + T*XI(J)
          S = S + XNEW(J)**2
   10 CONTINUE
      S = DSQRT(S)
      DO 20 J = 1,N
          XNEW(J) = XNEW(J)/S
   20 CONTINUE
      RETURN
      END
