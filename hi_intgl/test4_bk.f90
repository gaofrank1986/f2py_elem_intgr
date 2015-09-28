
    subroutine compute_coeff_GH(NDIM,NBDM,NODE,NPOWG,XP,XIP,XIQ,  &
     &                    COEFG,COEFH)
        ! set hi_beta to shared variable
        implicit none

        integer ::  NDIM,nBDM,NODE,NPOWG
        real(8) :: drdn

        real(8) ::XP(NDIM),XIP(NBDM),XIQ(NBDM),COSN(NDIM),        &
        &          GCD(NDIM,NBDM),XI(NDIM),RI(NDIM),SHAP(NODE),     &
        &          COEFG(0:NPOWG),COEFC(0:NPW),COEFH(0:NPW)

        integer :: i,n,ip,jp
       
        ! DETERMINE COEFFICIENTS Gn
       
        COEFG=0.D0
       
        CALL compute_coeff_G(NDIM,NBDM,NODE,NPOWG,XP,XIP,XIQ, DRDN,COSN,    &
        &            GCD,RI,SHAP,COEFG)
        
        ! DETERMINE COEFFICIENTS Cn USING Eq. (3-6-37)

        
        COEFC(0)=DSQRT(COEFG(0))
        DO N=1,NPW
            COEFC(N)=0.D0
            DO I=1,N-1
                COEFC(N)=COEFC(N)-COEFC(I)*COEFC(N-I)     
            ENDDO
            IF(N.LE.NPOWG)COEFC(N)=COEFG(N)+COEFC(N)   ! Eq. (3-6-37)
            COEFC(N)=COEFC(N)/(2.D0*COEFC(0))
        ENDDO
        ! DETERMINE COEFFICIENTS Hi USING Eq.(3-6-47a)
        COEFC(1:NPW)=COEFC(1:NPW)/COEFC(0)    ! Eq.(3-6-47b)
        COEFH(0)=1.D0/COEFC(0)
        COEFH(1)=COEFC(1)
        COEFH(2)=2.*COEFC(2)-COEFC(1)*COEFC(1)
        COEFH(3)=3.*COEFC(3)-3.*COEFC(1)*COEFC(2)+COEFC(1)**3
        COEFH(4)=4.*COEFC(4)+4.*COEFC(1)*COEFC(1)*COEFC(2)                &
        &        -4.*COEFC(1)*COEFC(3)-2.*COEFC(2)*COEFC(2)                &
        &        -COEFC(1)**4
   
    end subroutine
