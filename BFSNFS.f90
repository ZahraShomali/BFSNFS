!      Main program
       IMPLICIT NONE

!      Decleration
	   INTEGER           N,NRHS,N1
       PARAMETER         (N=16)
	   PARAMETER         (NRHS=1)

       INTEGER          LWORK,LRWORK,WORK

       INTEGER          NMAX
       PARAMETER        (NMAX=16)
       INTEGER          LDA
       PARAMETER        (LDA=NMAX)
!      .. Local Scalars ..
       INTEGER          INFO,INFO1,INFO2
!      .. Local Arrays ..
       COMPLEX *16     A(LDA,NMAX), B(NMAX,1)
       COMPLEX *16      C(2,2),D(2,1)
       INTEGER          IPIV(NMAX)

	   PARAMETER         (N1=16)
       INTEGER           LWORK2,LRWORK2
       PARAMETER         (LWORK2=2*N1,LRWORK2=2*N1)
       INTEGER           LWORKhat2
       PARAMETER         (LWORKhat2=4*2*N1)

       COMPLEX*16        VL(1,1),VR(N1,N1),W(N1),WORK2(LWORK2)
       COMPLEX*16        WORK3(LWORKhat2)
       DOUBLE PRECISION  RWORK2(LRWORK2)

!	   Decleration for Lapack library ZGSVD
	   INTEGER LDA1
       PARAMETER       (LDA1=NMAX)
       DOUBLE PRECISION S(NMAX)
       INTEGER LWORK1
       PARAMETER (LWORK1=2*N+N)
       COMPLEX*16 WORK1(LWORK1)
       DOUBLE PRECISION RWORK1(5*NMAX)

       INTEGER         LDU
       PARAMETER       (LDU=NMAX)
       COMPLEX*16 U(LDU,N)

       INTEGER         LDVT
       PARAMETER      (LDVT=NMAX)
       COMPLEX*16 VT(LDVT,N)


!      .. Local Arrays ..

      COMPLEX *16      B2(NMAX,1),A2(LDA1,NMAX),A3(LDA1,NMAX)
   
      EXTERNAL ZGESVD
	  EXTERNAL ZGESV
 	  EXTERNAL ZGETRF
 	  EXTERNAL ZGEEV

      
      COMPLEX*16 q1uph,q1upe,q1downh,q1downe,gamma1mupp,gamma1mupm,&
gamma1mdownp,gamma1mdownm,gamma1pupp,gamma1pupm,gamma1pdownp,&
gamma1pdownm,u1up,u1down,v1up,v1down,gamma1mbarupp,&
gamma1mbarupm,gamma1mbardownp,gamma1mbardownm,gamma1pbarupp,&
gamma1pbarupm,gamma1pbardownp,gamma1pbardownm,reup,redown,&
rhup,rhdown,teup,tedown,thup,thdown,c3p,c4p,Det
      REAL km2,q1upe2,q1uph2,q1downe2,q1downh2,q2upe2,q2uph2,q2downe2,&
q2downh2,kup,kdown
	
	  COMPLEX*16  kp,km,ky

	  COMPLEX*16 q2uph,q2upe,q2downh,q2downe,gamma2mupp,gamma2mupm,&
gamma2mdownp,gamma2mdownm,gamma2pupp,gamma2pupm,gamma2pdownp,&
gamma2pdownm,u2up,u2down,v2up,v2down,gamma2mbarupp,&
gamma2mbarupm,gamma2mbardownp,gamma2mbardownm,gamma2pbarupp,&
gamma2pbarupm,gamma2pbardownp,gamma2pbardownm,Del1up,Del1down,&
Del2up,Del2down,h1vEf,h2vEf
	

      COMPLEX*16 h1Ef,h2Ef,i,VL1(16),VL2(16),VE1(16,1),VE2(16,1),E1,E2
	  COMPLEX*16 VEN1(16,1),VEN2(16,1)
	  INTEGER j,h,k,kk,rr,t,tt,ff,jg,tg,hhh,kjj,ee,ISTEP,iii,s1,s2,NNN,&
ISTEP2,kjjj,uy,ww

 	  DOUBLE PRECISION phi,alpha,pi,abskp,abskm,L,Isx(50,50),Isy(50,50),&
Isz(50,50),Meanrho1(160,100),Meanrho2(160,100),Meanrho(160,100)
      DOUBLE PRECISION theta1up,theta1down,theta2up,theta2down
	  DOUBLE PRECISION errIsx,errIsy,errIsz,IEnsz(50,50),IEnsy(50,50),&
IEnsx(50,50),scale,errE,errEE,diffE,errorE,SSS(16),&
SS(100000,1),Isxlnfs1,xln
      DOUBLE PRECISION EvEf,tolInsz,tolInsx,Icharge(50,50),Ic,&
Ich(50,50),tolInsy,E0VEF,EFVEF,X0,DELTA,DELE(1000),E2vEf,STT,&
STTx,STTy,STTz,STTTn,STTy2,STTz2
      COMPLEX*16  FS2e1,FS2e2,FS2e3,FS2e4,sttn(1000),sttnx(1000),&
sttny(1000),sttnz(1000),sttny2(1000),sttnz2(1000)
	

      OPEN (UNIT=1,FILE='EvEf.txt',STATUS='unknown')
      OPEN (UNIT=2,FILE='E2vEf.txt',STATUS='unknown')
	  OPEN (UNIT=3,FILE='taux.txt',STATUS='unknown')
      OPEN (UNIT=4,FILE='tauy.txt',STATUS='unknown')
      OPEN (UNIT=5,FILE='tauz.txt',STATUS='unknown')

      i=(0.0,1.0)
	  pi=3.14159265358979D0

!     .. Errors

      errIsx=0.0001D0
	  errIsy=0.0001D0
      errIsz=0.0001D0
	  errE=0.0001D0
      errEE=0.0001D0
      errorE=0.00001D0
      do 8888 ww=1,800
        write (237,*) ((0.002D0*ww))
8888  end do
!
      DO 33 kjjj=1,800
 	  DO 20 t=2,2
	  alpha=0.75D0*(t-1)*pi
      stt=(0.0,0.0)
      DO 10 tt=2,2
	    write(1,*)0.05D0*(tt-1)
	    phi=0.25D0*(tt-1)*pi
	    DO 3 kjj=1,100,1
    	xln=0.2D0*kjjj

	    Isx(t,tt)=0.0D0
    	FS2e1=(0.0,0.0)
        FS2e2=(0.0,0.0)
        FS2e3=(0.0,0.0)
        FS2e4=(0.0,0.0)

	    IEnsx(t,tt)=0.0D0
        Isy(t,tt)=0.0D0
	    IEnsy(t,tt)=0.0D0
        Isz(t,tt)=0.0D0
	    Isxlnfs1=0.0D0
	    IEnsz(t,tt)=0.0D0
	    Ich(t,tt)=0.0D0
	    Icharge(t,tt)=0.0D0
	    rr=1
1       IF (rr.EQ.1) THEN

          ky=complex(0.0-0.01*(kjj-1),0.0)
	      h1vEf=complex(0.1,0.0)
	      h2vEf=complex(0.1,0.0)
	      Del1up=(0.01,0.0)
	      Del1down=(0.01,0.0)
	      Del2up=(0.01,0.0)
	      Del2down=(0.01,0.0)
	      theta1down=0.0D0*pi
	      theta1up=1.0D0*pi
	      theta2down=.0D0*pi
	      theta2up=1.0D0*pi

	      L=0.2D0

	      kup=(1.0,0.0)
	      kdown=(1.0,0.0)
	      ISTEP=0
	      SSS(1)=1.0D0
          E2vEf=0.0D0

	      DO WHILE (SSS(1).GT.errorE)
            ff=0
            11	DO 1234 k=1,1000,1
	
	          IF (ISTEP.EQ.0) THEN
                DELTA=1.0D0*0.00001D0*k
	          ELSE
	            DELTA=1.0D0*0.00001D0*(0.001D0**ISTEP)*k
	          END IF
              IF (ff.EQ.0) THEN
	            EvEf=E2vEf+DELTA
	          ELSE
	            EvEf=E2vEf-DELTA
	          ENDIF
        	  write(1,*)'E=',EvEf
              CALL EGV(i,EvEf,kp,km,ky,abskp,abskm,A,A2,A3,kup,kdown,&
Del1up,Del1down,Del2up,Del2down,u1up,u1down,v1up,v1down,&
u2up,u2down,v2up,v2down,alpha,h1vEf,h2vEf,Det,&
gamma1pbarupp,gamma1pbardownp,gamma2pbarupp,gamma2pbardownp,&
gamma1mbarupp,gamma1mbardownp,gamma2mbarupp,gamma2mbardownp,&
gamma1mbarupm,gamma1mbardownm,gamma2mbarupm,gamma2mbardownm,&
gamma1pbarupm,gamma1pbardownm,gamma2pbarupm,gamma2pbardownm,&
gamma1pupp,gamma1pdownp,gamma2pupp,gamma2pdownp,&
gamma1mupp,gamma1mdownp,gamma2mupp,gamma2mdownp,&
gamma1mupm,gamma1mdownm,gamma2mupm,gamma2mdownm,&
gamma1pupm,gamma1pdownm,gamma2pupm,gamma2pdownm,&
q1upe,q2upe,q1uph,q2uph,q1downe,q2downe,&
q1downh,q2downh,L,theta1up,theta1down,theta2up,theta2down,pi,phi)
     

              CALL ZGESVD('N','A',N,N,A2,&
LDA1,S,U,LDU,VT,LDVT,WORK1,LWORK1,RWORK1,INFO1)

	          SS(k,1)=S(16)
1234        CONTINUE

	        IF (ISTEP.NE.0) THEN
	          IF (SS(1,1).LT.SS(2,1)) THEN
	            ff=1
              GOTO 11
 	          ENDIF
	        ENDIF

            SSS(1)=SS(1,1)
            ISTEP2=0
      
            IF (ISTEP.EQ.0) THEN
              DO 4321 k=2,999,1
	            IF (ISTEP2.EQ.0) THEN
                  IF ((SS(k,1).LT.SS(k-1,1)).AND.(SS(k,1).LT.SS(k+1,1))) THEN
                    NNN=k
	                ISTEP2=ISTEP2+1
	                SSS(ISTEP2)=SS(k,1)
	              END IF
                END IF
4321          CONTINUE
              IF ((ISTEP2.EQ.0)) then
                GOTO 3
              ENDIF
            ELSE
              DO 3333 k=2,1000,1
                IF (SS(k,1).LT.SSS(1)) THEN
                  SSS(1)=SS(k,1)
                  NNN=k
                END IF
3333          CONTINUE
	        ENDIF

            ISTEP=ISTEP+1
            IF (ISTEP.EQ.1) THEN
              DELE(ISTEP)=1.0D0*(0.00001D0)*(NNN)
            ELSE
              DELE(ISTEP)=1.0D0*0.00001D0*(0.001D0**(ISTEP-1))*(NNN)*(-1.0D0**ff)
            END IF

            DO 5432 k=1,ISTEP,1
              E2vEf=E2vEf+DELE(ISTEP)
5432        CONTINUE
            write(2,*)'E=',E2vEf
	    END DO

        EvEf=E2vEf

        CALL EGV(i,EvEf,kp,km,ky,abskp,abskm,A,A2,A3,kup,kdown,&
Del1up,Del1down,Del2up,Del2down,u1up,u1down,v1up,v1down,&
u2up,u2down,v2up,v2down,alpha,h1vEf,h2vEf,Det,&
gamma1pbarupp,gamma1pbardownp,gamma2pbarupp,gamma2pbardownp,&
gamma1mbarupp,gamma1mbardownp,gamma2mbarupp,gamma2mbardownp,&
gamma1mbarupm,gamma1mbardownm,gamma2mbarupm,gamma2mbardownm,&
gamma1pbarupm,gamma1pbardownm,gamma2pbarupm,gamma2pbardownm,&
gamma1pupp,gamma1pdownp,gamma2pupp,gamma2pdownp,&
gamma1mupp,gamma1mdownp,gamma2mupp,gamma2mdownp,&
gamma1mupm,gamma1mdownm,gamma2mupm,gamma2mdownm,&
gamma1pupm,gamma1pdownm,gamma2pupm,gamma2pdownm,&
q1upe,q2upe,q1uph,q2uph,q1downe,q2downe,q1downh,q2downh,&
L,theta1up,theta1down,theta2up,theta2down,pi,phi)

        CALL ZGESVD('N','A',N,N,A2,&
LDA1,S,U,LDU,VT,LDVT,WORK1,LWORK1,RWORK1,INFO1)
        DO 234 kk=1,16,1
	     B(kk,1)=VT(16,kk)
234     CONTINUE

	    teup=B(13,1)
	    tedown=B(14,1)
	    thup=B(15,1)
        thdown=B(16,1)


	    FS2e1=teup*u2up*exp(i*q2upe*(xln))+thup*v2up*gamma2mbarupp*&
exp(-i*q2uph*(xln))

	    FS2e2=tedown*u2down*exp(i*q2downe*(xln))+thdown*v2down*&
gamma2mbardownp*exp(-i*q2downh*(xln))

	    FS2e3=teup*v2up*gamma2pupm*exp(i*q2upe*(xln))+thup*u2up*&
exp(-i*q2uph*(xln))

	    FS2e4=tedown*v2down*gamma2pdownm*exp(i*q2downe*(xln))+&
thdown*u2down*exp(-i*q2downh*(xln))

        sttnx(kjj)=h2vEf*sin(alpha)*(FS2e2*conjg(FS2e2)-FS2e1*conjg(FS2e1&
)+FS2e3*conjg(FS2e3)-FS2e4*conjg(FS2e4))+2*h2vEf*cos(alpha)*&
Imag(FS2e2*conjg(FS2e1)+conjg(FS2e3)*FS2e4)

       sttny(kjj)=h2vEf*sin(alpha)*(2*real(conjg(teup)*tedown*exp(i*(q2downe-q2upe)*xln)&
*((conjg(v2up)*v2down*gamma2pupp*gamma2mdownp)-(conjg(u2up)*u2down))&
+conjg(thup)*tedown*exp(i*(q2uph+q2downe)*xln)*((conjg(u2up)*v2down*&
gamma2mdownp)-(conjg(v2up)*u2down)*gamma2mbarupm)+&
conjg(teup)*thdown*exp(-i*(q2upe+q2downh)*xln)*(conjg(v2up)*&
u2down*gamma2pupp-conjg(u2up)*v2down*gamma2pbardownm)+&
conjg(thup)*thdown*exp(i*(q2uph-q2downh)*xln)*(conjg(u2up)*u2down-conjg(v2up)&
*v2down*gamma2pbardownm*gamma2mbarupm)))

       sttnz(kjj)=h2vEf*cos(alpha)*(2*real(conjg(teup)*tedown*exp(i*(q2downe-q2upe)*xln)&
*(conjg(u2up)*u2down-conjg(v2up)*v2down*gamma2pupp*gamma2pdownm)+&
conjg(teup)*thdown*exp(-i*(q2upe+q2downh)*xln)*(conjg(u2up)*v2down*&
gamma2mbardownp-conjg(v2up)*u2down*gamma2pupp)+conjg(thup)*tedown*&
exp(i*(q2uph+q2downe)*xln)*(conjg(v2up)*u2down*gamma2mbarupm-&
conjg(u2up)*v2down*gamma2pdownm)+conjg(thup)*thdown*exp(i*(q2uph-q2downh)&
*xln)*(conjg(v2up)*v2down*gamma2mbardownp*gamma2mbarupm-conjg(u2up)*u2down)))

       Isxlnfs1=imag(i*exp(i*(q2uph+q2downe)*xln)*tedown*conjg(thup)&
*q2downe*(u2down*conjg(v2up)*gamma2mbarupm+v2down*conjg(u2up)*&
gamma2pdownm)+i*exp(i*(q2downe-q2upe)*xln)*tedown*conjg(teup)*q2downe*&
(u2down*conjg(u2up)+v2down*conjg(v2up)*gamma2pupp*gamma2pdownm)&
-i*exp(-i*(q2downh+q2upe)*xln)*thdown*conjg(teup)*q2downh*&
(v2down*conjg(u2up)*gamma2mbardownp+u2down*conjg(v2up)&
*gamma2pupp)-i*exp(i*(-q2downh+q2uph)*xln)*thdown*conjg(thup)*&
q2downh*(v2down*conjg(v2up)*gamma2mbarupm+u2down*conjg(u2up))&
+i*exp(i*(q2upe-q2downe)*xln)*teup*conjg(tedown)*q2upe*&
(u2up*conjg(u2down)+v2up*conjg(v2down)*gamma2pdownp*&
gamma2pupm)-i*exp(-i*(q2downe+q2uph)*xln)&
*thup*conjg(tedown)*q2uph*(v2up*conjg(u2down)*&
gamma2mbarupp+u2up*conjg(v2down)*gamma2pdownp)&
+i*exp(i*(q2downh+q2upe)*xln)*teup*conjg(thdown)*&
q2upe*(u2up*conjg(v2down)*gamma2mbardownm+v2up*conjg(u2down)&
*gamma2pbarupm)-i*exp(i*(q2downh-q2uph)*xln)*thup*&
conjg(thdown)*q2uph*(u2up*conjg(u2down)+v2up*&
conjg(v2down)*gamma2mupp*gamma2mdownm))

      IEnsz(t,tt)=REAL(-kp*(B(1,1)*conjg(B(1,1))-B(2,1)*conjg(B(2,1))&
+B(4,1)&
*conjg(B(4,1))-B(3,1)*conjg(B(3,1)))-km*(-B(5,1)*conjg(B(5,1))&
+B(6,1)*conjg(B(6,1))-B(8,1)*conjg(B(8,1))+&
B(7,1)*conjg(B(7,1))))

      Isz(t,tt)=IEnsz(t,tt)+Isz(t,tt)

      IF (IEnsz(t,tt).NE.0.0D0) THEN
        tolInsz=abs(IEnsz(t,tt)/Isz(t,tt))
      ENDIF

      IEnsx(t,tt)=-kp*REAL(CONJG(B(1,1))*B(3,1)-CONJG(B(2,1))*B(4,1))-&
km*REAL(CONJG(B(6,1))*B(8,1)-CONJG(B(5,1))*B(7,1))

      Isx(t,tt)=IEnsx(t,tt)+Isx(t,tt)

	  IF (Isx(t,tt).NE.0.0D0) THEN
        tolInsx=abs(IEnsx(t,tt)/Isx(t,tt))
      ENDIF

      IEnsy(t,tt)=kp*AIMAG(CONJG(B(1,1))*B(3,1)-CONJG(B(2,1))*B(4,1))+&
km*AIMAG(CONJG(B(5,1))*B(7,1)-CONJG(B(6,1))*B(8,1))

	  Isy(t,tt)=IEnsy(t,tt)+Isy(t,tt)

      IF (IEnsy(t,tt).NE.0.0D0) THEN
         tolInsy=abs(IEnsy(t,tt))/abs(Isy(t,tt))
      ENDIF

      Icharge(t,tt)=-kp*(B(1,1)*conjg(B(1,1))+B(3,1)*conjg(B(3,1))-&
B(2,1)*conjg&
(B(2,1))-B(4,1)*conjg(B(4,1)))-km*(B(6,1)*conjg(B(6,1))+B(8,1)*&
conjg(B(8,1))-B(5,1)*conjg(B(5,1))-B(7,1)*conjg(B(7,1)))

      Ich(t,tt)=Icharge(t,tt)+Ich(t,tt)
	
      rr=0
      GOTO 1
    ENDIF

3  CONTINUE

      DO 41414 uy=1,97,1
       STTx=(3.0D0*0.01D0/8.0D0)*(real(sttnx(uy))+3.0D0*real(sttnx(uy+1))+&
       3.0D0*real(sttnx(uy+2))+real(sttnx(uy+3)))
41414 CONTINUE
      DO 41415 uy=1,97,1
       STTy=(3.0D0*0.01D0/8.0D0)*(real(sttny(uy))+3.0D0*real(sttny(uy+1))+&
       3.0D0*real(sttny(uy+2))+real(sttny(uy+3)))
41415 CONTINUE
      DO 41416 uy=1,97,1
       STTz=(3.0D0*0.01D0/8.0D0)*(real(sttnz(uy))+3.0D0*real(sttnz(uy+1))+&
       3.0D0*real(sttnz(uy+2))+real(sttnz(uy+3)))
41416 CONTINUE
      write(3,*)STTx
      write(4,*)STTy
      write(5,*)STTz

10    CONTINUE

!     Finding exchange coupling (not complete)
!
	  DO 4444 jg=1,21,1
       DO 5555 tg=1,1,1
         Meanrho1(jg,tg) = 0.0d0
         scale=1.0D0
         Meanrho2(jg,tg) = 0.0d0
         Do hhh=1, jg , 2
           Meanrho2(jg,tg)=Meanrho2(jg,tg)+ Isx(hhh-1,tg)/scale+&
4.0d0*Isz(hhh,tg)/scale+Isz(hhh+1,tg)/scale
         End do
         Meanrho2(jg,tg) = Meanrho2(jg,tg) / 3.d0
5555   continue
4444  continue

	 Ic=abs(Ich(6,1))
	 DO 412 k=1,20,1
      IF (abs(Ic).LT.abs(Ich(6,k+1))) Then
        Ic=(Ich(6,k+1))
	  ENDIF

412  CONTINUE

20	CONTINUE
33 CONTINUE

   END program

	
   SUBROUTINE EGV(i,EvEf,kp,km,ky,abskp,abskm,A,A2,A3,kup,kdown,&
Del1up,Del1down,Del2up,Del2down,u1up,u1down,v1up,v1down,&
u2up,u2down,v2up,v2down,alpha,h1vEf,h2vEf,Det,&
gamma1pbarupp,gamma1pbardownp,gamma2pbarupp,gamma2pbardownp,&
gamma1mbarupp,gamma1mbardownp,gamma2mbarupp,gamma2mbardownp,&
gamma1mbarupm,gamma1mbardownm,gamma2mbarupm,gamma2mbardownm,&
gamma1pbarupm,gamma1pbardownm,gamma2pbarupm,gamma2pbardownm,&
gamma1pupp,gamma1pdownp,gamma2pupp,gamma2pdownp,&
gamma1mupp,gamma1mdownp,gamma2mupp,gamma2mdownp,&
gamma1mupm,gamma1mdownm,gamma2mupm,gamma2mdownm,&
gamma1pupm,gamma1pdownm,gamma2pupm,gamma2pdownm,&
q1upe,q2upe,q1uph,q2uph,q1downe,q2downe,&
q1downh,q2downh,L,theta1up,theta1down,theta2up,theta2down,pi,phi)
	
	IMPLICIT NONE
 
	COMPLEX*16 q1uph,q1upe,q1downh,q1downe,gamma1mupp,gamma1mupm,&
gamma1mdownp,gamma1mdownm,gamma1pupp,gamma1pupm,gamma1pdownp,&
gamma1pdownm,u1up,u1down,v1up,v1down,gamma1mbarupp,Det,&
gamma1mbarupm,gamma1mbardownp,gamma1mbardownm,gamma1pbarupp,&
gamma1pbarupm,gamma1pbardownp,gamma1pbardownm

    REAL kup,kdown
    INTEGER s1,s2
	INTEGER N,N1,LDA,IPIV,INFO
    COMPLEX*16  i,kp,km,ky
    COMPLEX*16 q2uph,q2upe,q2downh,q2downe,gamma2mupp,gamma2mupm,&
gamma2mdownp,gamma2mdownm,gamma2pupp,gamma2pupm,gamma2pdownp,&
gamma2pdownm,u2up,u2down,v2up,v2down,gamma2mbarupp,&
gamma2mbarupm,gamma2mbardownp,gamma2mbardownm,gamma2pbarupp,&
gamma2pbarupm,gamma2pbardownp,gamma2pbardownm,Del1up,Del1down,&
Del2up,Del2down,h1vEf,h2vEf
    DOUBLE PRECISION phi,alpha,pi,abskp,abskm,L
    DOUBLE PRECISION theta1up,theta1down,theta2up,theta2down
    DOUBLE PRECISION EvEf
    COMPLEX*16 A
	DIMENSION A(16,16)
	COMPLEX*16 A2
	DIMENSION A2(16,16)
	COMPLEX*16 A3
	DIMENSION A3(16,16)

	kp=SQRT((1+EvEf)-(ky*ky))
	abskp=abs(kp+i*ky)
	km=SQRT((1-EvEf)-(ky*ky))

	abskm=abs(km+i*ky)
	q1upe=SQRT(1+h1vEf+SQRT((EvEf*EvEf)-(Del1up*Del1up))-(ky*ky))

	q1downe=SQRT(1-h1vEf+SQRT((EvEf*EvEf)-(Del1down*Del1down))-(ky*ky))
	q1uph=SQRT(1+h1vEf-SQRT((EvEf*EvEf)-(Del1up*Del1up))-(ky*ky))

	q1downh=SQRT(1-h1vEf-SQRT((EvEf*EvEf)-(Del1down*Del1down))-(ky*ky))

	q2upe=SQRT(1+h2vEf+SQRT((EvEf*EvEf)-(Del2up*Del2up))-(ky*ky))
	q2downe=SQRT(1-h2vEf+SQRT((EvEf*EvEf)-(Del2down*Del2down))&
-(ky*ky))

	q2uph=SQRT(1+h2vEf-SQRT((EvEf*EvEf)-(Del2up*Del2up))-(ky*ky))
	q2downh=SQRT(1-h2vEf-SQRT((EvEf*EvEf)-(Del2down*Del2down))-(ky*ky))

	gamma1pupp=((kp+i*ky)/(kup*abskp))*exp(i*theta1up)*exp(i*phi)
	gamma1pdownp=((kp+i*ky)/(kdown*abskp))*exp(i*theta1down)*exp(i*phi)
	gamma1pupm=((kp-i*ky)/(kup*abskp))*exp(-i*theta1up)*exp(-i*phi)
	gamma1pdownm=((kp-i*ky)/(kdown*abskp))*exp(-i*theta1down)*exp(-i*phi)

	gamma1mupp=((-kp+i*ky)/(kup*abskp))*exp(i*theta1up)*exp(i*phi)
	gamma1mdownp=((-kp+i*ky)/(kdown*abskp))*exp(i*theta1down)*exp(i*phi)
	gamma1mupm=((-kp-i*ky)/(kup*abskp))*exp(-i*theta1up)*exp(-i*phi)
	gamma1mdownm=((-kp-i*ky)/(kdown*abskp))*exp(-i*theta1down)*exp(-i*phi)

	gamma2pupp=((kp+i*ky)/(kup*abskp))*exp(i*theta2up)
	gamma2pdownp=((kp+i*ky)/(kdown*abskp))*exp(i*theta2down)
	gamma2pupm=((kp-i*ky)/(kup*abskp))*exp(-i*theta2up)
	gamma2pdownm=((kp-i*ky)/(kdown*abskp))*exp(-i*theta2down)

	gamma2mupp=((-kp+i*ky)/(kup*abskp))*exp(i*theta2up)
	gamma2mdownp=((-kp+i*ky)/(kdown*abskp))*exp(i*theta2down)
	gamma2mupm=((-kp-i*ky)/(kup*abskp))*exp(-i*theta2up)
	gamma2mdownm=((-kp-i*ky)/(kdown*abskp))*exp(-i*theta2down)

	gamma1pbarupp=((km+i*ky)/(kup*abskm))*exp(i*theta1up)*exp(i*phi)
	gamma1pbardownp=((km+i*ky)/(kdown*abskm))*exp(i*theta1down)*exp(i*&
phi)
	gamma1pbarupm=((km-i*ky)/(kup*abskm))*exp(-i*theta1up)*exp(-i*phi)
	gamma1pbardownm=((km-i*ky)/(kdown*abskm))*exp(-i*theta1down)&
*exp(-i*phi)

	gamma1mbarupp=((-km+i*ky)/(kup*abskm))*exp(i*theta1up)*exp(i*phi)
	gamma1mbardownp=((-km+i*ky)/(kdown*abskm))*exp(i*theta1down)&
*exp(i*phi)
	gamma1mbarupm=((-km-i*ky)/(kup*abskm))*exp(-i*theta1up)*exp(-i*phi)
	gamma1mbardownm=((-km-i*ky)/(kdown*abskm))*exp(-i*theta1down)&
*exp(-i*phi)

	gamma2pbarupp=((km+i*ky)/(kup*abskm))*exp(i*theta2up)
	gamma2pbardownp=((km+i*ky)/(kdown*abskm))*exp(i*theta2down)
	gamma2pbarupm=((km-i*ky)/(kup*abskm))*exp(-i*theta2up)
	gamma2pbardownm=((km-i*ky)/(kdown*abskm))*exp(-i*theta2down)

	gamma2mbarupp=((-km+i*ky)/(kup*abskm))*exp(i*theta2up)
	gamma2mbardownp=((-km+i*ky)/(kdown*abskm))*exp(i*theta2down)
	gamma2mbarupm=((-km-i*ky)/(kup*abskm))*exp(-i*theta2up)
	gamma2mbardownm=((-km-i*ky)/(kdown*abskm))*exp(-i*theta2down)

	u1up=SQRT((1+(SQRT((EvEf*EvEf)-(Del1up*Del1up))/EvEf))/2.0D0)
	u1down=SQRT((1+(SQRT((EvEf*EvEf)-(Del1down*Del1down))/EvEf))/2.0D0)

	u2up=SQRT((1+(SQRT((EvEf*EvEf)-(Del2up*Del2up))/EvEf))/2.0D0)
	u2down=SQRT((1+(SQRT((EvEf*EvEf)-(Del2down*Del2down))/EvEf))/2.0D0)

	v1up=SQRT((1-(SQRT((EvEf*EvEf)-(Del1up*Del1up))/EvEf))/2.0D0)
	v1down=SQRT((1-(SQRT((EvEf*EvEf)-(Del1down*Del1down))/EvEf))/2.0D0)
	
	v2up=SQRT((1-(SQRT((EvEf*EvEf)-(Del2up*Del2up))/EvEf))/2.0D0)
	v2down=SQRT((1-(SQRT((EvEf*EvEf)-(Del2down*Del2down))/EvEf))&
/2.0D0)

	A(1,1)=(-1.0,0.0);A(1,2)=(-1.0,0.0);A(1,3)=(0.0,0.0)
	A(1,4)=(0.0,0.0);A(1,5)=(0.0,0.0);A(1,6)=(0.0,0.0)
	A(1,7)=(0.0,0.0);A(1,8)=(0.0,0.0)
	A(1,9)=(cos(alpha/2.0D0))*u1up
	A(1,10)=i*(sin(alpha/2.0D0))*u1down
	A(1,11)=(cos(alpha/2.0D0))*v1up*gamma1pbarupp
	A(1,12)=i*(sin(alpha/2.0D0))*v1down*gamma1pbardownp
	A(1,13)=(0.0,0.0);A(1,14)=(0.0,0.0);A(1,15)=(0.0,0.0)
	A(1,16)=(0.0,0.0)
	
	A(2,1)=(0.0,0.0);A(2,2)=(0.0,0.0);A(2,3)=(-1.0,0.0)
	A(2,4)=(-1.0,0.0);A(2,5)=(0.0,0.0);A(2,6)=(0.0,0.0)
	A(2,7)=(0.0,0.0);A(2,8)=(0.0,0.0)
	A(2,9)=i*(sin(alpha/2.0D0))*u1up
	A(2,10)=(cos(alpha/2.0D0))*u1down
	A(2,11)=i*(sin(alpha/2.0D0))*v1up*gamma1pbarupp
	A(2,12)=(cos(alpha/2.0D0))*v1down*gamma1pbardownp
	A(2,13)=(0.0,0.0);A(2,14)=(0.0,0.0);A(2,15)=(0.0,0.0)
	A(2,16)=(0.0,0.0)

	A(3,1)=(0.0,0.0);A(3,2)=(0.0,0.0);A(3,3)=(0.0,0.0)
	A(3,4)=(0.0,0.0);A(3,5)=(-1.0,0.0);A(3,6)=(-1.0,0.0)
	A(3,7)=(0.0,0.0);A(3,8)=(0.0,0.0)
	A(3,9)=(cos(alpha/2.0D0))*v1up*gamma1mupm
	A(3,10)=-i*(sin(alpha/2.0D0))*v1down*gamma1mdownm
	A(3,11)=(cos(alpha/2.0D0))*u1up
	A(3,12)=-i*(sin(alpha/2.0D0))*u1down
	A(3,13)=(0.0,0.0);A(3,14)=(0.0,0.0);A(3,15)=(0.0,0.0)
	A(3,16)=(0.0,0.0)
	
	A(4,1)=(0.0,0.0);A(4,2)=(0.0,0.0);A(4,3)=(0.0,0.0)
	A(4,4)=(0.0,0.0);A(4,5)=(0.0,0.0);A(4,6)=(0.0,0.0)
	A(4,7)=(-1.0,0.0);A(4,8)=(-1.0,0.0)
	A(4,9)=-i*(sin(alpha/2.0D0))*v1up*gamma1mupm
	A(4,10)=(cos(alpha/2.0D0))*v1down*gamma1mdownm
	A(4,11)=-i*(sin(alpha/2.0D0))*u1up
	A(4,12)=(cos(alpha/2.0D0))*u1down
	A(4,13)=(0.0,0.0);A(4,14)=(0.0,0.0);A(4,15)=(0.0,0.0)
	A(4,16)=(0.0,0.0)

	A(5,1)=(-1.0,0.0);A(5,2)=(1.0,0.0);A(5,3)=(0.0,0.0)
	A(5,4)=(0.0,0.0);A(5,5)=(0.0,0.0);A(5,6)=(0.0,0.0)
	A(5,7)=(0.0,0.0);A(5,8)=(0.0,0.0)
	A(5,9)=-(cos(alpha/2.0D0))*u1up*(q1upe/kp)
	A(5,10)=-i*(sin(alpha/2.0D0))*u1down*(q1downe/kp)
	A(5,11)=(cos(alpha/2.0D0))*v1up*gamma1pbarupp*(q1uph/kp)
	A(5,12)=i*(sin(alpha/2.0D0))*v1down*gamma1pbardownp*(q1downh/kp)
	A(5,13)=(0.0,0.0);A(5,14)=(0.0,0.0);A(5,15)=(0.0,0.0)
	A(5,16)=(0.0,0.0)

	A(6,1)=(0.0,0.0);A(6,2)=(0.0,0.0);A(6,3)=(-1.0,0.0);
	A(6,4)=(1.0,0.0);A(6,5)=(0.0,0.0);A(6,6)=(0.0,0.0)
	A(6,7)=(0.0,0.0);A(6,8)=(0.0,0.0)
	A(6,9)=-i*(sin(alpha/2.0D0))*u1up*(q1upe/kp)
	A(6,10)=-(cos(alpha/2.0D0))*u1down*(q1downe/kp)
	A(6,11)=i*(sin(alpha/2.0D0))*v1up*gamma1pbarupp*(q1uph/kp)
	A(6,12)=(cos(alpha/2.0D0))*v1down*gamma1pbardownp*(q1downh/kp)
	A(6,13)=(0.0,0.0);A(6,14)=(0.0,0.0);A(6,15)=(0.0,0.0)
	A(6,16)=(0.0,0.0)
	
	A(7,1)=(0.0,0.0);A(7,2)=(0.0,0.0);A(7,3)=(0.0,0.0)
	A(7,4)=(0.0,0.0);A(7,5)=(1.0,0.0);A(7,6)=(-1.0,0.0)
	A(7,7)=(0.0,0.0);A(7,8)=(0.0,0.0)
	A(7,9)=-(cos(alpha/2.0D0))*v1up*gamma1mupm*(q1upe/km)
	A(7,10)=i*(sin(alpha/2.0D0))*v1down*gamma1mdownm*(q1downe/km)
	A(7,11)=(cos(alpha/2.0D0))*u1up*(q1uph/km)
	A(7,12)=-i*(sin(alpha/2.0D0))*u1down*(q1downh/km)
	A(7,13)=(0.0,0.0);A(7,14)=(0.0,0.0);A(7,15)=(0.0,0.0)
	A(7,16)=(0.0,0.0)

	A(8,1)=(0.0,0.0);A(8,2)=(0.0,0.0);A(8,3)=(0.0,0.0)
	A(8,4)=(0.0,0.0);A(8,5)=(0.0,0.0);A(8,6)=(0.0,0.0)
	A(8,7)=(1.0,0.0);A(8,8)=(-1.0,0.0)
	A(8,9)=i*(sin(alpha/2.0D0))*v1up*gamma1mupm*(q1upe/km)
	A(8,10)=-(cos(alpha/2.0D0))*v1down*gamma1mdownm*(q1downe/km)
	A(8,11)=-i*(sin(alpha/2.0D0))*u1up*(q1uph/km)
	A(8,12)=(cos(alpha/2.0D0))*u1down*(q1downh/km)
	A(8,13)=(0.0,0.0);A(8,14)=(0.0,0.0);A(8,15)=(0.0,0.0)
	A(8,16)=(0.0,0.0)

	A(9,1)=-exp(i*2*pi*kp*L);A(9,2)=-exp(-i*2*pi*kp*L)
	A(9,3)=(0.0,0.0);A(9,4)=(0.0,0.0);A(9,5)=(0.0,0.0)
	A(9,6)=(0.0,0.0);A(9,7)=(0.0,0.0);A(9,8)=(0.0,0.0)
	A(9,9)=(0.0,0.0);A(9,10)=(0.0,0.0);A(9,11)=(0.0,0.0)
	A(9,12)=(0.0,0.0);A(9,13)=u2up*exp(i*2*pi*q2upe*L)
	A(9,14)=(0.0,0.0);
	A(9,15)=v2up*exp(-i*2*pi*q2uph*L)*gamma2mbarupp
	A(9,16)=(0.0,0.0)

	A(10,1)=(0.0,0.0);A(10,2)=(0.0,0.0)
	A(10,3)=-exp(i*2*pi*kp*L);A(10,4)=-exp(-i*2*pi*kp*L)
	A(10,5)=(0.0,0.0);A(10,6)=(0.0,0.0);A(10,7)=(0.0,0.0)
	A(10,8)=(0.0,0.0);A(10,9)=(0.0,0.0);A(10,10)=(0.0,0.0)
	A(10,11)=(0.0,0.0);A(10,12)=(0.0,0.0);A(10,13)=(0.0,0.0)
	A(10,14)=u2down*exp(i*2*pi*q2downe*L);A(10,15)=(0.0,0.0)
	A(10,16)=v2down*exp(-i*2*pi*q2downh*L)*gamma2mbardownp

	A(11,1)=(0.0,0.0);A(11,2)=(0.0,0.0);A(11,3)=(0.0,0.0)
	A(11,4)=(0.0,0.0);A(11,5)=-exp(-i*2*pi*km*L)
	A(11,6)=-exp(i*2*pi*km*L);A(11,7)=(0.0,0.0)
	A(11,8)=(0.0,0.0);A(11,9)=(0.0,0.0);A(11,10)=(0.0,0.0)
	A(11,11)=(0.0,0.0);	A(11,12)=(0.0,0.0)
	A(11,13)=v2up*exp(i*2*pi*q2upe*L)*gamma2pupm
	A(11,14)=(0.0,0.0);A(11,15)=u2up*exp(-i*2*pi*q2uph*L)
	A(11,16)=(0.0,0.0)

	A(12,1)=(0.0,0.0);A(12,2)=(0.0,0.0);A(12,3)=(0.0,0.0)
	A(12,4)=(0.0,0.0);A(12,5)=(0.0,0.0);A(12,6)=(0.0,0.0)
	A(12,7)=-exp(-i*2*pi*km*L);A(12,8)=-exp(i*2*pi*km*L)
	A(12,9)=(0.0,0.0);A(12,10)=(0.0,0.0);A(12,11)=(0.0,0.0)
	A(12,12)=(0.0,0.0);A(12,13)=(0.0,0.0)
	A(12,14)=v2down*exp(i*2*pi*q2downe*L)*gamma2pdownm
	A(12,15)=(0.0,0.0);A(12,16)=u2down*exp(-i*2*pi*q2downh*L)

	A(13,1)=-exp(i*2*pi*kp*L);A(13,2)=exp(-i*2*pi*kp*L)
	A(13,3)=(0.0,0.0);A(13,4)=(0.0,0.0);A(13,5)=(0.0,0.0)
	A(13,6)=(0.0,0.0);A(13,7)=(0.0,0.0);A(13,8)=(0.0,0.0)
	A(13,9)=(0.0,0.0);A(13,10)=(0.0,0.0);A(13,11)=(0.0,0.0)
	A(13,12)=(0.0,0.0);
	A(13,13)=u2up*exp(i*2*pi*q2upe*L)*(q2upe/kp)
	A(13,14)=(0.0,0.0)
	A(13,15)=v2up*exp(-i*2*pi*q2uph*L)*gamma2mbarupp*(-q2uph/kp)
	A(13,16)=(0.0,0.0)

	A(14,1)=(0.0,0.0);A(14,2)=(0.0,0.0);A(14,3)=-exp(i*2*pi*kp*L)
	A(14,4)=exp(-i*2*pi*kp*L);A(14,5)=(0.0,0.0);A(14,6)=(0.0,0.0)
	A(14,7)=(0.0,0.0);A(14,8)=(0.0,0.0);A(14,9)=(0.0,0.0)
	A(14,10)=(0.0,0.0);A(14,11)=(0.0,0.0);A(14,12)=(0.0,0.0)
	A(14,13)=(0.0,0.0)
	A(14,14)=u2down*exp(i*2*pi*q2downe*L)*(q2downe/kp)
	A(14,15)=(0.0,0.0)
	A(14,16)=v2down*exp(-i*2*pi*q2downh*L)*gamma2mbardownp*(-q2downh/kp)

	A(15,1)=(0.0,0.0);A(15,2)=(0.0,0.0);A(15,3)=(0.0,0.0)
	A(15,4)=(0.0,0.0);A(15,5)=exp(-i*2*pi*km*L)
	A(15,6)=-exp(i*2*pi*km*L);A(15,7)=(0.0,0.0)
	A(15,8)=(0.0,0.0);A(15,9)=(0.0,0.0);A(15,10)=(0.0,0.0)
	A(15,11)=(0.0,0.0);A(15,12)=(0.0,0.0)
	A(15,13)=v2up*exp(i*2*pi*q2upe*L)*gamma2pupm*(q2upe/km)
	A(15,14)=(0.0,0.0)
	A(15,15)=u2up*exp(-i*2*pi*q2uph*L)*(-q2uph/km)
	A(15,16)=(0.0,0.0)

	A(16,1)=(0.0,0.0);A(16,2)=(0.0,0.0);A(16,3)=(0.0,0.0)
	A(16,4)=(0.0,0.0);A(16,5)=(0.0,0.0);A(16,6)=(0.0,0.0)
	A(16,7)=exp(-i*2*pi*km*L);A(16,8)=-exp(i*2*pi*km*L)
	A(16,9)=(0.0,0.0);A(16,10)=(0.0,0.0);A(16,11)=(0.0,0.0)
	A(16,12)=(0.0,0.0);A(16,13)=(0.0,0.0)
	A(16,14)=v2down*exp(i*2*pi*q2downe*L)*gamma2pdownm*(q2downe/km)
	A(16,15)=(0.0,0.0)
	A(16,16)=u2down*exp(-i*2*pi*q2downh*L)*(-q2downh/km)

	  
    DO 77 s1=1,16,1
      DO 777 s2=1,16,1
        A2(s1,s2)=A(s1,s2)
        A3(s1,s2)=A(s1,s2)
777   CONTINUE
77  CONTINUE

	Det=A(1,1)*A(2,2)*A(3,3)*A(4,4)*A(5,5)*A(6,6)*A(7,7)&
*A(8,8)*A(9,9)*A(10,10)*A(11,11)*A(12,12)*A(13,13)*A(14,14)&
*A(15,15)*A(16,16)

	return
	end subroutine
