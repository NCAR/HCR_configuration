        program liebe
	parameter (n=100)
	integer num_layers
	real freq,sfcpres,height(n),temp(n),humidity(n),
     .  vap(n),cloudlwc(n),extinct(n)
	open(unit=3,file='atmosn.dat',status='unknown')
	read(3,*) freq,sfcpres,num_layers
	write(6,*) freq,sfcpres,num_layers
	do 10 i=num_layers+1,1,-1
        read(3,*) height(i), temp(i),humidity(i),cloudlwc(i)
 10     continue
	CALL RDLINS
        call MWABS (FREQ, SFCPRES, NUM_LAYERS, HEIGHT, TEMP,
     .              HUMIDITY,vap, CLOUDLWC, EXTINCT)
	do 20 i=1,num_layers
20     write(16,*) height(i), temp(i),humidity(i),vap(i),
     .  cloudlwc(i), extinct(i)*4.343
	stop
	end


      SUBROUTINE MWABS (FREQ, SFCPRES, NUM_LAYERS, HEIGHT, TEMP,
     .              HUMIDITY, vap,CLOUDLWC, EXTINCT)
      INTEGER NUM_LAYERS
      REAL  FREQ, HEIGHT(*), TEMP(*), HUMIDITY(*), CLOUDLWC(*)
      REAL  EXTINCT(*), vap(*),SFCPRES,theta
C        MWABS CALCULATES THE EXTINCTION DUE TO ABSORPTION BY AIR,
C      WATER VAPOR AND CLOUD LIQUID WATER FOR MICROWAVE FREQUENCIES.
C      THE ATMOSPHERIC INFORMATION IS IN TERMS OF THE DISCRETE LAYERS
C      FOR THE RADIATIVE TRANSFER CALCULATION.  THE HEIGHT, TEMPERATURE,
C      AND HUMIDITY ARE GIVEN AT THE ITERFACES BETWEEN LAYERS AND ARE
C      LINEARLY INTERPOLATED ACROSS EACH LAYER. THE CLOUD LWC AND THE
C      RETURNED EXTINCTION ARE FOR THE WHOLE LAYER.
C
      INTEGER L, I, NUMZ, MAXN
      REAL    FRAC, PRES, TV, PATANG, CUMATT, PWCALC
      PARAMETER (MAXN=50)
      REAL    HT(MAXN), P(MAXN), T(MAXN), RH(MAXN), E(MAXN)
      REAL    W(MAXN), RR(MAXN), ATTEN(MAXN), DISPER(MAXN)

C     P=PRESSURE, T=TEMP., RH=RELATIVE HUMID, E=WATER VAP. PARTIAL PRESS.,
C     W=LIQUID WATER CONCENTRATION (G./CU.M), RR=RAIN RATE(MM/HR)
C     ATTEN=ATTENUATION (DB/KM),  DISPER=DISPERSION (PHASE/KM)


      PATANG = 0.0

      PRES = SFCPRES
      print *,'layers',num_layers
      DO 120 L = NUM_LAYERS, 1, -1
C              INTERPOLATE THE HEIGHT, TEMPERATURE, AND HUMIDITY
        NUMZ = MAX0( NINT((HEIGHT(L) - HEIGHT(L+1))/ 0.25) + 1, 2)
        NUMZ = MIN0 (NUMZ, MAXN)
        DO 100 I = 1, NUMZ
          FRAC = FLOAT(I-1)/(NUMZ-1)
          HT(I) = HEIGHT(L+1) + FRAC*(HEIGHT(L) - HEIGHT(L+1))
          T(I) = TEMP(L+1) + FRAC*(TEMP(L) - TEMP(L+1))
          RH(I) = HUMIDITY(L+1) + FRAC*(HUMIDITY(L) - HUMIDITY(L+1))
          W(I) = CLOUDLWC(L)
          RR(I) = 0.0
          theta=300/T(i)
	Vap(I) = 1.739e+11*RH(i)*theta**5*exp(-22.64*theta)
	write(*,*)ht(i),rh(i),T(i),theta, exp(-22.64*theta),vap(i)
100     CONTINUE
C              CALCULATE THE PRESSURE USING THE HYDROSTATIC RELATION,
C                INCLUDING THE EFFECTS OF WATER VAPOR ON DENSITY
        P(1) = PRES
        DO 110 I = 1, NUMZ-1
          E(I) = PWCALC (RH(I),P(I),T(I))
          TV = (T(I+1)+T(I))/2. *(1. + 0.61*E(I)/P(I))
          P(I+1) = P(I) - 9.8*P(I)/(287.*TV) *1000.* (HT(I+1)-HT(I))
110     CONTINUE
        E(NUMZ) = PWCALC (RH(NUMZ),P(NUMZ),T(NUMZ))
        PRES = P(NUMZ)

        CALL DOLEVS (NUMZ,P,T,RH,E,W,RR, FREQ, ATTEN,DISPER)

        CALL PTHINT (NUMZ,HT, PATANG, ATTEN, CUMATT)
        EXTINCT(L) = 0.1*DLOG(10.D0) *CUMATT/ (HEIGHT(L)-HEIGHT(L+1))
c	write(6,*)'extinction',l,extinct(l)
120   CONTINUE

      RETURN
      END

      SUBROUTINE OUTPUT (NUMLEV, FREQ, PATANG, HT,P,T,RH,W,RR,
     .                   ATTEN,DISPER, CUMATT,CUMDSP)
C       OUTPUTS THE ATTENUATION AND DISPERSION COEFFICIENTS FOR EACH
C         LEVEL AS WELL AS THE INTEGRATED VALUES.
      DIMENSION HT(100), P(100), T(100), RH(100), W(100), RR(100)
      DIMENSION ATTEN(100), DISPER(100)
      CHARACTER*64 FILENM

c      PRINT *,'OUTPUT FILE'
c      READ (*,'(A)') FILENM
      OPEN (UNIT=1, FILE='ext.dat',status='new')

      WRITE (1,'(A,F10.5,A)') 'ATTENUATION/DISPERSION FOR ',
     .                        FREQ, ' GHZ'
      WRITE (1,'(A,F6.2,A)') 'RAY PATH ANGLE: ',
     .                        PATANG, ' DEGREES FROM ZENITH'
      WRITE (1,*)
      WRITE (1,'(A,F10.3)') 'INTEGRATED ATTENUATION (DB) : ', CUMATT
      WRITE (1,'(A,F10.3)') 'INTEGRATED DISPERSION (PS) : ', CUMDSP
      WRITE (1,*)
      WRITE (1,'(A,A)') '   HEIGHT    ATTEN    DELAY   ',
     .               '    PRES     TEMP      RH       LWC      RR'
      WRITE (1,'(A,A)') '    (KM)    (DB/KM)  (PS/KM)  ',
     .               '    (KPA)     (K)             (G/M~)  (MM/HR)'
1     FORMAT (1X,F7.3,3X,F7.3,2X,F7.3,
     .        5X,F7.2,2X,F7.2,2X,F6.3,2X,F7.3,2X,F7.2)
      DO 100 I = 1, NUMLEV
          WRITE (1,1) HT(I), ATTEN(I), DISPER(I),
     .                P(I), T(I), RH(I), W(I), RR(I)
100   CONTINUE
      CLOSE (1)

      RETURN
      END






      SUBROUTINE PTHINT (NUMLEV,HT,PATANG, DATA, CUMDAT)
C       INTEGRATES THE DATA OVER THE RAY PATH BY SUMMING LEVELS.
C         THE LAYER VALUE IS TAKEN TO BE THE AVERAGE VALUES AT THE
C         BORDERING LEVELS.
      DIMENSION HT(100), DATA(100)

C.......... COMPUTE SLANT PATH LENGTH FACTOR.  (COSECANT OF RAY
C           PATH ANGLE)
      FACTOR = 1./(COS(3.1415926*PATANG/180.))

      CUMDAT = 0.0
      DO 100 L = 1, NUMLEV-1
          PATH = (HT(L+1) - HT(L)) * FACTOR
          CUMDAT = CUMDAT + PATH*(DATA(L)+DATA(L+1))/2.
100   CONTINUE

      RETURN
      END



      SUBROUTINE DOLEVS (NUMLEV,PZ,TZ,RHZ,EZ,WZ,RRZ, FREQ,ATTEN,DISPER)
C       CALLS THE BASIC MICROWAVE CALCULATION ROUTINE FOR EACH LEVEL

      DIMENSION PZ(100), TZ(100), RHZ(100), EZ(100), WZ(100), RRZ(100)
      DIMENSION ATTEN(100), DISPER(100)

      COMMON /COEFF/ C1,C2,X1,X2,X3,X4
      COMMON/CLIMDAT/P,T,THETA,RH,E,W,RR
      COMMON/LINES/F0(50),A1(50),A2(50),A3(50),A4(50),A5(50),A6(50),
     -F0W(35),B1(35),B2(35),B3(35)

      DO 100 L=1,NUMLEV
        P=PZ(L)
        T=TZ(L)
        THETA = 300./T
        RH=RHZ(L)
        E=EZ(L)
        W=WZ(L)
        RR=RRZ(L)
        CALL ATDCALC (FREQ,ATTEN(L),DISPER(L))
100   CONTINUE

      RETURN
      END

      FUNCTION PWCALC (RH,P,T)
C-------- THIS FUNCTION CALCULATES PARTIAL WATER VAPOR PRESSURE AT
C         A GIVEN TEMP. AND RELATIVE HUMIDITY.  PRESSURE IS NOT CURRENTLY
C         USED.
      ST = 300./T
      PWSAT = (2.4089*ST**5.)*(10.**(10.-9.834*ST))
      PWCALC = PWSAT * RH

      RETURN
      END



      SUBROUTINE RDLINS
C       READS IN THE O2 AND H2O LINE DATABASE IN FILE 'LNDAT.DAT'.
C       ABSORPTION LINE DATA PASSED THRU /LINES/ COMMON.
C       ALSO SETS UP THE NON-RESONANT H20 (CONTINUUM) COEFFICIENTS
C         WHICH ARE PASSED THRU /COEFF/ COMMON.
C       DATABASE FILE CONTAINS FREQ AND COEFS FOR LINES FOLLOWED
C       BY -1.0 FLAG.
C
      COMMON /COEFF/ C1,C2,X1,X2,X3,X4
      COMMON/LINES/F0(50),A1(50),A2(50),A3(50),A4(50),A5(50),A6(50),
     -F0W(35),B1(35),B2(35),B3(35)

C................................................................
C          READ RESONANCE LINE DATA FROM FILE AND STORE IN
C          ARRAYS PASSED THROUGH COMMON /LINES/
C................................................................
      OPEN (UNIT=2, FILE='line.dat',STATUS='OLD')
      N = 1
100   CONTINUE
C          READ THE O2 LINES AND THEIR COEFFICIENTS
      IF (N .GT. 50) THEN
        PRINT *,'ERROR: MORE THAN 50 O2 LINES'
        STOP
      ENDIF
      READ (2,1)F0(N),A1(N),A2(N),A3(N),A4(N),A5(N),A6(N)
    1 FORMAT (F10.6,F10.4,F9.3,4F9.2)
      IF(F0(N).LT.0.)GO TO 200
      N = N + 1
      GOTO 100

200   N = 1
210   CONTINUE
C          READ THE H2O LINES AND THEIR COEFFICIENTS
      IF (N .GT. 35) THEN
        PRINT *,'ERROR: MORE THAN 35 H2O LINES'
        STOP
      ENDIF
      READ (2,1) F0W(N),B1(N),B2(N),B3(N)
      IF(F0W(N).LT.0.)GO TO 300
      N = N + 1
      GOTO 210

300   CONTINUE
      CLOSE (2)


C---------  SPECIFY NON-RESONANT H2O COEFFICIENTS
C           LAST MODIFIED 7/24/84
      C1= 1.40
      C2= 54.1
      X1= 1.0
      X2= 2.5
      X3= 1.0
      X4=  5.5

      RETURN
      END






      SUBROUTINE ATDCALC (F,FNPP,D)

C..........THIS ROUTINE IS A NUMERICAL MODEL OF ATTENUATION AND PHASE
C          OF RADIO WAVES IN AIR FOR FREQUENCIES 1-1000 GHZ.
C          THE MODEL IS DESCRIBED COMPLETELY IN THE FOLLOWING PAPER: 
C                MODELING ATTENUATION AND PHASE OF RADIO WAVES IN
C                AIR AT FREQUENCIES BELOW 1000 GHZ
C
C                BY HANS J. LIEBE, NTIA/ITS, BOULDER, CO  80303
C          ALL VARIABLE NAMES AND EQUATION NUMBERS ARE AS THEY APPEAR
C          IN THIS PAPER IN RADIO SCIENCE 16 (1981).
C
C          ROUTINE WRITTEN BY BEN SHAW, NTIA/ITS,  7/81

C...............INPUT VARIABLES: 
C               1) /LINES/ -- COMMON BLOCK WITH ARRAYS CONTAINING THE
C                             48 OXYGEN LINES, 30 WATER VAPOR LINES, AND
C                             THEIR ASSOCIATED COEFFICIENTS.
C               2) /CLIMDAT/ -COMMON BLOCK WITH CLIMATALOGICAL DATA FOR
C                             CALCULATIONS: 
C                                   P -- DRY AIR PRESSURE, (KPA)
C                                   T -- TEMPERATURE (KELVINS)
C                                   THETA -- RELATIVE INVERSE TEMP,
C                                            = 300/T
C                                   RH -- RELATIVE HUMIDITY (0 < RH < 1)
C                                   E -- WATER VAPOR PARTIAL PRESS (KPA)
C                                   W -- LIQUID WATER CONCENTRATION,
C                                          (GM./CU.MTR.)
C                                   RR -- RAIN RATE (MM/HR)
C               3) /COEFF/ -- COMMON BLOCK WITH COEFFICIENTS FOR FAR-WING
C                             NON-RESONANT H2O TERMS.
C               4) F       -- FREQUENCY (GHZ) AT WHICH ATTENUATION AND
C                             DISPERSION ARE TO BE EVALUATED.

C...............OUTPUT VARIABLES: 
C              1)  FNPP   -- ATTENUATION, (DB/KM)
C              2)     D   -- REFRACTIVE DISPERSION, (*1.0E-6)
C.......................................................................


      COMMON /CLIMDAT/ P,T,THETA,RH,E,W,RR
      COMMON /COEFF/ C1,C2,X1,X2,X3,X4
      COMMON /LINES/ F0(50),A1(50),A2(50),A3(50),A4(50),A5(50),A6(50),
     -F0W(35),B1(35),B2(35),B3(35)


C......................................................................
C       IMPORTANT INTERNAL VARIABLES: 
C            -- SUMSFPP = SUMATION OF PRODUCT OF S AND F DOUBLE PRIME
C                         (LINE BY LINE ATTENUATION CONTRIBUTION)
C            -- SUMSFP  = SUMATION OF PRODUCT OF S AND F PRIME
C                         (LINE BY LINE DISPERSION CONTRIBUTION)
C            -- FNPPV   = FLOATING POINT N DOUBLE PRIME SUB V
C                         (FAR WING NON-RES H2O ATTENUATION CONTRIBUTION.)
C            -- FNPPD   = FLOATING POINT N DOUBLE PRIME SUB D
C                         (NON-RES DRY AIR ATTENUATION CONTRIBUTION)
C            -- FNPPW   = FLOATING POINT N DOUBLE PRIME SUB W
C                         (LIQUID WATER ATTENUATION CONTRIBUTION)
C.....................................................................

C.....................................................................
C        EVALUATE LINE BY LINE CONTRIBUTION TO ATTENUATION AND DISPERSION
C        FIRST FOR O2 LINES, THEN FOR H2O LINES.  THE SENTINEL FOR THE END
C        OF A GROUP OF LINES IS A NEGATIVE LINE CENTER FREQUENCY (F0 OR F0W)
C......................................................................

C    IF FREQUENCY <= 0.0, SET OUTPUT VARIABLES TO 0 AND RETURN
      IF (F.GT.0.0) GO TO 100
C        SINCE ATTENUATION IS PLOTTED ON LOG AXIS, THE VALUE CANNOT BE
C        LESS THAN OR EQUAL TO 0.0
      FNPP =0.00000001
      D = 0.0
      RETURN

  100 SUMSFPP = 0.0
      SUMSFP = 0.0
      ONEMTTA = 1. - THETA

C..........BEGIN CALCULATION OF O2 RESONANCE LINE CONTRIBUTION.
C          REPEAT FOR EACH LINE UNTIL NEGATIVE FREQUENCY SENTINEL.
C****         NOTE:  SLIGHT MODIFICATION OF MODEL ON 8/31/82.
C                     LINES THAT WERE CHANGED REMAIN AS COMMENTS.
C
      DO 101 N=1,1000
      IF (F0(N).LT.0) GO TO 199
         S = A1(N)*P*THETA**3.*EXP(A2(N)*ONEMTTA)*1.E-6
C*       GAMMA = A3(N)*(P+1.3*E)*THETA**0.9*1.E-3
         GAMMA = A3(N)*(P*THETA**(0.8-A6(N))+1.1*E*THETA)*1.E-3
         PSI = A4(N)*P*THETA**A5(N)*1.E-3

         XNU = F0(N)
         XNUMF = XNU - F
         XNUPF = XNU + F
         XNUMFSQ = XNUMF*XNUMF
         XNUPFSQ = XNUPF * XNUPF
         GAMMASQ = GAMMA * GAMMA
C*       GAMPSI = GAMMA * PSI
         PSIF = PSI * F


         SUMSFPP = SUMSFPP + S * F/XNU*(((GAMMA - XNUMF*PSI)/(XNUMFSQ+
     -GAMMASQ))+((GAMMA-XNUPF*PSI)/(XNUPFSQ+GAMMASQ)))
C*       SUMSFP=SUMSFP+S*(((XNUMF+GAMPSI)/(XNUMFSQ+GAMMASQ))+((XNUPF+
C*   -GAMPSI)/(XNUPFSQ+GAMMASQ))-(2./XNU))
C*       SUMSFP=SUMSFP+S*(((XNUMF+GAMMA*(GAMMA+PSIF)/XNU)/(XNUMFSQ+
C*   - GAMMASQ))+((XNUPF+GAMMA*(GAMMA-PSIF)/XNU)/(XNUPFSQ+
C*   - GAMMASQ))-(2./XNU))
C*       SUMSFPP = SUMSFPP + S * ((F*GAMMA/XNU-XNUMF*PSI)/(XNUMFSQ+
C*   -GAMMASQ)+(F*GAMMA/XNU-XNUPF*PSI)/(XNUPFSQ+GAMMASQ))
         SUMSFP = SUMSFP + S * (((XNUMF+GAMMA*((GAMMA/XNU)+PSI))/
     -(XNUMFSQ+GAMMASQ))+((XNUPF+GAMMA*((GAMMA/XNU)-PSI))/
     -(XNUPFSQ+GAMMASQ))-(2./XNU))

  101 CONTINUE


C..........BEGIN CALCULATION OF H2O RESONANCE LINE CONTRIBUTION.
C          REPEAT FOR EACH LINE UNTIL NEGATIVE FREQUENCY SENTINEL.
C****         NOTE:  SLIGHT MODIFICATION OF MODEL ON 8/31/82.
C                    LINES THAT WERE CHANGED REMAIN AS COMMENTS.
C
  199 DO 200 N=1,1000
      IF (F0W(N).LT.0) GO TO 299

         S = B1(N)*E*THETA**3.5*EXP(B2(N)*ONEMTTA)
C*       GAMMA = B3(N)*(4.8*E+P)*THETA**.6*1.E-3
         GAMMA = B3(N)*(4.80*E*THETA**0.2+P)*THETA**0.8*1.E-3

         XNU=F0W(N)
         XNUMF = XNU -F
         XNUPF = XNU+ F
         XNUMFSQ = XNUMF*XNUMF
         XNUPFSQ = XNUPF*XNUPF
         GAMMASQ = GAMMA*GAMMA
         GASQONU = GAMMASQ/XNU

C.......... NOTE: LINE SHAPE FUNCTIONS SIMPLIFIED FOR H2O BECAUSE
C           PSI = 0.
         SUMSFPP = SUMSFPP + S*F*GAMMA/XNU*(1./(XNUMFSQ+GAMMASQ)+1./(
     -XNUPFSQ+GAMMASQ))
C*       SUMSFP = SUMSFP + S*((XNUMF)/(XNUMFSQ+GAMMASQ)+(XNUPF)/
C*   -(XNUPFSQ+GAMMASQ)-(2./XNU))
         SUMSFP = SUMSFP + S*((XNUMF+GASQONU)/(XNUMFSQ+GAMMASQ)+(XNUPF+
     - GASQONU)/(XNUPFSQ+GAMMASQ)-(2./XNU))

  200 CONTINUE


C.................................................................
C        EVALUATE FAR-WING H20 CONTRIBUTION
C.................................................................
  299 FNPPV = (C1*F**X1*E*P*THETA**X2+C2*F**X3*E*E*THETA**X4)*1.0E-6


C................................................................
C       EVALUATE NON-RESONANT DRY AIR CONTRIBUTION
C.................................................................
C*    GAMMA0 = .012*(P+1.3*E)*THETA**.9
      GAMMA0 = .0056*(P+1.1*E)*THETA**.75
      GAM0SQ = GAMMA0*GAMMA0
      FNPPD = 3.1E-4*P*THETA*THETA*(2.*GAMMA0*F/(F*F+GAM0SQ)+4.2E-7*
     -F*P*THETA**.5)


C.................................................................
C       EVALUATE LIQUID WATER EXTINCTION TERM
C.................................................................
      IF (F.GT.300.) GO TO 500

C           DEBYE MODEL APPLIES FOR F < 300. GHZ
      XNUE = 2.4E4/THETA*EXP(-7.13*THETA)
      FONUE = F/XNUE
      EPSPP = (185.1 - 113./THETA)*FONUE/(1.+FONUE*FONUE)
      EPSP = 4.9 + EPSPP/FONUE
      FNPPW = 4.49 * W * EPSPP/((EPSP+2.)**2.+EPSPP*EPSPP)
      GO TO 604


C            SIMPLIFIED APPROXIMATION APPLIES FOR F > 300. GHZ
  500 FNPPW = W*.549451*F**(-.1)*THETA**(-6.)
C
C....................................................................
C------EVALUATE RAIN ATTENUATION TERM---------------------------------
C        THE MODEL OF OLSEN, ROGERS AND HODGE (AS REPORTED IN "IEEE
C        TRANSACTIONS", MARCH 1978 NUMBER 2) WAS UTILIZED TO CALCULATE
C        THE RAIN ATTENUATION TERM
C....................................................................
C
  604 CONTINUE
      ARAIN=0.
      BRAIN=0.
      ATRAN=0.
      IF (RR.EQ.0.) GO TO 507
C-------------------------ALPHA CALCULATION-------------------------
      IF (F.GE.2.9) GO TO 300
      GA=6.39E-5
      EA=2.03
      GO TO 330
  300 CONTINUE
      IF (F.GE.54.) GO TO 310
      GA=4.21E-5
      EA=2.42
      GO TO 330
  310 IF (F.GE.180.) GO TO 320
      GA=4.09E-2
      EA=0.699
      GO TO 330
  320 GA=3.38
      EA=-0.151
  330 ARAIN=GA*(F**(EA))
C
C-------------------------BETA CALCULATION-------------------------
      IF (F.GE.8.5) GO TO 340
      GB=0.851
      EB=0.158
      GO TO 370
  340 IF (F.GE.25.) GO TO 350
      GB=1.41
      EB=-0.0779
      GO TO 370
  350 IF (F.GE.164.) GO TO 360
      GB=2.63
      EB=-0.272
      GO TO 370
  360 GB=0.616
      EB=0.0126
  370 BRAIN=GB*(F**(EB))
      ATRAN=ARAIN*RR**(BRAIN)

C......................................................................
C------EFFECTIVE PATHLENGTH CALCULATIONS-------------------------------
C        THE MODEL OF STUTZMAN AND DISHMAN (AS REPORTED IN "RADIO SCIENCE,
C        VOL. 17, NUMBER 6, NOV.-DEC. 1982) WAS UTILIZED TO CALCULATE THE
C        RAIN ATTENUATION OVER A HORIZONTAL PATH.

      AMG=0.
      GMA=0.
      ZMAG=0.
C        EFFECTIVE HORIZONTAL PATHLENGTH MODE NOT USED
       GOTO 507
C      IF (ENGT.EQ.0.) GO TO 507
      IF (RR.LE.10.) GO TO 111
      AMG=ALOG(RR/10.)
      GMA=1./22.
      ZMAG=(1.-EXP(-GMA*AMG*ENGT))/(GMA*AMG*ENGT)
      ATRAN=ATRAN*ZMAG
      GO TO 222
  111 CONTINUE
      ATRAN=ATRAN*ENGT
  222 CONTINUE
      ZKE=(SUMSFPP+FNPPV+FNPPD+FNPPW)*.182*F
      ZEK=ZKE*ENGT
      FNPP=ZEK+ATRAN
      GO TO 509
C....................................................................
C       SUM THE TERMS AND RETURN RESULTS
C...................................................................
C       ATTENUATION
  507 CONTINUE
  501 FNPP = (SUMSFPP + FNPPV + FNPPD + FNPPW) * .182 * F+ ATRAN
  509 CONTINUE
C       DISPERSION, INCLUDING NONRESONANT O2
      D =(SUMSFP+3.07E-4*P*THETA*THETA*(GAMOSQ/(F*F+GAMOSQ)-1))*3.336

C      IF ((IPTHTYP.EQ.1).AND.(ENGT.NE.0.)) D=D*ENGT

      RETURN
      END
