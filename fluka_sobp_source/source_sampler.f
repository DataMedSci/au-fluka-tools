*$ CREATE SOURCE.FOR
*COPY SOURCE
*
* Useful materials on how to prepare user source can be found here:
* http://www.fluka.org/content/course/NEA/lectures/UserRoutines.pdf
*
* To compile this code, run:
* $FLUPRO/flutil/ldpm3qmd source_sampler.f -o flukadpm3_sobp
*
* In order to use it, pass it as a parameter to rfluka executable:
* rfluka -N0 -M1 -e flukadpm3_sobp your_input_file
*
* Documentation of the SOURCE is here http://www.fluka.org/fluka.php?id=man_onl&sub=71
* example source card:
* SOURCE           1.0      12.0      13.0      14.0      15.0      16.0sobp.dat
*===  source ===========================================================*
*

* TODO include description of 3 sobp.dat file formats (5,6,7 columns)
* TODO 7 columns format implementation
* TODO rescaling of sigma using mutliple scattering theory
* TODO energy reduction (nozzle exit -> virtual source)

!> @brief
!! Returns number of text blocks separated with white characters (spaces)
!! It corresponds to number of columns in CSV file.
!! @param[in] STRING
!! @retval LENGTH number of columns
      INTEGER FUNCTION NOCOLS(ASTRING)
      IMPLICIT NONE
      CHARACTER*(*) ASTRING   ! input - string of arbitrary length
      INTEGER I               !
      NOCOLS = 0

! make a loop and check for non-space followed by space
! each occurence will mark end of a column
      DO I = 1, LEN(ASTRING)-1
         IF((ASTRING(I:I) .NE. ' ') .AND.
     & (ASTRING(I+1:I+1) .EQ. ' ')) THEN
            NOCOLS = NOCOLS + 1
         END IF

! in a special case when last character isn't a space, add another column
      IF ( ASTRING(LEN(ASTRING):LEN(ASTRING)) .NE. ' ' ) THEN
        NOCOLS = NOCOLS + 1
      END IF

      END DO

      RETURN
      END


      SUBROUTINE SOURCE ( NOMORE )

      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 1990-2010      by    Alfredo Ferrari & Paola Sala  *
*     All Rights Reserved.                                             *
*                                                                      *
*                                                                      *
*     New source for FLUKA9x-FLUKA20xy:                                *
*                                                                      *
*     Created on 07 January 1990   by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on  17-Oct-10    by    Alfredo Ferrari               *
*                                                                      *
*  This is just an example of a possible user written source routine.  *
*  note that the beam card still has some meaning - in the scoring the *
*  maximum momentum used in deciding the binning is taken from the     *
*  beam momentum.  Other beam card parameters are obsolete.            *
*                                                                      *
*       Output variables:                                              *
*                                                                      *
*              Nomore = if > 0 the run will be terminated              *
*                                                                      *
*----------------------------------------------------------------------*
*
      INCLUDE '(BEAMCM)'
      INCLUDE '(FHEAVY)'
      INCLUDE '(FLKSTK)'
      INCLUDE '(IOIOCM)'
      INCLUDE '(LTCLCM)'
      INCLUDE '(PAPROP)'
      INCLUDE '(SOURCM)'
      INCLUDE '(SUMCOU)'
*
      INCLUDE '(CASLIM)'
*
*
*
      DOUBLE PRECISION ENERGY(65000), DELTAE(65000)
      DOUBLE PRECISION XPOS(65000), YPOS(65000)
      DOUBLE PRECISION FWHMX(65000), FWHMY(65000)
      DOUBLE PRECISION PART(65000)
      INTEGER NWEIGHT, NOCOLUMNS
      LOGICAL LEXISTS, LPNTSRC
      DOUBLE PRECISION SRC2SPOT
      CHARACTER(8192) LINE

      SAVE ENERGY, DELTAE, XPOS, YPOS
      SAVE FWHMX, FWHMY, PART
      SAVE NWEIGHT

      LOGICAL LFIRST
*
      SAVE LFIRST
      DATA LFIRST / .TRUE. /
*======================================================================*
*                                                                      *
*                 BASIC VERSION                                        *
*                                                                      *
*======================================================================*
      NOMORE = 0
*  +-------------------------------------------------------------------*
*  |  First call initializations:
      IF ( LFIRST ) THEN
*  |  *** The following 3 cards are mandatory ***
	 WRITE(LUNOUT,*) 'SOBP SOURCE INVOKED'
         TKESUM = ZERZER
         LFIRST = .FALSE.
         LUSSRC = .TRUE.
*  |  *** User initialization ***
* Fluka run happens in a temporary directory, created in same level as input file
* sobp.dat is not copied there, so reach one level up to get it

* Warn if sobp.dat file is missing
         inquire(file='../sobp.dat',exist=lexists)
         IF(.NOT. LEXISTS) THEN
           WRITE(LUNOUT,*) 'SOBP FILE sobp.dat missing'
         END IF

         OPEN(44, FILE = '../sobp.dat',
     $        STATUS = 'OLD')


         LINE(1:1) = '*'
* skip comment lines, first non-comment line will be saved to LINE
         DO WHILE( LINE(1:1) .EQ. '*' )
            READ(44,'(A)',END=10) LINE
         END DO

         WRITE(LUNOUT,*) 'SOBP FIRST LINE', TRIM(LINE)
         NOCOLUMNS = NOCOLS(LINE)
         WRITE(LUNOUT,*) 'SOBP NUMBER OF COLUMNS', NOCOLUMNS

         REWIND(44)

         WRITE(LUNOUT,*) 'SOBP SOURCE ZPOS fixed to', ZBEAM
         WRITE(LUNOUT,*) 'SOBP USER PARAM LIST', WHASOU

* First parameter in SOURCE indicated whether point-like or parallel source is used
* Positive number indicates point-like source
         IF ( WHASOU(1) .GT. 0.0D0 ) THEN
            WRITE(LUNOUT,*) 'SOBP POINT-LIKE VIRTUAL SOURCE'
            LPNTSRC = .TRUE.
         ELSE
            WRITE(LUNOUT,*) 'SOBP PARALLEL VIRTUAL SOURCE'
            LPNTSRC = .FALSE.
         END IF

         IF ( IJBEAM .EQ. -2 .AND. LRDBEA ) THEN
          IJHION = IPROZ  * 1000 + IPROA
          IJHION = IJHION * 100 + KXHEAV
          IONID  = IJHION
*  +-------------------------------------------------------------------*
*  |  Heavy ion:
         ELSE IF ( IJBEAM .EQ. -2 ) THEN
          IJHION = IPROZ  * 1000 + IPROA
          IJHION = IJHION * 100 + KXHEAV
          IONID  = IJHION
*  +-------------------------------------------------------------------*
*  |  Normal hadron:
         ELSE
          IONID = IJBEAM
         END IF

         NWEIGHT = 0
         WSUM = 0.0

         DO
*  fortran arrays start with 1

            NWEIGHT = NWEIGHT + 1
            IF (NWEIGHT .GT. 65000) THEN
               WRITE(LUNOUT,*) 'SOBP SOURCE ERROR: too many beamlets'
            ENDIF

            LINE(1:1) = '*'
* skip comment lines, first non-comment line will be saved to LINE
            DO WHILE( LINE(1:1) .EQ. '*' )
              READ(44,'(A)',END=10) LINE
            END DO


            IF( NOCOLUMNS .EQ. 5 ) THEN

              READ(LINE,*,END=10) ENERGY(NWEIGHT),XPOS(NWEIGHT),
     $            YPOS(NWEIGHT),FWHMX(NWEIGHT),PART(NWEIGHT)
              FWHMY = FWHMX
              DELTAE = 0.0D0

            ELSE IF ( NOCOLUMNS .EQ. 6 ) THEN

              READ(LINE,*,END=10) ENERGY(NWEIGHT),XPOS(NWEIGHT),
     $            YPOS(NWEIGHT),FWHMX(NWEIGHT),FWHMY(NWEIGHT),
     $            PART(NWEIGHT)
              DELTAE = 0.0D0

            ELSE IF ( NOCOLUMNS .EQ. 7 ) THEN

              READ(LINE,*,END=10) ENERGY(NWEIGHT), DELTAE(NWEIGHT),
     $            XPOS(NWEIGHT), YPOS(NWEIGHT),
     $            FWHMX(NWEIGHT),FWHMY(NWEIGHT),
     $            PART(NWEIGHT)

            ELSE
                WRITE(LUNOUT,*) 'INCOMPATIBLE NO OF COLUMNS ',NOCOLUMNS
            END IF

            WSUM = WSUM + PART(NWEIGHT)
*
* In sobp.dat file energy is saved in MeV/amu, while
* in Fluka we need it in not per amu, but simply in MeV
* To achieve this we multiply energy by mass number A
* Fluka doesn't have any consistent method of calculating
* mass number (number of nucleons) for simple particles
* (such as protons and alphas) and heavy ions
* we use a method IBARCH to get baryonic charge which
* reduces to mass number for particles of interest
           ENERGY(NWEIGHT) = ENERGY(NWEIGHT) * IBARCH(IONID)

         ENDDO
 10      CONTINUE
*        fix index
	 NWEIGHT = NWEIGHT - 1
         WRITE(LUNOUT,*) 'SOBP SOURCE beamlets found:', NWEIGHT
         WRITE(LUNOUT,*) 'SOBP SOURCE Particle sum (float) :', WSUM
         WRITE(LUNOUT,*) 'SOBP SOURCE TODO: particle sum is not exact.'

* check for gaussian, for future implementation
         IF ((Ldygss) .AND. (Ldxgss)) THEN
            WRITE(LUNOUT,*) 'SOBP SOURCE GAUSSIAN: TRUE'
         ELSE
            WRITE(LUNOUT,*) 'SOBP SOURCE GAUSSIAN: FALSE'
         ENDIF

      END IF

*** Sample a beamlet ****************************

      RAN = FLRNDM(111)

*     http://infohost.nmt.edu/tcc/help/lang/fortran/scaling.html
*     If you want an integer between i and j inclusive
*     use int(rand(0)*(j+1-i))+i
*     i hope hope FLRNDM [0,1[ ??
      NRAN = INT(RAN * NWEIGHT) + 1

*     If you want a real number in the interval [x,y),
*      use this expression:
*     (rand(0)*(y-x))+x

      IF ((NRAN .GT. NWEIGHT) .OR. (NRAN .LT. 1)) THEN
         WRITE(LUNOUT,*) 'SOBP SOURCE ERROR. NRAN, RAN:', NRAN, RAN
      END IF

*     No energy spread at the moment
      ENK = ENERGY(NRAN)

*     First goes point like source. Source is located at (0,0,ZBEAM)
*     whatever user provides as XBEAM and YBEAM is overriden with zeros
      IF( LPNTSRC ) THEN
        XBEAM = 0.0D0
        YBEAM = 0.0D0
*     Second goes parallel source. Whatever user provides as XBEAM and YBEAM is overriden.
*     For each spot center os spot positions goes to XBEAM and YBEAM.
      ELSE
        XBEAM = XPOS(NRAN)
        YBEAM = YPOS(NRAN)
      END IF

*     Now we set initial displacement (relative to the spot center)
*     Starting value -> sigma
      XSPOT = FWHMX(NRAN)/2.35482
      YSPOT = FWHMY(NRAN)/2.35482


*** End of beamlet sample ********************************************


*  +-------------------------------------------------------------------*
*  Push one source particle to the stack. Note that you could as well
*  push many but this way we reserve a maximum amount of space in the
*  stack for the secondaries to be generated
* Npflka is the stack counter: of course any time source is called it
* must be =0
      NPFLKA = NPFLKA + 1
* Wt is the weight of the particle
**      WTFLK  (NPFLKA) = ONEONE         set new weight
      WTFLK  (NPFLKA) = PART(NRAN)
      WEIPRI = WEIPRI + WTFLK (NPFLKA)
* Particle type (1=proton.....). Ijbeam is the type set by the BEAM
* card
*  +-------------------------------------------------------------------*
*  |  (Radioactive) isotope:
      IF ( IJBEAM .EQ. -2 .AND. LRDBEA ) THEN
         IARES  = IPROA
         IZRES  = IPROZ
         IISRES = IPROM
         CALL STISBM ( IARES, IZRES, IISRES )
         IJHION = IPROZ  * 1000 + IPROA
         IJHION = IJHION * 100 + KXHEAV
         IONID  = IJHION
         CALL DCDION ( IONID )
         CALL SETION ( IONID )
*  |
*  +-------------------------------------------------------------------*
*  |  Heavy ion:
      ELSE IF ( IJBEAM .EQ. -2 ) THEN
         IJHION = IPROZ  * 1000 + IPROA
         IJHION = IJHION * 100 + KXHEAV
         IONID  = IJHION
         CALL DCDION ( IONID )
         CALL SETION ( IONID )
         ILOFLK (NPFLKA) = IJHION
*  |  Flag this is prompt radiation
         LRADDC (NPFLKA) = .FALSE.
*  |  Group number for "low" energy neutrons, set to 0 anyway
         IGROUP (NPFLKA) = 0
*  |
*  +-------------------------------------------------------------------*
*  |  Normal hadron:
      ELSE
         IONID = IJBEAM
         ILOFLK (NPFLKA) = IJBEAM
*  |  Flag this is prompt radiation
         LRADDC (NPFLKA) = .FALSE.
*  |  Group number for "low" energy neutrons, set to 0 anyway
         IGROUP (NPFLKA) = 0
      END IF
*  |
*  +-------------------------------------------------------------------*
* From this point .....
* Particle generation (1 for primaries)
      LOFLK  (NPFLKA) = 1
* User dependent flag:
      LOUSE  (NPFLKA) = 0
*  No channeling:
      LCHFLK (NPFLKA) = .FALSE.
      DCHFLK (NPFLKA) = ZERZER
* User dependent spare variables:
      DO 100 ISPR = 1, MKBMX1
         SPAREK (ISPR,NPFLKA) = ZERZER
 100  CONTINUE
* User dependent spare flags:
      DO 200 ISPR = 1, MKBMX2
         ISPARK (ISPR,NPFLKA) = 0
 200  CONTINUE
* Save the track number of the stack particle:
      ISPARK (MKBMX2,NPFLKA) = NPFLKA
      NPARMA = NPARMA + 1
      NUMPAR (NPFLKA) = NPARMA
      NEVENT (NPFLKA) = 0
      DFNEAR (NPFLKA) = +ZERZER
* ... to this point: don't change anything
* Particle age (s)
      AGESTK (NPFLKA) = +ZERZER
      AKNSHR (NPFLKA) = -TWOTWO
* Group number for "low" energy neutrons, set to 0 anyway
      IGROUP (NPFLKA) = 0
****************************************************************

*sample a gaussian position
      CALL FLNRR2 (RGAUS1, RGAUS2)

* use center position and displacement multiplied by a
* random number sampled from gaussian distribution N(0,1)
      XFLK(NPFLKA) = XBEAM + XSPOT * RGAUS1
      YFLK(NPFLKA) = YBEAM + YSPOT * RGAUS2
      ZFLK(NPFLKA) = ZBEAM

* calculate distance from virtual beam source to the center of spot position
      SRC2SPOT = SQRT(ZBEAM*ZBEAM+XPOS(NRAN)*XPOS(NRAN)
     & + YPOS(NRAN)*YPOS(NRAN))

*      WRITE(LUNOUT,*) 'SOBP SOURCE gaussian sampled'

* Direction cosines (tx,ty,tz)
      IF( LPNTSRC ) THEN

* we use a vector from virtual source to the spot center (at isocenter plane)
* to calculate direction cosines
        TXFLK  (NPFLKA) = XPOS(NRAN) / SRC2SPOT
        TYFLK  (NPFLKA) = YPOS(NRAN) / SRC2SPOT
      ELSE

* parallel beam, direction cosines are (0,0,1)
        TXFLK  (NPFLKA) = ZERZER
        TYFLK  (NPFLKA) = ZERZER
      END IF

      WRITE(LUNOUT,*) 'SOBP SOURCE cosine X', TXFLK  (NPFLKA)
      WRITE(LUNOUT,*) 'SOBP SOURCE cosine Y', TYFLK  (NPFLKA)

* to ensure proper normalization last cosine (tz) is calculated
* from the first two (tx, ty)
      TZFLK (NPFLKA) = SQRT ( ONEONE - TXFLK(NPFLKA)*TXFLK(NPFLKA)
     & - TYFLK(NPFLKA)*TYFLK(NPFLKA))

      WRITE(LUNOUT,*) 'SOBP SOURCE cosine Z', TZFLK  (NPFLKA)

      WRITE(LUNOUT,*) 'SOBP SOURCE cosines set'
*********************************************************************
* Particle momentum
*      PMOFLK (NPFLKA) = PBEAM
*      WRITE(LUNOUT,*) 'SOBP SOURCE rest mass',AM (IONID)
       CALL FLNRRN(RGAUSS)
       PMOFLK (NPFLKA) = SQRT ( ENK* ( ENK
     &     + TWOTWO * AM (IONID) ))
     &     +DPBEAM*RGAUSS/2.35482


* Kinetic energy of the particle (GeV)
* set energy
      TKEFLK (NPFLKA) = SQRT(PMOFLK(NPFLKA)**2 + AM(IONID)**2)
     &      -AM(IONID)

*      WRITE(LUNOUT,*) 'SOBP SOURCE set ekin'


* Polarization cosines (TXPOL=-2 flag for "no polarization"):
      TXPOL  (NPFLKA) = -TWOTWO
      TYPOL  (NPFLKA) = +ZERZER
      TZPOL  (NPFLKA) = +ZERZER
*      WRITE(LUNOUT,*) 'SOBP SOURCE pol set'
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*  Calculate the total kinetic energy of the primaries: don't change
      IF ( ILOFLK (NPFLKA) .EQ. -2 .OR. ILOFLK (NPFLKA) .GT. 100000 )
     &   THEN
         TKESUM = TKESUM + TKEFLK (NPFLKA) * WTFLK (NPFLKA)
      ELSE IF ( ILOFLK (NPFLKA) .NE. 0 ) THEN
         TKESUM = TKESUM + ( TKEFLK (NPFLKA) + AMDISC (ILOFLK(NPFLKA)) )
     &          * WTFLK (NPFLKA)
      ELSE
         TKESUM = TKESUM + TKEFLK (NPFLKA) * WTFLK (NPFLKA)
      END IF
      RADDLY (NPFLKA) = ZERZER

*  Here we ask for the region number of the hitting point.
*     NREG (NPFLKA) = ...
*  The following line makes the starting region search much more
*  robust if particles are starting very close to a boundary:
      CALL GEOCRS ( TXFLK (NPFLKA), TYFLK (NPFLKA), TZFLK (NPFLKA) )
      CALL GEOREG ( XFLK  (NPFLKA), YFLK  (NPFLKA), ZFLK  (NPFLKA),
     &              NRGFLK(NPFLKA), IDISC )
*  Do not change these cards:
      CALL GEOHSM ( NHSPNT (NPFLKA), 1, -11, MLATTC )
      NLATTC (NPFLKA) = MLATTC
      CMPATH (NPFLKA) = ZERZER
      CALL SOEVSV


*      WRITE(LUNOUT,*) 'SOBP SOURCE END'
      CLOSE(44)
      RETURN
*=== End of subroutine Source =========================================*
      END