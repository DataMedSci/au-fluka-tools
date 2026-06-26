*$ CREATE SOURCE.FOR
*COPY SOURCE

!> @brief
!! Particle source for pencil beam scanning in hadrontherapy.
!! It mimics the beam coming out from gantry nozzle and aims
!! at creating spread-out Bragg peak shape in depth
!! Simulated source consists of several beamlets,
!! each of them characterized by relative weight,
!! energy and size (both defined at source location).
!! Additionally energy spread may be specified by the user.
!! User can also decide whether divergent beam is emitted from
!! a nozzle plane or point-like virtual source is used.
!! This source needs an additional file (typically sobp.dat)
!! with description of pencil beam geometry and kinematics.
!!
!! This implementation is based on a template from $FLUPRO/usermvax/source.f
!!
!! In order to use the source, first compile this file using
!! following command or Flair GUI:
!!
!!  ldpmqmd -oflukadpm_sobp source_sampler.f
!!
!! Then get a file called sobp.dat and put it in the same directory as
!! your Fluka input file. In the input file add a card called SOURCE
!! to activate this custom source. To run it, call (or use Flair):
!!
!! rfluka -N0 -M1 -e flukadpm_sobp your_input_file
!!
!!
!!
!! @details
!!
!! -------------------- SOBP CONFIG FILE -------------------------------------
!! Input file (typically sobp.dat) is a text file.
!! Comment lines start with #. It may have 5,6,7,9 or 11 columns.
!! These columns can be:
!!
!!  - 5 columns: E, X, Y, FWHM, W
!!  - 6 columns: E, X, Y, FWHM_X, FWHM_Y, W
!!  - 7 columns: E, DE, X, Y, FWHM_X, FWHM_Y, W
!!  - 9 columns: E, DE, X, Y, FWHM_X, FWHM_Y, DIVX, DIVY, W
!!  - 11 columns: E, DE, X, Y, FWHM_X, FWHM_Y, DIVX, DIVY, CORX, CORY, W
!!
!! where:
!!
!!  - E : kinetic energy in GeV/nucleon
!!  - DE : energy spread (sigma) in GeV/nucleon
!!  - X, Y : position (in cm) of the beamlet/spot center
!!  - FWHM_X, FWHM_Y, FWHM : spot size (in cm)
!!  - DIVX, DIVY : angular divergence (in mrad)
!!  - CORX, CORY : correlation coefficient rho (dimensionless)
!!
!! All above-mentioned quantities are defined at beam source location (typically nozzle exit)
!!
!! For more details, see SHIELD-HIT12A manual http://shieldhit.org/index.php?id=documentation
!!
!!
!!
!! -------------------- Beam configuration in the INPUT FILE -------------------------------------
!! Beam momentum/energy, its spread and divergence is usually specified in BEAM input card.
!! Also beam width is specified there.
!! All these values except momentum spread will be ignored and overridden.
!! In case 5- or 6-columns sobp.dat file format (the case where DE is missing) is used,
!! then beam momentum will be taken from the input card, otherwise is will be ignored and overriden.
!!
!! Center of the beam spot is usually defined in BEAMPOS card in the input file,
!! which contains X,Y,Z positions and direction cosines.
!! Whatever is specified as X,Y and direction cosines will be ignored and overridden
!! by this custom source. SDUM value in BEAMPOS card will also be ignored.
!!
!! In order to activate this custom source, please add SOURCE card to the input file.
!! WHAT(1), if nonzero, enables virtual-source geometry using SADx and SADy
!! from WHAT(3) and WHAT(4). If WHAT(1) is zero or omitted, this mode is off.
!! The SOURCE card SDUM is the filename containing the spotlist. If SDUM is
!! omitted, this routine reads sobp.dat.
!! An example SOURCE card looks like that:
!!
!! SOURCE             1.0       0.0     205.0     205.0                  sobp.dat
!!
!!
!! For more details, see FLUKA documentation:
!!  - BEAM    card http://www.fluka.org/fluka.php?id=man_onl&sub=12
!!  - BEAMPOS card http://www.fluka.org/fluka.php?id=man_onl&sub=14
!!  - SOURCE  card http://www.fluka.org/fluka.php?id=man_onl&sub=71
!!
!! -------------------------------------------------------------------------------------------------
!! ==================================================================================================
!! ==================================================================================================

!> @brief
!! Returns number of text blocks separated with white characters (spaces)
!! It corresponds to number of columns in CSV file.
!! @param[in] STRING
!! @retval LENGTH number of columns
      INTEGER FUNCTION NCOLS( ASTRING )
      IMPLICIT NONE
      CHARACTER*(*) ASTRING   ! input - string of arbitrary length
      INTEGER I               !
      NCOLS = 0

*     make a loop and check for non-space followed by space
*     each occurence will mark end of a column
      DO I = 1, LEN (ASTRING) - 1
         IF ((ASTRING (I:I) .NE. ' ') .AND.
     & (ASTRING(I+1:I+1) .EQ. ' ')) THEN
            NCOLS = NCOLS + 1
         END IF

*     in a special case when last character isn't a space, add another column
      IF ( ASTRING( LEN(ASTRING):LEN(ASTRING) ) .NE. ' ' ) THEN
        NCOLS = NCOLS + 1
      END IF

      END DO

      RETURN
      END


!! ==================================================================================================
!! ==================================================================================================

!> @brief
!! Reads beam configuration file
!! @param[in] FILEPATH  path to the beam configuration file
!! @param[out] ENERGY  particle energy in GeV/nucleon
!! @param[out] DELTAE  particle energy spread (sigma) in GeV/nucleon
!! @param[out] XPOS beam spot center (X coordinate), in cm
!! @param[out] YPOS beam spot center (Y coordinate), in cm
!! @param[out] FWHMX beam spot size in X axis, in cm
!! @param[out] FWHMY beam spot size in Y axis, in cm
!! @param[out] DIVX beam spot angular divergence in X axis, in mrad
!! @param[out] DIVY beam spot angular divergence in Y axis, in mrad
!! @param[out] CORX correlation coefficient rho(x,tx) (dimensionless)
!! @param[out] CORY correlation coefficient rho(y,ty) (dimensionless)
!! @param[out] PART beamlet weight
!! @param[out] NCOLUMNS number of columns in the file, negative if file missing or corrupted
!! @param[out] NWEIGHT number of data rows in the file, negative if file missing or corrupted
      SUBROUTINE READSOBP ( FILEPATH,
     $   ENERGY, DELTAE, XPOS, YPOS, FWHMX, FWHMY,
     $   DIVX, DIVY, CORX, CORY, PART,
     $   NCOLUMNS, NWEIGHT )

      INCLUDE 'dblprc.inc'
      INCLUDE 'dimpar.inc'
      INCLUDE 'iounit.inc'

      CHARACTER*(*) FILEPATH   ! path to the sobp file
      CHARACTER(8192) LINE
      DOUBLE PRECISION ENERGY(65000), DELTAE(65000)
      DOUBLE PRECISION XPOS(65000), YPOS(65000)
      DOUBLE PRECISION FWHMX(65000), FWHMY(65000)
      DOUBLE PRECISION DIVX(65000), DIVY(65000)
      DOUBLE PRECISION CORX(65000), CORY(65000)
      DOUBLE PRECISION PART(65000)
      INTEGER NCOLUMNS
      INTEGER NWEIGHT
      LOGICAL LEXISTS

      NWEIGHT = 0
      NCOLUMNS = 0

      WRITE(LUNOUT,*) 'SOBP SOURCE READING ', FILEPATH

*     warn if the spotlist file is missing, stop calculation
      INQUIRE( FILE=FILEPATH, EXIST=LEXISTS )
      IF ( .NOT. LEXISTS ) THEN
         WRITE(LUNOUT,*) 'SOBP FILE missing: ', FILEPATH
         RETURN
      END IF
*
*     open spotlist file for reading
      OPEN( 44, FILE=FILEPATH, STATUS='OLD' )
*
*     we will now probe the file to get the number of columns with numbers
*     we skip comment lines, first non-comment line will be saved to LINE
*     later we rewind the file to the beginning and read it with proper column format
      LINE(1:1) = '#'
      DO WHILE( LINE(1:1) .EQ. '#' )
         READ(44,'(A)',END=10) LINE
      END DO
      REWIND(44)
*     WRITE(LUNOUT,*) 'SOBP FIRST LINE', TRIM(LINE)
*
*     get number of columns from the first non-comment line
      NCOLUMNS = NCOLS(LINE)
      WRITE(LUNOUT,*) 'SOBP NUMBER OF COLUMNS', NCOLUMNS


      DO
*        fortran arrays start with 1, so we increase the counter as the loop starts
         NWEIGHT = NWEIGHT + 1
         IF (NWEIGHT .GT. 65000) THEN
            WRITE(LUNOUT,*) 'SOBP SOURCE ERROR: too many beamlets'
            RETURN
         ENDIF

*        skip comment lines, first non-comment line will be saved to LINE
         LINE(1:1) = '#'
         DO WHILE( LINE(1:1) .EQ. '#' )
            READ(44,'(A)',END=10) LINE
         END DO

*        read the line, and guess the format from number of columns
         IF ( NCOLUMNS .EQ. 5 ) THEN
            READ(LINE,*,END=10) ENERGY(NWEIGHT), XPOS(NWEIGHT),
     $         YPOS(NWEIGHT), FWHMX(NWEIGHT), PART(NWEIGHT)
            FWHMY(NWEIGHT)  = FWHMX(NWEIGHT)
            DELTAE(NWEIGHT) = 0.0D0
            DIVX(NWEIGHT)    = 0.0D0
            DIVY(NWEIGHT)    = 0.0D0
            CORX(NWEIGHT)   = 0.0D0
            CORY(NWEIGHT)   = 0.0D0

         ELSE IF ( NCOLUMNS .EQ. 6 ) THEN
            READ(LINE,*,END=10) ENERGY(NWEIGHT), XPOS(NWEIGHT),
     $         YPOS(NWEIGHT), FWHMX(NWEIGHT), FWHMY(NWEIGHT),
     $         PART(NWEIGHT)
            DELTAE(NWEIGHT) = 0.0D0
            DIVX(NWEIGHT)    = 0.0D0
            DIVY(NWEIGHT)    = 0.0D0
            CORX(NWEIGHT)   = 0.0D0
            CORY(NWEIGHT)   = 0.0D0

         ELSE IF ( NCOLUMNS .EQ. 7 ) THEN
            READ(LINE,*,END=10) ENERGY(NWEIGHT), DELTAE(NWEIGHT),
     $         XPOS(NWEIGHT), YPOS(NWEIGHT),
     $         FWHMX(NWEIGHT), FWHMY(NWEIGHT),
     $         PART(NWEIGHT)
            DIVX(NWEIGHT)    = 0.0D0
            DIVY(NWEIGHT)    = 0.0D0
            CORX(NWEIGHT)   = 0.0D0
            CORY(NWEIGHT)   = 0.0D0

         ELSE IF ( NCOLUMNS .EQ. 9 ) THEN
            READ(LINE,*,END=10) ENERGY(NWEIGHT), DELTAE(NWEIGHT),
     $         XPOS(NWEIGHT), YPOS(NWEIGHT),
     $         FWHMX(NWEIGHT), FWHMY(NWEIGHT),
     $         DIVX(NWEIGHT), DIVY(NWEIGHT),
     $         PART(NWEIGHT)
            CORX(NWEIGHT)   = 0.0D0
            CORY(NWEIGHT)   = 0.0D0

         ELSE IF ( NCOLUMNS .EQ. 11 ) THEN
            READ(LINE,*,END=10) ENERGY(NWEIGHT), DELTAE(NWEIGHT),
     $         XPOS(NWEIGHT), YPOS(NWEIGHT),
     $         FWHMX(NWEIGHT), FWHMY(NWEIGHT),
     $         DIVX(NWEIGHT), DIVY(NWEIGHT),
     $         CORX(NWEIGHT), CORY(NWEIGHT),
     $         PART(NWEIGHT)

         ELSE
            WRITE(LUNOUT,*) 'SOBP WRONG NO OF COLUMNS',NCOLUMNS
            RETURN
         END IF

      ENDDO


 10   CONTINUE
*     fix index
	  NWEIGHT = NWEIGHT - 1
      WRITE(LUNOUT,*) 'SOBP SOURCE beamlets found:', NWEIGHT

      RETURN
      END

!! ==================================================================================================
!! ==================================================================================================

!> @brief
!! Main user source subroutine
      SUBROUTINE SOURCE ( NOMORE )

      INCLUDE 'dblprc.inc'
      INCLUDE 'dimpar.inc'
      INCLUDE 'iounit.inc'
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
*  Useful materials on how to prepare user source can be found here:   *
*  http://www.fluka.org/content/course/NEA/lectures/UserRoutines.pdf   *
*                                                                      *
*       Output variables:                                              *
*                                                                      *
*              Nomore = if > 0 the run will be terminated              *
*                                                                      *
*----------------------------------------------------------------------*

      INCLUDE 'beamcm.inc'
      INCLUDE 'fheavy.inc'
      INCLUDE 'flkstk.inc'
      INCLUDE 'ioiocm.inc'
      INCLUDE 'ltclcm.inc'
      INCLUDE 'paprop.inc'
      INCLUDE 'sourcm.inc'
      INCLUDE 'sumcou.inc'
      INCLUDE 'caslim.inc'
*
*
*     containers to store data from sobp.dat file
      DOUBLE PRECISION ENERGY(65000), DELTAE(65000)
      DOUBLE PRECISION XPOS(65000), YPOS(65000)
      DOUBLE PRECISION FWHMX(65000), FWHMY(65000)
      DOUBLE PRECISION DIVX(65000), DIVY(65000)
      DOUBLE PRECISION CORX(65000), CORY(65000)
      DOUBLE PRECISION PART(65000)
*
      DOUBLE PRECISION CUMW(65000), TOTW
      SAVE CUMW, TOTW

      INTEGER NWEIGHT, NCOLUMNS
      INTEGER CDF_BINSEARCH
      LOGICAL LPNTSRC
      LOGICAL LENMOMPOS
      DOUBLE PRECISION FWHM2SIGMA
      DOUBLE PRECISION DES
      DOUBLE PRECISION XSPOT, YSPOT
      DOUBLE PRECISION SADX, SADY

      INTEGER IPOS
      INTEGER I, NRAN
      INTEGER IONA
      DOUBLE PRECISION RW, ES, RAN

      CHARACTER*256 FNAME, FPATH

      SAVE LPNTSRC
      SAVE ENERGY, DELTAE, XPOS, YPOS
      SAVE FWHMX, FWHMY, DIVX, DIVY, CORX, CORY, PART
      SAVE NWEIGHT

      LOGICAL LFIRST

      SAVE LFIRST
      DATA LFIRST / .TRUE. /
*
      NOMORE = 0
*  +-------------------------------------------------------------------*
*  |  First call initializations:
      IF ( LFIRST ) THEN
*  |  *** The following 3 cards are mandatory ***
         TKESUM = ZERZER
         LFIRST = .FALSE.
         LUSSRC = .TRUE.
*  |  *** User initialization ***

* set default file name for sobp.dat if not provided by user in SOURCE card
         FNAME = TRIM(SDUSOU)
         IF ( FNAME.EQ. '' ) THEN
             FNAME = 'sobp.dat'
         END IF

         FWHM2SIGMA = 2.0D0 * SQRT( 2.0D0 * LOG(2.0D0))
         WRITE(LUNOUT,*) 'SPOTLIST SOURCE ZPOS fixed to', ZBEAM
         WRITE(LUNOUT,*) 'SPOTLIST USER PARAM LIST', WHASOU
         WRITE(LUNOUT,*) 'SPOTLIST FNAME', FNAME

         FPATH = '../' // TRIM(FNAME)

*        Fluka run happens in a temporary directory,
*        created in same level as input file
*        sobp.dat is not copied there,
*        We reach one level up to get it via ../sobp.dat
         CALL READSOBP ( TRIM(FPATH), ENERGY, DELTAE,
     $            XPOS, YPOS, FWHMX, FWHMY,
     $            DIVX, DIVY, CORX, CORY, PART, NCOLUMNS, NWEIGHT )

*        In case of problem with reading sobp.dat file
         IF ( (NCOLUMNS .LE. ZERZER) .OR. (NWEIGHT .LE. ZERZER)) THEN
            NOMORE = 1
            RETURN
         ENDIF

*        Build cumulative weights only over strictly positive weights
*        and compact beamlet arrays in-place to exclude zero-weight rows.
         IPOS = 0
         TOTW = 0.0D0
         DO I = 1, NWEIGHT
            IF (PART(I) .LT. 0.0D0) PART(I) = 0.0D0
            IF (PART(I) .GT. 0.0D0) THEN
               IPOS = IPOS + 1
*              Compact all parameter arrays to keep only positive-weight rows
               ENERGY(IPOS) = ENERGY(I)
               DELTAE(IPOS) = DELTAE(I)
               XPOS(IPOS)   = XPOS(I)
               YPOS(IPOS)   = YPOS(I)
               FWHMX(IPOS)  = FWHMX(I)
               FWHMY(IPOS)  = FWHMY(I)
               DIVX(IPOS)    = DIVX(I)
               DIVY(IPOS)    = DIVY(I)
               CORX(IPOS)   = CORX(I)
               CORY(IPOS)   = CORY(I)
               PART(IPOS)   = PART(I)
               TOTW = TOTW + PART(IPOS)
               CUMW(IPOS) = TOTW
            END IF
         END DO
*        Update NWEIGHT to the number of strictly positive-weight rows
         NWEIGHT = IPOS

         IF (TOTW .LE. 0.0D0) THEN
            WRITE(LUNOUT,*) 'SOBP SOURCE ERROR: total weight <= 0'
            NOMORE = 1
            RETURN
         END IF

*        First parameter in SOURCE enables virtual-source geometry.
*        If the number is present and non-zero we expect point-like source
         SADX = WHASOU(3)
         SADY = WHASOU(4)
         IF ( WHASOU(1) .NE. 0.0D0 ) THEN
            WRITE(LUNOUT,*) 'SOBP POINT-LIKE VIRTUAL SOURCE'
            IF ( SADX .LE. 0.0D0 .OR. SADY .LE. 0.0D0 ) THEN
               WRITE(LUNOUT,*) 'SOBP SOURCE ERROR: point-like source'
               WRITE(LUNOUT,*) 'requires positive SADx and SADy, got',
     &            SADX, SADY
               NOMORE = 6
               RETURN
            END IF
            LPNTSRC = .TRUE.
         ELSE
            WRITE(LUNOUT,*) 'SOBP PARALLEL VIRTUAL SOURCE'
            LPNTSRC = .FALSE.
         END IF

      END IF
*  |  End of first call initializations.
*  +-------------------------------------------------------------------*

*  +-------------------------------------------------------------------*
*  |  sample beamlet collection :


*     Get a 64-bit pseudo random number in the interval [0.D+00,1.D+00),
*     1 being not included
      RAN = FLRNDM(111)
      RW  = RAN * TOTW

      NRAN = CDF_BINSEARCH( CUMW, NWEIGHT, RW )

*     check if random number selected properly
      IF ((NRAN .GT. NWEIGHT) .OR. (NRAN .LT. 1)) THEN
         WRITE(LUNOUT,*) 'SOBP SOURCE ERROR. NRAN, RAN:', NRAN, RAN
         NOMORE = 3
         RETURN
      END IF

*  +-------------------------------------------------------------------*
*  Push one source particle to the stack. Note that you could as well
*  push many but this way we reserve a maximum amount of space in the
*  stack for the secondaries to be generated
* Npflka is the stack counter: of course any time source is called it
* must be =0
      NPFLKA = NPFLKA + 1
* Wt is the weight of the particle
      WTFLK (NPFLKA) = 1.0D0   ! particles were already sampled per weight.
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
         IONA   = IPROA
         CALL DCDION ( IONID )
         CALL SETION ( IONID )
*  |
*  +-------------------------------------------------------------------*
*  |  Heavy ion:
      ELSE IF ( IJBEAM .EQ. -2 ) THEN
         IJHION = IPROZ  * 1000 + IPROA
         IJHION = IJHION * 100 + KXHEAV
         IONID  = IJHION
         IONA   = IPROA
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
         IONA  = IBARCH(IONID)
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
****************************************************************



*  +-------------------------------------------------------------------*
*  |  Particle momentum and energy

*      ENERGY and DELTAE from sobp.dat are kinetic energy and sigma
*      in GeV/nucleon. Convert them to total kinetic energy in GeV
*      by multiplying by the integer mass number A, not by the
*      physical ion mass in amu.
        IF (IONA .LE. 0) THEN
           WRITE(LUNOUT,*) 'SOBP SOURCE ERROR: invalid mass number',
     &        IONA, ' for particle ', IONID
           NOMORE = 5
           RETURN
        END IF
        ES   = ENERGY(NRAN)  * DBLE(IONA)
        DES  = DELTAE(NRAN)  * DBLE(IONA)
        IF (DES .LT. 0.0D0) DES = 0.0D0

*  ....................................................................................
*      Basic equation which relates particle energy and momentum is:
*
*         p^2 c^2 + m0^2 c^4 = E^2                 (1)
*
*      Total energy E can be written as the sum of kinetic energy T and rest energy:
*
*         p^2 c^2 + m0^2 c^4 = (T + m0 c^2)^2      (2)
*
*      By doing simple math we can extract from this formula following relations:
*
*         T = sqrt(p^2 c^2 + m0^2 c^4) - m0 c^2    (3)
*         p c = sqrt( E (E + 2 m0 c^2) )           (4)
*
*      We use following variables below:
*        @   PMOFLK - particle momentum, in GeV/c
*        @   AM(IONID) - particle rest energy, in GeV
*        @   ES - total particle kinetic energy, in GeV
*  ....................................................................................
*
*      There is always a small chance that randomly sampled energy or momentum
*      will get negative (i.e. mean energy 10 MeV with energy spread 8 MeV).
*      We will sample gaussian distribution in a loop until positive value is obtained
*
       LENMOMPOS = .FALSE.
       DO WHILE ( LENMOMPOS .EQV. .FALSE. )

          IF ( ES .LE. ZERZER ) THEN
             WRITE(LUNOUT,*) 'SOBP NEGATIVE EN:', ES
             NOMORE = 4
             RETURN
          END IF

*         Lets get a normally distributed random number RGAUSS
          CALL FLNRRN(RGAUSS)

          IF ( NCOLUMNS .GE. 7 ) THEN

*            In case sobp.dat file has 7 columns, we expect than
*            energy spread will be there, set as standard deviation
*            First lets sample gaussian energy distribution.
*            Mean energy is set to ES, value from sobp.dat file
*            Standard deviation is also taken from sobp.dat file
             TKEFLK (NPFLKA) = ES + DES * RGAUSS

*            Exit the loop if kinetic energy is positive (momentum will also be positive)
             IF ( TKEFLK (NPFLKA) .GT. ZERZER ) THEN
                LENMOMPOS = .TRUE.
             ELSE
                WRITE(LUNOUT,*) 'KINETIC ENERGY:', TKEFLK (NPFLKA)
                WRITE(LUNOUT,*) 'RESAMPLING TO GET POSITIVE NUMBER'
             END IF

*            Momentum of the particle, according to eq (4)
             PMOFLK (NPFLKA) = SQRT ( TKEFLK (NPFLKA)*
     &        ( TKEFLK (NPFLKA) + TWOTWO * AM (IONID) ))

          ELSE

*            In case energy spread is not provided in sobp.dat file,
*            we can use a momentum spread from BEAM input card (which by default is set to 0).
*            Mean momentum is calculated from kinetic energy of the spot, using eq (4)
*            Standard deviation is obtained from FWHM value set by user in the input card (DPBEAM)
             PMOFLK (NPFLKA) = SQRT ( ES *
     &          ( ES + TWOTWO * AM (IONID) ))
     &          + DPBEAM*RGAUSS/FWHM2SIGMA

*
*            Exit the loop if momentum is positive (kinetic energy will also be positive)
             IF ( PMOFLK (NPFLKA) .GT. ZERZER ) THEN
                LENMOMPOS = .TRUE.
             ELSE
                WRITE(LUNOUT,*) 'MOMENTUM:', PMOFLK (NPFLKA)
                WRITE(LUNOUT,*) 'RESAMPLING TO GET POSITIVE NUMBER'
             END IF

*            Kinetic energy of the particle (GeV), according to eq (3)
             TKEFLK (NPFLKA) = SQRT(PMOFLK(NPFLKA)**2 + AM(IONID)**2)
     &         -AM(IONID)
*

          END IF

       END DO

*      debugging printouts only if requested by user
       IF ( WHASOU(2) .NE. 0.0D0 ) THEN
          WRITE(LUNOUT,*) 'SOBP E:', 1.D3*ES / DBLE(IONA),
     &       'MeV/nucleon'
          WRITE(LUNOUT,*) 'SOBP Z:', ZBEAM, 'cm'
          WRITE(LUNOUT,*) 'SOBP SOURCE EKIN', TKEFLK(NPFLKA)
          WRITE(LUNOUT,*) 'SOBP SOURCE MOMENTUM', PMOFLK(NPFLKA)
       ENDIF

*  +-------------------------------------------------------------------*
*  |  Polarization cosines (TXPOL=-2 flag for "no polarization"):

      TXPOL  (NPFLKA) = -TWOTWO
      TYPOL  (NPFLKA) = +ZERZER
      TZPOL  (NPFLKA) = +ZERZER

*  *  +-------------------------------------------------------------------*
*  |  Particle coordinates + phase space (local beam frame)
      XSPOT = XPOS(NRAN)
      YSPOT = YPOS(NRAN)

*     Interpret XSPOT/YSPOT as SHIELD-HIT reference-plane spot
*     coordinates [cm], not as particle birth coordinates at BEAMPOS.
*     For point-like virtual source/SAD mode, back-project the birth
*     point to the FLUKA source plane ZBEAM [cm]. This makes the
*     SAD-steered central ray pass through the requested spot coordinate
*     near the SHIELD-HIT reference plane z = 0 cm.
      XBEAM = XSPOT
      YBEAM = YSPOT
      IF ( LPNTSRC .AND. SADX .GT. 0.0D0 ) THEN
         XBEAM = XSPOT - (XSPOT / SADX) * (0.0D0 - ZBEAM)
      ENDIF
      IF ( LPNTSRC .AND. SADY .GT. 0.0D0 ) THEN
         YBEAM = YSPOT - (YSPOT / SADY) * (0.0D0 - ZBEAM)
      ENDIF

*     Sample (X,Y,TX,TY) in the local frame (mean angles = 0 here)
      CALL SAMPLE_PHASESPACE(
     &   XBEAM, YBEAM, FWHMX(NRAN), FWHMY(NRAN),
     &   DIVX(NRAN),  DIVY(NRAN),  CORX(NRAN), CORY(NRAN),
     &   FWHM2SIGMA,
     &   XFLK(NPFLKA), YFLK(NPFLKA),
     &   TXFLK(NPFLKA), TYFLK(NPFLKA) )


*     Optional scanning steering (flag-controlled)
      IF ( LPNTSRC ) THEN
         CALL APPLY_SAD_TILT(
     &      XSPOT, YSPOT,
     &      SADX, SADY,
     &      TXFLK(NPFLKA), TYFLK(NPFLKA) )
      ENDIF

*     Renormalize direction
      IF (TXFLK(NPFLKA)**2 + TYFLK(NPFLKA)**2 .GE. 1.0D0) THEN
         TXFLK(NPFLKA) = 0.0D0
         TYFLK(NPFLKA) = 0.0D0
      ENDIF
      TZFLK(NPFLKA) = SQRT(1.0D0
     &   - TXFLK(NPFLKA)**2 - TYFLK(NPFLKA)**2 )

      ZFLK(NPFLKA) = ZBEAM


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
      RETURN
*=== End of subroutine Source =========================================*
      END


      SUBROUTINE SAMPLE_PHASESPACE(
     &   X0, Y0, FWHMX, FWHMY, DIVX, DIVY, CORX, CORY,
     &   FWHM2SIGMA,
     &   X, Y, TX, TY )

      IMPLICIT NONE
      DOUBLE PRECISION X0, Y0, FWHMX, FWHMY, DIVX, DIVY, CORX, CORY
      DOUBLE PRECISION FWHM2SIGMA
      DOUBLE PRECISION X, Y, TX, TY

      DOUBLE PRECISION SIGX, SIGY, SIGTX, SIGTY
      DOUBLE PRECISION RHOX, RHOY, COVX, COVY
      DOUBLE PRECISION G1, G2, G3, G4
      DOUBLE PRECISION DX, DY, DTX, DTY
      DOUBLE PRECISION VARX, VARY

*     Convert spot size FWHM -> sigma (cm)
      SIGX = 0.0D0
      SIGY = 0.0D0
      IF (FWHM2SIGMA .GT. 0.0D0) THEN
         SIGX = FWHMX / FWHM2SIGMA
         SIGY = FWHMY / FWHM2SIGMA
      ENDIF
      IF (SIGX .LT. 0.0D0) SIGX = -SIGX
      IF (SIGY .LT. 0.0D0) SIGY = -SIGY

*     Angular spreads: DIVX/DIVY are in mrad -> convert to rad
      SIGTX = DIVX * 1.0D-3
      SIGTY = DIVY * 1.0D-3
      IF (SIGTX .LT. 0.0D0) SIGTX = -SIGTX
      IF (SIGTY .LT. 0.0D0) SIGTY = -SIGTY

*     CORX/CORY are correlation coefficients rho in [-1,1]
      RHOX = CORX
      RHOY = CORY

*     Clamp rho for numerical safety
      IF (RHOX .GT.  0.999999D0) RHOX =  0.999999D0
      IF (RHOX .LT. -0.999999D0) RHOX = -0.999999D0
      IF (RHOY .GT.  0.999999D0) RHOY =  0.999999D0
      IF (RHOY .LT. -0.999999D0) RHOY = -0.999999D0

*     Convert correlation -> covariance (cm*rad)
      COVX = RHOX * SIGX * SIGTX
      COVY = RHOY * SIGY * SIGTY

*     Draw 4 independent standard normals
      CALL FLNRR2(G1, G2)
      CALL FLNRR2(G3, G4)

*     X-plane: correlated (DX, DTX)
      DX = SIGX * G1

      IF (SIGX .GT. 0.0D0) THEN
         VARX = SIGTX*SIGTX - (COVX*COVX)/(SIGX*SIGX)
         IF (VARX .LT. 0.0D0) VARX = 0.0D0
         DTX = (COVX / SIGX) * G1 + SQRT(VARX) * G2
      ELSE
*        if SIGX=0, correlation is meaningless -> just angle spread
         DTX = SIGTX * G2
      ENDIF

*     Y-plane: correlated (DY, DTY)
      DY = SIGY * G3

      IF (SIGY .GT. 0.0D0) THEN
         VARY = SIGTY*SIGTY - (COVY*COVY)/(SIGY*SIGY)
         IF (VARY .LT. 0.0D0) VARY = 0.0D0
         DTY = (COVY / SIGY) * G3 + SQRT(VARY) * G4
      ELSE
         DTY = SIGTY * G4
      ENDIF

*     Apply to particle (local frame, mean angles = 0)
      X  = X0 + DX
      Y  = Y0 + DY
      TX = DTX
      TY = DTY

      RETURN
      END


      INTEGER FUNCTION CDF_BINSEARCH( CUMW, N, X )
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION CUMW(N), X
      INTEGER LO, HI, MID

*     For N <= 1, always return the first (and only) bin
      IF (N .LE. 1) THEN
         CDF_BINSEARCH = 1
         RETURN
      ENDIF

*     Find the first index with CUMW(i) > X (strictly greater).
*     This is a binary search over the monotonic cumulative weights.
      LO = 1
      HI = N
      DO WHILE (LO .LT. HI)
         MID = (LO + HI) / 2
         IF (CUMW(MID) .GT. X) THEN
            HI = MID
         ELSE
            LO = MID + 1
         ENDIF
      END DO

      CDF_BINSEARCH = LO
      RETURN
      END


      SUBROUTINE APPLY_SAD_TILT(XC, YC, SADX, SADY, TX, TY)
      IMPLICIT NONE
      DOUBLE PRECISION XC, YC, SADX, SADY, TX, TY

*     Apply point-like virtual-source angular correction.
*     SADX and SADY are positive downstream source-axis distances [cm].
*     XSPOT/YSPOT are reference-plane coordinates, so positive X/Y spots
*     require positive TX/TY contributions for a +z travelling beam.
      IF (SADX .GT. 0.0D0) TX = TX + XC / SADX
      IF (SADY .GT. 0.0D0) TY = TY + YC / SADY

      RETURN
      END
