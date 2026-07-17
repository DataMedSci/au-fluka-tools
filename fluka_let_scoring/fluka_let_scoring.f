* FLUKA LET scoring routines for averaged LET-moment scoring
*======================================================================
*
* Reference:
* Kalholm F, Grzanka L, Traneus E, Bassler N. A systematic review on
* the usage of averaged LET in radiation biology for particle therapy.
* Radiotherapy and Oncology. 2021 Aug 1;161:211-21.
*
* This file implements FLUSCW and COMSCW scoring weights for LET-moment,
* fluence-filter, and dose-filter scorers.
*

*                                                                      *
*=== fluscw ===========================================================*
*                                                                      *
      DOUBLE PRECISION FUNCTION FLUSCW ( IJ    , PLA   , TXX   , TYY   ,
     &                                   TZZ   , WEE   , XX    , YY    ,
     &                                   ZZ    , NREG  , IOLREG, LLO   ,
     &                                   NSURF )

      INCLUDE 'dblprc.inc'
      INCLUDE 'dimpar.inc'
      INCLUDE 'iounit.inc'
      INCLUDE 'scohlp.inc'
      INCLUDE 'usrbin.inc'
      INCLUDE 'flkmat.inc'
      INCLUDE 'trackr.inc'
      INCLUDE 'paprop.inc'
      INCLUDE 'fheavy.inc'
      

      DOUBLE PRECISION GETLET
      DOUBLE PRECISION EKIN, LETW, LETLIN, SUMT, SUMD
      INTEGER MATLET, IHEAV, II, MWATER
      LOGICAL LETFRS
      CHARACTER*8 SCONAM

C
C     Water-equivalent material index for the water-reference LET
C     scorers (PAW1/PAW2, P1W1/P1W2, ALW1/ALW2). Resolved once, on the
C     first scoring call, into MWATER. There are no hardcoded material
C     numbers anywhere in this routine.
C
C     Local-material scorers use MEDFLK(NREG,1) directly, i.e. LET is
C     evaluated in whatever material the particle is currently in.
C     Water-reference scorers use MWATER instead.
C
C     MWATER is resolved as:
C        1. MATQLT -- FLUKA's built-in "extra water material for Q(L)
C           calculations" (flkmat.inc). This exists for dose-equivalent
C           / quality-factor scoring and is available even when the
C           input deck defines no explicit WATER material.
C        2. otherwise, the first material named 'WATER' in MATNAM.
C        3. otherwise MWATER stays .LE. 0 and the water-reference
C           scorers return zero (a warning is written to LUNOUT).
C
C     The resolved value is written to the FLUKA output (LUNOUT) on the
C     first scoring call so it can be verified.
C
      SAVE MWATER, LETFRS
      DATA MWATER / -1 /
      DATA LETFRS / .TRUE. /


      FLUSCW = ONEONE
      LSCZER = .FALSE.
      SCONAM = TRIM(ADJUSTL(TITUSB(JSCRNG)))

      IF ( LETFRS ) THEN
         MWATER = MATQLT
         IF ( MWATER .LE. 0 ) THEN
            DO II = 1, NMAT
               IF ( MWATER .LE. 0 .AND.
     &              MATNAM(II)(1:5) .EQ. 'WATER' ) MWATER = II
            END DO
         END IF
         IF ( MWATER .LE. 0 ) THEN
            WRITE(LUNOUT,*) ' fluka_let_scoring: WARNING no WATER',
     &                      ' material found; water scorers return 0'
         ELSE
            WRITE(LUNOUT,*) ' fluka_let_scoring: water MWATER =', MWATER
         END IF
         LETFRS = .FALSE.
      END IF

C
C     Scorer-key naming convention:
C        The scorer identifiers are kept to four characters because the
C        FLUKA USRBIN/AUXSCORE workflow used here relies on four-character
C        scorer keys.
C
C     ------------------------------------------------------------------
C     Light-fragment LET weighting branches for FLUSCW.
C
C     These branches score LET-weighted fluence contributions for light
C     charged fragments transported by FLUKA.
C
C     IJ is the FLUKA particle identifier passed to FLUSCW/COMSCW. The
C     light-fragment IJ values used here are:
C
C        IJ = -3    deuteron, 2H
C        IJ = -4    triton, 3H
C        IJ = -5    helium-3, 3He
C        IJ = -6    helium-4 / alpha, 4He
C
C     Scorer-key convention:
C
C        D2L1 / D2L2    deuteron LET / LET^2
C        T3L1 / T3L2    triton LET / LET^2
C        H3L1 / H3L2    helium-3 LET / LET^2
C        H4L1 / H4L2    helium-4 LET / LET^2
C
C     The L1 scorers return the first raw LET moment: a contribution
C     weighted by LET [keV/um].
C     The L2 scorers return the second raw LET moment: a contribution
C     weighted by LET^2 [(keV/um)^2].
C
C     Unlike the Li-6/Li-7 branches below, these light fragments use
C     GETLET directly. Lithium is handled separately because the tested
C     GETLET calls returned zero for transported lithium ions in this
C     implementation; Li-6/Li-7 are therefore handled through FLUKA's
C     heavy-fragment bookkeeping and reconstructed from TRACKR quantities.
C
C     MATLET is taken from MEDFLK(NREG,1), i.e. the material assigned to
C     the current FLUKA region, so LET is evaluated in whatever material
C     the particle is currently in. The only guard applied is to skip
C     vacuum / non-material regions (MATLET .LE. 0 or RHO .LE. 0), where
C     GETLET/RHO would be meaningless. Spatial restriction of the scoring
C     is the job of the USRBIN geometry, not of a hardcoded material list.
C
C     These branches classify the particle currently being transported.
C     They do not record where the fragment was produced or which parent
C     particle produced it. That would require STUPRF or MDSTCK ancestry
C     tagging.
C     ------------------------------------------------------------------

C     Deuteron LET weighting branch for fluence-type USRBIN.
C     Scorer name first four characters: D2L1.

      IF ( SCONAM(1:4) .EQ. 'D2L1' ) THEN
         FLUSCW = ZERZER

         IF ( IJ .NE. -3 ) THEN
            RETURN
         END IF
         EKIN = -PLA
         IF ( EKIN .LE. 1.0D-09 ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF

         MATLET = MEDFLK(NREG,1)
         IF ( MATLET .LE. 0 .OR. RHO(MATLET) .LE. ZERZER ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF

         LETW = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)
         LETLIN = RHO(MATLET) * LETW
         FLUSCW = LETLIN
         RETURN
      END IF
C     Deuteron LET^2 weighting branch for DLET numerator.
C     Scorer name first four characters: D2L2.

      IF ( SCONAM(1:4) .EQ. 'D2L2' ) THEN
         FLUSCW = ZERZER

         IF ( IJ .NE. -3 ) THEN
            RETURN
         END IF

         EKIN = -PLA
         IF ( EKIN .LE. 1.0D-09 ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF

         MATLET = MEDFLK(NREG,1)
         IF ( MATLET .LE. 0 .OR. RHO(MATLET) .LE. ZERZER ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF

         LETW = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)
         LETLIN = RHO(MATLET) * LETW
         FLUSCW = LETLIN * LETLIN
         RETURN
      END IF

C     Triton LET weighting branch for fluence-type USRBIN.
C     Scorer name first four characters: T3L1.

      IF ( SCONAM(1:4) .EQ. 'T3L1' ) THEN
         FLUSCW = ZERZER

         IF ( IJ .NE. -4 ) THEN
            RETURN
         END IF
         EKIN = -PLA
         IF ( EKIN .LE. 1.0D-09 ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF

         MATLET = MEDFLK(NREG,1)
         IF ( MATLET .LE. 0 .OR. RHO(MATLET) .LE. ZERZER ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF

         LETW = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)
         LETLIN = RHO(MATLET) * LETW
         FLUSCW = LETLIN
         RETURN
      END IF
C     Triton LET^2 weighting branch for DLET numerator.
C     Scorer name first four characters: T3L2.

      IF ( SCONAM(1:4) .EQ. 'T3L2' ) THEN
         FLUSCW = ZERZER

         IF ( IJ .NE. -4 ) THEN
            RETURN
         END IF

         EKIN = -PLA
         IF ( EKIN .LE. 1.0D-09 ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF

         MATLET = MEDFLK(NREG,1)
         IF ( MATLET .LE. 0 .OR. RHO(MATLET) .LE. ZERZER ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF

         LETW = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)
         LETLIN = RHO(MATLET) * LETW
         FLUSCW = LETLIN * LETLIN
         RETURN
      END IF

C     Helium-3 LET weighting branch for fluence-type USRBIN.
C     Scorer name first four characters: H3L1.

      IF ( SCONAM(1:4) .EQ. 'H3L1' ) THEN
         FLUSCW = ZERZER

         IF ( IJ .NE. -5 ) THEN
            RETURN
         END IF
         EKIN = -PLA
         IF ( EKIN .LE. 1.0D-09 ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF

         MATLET = MEDFLK(NREG,1)
         IF ( MATLET .LE. 0 .OR. RHO(MATLET) .LE. ZERZER ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF

         LETW = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)
         LETLIN = RHO(MATLET) * LETW
         FLUSCW = LETLIN
         RETURN
      END IF
C     Helium-3 LET^2 weighting branch for H3LET numerator.
C     Scorer name first four characters: H3L2.

      IF ( SCONAM(1:4) .EQ. 'H3L2' ) THEN
         FLUSCW = ZERZER

         IF ( IJ .NE. -5 ) THEN
            RETURN
         END IF

         EKIN = -PLA
         IF ( EKIN .LE. 1.0D-09 ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF

         MATLET = MEDFLK(NREG,1)
         IF ( MATLET .LE. 0 .OR. RHO(MATLET) .LE. ZERZER ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF

         LETW = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)
         LETLIN = RHO(MATLET) * LETW
         FLUSCW = LETLIN * LETLIN
         RETURN
      END IF

C     Helium-4 / alpha LET weighting branch for fluence-type USRBIN.
C     Scorer name first four characters: H4L1.

      IF ( SCONAM(1:4) .EQ. 'H4L1' ) THEN
         FLUSCW = ZERZER

         IF ( IJ .NE. -6 ) THEN
            RETURN
         END IF
         EKIN = -PLA
         IF ( EKIN .LE. 1.0D-09 ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF
         MATLET = MEDFLK(NREG,1)
         IF ( MATLET .LE. 0 .OR. RHO(MATLET) .LE. ZERZER ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF

         LETW = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)
         LETLIN = RHO(MATLET) * LETW
         FLUSCW = LETLIN
         RETURN
      END IF
C     Helium-4 / alpha LET^2 weighting branch for H4LET numerator.
C     Scorer name first four characters: H4L2.

      IF ( SCONAM(1:4) .EQ. 'H4L2' ) THEN
         FLUSCW = ZERZER

         IF ( IJ .NE. -6 ) THEN
            RETURN
         END IF

         EKIN = -PLA
         IF ( EKIN .LE. 1.0D-09 ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF

         MATLET = MEDFLK(NREG,1)
         IF ( MATLET .LE. 0 .OR. RHO(MATLET) .LE. ZERZER ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF

         LETW = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)
         LETLIN = RHO(MATLET) * LETW
         FLUSCW = LETLIN * LETLIN
         RETURN
      END IF

C     ------------------------------------------------------------------
C     Li-6 LET weighting branch for FLUSCW track-length scoring.
C
C     Scorer key:
C        SCONAM = 'L6L1'
C
C     Physics meaning:
C        Keep only transported lithium-6 fragments, identified as
C        charge Z = 3 and mass A = 6, and return the first raw LET moment.
C        This contributes to a track-length-weighted LET numerator.
C
C     FLUKA bookkeeping:
C        JTRACK .LT. -6       means current particle is transported as a
C                             heavy ion / nuclear fragment.
C        NPHEAV .GT. 0        means an entry in FHEAVY is available.
C        KHEAVY(NPHEAV)       maps the current heavy fragment to an
C                             isotope table index IHEAV.
C        ICHEAV(IHEAV)        is fragment charge number Z.
C        IBHEAV(IHEAV)        is fragment mass number A.
C
C     LET reconstruction:
C        GETLET is not used for Li here, because the GETLET call returns
C        zero for these transported lithium ions in this implementation.
C        Instead, local LET is reconstructed from TRACKR step data:
C
C           LET [keV/um] = 100 * SUMD [GeV] / SUMT [cm]
C
C        because 1 GeV/cm = 100 keV/um.
C     ------------------------------------------------------------------


      IF ( SCONAM(1:4) .EQ. 'L6L1' ) THEN
         FLUSCW = ZERZER

         IF ( JTRACK .LT. -6 .AND. NPHEAV .GT. 0 ) THEN
            IHEAV = KHEAVY(NPHEAV)

            IF ( IHEAV .GE. 1 .AND. IHEAV .LE. KXHEAV ) THEN
               IF ( ICHEAV(IHEAV) .EQ. 3 .AND.
     &              IBHEAV(IHEAV) .EQ. 6 ) THEN

                  SUMT = ZERZER
                  DO II = 1, NTRACK
                     SUMT = SUMT + TTRACK(II)
                  END DO

                  SUMD = ZERZER
                  DO II = 1, MTRACK
                     SUMD = SUMD + DTRACK(II)
                  END DO

                  IF ( SUMT .GT. ZERZER ) THEN
                     LETW = 100.0D0 * SUMD / SUMT
                     FLUSCW = LETW
                  END IF
               END IF
            END IF
         END IF

         RETURN
      END IF
C     ------------------------------------------------------------------
C     Li-6 LET^2 weighting branch for FLUSCW track-length scoring.
C
C     Scorer key:
C        SCONAM = 'L6L2'
C
C     Physics meaning:
C        Keep only transported lithium-6 fragments, identified as
C        charge Z = 3 and mass A = 6, and return LET^2. Together with
C        the L6L1 scorer, this allows reconstruction of a dose-like
C        or LET-weighted mean for Li-6, depending on the post-processing
C        denominator used.
C
C     Unit reconstruction is identical to L6L1:
C        LET [keV/um] = 100 * SUMD [GeV] / SUMT [cm].
C     ------------------------------------------------------------------

      IF ( SCONAM(1:4) .EQ. 'L6L2' ) THEN
         FLUSCW = ZERZER

         IF ( JTRACK .LT. -6 .AND. NPHEAV .GT. 0 ) THEN
            IHEAV = KHEAVY(NPHEAV)

            IF ( IHEAV .GE. 1 .AND. IHEAV .LE. KXHEAV ) THEN
               IF ( ICHEAV(IHEAV) .EQ. 3 .AND.
     &              IBHEAV(IHEAV) .EQ. 6 ) THEN

                  SUMT = ZERZER
                  DO II = 1, NTRACK
                     SUMT = SUMT + TTRACK(II)
                  END DO

                  SUMD = ZERZER
                  DO II = 1, MTRACK
                     SUMD = SUMD + DTRACK(II)
                  END DO

                  IF ( SUMT .GT. ZERZER ) THEN
                     LETW = 100.0D0 * SUMD / SUMT
                     FLUSCW = LETW * LETW
                  END IF
               END IF
            END IF
         END IF

         RETURN
      END IF


C     ------------------------------------------------------------------
C     Li-7 LET weighting branch for FLUSCW track-length scoring.
C
C     Scorer key:
C        SCONAM = 'L7L1'
C
C     Physics meaning:
C        Keep only transported lithium-7 fragments, identified as
C        charge Z = 3 and mass A = 7, and return the first raw LET moment.
C        This is the Li-7 analogue of the L6L1 branch.
C
C     FLUKA isotope identification:
C        JTRACK .LT. -6       selects transported heavy ions/fragments.
C        NPHEAV .GT. 0        requires a valid FHEAVY fragment entry.
C        KHEAVY(NPHEAV)       gives the fragment table index IHEAV.
C        ICHEAV(IHEAV) = 3    requires lithium charge Z = 3.
C        IBHEAV(IHEAV) = 7    requires lithium-7 mass A = 7.
C
C     LET reconstruction:
C        LET [keV/um] = 100 * SUMD [GeV] / SUMT [cm].
C     ------------------------------------------------------------------


      IF ( SCONAM(1:4) .EQ. 'L7L1' ) THEN
         FLUSCW = ZERZER

         IF ( JTRACK .LT. -6 .AND. NPHEAV .GT. 0 ) THEN
            IHEAV = KHEAVY(NPHEAV)

            IF ( IHEAV .GE. 1 .AND. IHEAV .LE. KXHEAV ) THEN
               IF ( ICHEAV(IHEAV) .EQ. 3 .AND.
     &              IBHEAV(IHEAV) .EQ. 7 ) THEN

                  SUMT = ZERZER
                  DO II = 1, NTRACK
                     SUMT = SUMT + TTRACK(II)
                  END DO

                  SUMD = ZERZER
                  DO II = 1, MTRACK
                     SUMD = SUMD + DTRACK(II)
                  END DO

                  IF ( SUMT .GT. ZERZER ) THEN
                     LETW = 100.0D0 * SUMD / SUMT
                     FLUSCW = LETW
                  END IF
               END IF
            END IF
         END IF

         RETURN
      END IF

C     ------------------------------------------------------------------
C     Li-7 LET^2 weighting branch for FLUSCW track-length scoring.
C
C     Scorer key:
C        SCONAM = 'L7L2'
C
C     Physics meaning:
C        Keep only transported lithium-7 fragments, identified as
C        charge Z = 3 and mass A = 7, and return LET^2. Together with
C        the L7L1 scorer, this allows reconstruction of a Li-7
C        LET-weighted quantity in post-processing.
C
C     Unit reconstruction is identical to L7L1:
C        LET [keV/um] = 100 * SUMD [GeV] / SUMT [cm].
C     ------------------------------------------------------------------

      IF ( SCONAM(1:4) .EQ. 'L7L2' ) THEN
         FLUSCW = ZERZER

         IF ( JTRACK .LT. -6 .AND. NPHEAV .GT. 0 ) THEN
            IHEAV = KHEAVY(NPHEAV)

            IF ( IHEAV .GE. 1 .AND. IHEAV .LE. KXHEAV ) THEN
               IF ( ICHEAV(IHEAV) .EQ. 3 .AND.
     &              IBHEAV(IHEAV) .EQ. 7 ) THEN

                  SUMT = ZERZER
                  DO II = 1, NTRACK
                     SUMT = SUMT + TTRACK(II)
                  END DO

                  SUMD = ZERZER
                  DO II = 1, MTRACK
                     SUMD = SUMD + DTRACK(II)
                  END DO

                  IF ( SUMT .GT. ZERZER ) THEN
                     LETW = 100.0D0 * SUMD / SUMT
                     FLUSCW = LETW * LETW
                  END IF
               END IF
            END IF
         END IF

         RETURN
      END IF
C     ------------------------------------------------------------------
C     Li-6 fluence/filter branch for FLUSCW.
C
C     Scorer key:
C        SCONAM starts with 'L6FL'
C
C     Physics meaning:
C        Keep only transported lithium-6 fragments, identified as
C        charge Z = 3 and mass A = 6. For matching Li-6 tracks this
C        branch returns ONEONE, so the underlying estimator is scored
C        without additional LET weighting.
C
C     Interpretation:
C        This is an isotope-selection filter. It answers "is the current
C        transported heavy fragment Li-6?" It does not determine where
C        the Li-6 fragment was produced or which parent particle produced
C        it. That ancestry information would require production-time
C        tagging with STUPRF or MDSTCK.
C     ------------------------------------------------------------------

      IF ( SCONAM(1:4) .EQ. 'L6FL' ) THEN
         FLUSCW = ZERZER

         IF ( JTRACK .LT. -6 .AND. NPHEAV .GT. 0 ) THEN
            IHEAV = KHEAVY(NPHEAV)

            IF ( IHEAV .GE. 1 .AND. IHEAV .LE. KXHEAV ) THEN
               IF ( ICHEAV(IHEAV) .EQ. 3 .AND.
     &              IBHEAV(IHEAV) .EQ. 6 ) THEN
                  FLUSCW = ONEONE
               END IF
            END IF
         END IF

         RETURN
      END IF
C     ------------------------------------------------------------------
C     Li-7 fluence/filter branch for FLUSCW.
C
C     Scorer key:
C        SCONAM starts with 'L7FL'
C
C     Physics meaning:
C        Keep only transported lithium-7 fragments, identified as
C        charge Z = 3 and mass A = 7. For matching Li-7 tracks this
C        branch returns ONEONE, so the underlying estimator is scored
C        without additional LET weighting.
C
C     Interpretation:
C        This is an isotope-selection filter during particle transport.
C        It separates Li-7 from other transported heavy fragments, but
C        it does not record the production vertex, parent particle, or
C        nuclear reaction channel. That would require production-time
C        ancestry tagging with STUPRF or MDSTCK.
C     ------------------------------------------------------------------

      IF ( SCONAM(1:4) .EQ. 'L7FL' ) THEN
         FLUSCW = ZERZER

         IF ( JTRACK .LT. -6 .AND. NPHEAV .GT. 0 ) THEN
            IHEAV = KHEAVY(NPHEAV)

            IF ( IHEAV .GE. 1 .AND. IHEAV .LE. KXHEAV ) THEN
               IF ( ICHEAV(IHEAV) .EQ. 3 .AND.
     &              IBHEAV(IHEAV) .EQ. 7 ) THEN
                  FLUSCW = ONEONE
               END IF
            END IF
         END IF

         RETURN
      END IF

C     ==================================================================
C     All-charged-particle LET moments -- the primary quantity.
C
C     Scorer keys:
C        ALL1   first  LET moment  [keV/um]
C        ALL2   second LET moment  [(keV/um)^2]
C        ALFL   unweighted fluence over the SAME particle set
C
C     LET here is the true LOCAL energy-deposition LET of the currently
C     transported charged particle, reconstructed from TRACKR step data:
C
C        LET [keV/um] = 100 * SUMD [GeV] / SUMT [cm]
C
C     (1 GeV/cm = 100 keV/um). This definition is universal: it is valid
C     for every charged hadron and ion FLUKA transports -- protons, light
C     ions, and heavy fragments -- independent of GETLET or an explicit
C     material-index lookup. The quantity still depends physically on the
C     local transported material, through the energy deposited along the
C     step (DTRACK) per unit track length (SUMT).
C
C     Particle selection:
C        Neutral particles carry no LET and are skipped.
C        Electrons and positrons (JTRACK = 3, 4) are deliberately
C        EXCLUDED even when EMF transport is active. Averaged-LET
C        reporting for particle therapy conventionally covers the
C        hadron/ion field only; delta electrons are the energy-transfer
C        mechanism rather than separate LET carriers. Including them
C        would also make track-averaged LET depend on the EMF transport
C        threshold instead of on the physics.
C
C     Post-processing:
C
C        dose-averaged LET  = ALL2 / ALL1
C        track-averaged LET = ALL1 / ALFL
C
C     ALFL is the correct denominator for the track average and returns
C     ONEONE for exactly the particles ALL1/ALL2 accept, so numerator and
C     denominator cover the same particle set by construction.
C
C     Do NOT use a plain ALL-PART track-length USRBIN as the denominator.
C     The ALL-PART generalized particle counts neutrons, photons and
C     electrons, none of which contribute to the ALL1 numerator, so that
C     ratio underestimates the track-averaged LET. (ALL2/ALL1 is immune,
C     since both moments share this particle selection.)
C     ==================================================================

      IF ( SCONAM(1:4) .EQ. 'ALL1' .OR. SCONAM(1:4) .EQ. 'ALL2' .OR.
     &     SCONAM(1:4) .EQ. 'ALFL' ) THEN
         FLUSCW = ZERZER

C        Skip neutral particles (regular particles with zero charge) and
C        electrons/positrons (JTRACK = 3, 4). Ions and nuclear fragments
C        are transported with JTRACK .LT. 0 and are always charged.
         IF ( JTRACK .GT. 0 ) THEN
            IF ( ICHRGE(JTRACK) .EQ. 0 ) RETURN
            IF ( JTRACK .EQ. 3 .OR. JTRACK .EQ. 4 ) RETURN
         END IF

C        Unweighted fluence over the accepted particle set: the matching
C        denominator for the track average. No LET weight is applied.
         IF ( SCONAM(1:4) .EQ. 'ALFL' ) THEN
            FLUSCW = ONEONE
            RETURN
         END IF

         SUMT = ZERZER
         DO II = 1, NTRACK
            SUMT = SUMT + TTRACK(II)
         END DO

         SUMD = ZERZER
         DO II = 1, MTRACK
            SUMD = SUMD + DTRACK(II)
         END DO

         IF ( SUMT .GT. ZERZER ) THEN
            LETW = 100.0D0 * SUMD / SUMT
            IF ( SCONAM(1:4) .EQ. 'ALL1' ) THEN
               FLUSCW = LETW
            ELSE
               FLUSCW = LETW * LETW
            END IF
         END IF

         RETURN
      END IF

C     ==================================================================
C     All-particle water-reference LET moments (GETLET evaluated in
C     water, MATLET = MWATER).
C
C     Scorer keys:
C        ALW1   first  LET moment  [keV/um], LET evaluated in water
C        ALW2   second LET moment  [(keV/um)^2], LET evaluated in water
C        ALWF   unweighted fluence over the SAME particle set
C
C     Post-processing, exactly as for the ALL family:
C        dose-averaged LET  = ALW2 / ALW1
C        track-averaged LET = ALW1 / ALWF
C
C     ALWF, not a plain ALL-PART bin, is the denominator for the track
C     average: the ALW particle set is narrower still (light particles
C     only), so an ALL-PART denominator would be wrong by a larger margin
C     here than for the ALL family.
C
C     Water-reference LET is available through GETLET only for the light
C     particles it supports: proton, deuteron, triton, He-3, He-4
C     (IJ = 1, -3, -4, -5, -6). Heavier fragments (Li and above) cannot
C     be evaluated in water with this GETLET build and are NOT included
C     in the water-reference variant. Use ALL1/ALL2 for the complete
C     all-particle (local) quantity.
C     ==================================================================

      IF ( SCONAM(1:4) .EQ. 'ALW1' .OR. SCONAM(1:4) .EQ. 'ALW2' .OR.
     &     SCONAM(1:4) .EQ. 'ALWF' ) THEN
         FLUSCW = ZERZER

         IF ( MWATER .LE. 0 ) RETURN

         IF ( IJ .NE. 1  .AND. IJ .NE. -3 .AND. IJ .NE. -4 .AND.
     &        IJ .NE. -5 .AND. IJ .NE. -6 ) THEN
            RETURN
         END IF

         EKIN = -PLA
         IF ( EKIN .LE. 1.0D-09 ) THEN
            RETURN
         END IF

C        Unweighted fluence over the accepted particle set.
         IF ( SCONAM(1:4) .EQ. 'ALWF' ) THEN
            FLUSCW = ONEONE
            RETURN
         END IF

         LETW = GETLET(IJ, EKIN, PLA, ZERZER, MWATER)
         LETLIN = RHO(MWATER) * LETW
         IF ( SCONAM(1:4) .EQ. 'ALW1' ) THEN
            FLUSCW = LETLIN
         ELSE
            FLUSCW = LETLIN * LETLIN
         END IF

         RETURN
      END IF
C     ------------------------------------------------------------------
C     Proton FLUSCW branch for LET and primary-proton fluence scoring.
C
C     Trigger condition:
C        IJ .EQ. 1          current scored particle is a proton.
C        ISCRNG .EQ. 2      FLUSCW is being called for fluence-like
C                           estimators, including track-length USRBIN.
C
C     Generation convention:
C        LTRACK .EQ. 1      source-generation proton, i.e. one of the
C                           original protons sampled by SOURCE.
C        LTRACK .GT. 1      non-primary proton, i.e. a secondary or later
C                           proton created by a discrete interaction.
C
C     Therefore:
C        PAL1, PAL2, PAW1, PAW2
C           score all transported protons, including source protons and
C           secondary/later-generation protons.
C
C        P1FL, P1L1, P1L2, P1W1, P1W2
C           score only source-generation protons because they require
C           LTRACK .EQ. 1.
C
C     Limitation:
C        LTRACK gives generation number, but not the production vertex,
C        parent particle, target nucleus, or reaction channel. Those
C        ancestry details would require production-time tagging with
C        STUPRF or MDSTCK.
C     ------------------------------------------------------------------

      IF ( ISCRNG .EQ. 2 .AND.
     &     ( SCONAM(1:4) .EQ. 'PAL1' .OR.
     &       SCONAM(1:4) .EQ. 'PAL2' .OR.
     &       SCONAM(1:4) .EQ. 'PAW1' .OR.
     &       SCONAM(1:4) .EQ. 'PAW2' .OR.
     &       SCONAM(1:4) .EQ. 'P1FL' .OR.
     &       SCONAM(1:4) .EQ. 'P1L1' .OR.
     &       SCONAM(1:4) .EQ. 'P1L2' .OR.
     &       SCONAM(1:4) .EQ. 'P1W1' .OR.
     &       SCONAM(1:4) .EQ. 'P1W2' ) ) THEN
         IF ( IJ .NE. 1 ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF
      END IF

      IF ( IJ .EQ. 1 .AND. ISCRNG .EQ. 2 ) THEN
         EKIN = -PLA
         IF ( EKIN .LE. 1.0D-09 ) THEN
            FLUSCW = ZERZER
            RETURN
         END IF
C        Proton scorer-key map inside this branch:
C
C        PAL1:
C           All-proton LET in the local transport material.
C           Returns LET [keV/um].
C
C        PAL2:
C           All-proton LET^2 in the local transport material.
C           Returns LET^2 [(keV/um)^2].
C
C        PAW1:
C           All-proton LET evaluated in water, independent of local
C           material. Uses MATLET = MWATER.
C
C        PAW2:
C           All-proton LET^2 evaluated in water. Uses MATLET = MWATER.
C
C        P1FL:
C           Primary/source-generation proton fluence filter.
C           Returns ONEONE only when LTRACK .EQ. 1.
C
C        P1L1, P1L2:
C           Primary/source-generation proton LET and LET^2 in the local
C           transport material. Require LTRACK .EQ. 1.
C
C        P1W1, P1W2:
C           Primary/source-generation proton LET and LET^2 evaluated in
C           water. Require LTRACK .EQ. 1 and use MATLET = MWATER.

         IF (SCONAM(1:4) .EQ. 'PAL1') THEN
            MATLET = MEDFLK(NREG,1)
            IF ( MATLET .LE. 0 .OR. RHO(MATLET) .LE. ZERZER ) THEN
               FLUSCW = ZERZER
               RETURN
            END IF
            LETW = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)
            LETLIN = RHO(MATLET) * LETW
            FLUSCW = LETLIN

         ELSE IF (SCONAM(1:4) .EQ. 'PAL2') THEN
            MATLET = MEDFLK(NREG,1)
            IF ( MATLET .LE. 0 .OR. RHO(MATLET) .LE. ZERZER ) THEN
               FLUSCW = ZERZER
               RETURN
            END IF
            LETW = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)
            LETLIN = RHO(MATLET) * LETW
            FLUSCW = LETLIN * LETLIN

         ELSE IF (SCONAM(1:4) .EQ. 'PAW1') THEN
            IF ( MWATER .LE. 0 ) THEN
               FLUSCW = ZERZER
               RETURN
            END IF
            MATLET = MWATER
            LETW = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)
            LETLIN = RHO(MATLET) * LETW
            FLUSCW = LETLIN

         ELSE IF (SCONAM(1:4) .EQ. 'PAW2') THEN
            IF ( MWATER .LE. 0 ) THEN
               FLUSCW = ZERZER
               RETURN
            END IF
            MATLET = MWATER
            LETW = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)
            LETLIN = RHO(MATLET) * LETW
            FLUSCW = LETLIN * LETLIN

         ELSE IF (SCONAM(1:4) .EQ. 'P1FL') THEN
            IF ( LTRACK .EQ. 1 ) THEN
               FLUSCW = ONEONE
            ELSE
               FLUSCW = ZERZER
            END IF

         ELSE IF (SCONAM(1:4) .EQ. 'P1L1') THEN
            IF ( LTRACK .EQ. 1 ) THEN
               MATLET = MEDFLK(NREG,1)
               IF ( MATLET .LE. 0 .OR. RHO(MATLET) .LE. ZERZER ) THEN
                  FLUSCW = ZERZER
                  RETURN
               END IF
               LETW = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)
               LETLIN = RHO(MATLET) * LETW
               FLUSCW = LETLIN
            ELSE
               FLUSCW = ZERZER
            END IF

         ELSE IF (SCONAM(1:4) .EQ. 'P1L2') THEN
            IF ( LTRACK .EQ. 1 ) THEN
               MATLET = MEDFLK(NREG,1)
               IF ( MATLET .LE. 0 .OR. RHO(MATLET) .LE. ZERZER ) THEN
                  FLUSCW = ZERZER
                  RETURN
               END IF
               LETW = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)
               LETLIN = RHO(MATLET) * LETW
               FLUSCW = LETLIN * LETLIN
            ELSE
               FLUSCW = ZERZER
            END IF

         ELSE IF (SCONAM(1:4) .EQ. 'P1W1') THEN
            IF ( LTRACK .EQ. 1 ) THEN
               IF ( MWATER .LE. 0 ) THEN
                  FLUSCW = ZERZER
                  RETURN
               END IF
               MATLET = MWATER
               LETW = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)
               LETLIN = RHO(MATLET) * LETW
               FLUSCW = LETLIN
            ELSE
               FLUSCW = ZERZER
            END IF

         ELSE IF (SCONAM(1:4) .EQ. 'P1W2') THEN
            IF ( LTRACK .EQ. 1 ) THEN
               IF ( MWATER .LE. 0 ) THEN
                  FLUSCW = ZERZER
                  RETURN
               END IF
               MATLET = MWATER
               LETW = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)
               LETLIN = RHO(MATLET) * LETW
               FLUSCW = LETLIN * LETLIN
            ELSE
               FLUSCW = ZERZER
            END IF

         END IF

      END IF

      RETURN
*=== End of function Fluscw ===========================================*
      END
C=======================================================================
C COMSCW dose-scoring isotope filter.
C=======================================================================
C
C     COMSCW is the user weighting function used by dose-like estimators.
C     This is separate from FLUSCW, which handles fluence-like and
C     track-length estimators.
C
C     Therefore Li-6 and Li-7 isotope filtering for DOSE-like USRBIN
C     scorers must be done here. If this filtering were only present in
C     FLUSCW, the Li-specific dose scorers would not be restricted to the
C     intended lithium isotope.
C
C     Current scope:
C
C        L6DO:
C           keep only transported Li-6 fragments with Z = 3 and A = 6.
C
C        L7DO:
C           keep only transported Li-7 fragments with Z = 3 and A = 7.
C
C     This is transport-time isotope filtering. It does not identify the
C     production vertex, parent particle, target nucleus, or reaction
C     channel. Those ancestry details would require STUPRF or MDSTCK.
C=======================================================================
      DOUBLE PRECISION FUNCTION COMSCW ( IJ    , XA    , YA    , ZA    ,
     &                                   MREG  , RULL  , LLO   , ICALL )

      INCLUDE 'dblprc.inc'
      INCLUDE 'dimpar.inc'
      INCLUDE 'iounit.inc'
      INCLUDE 'scohlp.inc'
      INCLUDE 'usrbin.inc'
      INCLUDE 'trackr.inc'
      INCLUDE 'fheavy.inc'

      INTEGER IHEAV
      CHARACTER*8 SCONAM

      LSCZER = .FALSE.
      COMSCW = ONEONE
      SCONAM = TRIM(ADJUSTL(TITUSB(JSCRNG)))

C     ------------------------------------------------------------------
C     Primary-proton DOSE filter for COMSCW.
C
C     Scorer key:
C        SCONAM = 'P1DO'
C
C     Physics meaning:
C        For primary-proton dose-like USRBIN scorers, reject every
C        energy-deposition contribution except contributions from
C        source-generation protons.
C
C        IJ .EQ. 1      selects protons.
C        LTRACK .EQ. 1  selects source-generation / primary protons.
C
C     COMSCW return value:
C        COMSCW = ONEONE   keep this dose contribution.
C        COMSCW = ZERZER   reject this contribution for this scorer.
C
C     This is different from PDOSE_ZN, which uses AUXSCORE PROTON and
C     therefore includes both primary and secondary protons.
C     ------------------------------------------------------------------

      IF ( ISCRNG .EQ. 1 .AND. SCONAM(1:4) .EQ. 'P1DO' ) THEN
         IF ( IJ .EQ. 1 .AND. LTRACK .EQ. 1 ) THEN
            COMSCW = ONEONE
         ELSE
            COMSCW = ZERZER
         END IF

         RETURN
      END IF
C     ------------------------------------------------------------------
C     Li-6 DOSE filter for COMSCW.
C
C     Scorer key:
C        SCONAM starts with 'L6DO'
C
C     Physics meaning:
C        For Li-6 dose-like USRBIN scorers, reject every transported
C        particle except lithium-6 fragments. A matching Li-6 fragment is
C        identified by charge Z = 3 and mass A = 6.
C
C     COMSCW return value:
C        COMSCW = ONEONE   keep this energy-deposition contribution.
C        COMSCW = ZERZER   reject this contribution for this scorer.
C
C     Note:
C        This filters dose contributions during transport. It does not
C        identify the production site or parent particle of the Li-6.
C     ------------------------------------------------------------------

      IF ( ISCRNG .EQ. 1 .AND.
     &     SCONAM(1:4) .EQ. 'L6DO' ) THEN
         COMSCW = ZERZER

         IF ( JTRACK .LT. -6 .AND. NPHEAV .GT. 0 ) THEN
            IHEAV = KHEAVY(NPHEAV)

            IF ( IHEAV .GE. 1 .AND. IHEAV .LE. KXHEAV ) THEN
               IF ( ICHEAV(IHEAV) .EQ. 3 .AND.
     &              IBHEAV(IHEAV) .EQ. 6 ) THEN
                  COMSCW = ONEONE
               END IF
            END IF
         END IF

         RETURN
      END IF

C     ------------------------------------------------------------------
C     Li-7 DOSE filter for COMSCW.
C
C     Scorer key:
C        SCONAM starts with 'L7DO'
C
C     Physics meaning:
C        For Li-7 dose-like USRBIN scorers, reject every transported
C        particle except lithium-7 fragments. A matching Li-7 fragment is
C        identified by charge Z = 3 and mass A = 7.
C
C     COMSCW return value:
C        COMSCW = ONEONE   keep this energy-deposition contribution.
C        COMSCW = ZERZER   reject this contribution for this scorer.
C
C     Note:
C        This filters dose contributions during transport. It does not
C        identify the production site or parent particle of the Li-7.
C     ------------------------------------------------------------------

      IF ( ISCRNG .EQ. 1 .AND.
     &     SCONAM(1:4) .EQ. 'L7DO' ) THEN
         COMSCW = ZERZER

         IF ( JTRACK .LT. -6 .AND. NPHEAV .GT. 0 ) THEN
            IHEAV = KHEAVY(NPHEAV)

            IF ( IHEAV .GE. 1 .AND. IHEAV .LE. KXHEAV ) THEN
               IF ( ICHEAV(IHEAV) .EQ. 3 .AND.
     &              IBHEAV(IHEAV) .EQ. 7 ) THEN
                  COMSCW = ONEONE
               END IF
            END IF
         END IF

         RETURN
      END IF

      RETURN
*=== End of function Comscw ===========================================*
      END
