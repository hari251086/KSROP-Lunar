c     ============================================================
      subroutine moon_rotation_angle_deg(djd, theta)
c     Moon prime-meridian rotation angle W(d), degrees, normalized
c     to [0,360) -- the lunar analogue of KSROP's own gmst_deg and
c     KSROP-Mars's mars_rotation_angle_deg, for use as the body-fixed
c     reference angle in tesseral-harmonic longitude (KSROP#30
c     rollout, KSROP-Lunar issue #26).
c
c     W(d) = 38.3213 + 13.17635815*d - 1.4e-12*d^2, d = days (TDB)
c     from J2000.0 (JD 2451545.0) -- the SECULAR part only of the
c     IAU 2009 Moon rotation model (Archinal et al. 2011). Source:
c     NASA NAIF generic PCK kernel pck00011.tpc, BODY301_PM =
c     (38.3213, 13.17635815, -1.4e-12), fetched directly from
c     https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/
c     pck00011.tpc (2026-08-08). The rate (13.17635815 deg/day)
c     matches the Moon's known synchronous rotation: 360 deg /
c     27.321661 days (sidereal month) = 13.17640 deg/day, agreeing
c     to 5 significant figures -- cross-confirms the constant.
c
c     IMPORTANT LIMITATION (worth stating plainly, more so than the
c     analogous Earth/Mars gap): the full IAU lunar rotation model
c     also includes 13 periodic "physical libration" correction
c     terms (BODY301_NUT_PREC_PM in the same kernel), NOT included
c     here. Real lunar physical libration has amplitude of order a
c     few degrees in longitude -- much larger, relative to this
c     model's own precision, than the nutation/UT1 corrections
c     already-documented as missing from gmst_deg (Earth) or
c     mars_rotation_angle_deg (Mars). This mean-only formula is
c     adequate for the same tier of use those two are (multi-day to
c     multi-month secular/resonance studies), but genuinely less
c     precise, in relative terms, for anything needing sub-degree
c     selenographic longitude accuracy. Not implemented now; a real
c     future upgrade if that precision is ever needed, not a silent
c     gap -- see README/ALGORITHM.md.
c     ============================================================
      implicit double precision (a-h,o-z)

      dtd = djd - 2451545.0d0

      theta = 38.3213d0 + 13.17635815d0*dtd - 1.4d-12*dtd*dtd

      theta = dmod(theta, 360.0d0)
      if (theta .lt. 0.0d0) theta = theta + 360.0d0

      return
      end
