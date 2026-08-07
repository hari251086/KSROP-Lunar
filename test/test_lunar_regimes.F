c     test_lunar_regimes.F  -  Unit tests for KSROP-Lunar
c     -------------------------------------------------------
c     Compile:  fpm test --compiler ifx   (or gfortran on Linux/CI --
c               see README Known Issues for the Windows/gfortran
c               caveat this repo works around)
c     Run:      exit 0 = all pass, exit 1 = any fail
c     -------------------------------------------------------
c     Tests covered:
c       geo_coeff_body   lunar J2 coefficient file loads correctly
c       car2ks/ks2car    KS transform round-trip at lunar orbital
c                        scale (r~1800 km, v~1.6 km/s -- very
c                        different magnitude than KSROP's usual
c                        Earth GTO/LEO scale, worth its own check)
c     -------------------------------------------------------
      program test_lunar_regimes

      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      dimension c_j(4)
      dimension xr(3), xdr(3), ur(4), usr(4), xout(3), xdout(3)
      integer npass, nfail

      amue = 4902.8001d0
      R_Earth = 1737.4d0
      pi = 4.d0*datan(1.d0)

      npass = 0
      nfail = 0
      eps8 = 1.d-8

      write(*,*) '=============================================='
      write(*,*) ' KSROP-Lunar Unit Test Suite'
      write(*,*) '=============================================='

c     ==============================================================
c     GC1: lunar J2 coefficient file loads via geo_coeff_body
c     input/GRAIL_lowdeg_zonal.dat row: 2 0 -0.906501958078415D-04 ...
c     Expected c_j(2) = -sqrt(5)*(-0.906501958078415D-04)
c                      = J2_moon = 202.7d-6 (JPL/GRAIL, see README)
c     ==============================================================
      call geo_coeff_body(2, c_j, 'input/GRAIL_lowdeg_zonal.dat')
      call assert('GC1 lunar J2 from coefficient file',
     &            c_j(2), 202.7d-6, 1.d-9, npass, nfail)

c     ==============================================================
c     KS1: car2ks -> ks2car round-trip at lunar orbital scale.
c     A representative low circular lunar orbit: r=1837.4 km
c     (100 km altitude), v=1.6335 km/s (circular).
c     ==============================================================
      xr(1) = 1837.4d0
      xr(2) = 0.d0
      xr(3) = 0.d0
      xdr(1) = 0.d0
      xdr(2) = dsqrt(amue/1837.4d0)
      xdr(3) = 0.d0

      w = dsqrt(0.5d0*(amue/1837.4d0 - dotp3(xdr,xdr)/2.d0))
      call car2ks(xr, xdr, ur, usr, w)
      call ks2car(ur, usr, xout, xdout, w)

      call assert('KS1 round-trip X', xout(1), xr(1), eps8,
     &            npass, nfail)
      call assert('KS1 round-trip Y', xout(2), xr(2), eps8,
     &            npass, nfail)
      call assert('KS1 round-trip Z', xout(3), xr(3), eps8,
     &            npass, nfail)
      call assert('KS1 round-trip XDOT', xdout(1), xdr(1), eps8,
     &            npass, nfail)
      call assert('KS1 round-trip YDOT', xdout(2), xdr(2), eps8,
     &            npass, nfail)

c     ==============================================================
c     Summary
c     ==============================================================
      write(*,*) '=============================================='
      write(*,200) npass + nfail, npass, nfail
 200  format(' Total: ',i3,'   Passed: ',i3,'   Failed: ',i3)
      write(*,*) '=============================================='

      if (nfail .gt. 0) stop 1
      stop
      end

c     ============================================================
      subroutine assert(label, got, expected, tol, npass, nfail)
c     Absolute tolerance check: |got - expected| <= tol
c     ============================================================
      implicit double precision (a-h, o-z)
      character*(*) label
      integer npass, nfail

      diff = dabs(got - expected)
      if (diff .le. tol) then
         write(*,10) 'PASS', label
         npass = npass + 1
      else
         write(*,20) 'FAIL', label, got, expected, diff
         nfail = nfail + 1
      end if
 10   format('  [',a4,']  ',a)
 20   format('  [',a4,']  ',a,/,
     &       '           got=',e20.13,'  exp=',e20.13,
     &       '  diff=',e10.3)
      return
      end
