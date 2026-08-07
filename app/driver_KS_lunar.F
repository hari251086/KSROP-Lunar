c     COPYRIGHT 2026, HARISHKUMAR SELLAMUTHU, ALL RIGHTS RESERVED
C     -------------------------------------------------------------
c     KS - Numerical driver (Moon-centered)
c       Description:
c     Orbit propagation using Kustaanheimo-Stiefel regular elements
c     with numerical integrator RKG4, for a Moon-centered orbiter.
c     Adapted from KSROP's app/driver_KS.F (Earth-centered) --
c     see KSROP issue #26. The KS engine, force-model math
c     (third_body_aux/qsun/qmoon convention, oblateness sum, RKG4
c     integrator) is entirely reused unmodified from the KSROP
c     library dependency; only this driver's central body and the
c     physical meaning of the two third-body "slots" change.
c
c     Perturbations included:
c     - Lunar oblateness (zonal, degree configurable; ships with a
c       J2-only coefficient file -- see input/GRAIL_lowdeg_zonal.dat
c       and README Known Limitations)
c     - Earth + Sun third-body gravity (see slot repurposing below)
c     - Solar radiation pressure (cannon-ball, shadow-by-Moon)
c     - NO atmospheric drag: the Moon is airless. This driver never
c       reads an atmosphere table and never applies a drag force.
c
c     -------------------------------------------------------------
c     THIRD-BODY SLOT REPURPOSING (read this before touching the
c     third-body code below):
c
c     KSROP's Earth-centered driver has two third-body "slots",
c     historically named for the Sun and the Moon (variables ts/tm,
c     amuS/amuM, nsun_deg/nmoon_deg, qsun/qmoon, ...). The underlying
c     force math (third_body_aux, the qsun/qmoon KS-EOM convention)
c     only ever consumes a generic position vector + GM -- nothing
c     in it is Sun- or Moon-specific. For a Moon-centered orbiter,
c     the two dynamically relevant third bodies are EARTH and SUN,
c     so this driver keeps the original variable names (to reuse the
c     force-model code verbatim, unmodified from KSROP) but assigns:
c
c       ts / amuS / nsun_deg / qsun / RbyRs / Rs / cn_s / Tau_3body1
c         -> physically the EARTH slot
c       tm / amuM / nmoon_deg / qmoon / RbyRm / Rm / cn_m / Tau_3body2
c         -> physically the SUN slot
c
c     Positions are obtained by vector arithmetic on KSROP's existing
c     geocentric ephemerides -- no new ephemeris series needed:
c       Earth position relative to Moon  = -lunarpv(dj)
c       Sun position relative to Moon    = solarnpv(dj) - lunarpv(dj)
c     (lunarpv/solarnpv return Moon/Sun position AS SEEN FROM EARTH;
c     negating lunarpv gives Earth as seen from the Moon exactly,
c     since it's the same vector reversed. The Sun-relative-to-Moon
c     combination is the usual A-relative-to-C = A-rel-B - C-rel-B
c     construction.) SRP (which needs the true Sun direction) uses
c     `tm`, not `ts` -- the opposite of the Earth-centered driver.
c     -------------------------------------------------------------
c     Author: Harishkumar Sellamuthu
c     hari251086@gmail.com
c     -------------------------------------------------------------
c     Input file parameters <input.dat> by line (lunar-specific,
c     shorter than KSROP's Earth input.dat -- no drag line):
c       nrev  istep  tole
c       n_force(1)=lunar geo, n_force(2)=Earth 3rd-body,
c         n_force(3)=Sun 3rd-body
c       CR  AM(m2/kg)  IPSR(0/1)  ISHAD(0/1/2)
c     -------------------------------------------------------------
c     Subroutines used (all from the KSROP library dependency):
c     geo_coeff_body, car2ks, ks2car, aLegP, car2oe, ks2ksr
c     -------------------------------------------------------------

      implicit double precision (a-h, o-z)
      integer cn0
      double precision iE
      parameter(n = 10, cn0 = 2190)
c     Value of pi, degree to radians and radians to degree.
c     R_Earth (from the shared /xy/ common block, KSROP naming) holds
c     the MOON's radius here -- it is populated from whatever radius
c     const_moon.dat's first line supplies, and is consumed as-is by
c     the shadow-model subroutines (shadfncyl/shadfncone), which are
c     already body-agnostic: they just need "radius of the body the
c     satellite orbits," which for a lunar orbiter's eclipse geometry
c     is correctly the Moon's own radius.
      common /xy/ pi,d2r,r2d,amue,AU,R_Earth
      dimension pek(7),x(3),xd(3),q(4),u(4),us(4),y(10),z(10)
      dimension qu(3),q_til(4),u_til(4),uu(4)
      dimension ZbyR(cn0*2),ObyR(cn0*2),Re(cn0)
      dimension cn_j(cn0),c_j(cn0),V_p(cn0),aleg(cn0)
      dimension sleg(cn0),oleg(cn0),R(n*2)
      dimension cal(6),ts(3),tm(3),opt(3),n_force(3)
      dimension RbyRs(cn0),RbyRm(cn0),cn_m(cn0),cn_s(cn0),Rm(n*2)
      dimension qsun(4),qmoon(4),slambda(cn0),olambda(cn0),Rs(n*2)
c     Ephemeris scratch (Earth/Sun-relative-to-Earth, before the
c     Moon-centered vector combination described above).
c     DELIBERATELY two separate DIMENSION statements, not one
c     combined "dimension moon_wrt_earth(3), sun_wrt_earth(3)":
c     that form, followed by a scalar declaration and a call
c     passing one of the two arrays to an externally-compiled
c     subroutine (lunarpv/solarnpv, compiled in the KSROP
c     dependency, a separate translation unit), reproducibly
c     corrupts the scalar via a gfortran stack-layout/ABI
c     interaction confirmed independent of -fcheck and of
c     optimization level (isolated in a 15-line repro outside
c     this file). Splitting the declaration is a complete,
c     verified fix -- not a superstition.
      dimension moon_wrt_earth(3)
      dimension sun_wrt_earth(3)
c     No drag: this driver has no atmosphere table, no ALT_atm/
c     DEN_atm/SCH_atm arrays, and no qdrag/p_drag -- q(j)drag terms
c     below are permanently zero.
c     SRP: cannonball SRP perturbation (reused from KSROP unchanged)
      dimension p_srp(3), qsrp(4), r_srp(3)
      integer IPSR, ISHAD
c     OEM/OPM epoch strings (YYYY-MM-DDTHH:MM:SS.sss = 23 chars)
      character*23 epochstr, epoch_stop, oem_start
c     Output filenames: output/KSROP_LUNAR_YYYYMMDDTHHMMSS.*
      character*46 oem_fname, reg_fname
      character*15 compact_epoch
c     Trajectory buffer: allocated once nrev is known; written after run
      double precision, allocatable :: traj_jd(:)
      double precision, allocatable :: traj_x(:,:)
      double precision, allocatable :: traj_xd(:,:)

      pi=4.d0*datan(1.d0)
      d2r=pi/180.d0
      r2d=180.d0/pi

c     Get UTC timestamp for output filenames
      call utc_now_epoch(epochstr, compact_epoch)
      reg_fname = 'output/KSROP_LUNAR_'//compact_epoch//'_Regular.out'

c     Input Files
      open (unit = 1, file = 'input/const_moon.dat')  !Constants
      open (unit = 2, file = 'input/input.dat')       !Sim params
      open (unit = 8, file = 'input/input.opm')       !OPM input

c     Output files (unit 3 opened later with timestamped filename)
      open (unit = 4, file = reg_fname)                !KS elements (debug)
      open (unit = 5, file = 'output/ksrop_lunar.opm')  !OPM output

c     const_moon.dat line 1: mu_Moon R_Moon AU mu_Earth mu_Sun
c     (see THIRD-BODY SLOT REPURPOSING above for why amuS=Earth,
c     amuM=Sun here rather than the Earth-driver's Sun/Moon)
      read(1,*) amue,R_Earth,AU,amuS,amuM
      read(1,*) ngeo_deg,nsun_deg,nmoon_deg
      read(1,*) PSR_srp

c     Read initial state (position, velocity, epoch) from OPM file --
c     Moon-centered inertial frame, EME2000 axes (parallel-axis
c     translation from the Earth-centered frame, not rotated to the
c     Moon's body-fixed frame -- see README).
      call read_opm(8, x, xd, cal)
      close(8)

c     Read simulation parameters from input.dat (lunar-specific
c     format: no drag line)
      read(2,*,end=99) nrev,istep,tole !Revolutions, steps, tolerance
      read(2,*) n_force(1), n_force(2), n_force(3)
c     SRP parameters: reflectivity coeff, area-to-mass (m^2/kg),
c     SRP on/off, shadow model (0=sunlit,1=cyl,2=cone)
      read(2,*) CR_srp, AM_srp, IPSR, ISHAD

      allocate(traj_jd(nrev*istep+1),
     &         traj_x(3,nrev*istep+1),
     &         traj_xd(3,nrev*istep+1))
      idump = 0

      call car2oe(x,xd,pek)
      call oe2car(pek,x,xd,tole)
      call cal2jd(cal,djulian)
      call force_models(n_force,ngeo_deg,amuS,amuM)
      if (IPSR .eq. 0) PSR_srp = 0.d0
      SRPC_srp = CR_srp * AM_srp

      print *, x, xd

      write(*,*) '=============================================='
      write(*,*) 'KSROP-Lunar: Moon-centered propagation'

      if(ngeo_deg.eq.0.or.ngeo_deg.eq.1) then
         do i = 1,cn0
            c_j(i) = 0.d0
         end do
      else
!     compute j_n values from the lunar gravity coefficient file
!     (geo_coeff_body: KSROP v2.3.0+, see KSROP issue #26)
         call geo_coeff_body(ngeo_deg,c_j,
     &        'input/GRAIL_lowdeg_zonal.dat')
      end if

!     Initialization
      E_0 = 0.d0                !Gen. Eccentric Anomaly in radians
      E = E_0
      T = 0.d0                  !Physical time in seconds
      nt = 0                    !Counter in the integrator

      dE_0 = 2.d0*pi/istep      !Step-size in radians

      num_steps = nrev*istep

c     Earth/Sun position relative to the Moon (see THIRD-BODY SLOT
c     REPURPOSING above): ts=Earth, tm=Sun.
      call solarnpv(djulian, sun_wrt_earth)
      call lunarpv(djulian, moon_wrt_earth)
      do i = 1,3
         ts(i) = -moon_wrt_earth(i)
         tm(i) = sun_wrt_earth(i) - moon_wrt_earth(i)
      end do

      R(1) = vmn(x)
      V = vmn(xd)

      R(2) = R(1)*R(1)
      V2 = V*V

      do i = 2, 20
         R(i) = R(i-1) * R(1)
      end do

!     Perturbing potential, V_pot

!     Oblateness terms (Moon)
      Re(1) = R_Earth
      ZbyR(1) = x(3)/R(1)
      ObyR(1) = 1.d0/R(1)

!     Third-body terms: Earth (slot "s") + Sun (slot "m")
      call third_body_aux(nsun_deg,  x, ts, R(1), amuS,
     &                    Rs(1), RbyRs, first_s, cphi_s, cphi_s2)
      call third_body_aux(nmoon_deg, x, tm, R(1), amuM,
     &                    Rm(1), RbyRm, first_m, cphi_m, cphi_m2)

      do i = 1,ngeo_deg
         ObyR(i+1) = ObyR(i)*ObyR(1)
         Re(i+1) = Re(i)*Re(1)
      end do

      V_pot = 0.d0

      call aLegP(ngeo_deg+1,ZbyR(1),aleg)
      call aLegP(nsun_deg+1,cphi_s,sleg)
      call aLegP(nmoon_deg+1,cphi_m,oleg)

      do i = 2,ngeo_deg
         cn_j(i) = amue*Re(i)*c_j(i)*ObyR(i+1)
         V_p(i) = cn_j(i)*aleg(i+1)
         V_pot = V_pot+V_p(i)
      end do

      do i = 2,nsun_deg
         cn_s(i) = first_s * RbyRs(i)
         V_p(i) = cn_s(i) * sleg(i+1)
         V_pot = V_pot + V_p(i)
      end do

      do i = 2,nmoon_deg
         cn_m(i) = first_m * RbyRm(i)
         V_p(i) = cn_m(i) * oleg(i+1)
         V_pot = V_pot + V_p(i)
      end do

!     Angular frequency
      w = dsqrt(0.5d0*(amue/R(1)-V2/2.d0-V_pot))
      w_kep = dsqrt(0.5d0*(amue/R(1)-V2/2.d0))
      Gam = w/w_kep

      dE = dE_0*Gam

      call car2ks(x,xd,u,us,w)
      call ks2car(u,us,x,xd,w)

      Tow = T+dotp4(u,us)/w     !Time-element

!     Initial KS regular elements
      y(1) = Tow
      y(10) = w

      do i = 1,4
         y(i+1) = u(i)
         y(i+5) = 2.d0*us(i)
      end do

      sie = dsin(E/2.d0)
      cse = dcos(E/2.d0)
      iE = dint(E)

      call car2oe(x,xd,pek)
      call jd2epoch(djulian, epochstr)

c     --- CCSDS OPM v2.0 (initial epoch) — written immediately ---
      call write_opm(5, epochstr, x, xd, pek)

c     --- Buffer initial state point (OEM written after run) ---
      idump = idump + 1
      traj_jd(idump) = djulian + T/86400.d0
      do i = 1,3
         traj_x(i,idump)  = x(i)
         traj_xd(i,idump) = xd(i)
      end do
      write(4,101) iE,y

c     Integration
      do ik = 1,nrev            !Revolutions

         do ki = 1,istep        !Steps
 60         continue
            call u2uu(u,uu)

c     -----------------------------------------------------------
c     SRP: cannonball solar radiation pressure perturbation.
c     Uses `tm` (Sun-relative-to-Moon, see slot repurposing note),
c     NOT `ts` -- the opposite of KSROP's Earth-centered driver,
c     where the Sun IS the ts slot.
c     Shadow factor ROHP: 1=sunlit, 0=eclipse (by the Moon itself).
c     -----------------------------------------------------------
            if (PSR_srp .eq. 0.d0) then
               do i = 1,4
                  qsrp(i) = 0.d0
               end do
            else
               if (ISHAD .eq. 1) then
                  ROHP = shadfncyl(x(1),x(2),x(3),tm(1),tm(2),tm(3))
               else if (ISHAD .eq. 2) then
                  ROHP = shadfncone(x(1),x(2),x(3),tm(1),tm(2),tm(3))
               else
                  ROHP = 1.d0
               end if
               r_srp(1) = tm(1) - x(1)
               r_srp(2) = tm(2) - x(2)
               r_srp(3) = tm(3) - x(3)
               xsrr = vmn(r_srp)
               do i = 1,3
                  p_srp(i) = -ROHP*PSR_srp*SRPC_srp*r_srp(i)/xsrr*1.d-3
               end do
               qsrp(1) =  u(1)*p_srp(1)+u(2)*p_srp(2)+u(3)*p_srp(3)
               qsrp(2) = -u(2)*p_srp(1)+u(1)*p_srp(2)+u(4)*p_srp(3)
               qsrp(3) = -u(3)*p_srp(1)-u(4)*p_srp(2)+u(1)*p_srp(3)
               qsrp(4) =  u(4)*p_srp(1)-u(3)*p_srp(2)+u(2)*p_srp(3)
            end if

c     Time-element rate (KS regular elements eq., Sellamuthu 2018 thesis
c     eq. 2.56), reused unmodified from KSROP's Earth driver.
            Tau_geo = 0.d0
            do i = 2,ngeo_deg
               Tau_geo = Tau_geo+amue*(i-1)*c_j(i)*Re(i)*aleg(i+1)
     &                   *ObyR(i)
            end do

            Tau_3body1 = 0.d0
            do i = 2,nsun_deg
               Tau_3body1=Tau_3body1
     &              +(i+2)*sleg(i+1)*RbyRs(i)*RbyRs(1)*amuS
            end do
            Tau_3body2 = 0.d0
            do i = 2,nmoon_deg
               Tau_3body2=Tau_3body2
     &              +(i+2)*oleg(i+1)*RbyRm(i)*RbyRm(1)*amuM
            end do

c     -----------------------------------------------------------
c     Energy element: non-conservative SRP contribution only
c     (no drag term -- airless Moon)
c     dw/ds = -0.5 * (us . qsrp) / w
c     -----------------------------------------------------------
            z(10) = -0.5d0*(us(1)*qsrp(1)
     &                    + us(2)*qsrp(2)
     &                    + us(3)*qsrp(3)
     &                    + us(4)*qsrp(4))/w

c     -----------------------------------------------------------
c     Time element: conservative + SRP contributions
c     -----------------------------------------------------------
            z(1) = (amue + Tau_geo + Tau_3body1 + Tau_3body2
     &             + R(1)*dotp4(u,qsrp))/(8.d0*w*w*w)
     &             - z(10)*(2.d0/(w*w))*dotp4(u,us)

            if (ngeo_deg  .gt. 1) call aLegP(ngeo_deg+1, ZbyR(1), aleg)
            if (nsun_deg  .gt. 1) call aLegP(nsun_deg+1, cphi_s,  sleg)
            if (nmoon_deg .gt. 1) call aLegP(nmoon_deg+1,cphi_m,  oleg)

            do j = 1,4
               q(j) = 0.d0
               qsun(j) = 0.d0
               qmoon(j) = 0.d0
            end do

            do j = 1,4

!     oblateness terms (lunar)
               z2mR2=x(3)**2-R(1)**2
               qj=0.d0
               do i = 2,ngeo_deg
                  qj = dble(ObyR(i+1)*Re(i)*c_j(i)*
     &                 (R(1)*(i+1)*aleg(i+2)*(R(1)*uu(j)-x(3)*u(j))+
     &                 (aleg(i+1)*(u(j)*(i*R(2)+x(3)**2)
     &                 -uu(j)*(i+1)*R(1)*x(3))))/z2mR2)
                  q(j)=qj+q(j)
               end do

!     third-body terms: Earth (slot "s") + Sun (slot "m")
               call u2qu(u,qu,j)

               den_s = cphi_s2 - 1.d0
               den_m = cphi_m2 - 1.d0

               sq = ts(1)*qu(1)+ts(2)*qu(2)+ts(3)*qu(3)

               do i = 2,nsun_deg
                  qs = RbyRs(i)*(sleg(i+1)*u(j)
     &                 -(dble(i)*(sleg(i+1)-cphi_s*sleg(i))*u(j)
     &                 +dble(i+1)*(cphi_s*sleg(i+1)-sleg(i+2))*sq
     &                 /Rs(1))/den_s)
                  qsun(j) = qs + qsun(j)
               end do

               oq = tm(1)*qu(1)+tm(2)*qu(2)+tm(3)*qu(3)

               do i = 2,nmoon_deg
                  qm = RbyRm(i)*(oleg(i+1)*u(j)
     &                 -(dble(i)*(oleg(i+1)-cphi_m*oleg(i))*u(j)
     &                 +dble(i+1)*(cphi_m*oleg(i+1)-oleg(i+2))*oq
     &                 /Rm(1))/den_m)
                  qmoon(j) = qm + qmoon(j)
               end do

c     -----------------------------------------------------------
c     State: conservative + SRP, energy corr (no drag term)
c     -----------------------------------------------------------
               aqt = sie * 0.25d0
               bqt = cse * 0.25d0

               z(j+1) = (q(j)*amue + first_s*qsun(j)
     &              + first_m*qmoon(j)
     &              - R(1)*qsrp(j))
     &              * aqt/(w*w)
     &              + 2.d0*z(10)*us(j)*sie/w

               z(j+5) = -(q(j)*amue + first_s*qsun(j)
     &              + first_m*qmoon(j)
     &              - R(1)*qsrp(j))
     &              * bqt/(w*w)
     &              - 2.d0*z(10)*us(j)*cse/w

            end do

            I_rkg = rkgil(n,y,z,E,dE,nt)
            call ks2ksr(y,u,us,E,cse,sie)
            call ks2car(u,us,x,xd,w)
            call car2oe(x,xd,pek)

            w = y(10)
            T = y(1) - dotp4(u,us)/w

            R(1) = dotp4(u,u)

            do i = 2, 20
               R(i) = R(i-1) * R(1)
            end do

            V = vmn(xd)
            ZbyR(1) = x(3)/R(1)
            ObyR(1) = 1.d0/R(1)
            do i = 1,ngeo_deg
               ZbyR(i+1) = ZbyR(i)*ZbyR(1)
               ObyR(i+1) = ObyR(i)*ObyR(1)
            end do
            FracJd = T/86400.d0
            djulian1 = djulian + FracJd
            call solarnpv(djulian1, sun_wrt_earth)
            call lunarpv(djulian1, moon_wrt_earth)
            do i = 1,3
               ts(i) = -moon_wrt_earth(i)
               tm(i) = sun_wrt_earth(i) - moon_wrt_earth(i)
            end do

            call third_body_aux(nsun_deg,  x, ts, R(1), amuS,
     &                          Rs(1), RbyRs, first_s, cphi_s, cphi_s2)
            call third_body_aux(nmoon_deg, x, tm, R(1), amuM,
     &                          Rm(1), RbyRm, first_m, cphi_m, cphi_m2)

            if (I_rkg .le. 1) goto 60

            E_0 = E_0+dE_0
            iE = dint(E_0*r2d)

!     Buffer state at every completed integration step
            idump = idump + 1
            traj_jd(idump) = djulian + T/86400.d0
            do i = 1,3
               traj_x(i,idump)  = x(i)
               traj_xd(i,idump) = xd(i)
            end do

c     -----------------------------------------------------------
c     Impact / divergence check (no re-entry concept for an
c     airless body -- "impact" is surface contact, R_Earth here
c     holding the Moon's own radius per the common-block note
c     above). Same divergence handling as KSROP's Earth driver.
c     -----------------------------------------------------------
            h_alt = R(1) - R_Earth
            if (h_alt .ne. h_alt) then
               idump = idump - 1
               call jd2epoch(traj_jd(idump), epochstr)
               write(*,'(A,A,A)')
     &            'Integration diverged (non-finite state); ',
     &            'last valid epoch was ', epochstr
               write(*,'(A)') 'Propagation stopped; last valid '//
     &            'state saved as the last OEM data point.'
               goto 90
            else if (h_alt .lt. 0.d0) then
               call jd2epoch(traj_jd(idump), epochstr)
               write(*,'(A,F10.3,A,A)')
     &            'Lunar surface impact: altitude = ', h_alt,
     &            ' km (< 0 km) at epoch ', epochstr
               write(*,'(A)') 'Propagation stopped; final state '//
     &            'saved as the last OEM data point.'
               goto 90
            end if

         end do

         write(4,101) iE,y
      end do

  90  continue

c     --- Write CCSDS OEM v2.0 after run completes ---
      oem_fname = 'output/KSROP_LUNAR_'//compact_epoch//'.oem'
      open (unit = 3, file = oem_fname)

      call jd2epoch(traj_jd(1),     oem_start)
      call jd2epoch(traj_jd(idump), epoch_stop)

      call write_oem(3, epochstr, oem_start, epoch_stop,
     &               traj_jd, traj_x, traj_xd, idump)
      close(3)
      write(*,'(A,A)') 'OEM written: ', oem_fname
      deallocate(traj_jd, traj_x, traj_xd)
      write(*,*) x

c     format specifiers
 101  format(10f20.9)           !KS-regular elements (regular.out)

 99   do i = 1,5
         close(i)
      end do
      stop
      end
