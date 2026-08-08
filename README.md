# KSROP-Lunar — Moon-Centered KS Regular Orbit Propagator

Orbit propagation using **Kustaanheimo–Stiefel (KS) regular elements** with a
Runge–Kutta–Gill 4th-order integrator, for a **Moon-centered** orbiter. Built
on top of [KSROP](https://github.com/hari251086/KSROP) (the Earth-orbiting
propagator) as an fpm library dependency — the KS engine, integrator, and
force-model math are reused unmodified; only the central body and the
physical meaning of the two third-body force "slots" change.

**Author:** Harishkumar Sellamuthu · hari251086@gmail.com
**Copyright:** 2026, Harishkumar Sellamuthu, All Rights Reserved

---

## 1. Overview

KSROP-Lunar propagates a spacecraft orbiting the Moon under lunar
oblateness, Earth + Sun third-body gravity, and solar radiation pressure.
It exists to answer a specific question raised in
[KSROP issue #26](https://github.com/hari251086/KSROP/issues/26): does
KSROP's "one dominant central body + weaker third-body perturbations"
architecture generalize to a non-Earth central body with only configuration
changes, or does it need new dynamics? The answer, worked out here: **it
generalizes almost entirely as a config/data change** — the KS transform,
RKG4 integrator, and third-body force formulas (`third_body_aux`,
`qsun`/`qmoon` KS-EOM convention) are all central-body-agnostic in KSROP
itself, and Earth's position relative to the Moon and the Sun's position
relative to the Moon are both obtainable by vector arithmetic on KSROP's
*existing* geocentric `solarnpv`/`lunarpv` ephemerides — no new ephemeris
series needed. See `app/driver_KS_lunar.F`'s header comment (THIRD-BODY SLOT
REPURPOSING) for the exact construction.

This repo depends on `KSROP` (git tag `v2.8.0`) and is independent of every
other repo under `GitHub\` otherwise. It does **not** cover getting a
spacecraft *to* the Moon — a genuine translunar transfer passes through a
regime where Earth and Moon gravity are locally comparable, which breaks
this single-dominant-center architecture; that is tracked separately as
[KSROP issue #28](https://github.com/hari251086/KSROP/issues/28) and is
explicitly blocked on this repo landing first.

---

## 2. Project Structure

```
KSROP-Lunar/
├── fpm.toml                        fpm manifest; depends on KSROP v2.8.0
├── src/
│   └── moon_rotation.F             Moon prime-meridian rotation angle (issue #30)
├── app/
│   └── driver_KS_lunar.F           Moon-centered propagator (see header
│                                    comment for the third-body slot
│                                    repurposing convention)
├── input/
│   ├── const_moon.dat              mu_Moon, R_Moon, AU, mu_Earth, mu_Sun
│   ├── GRAIL_lowdeg_zonal.dat       Lunar J2 zonal coefficient (see §8)
│   ├── input.opm                   Initial state (CCSDS OPM v2.0)
│   └── input.dat                   Sim params (nrev/istep/tole, force
│                                    flags, SRP -- no drag line, see §7)
├── output/                         Runtime-generated OEM/OPM/debug files
├── test/
│   └── test_lunar_regimes.F        Fortran unit tests (6 checks)
├── test_lunar_regimes.py           Multi-regime integration test (110
│                                    checks, 10 lunar orbit regimes)
└── .github/workflows/ci.yml        CI (gfortran on Ubuntu, matching
                                     KSROP's own CI pattern)
```

---

## 3. Quick Start

```bash
fpm build --compiler gfortran
fpm test  --compiler gfortran   # 6 unit checks
```

Then run the propagator directly (reads the tracked `input/` files as-is —
a 100 km altitude circular low lunar orbit):

```bash
fpm run driver_KS_lunar --compiler gfortran
```

Writes `output/KSROP_LUNAR_<timestamp>.oem` (trajectory) and
`output/ksrop_lunar.opm` (initial elements).

---

## 4. Building

Requires **Intel oneAPI Fortran (`ifx`)** or **GNU Fortran (`gfortran`)**,
same as KSROP itself.

```bash
fpm build --compiler gfortran   # works on Windows and Linux -- see note below
fpm build --compiler ifx        # also works
```

### Resolved compiler issue (both toolchains work now)

During initial development, a `gfortran`-built `driver_KS_lunar.exe` crashed
with a memory corruption (a loop-control integer got clobbered during an
ephemeris call) — reproduced on **both** Windows/MinGW gfortran and Linux
gfortran (caught by this repo's own CI), so it was a genuine gfortran
stack-layout/ABI quirk, not a Windows-only issue, and not a bug in KSROP's
physics. Root-caused via an isolated repro to a specific pattern (two
same-size local arrays declared in one combined `dimension` statement,
followed by a scalar, with one array passed to an externally-compiled
subroutine) and **fixed** by moving the vulnerable loop-control scalars
(`nrev`, `istep`, `ik`, `ki`, `idump`) into their own `common /loopctl/`
block — static storage, not stack-allocated, so immune to the mechanism
regardless of its exact cause. See `app/driver_KS_lunar.F`'s header comment
for the full writeup. Confirmed clean on `gfortran` (Windows and Linux CI)
and `ifx` alike after the fix.

`fpm.toml` sets `[fortran] source-form = "fixed"` /
`implicit-typing = true` / `implicit-external = true`, same reasons as
KSROP itself (F77-style fixed-form, implicit typing).

---

## 5. Running

```bash
fpm run driver_KS_lunar --compiler gfortran
```

Reads `input/input.opm` (initial state, Moon-centered), `input/input.dat`
(sim params), `input/const_moon.dat` (constants + force-model degrees).
Writes `output/KSROP_LUNAR_<timestamp>.oem` (CCSDS OEM v2.0 trajectory),
`output/KSROP_LUNAR_<timestamp>_Regular.out` (KS-elements debug dump), and
`output/ksrop_lunar.opm` (initial osculating elements).

---

## 6. Testing

```bash
fpm test --compiler gfortran                          # 6 unit checks
python test_lunar_regimes.py <path-to-driver_KS_lunar.exe>   # 110 checks
```

**Total: 116 automated checks**, all passing (2026-08-07, `gfortran` on
Windows and Linux CI, and `ifx`). The Python integration test sweeps 10 lunar orbit regimes — low
lunar orbit (equatorial and polar/LRO-like), an eccentric case, a
near-equatorial case, the four inclinations at which lunar orbits are known
to be long-term "frozen" (27°, 50°, 76°, 86° — Elipe & Lara 2003; Nie &
Gurfil 2018), a higher circular case, and a retrograde case — each run
twice: pure two-body (exact conservation checks) and full dynamics (lunar
J2 + Earth/Sun third-body + SRP, bounded-drift checks). This mirrors
KSROP's own `test_initial_conditions.py` pattern and orbit-regime coverage
philosophy, adapted to lunar orbital scale (mu ≈ 4903 km³/s², r ≈
1800–2300 km — three orders of magnitude smaller than KSROP's usual Earth
GTO/GEO range).

CI (`.github/workflows/ci.yml`) runs `fpm build`/`fpm test` with `gfortran`
on `ubuntu-latest`, matching KSROP's own CI job.

---

## 7. Inputs & Outputs

### `const_moon.dat` — Physical constants

```
mu_Moon  R_Moon  AU  mu_Earth  mu_Sun
ngeo_deg  n_earth_deg  n_sun_deg
PSR_srp
```

| Parameter | Value | Source |
|---|---|---|
| `mu_Moon` | 4902.8001 km³/s² | JPL Horizons / GRAIL |
| `R_Moon` | 1737.4 km | Volumetric mean radius (JPL) |
| `mu_Earth` | 398600.4415 km³/s² | Same constant KSROP uses for Earth |
| `mu_Sun` | 1.32712440018×10¹¹ km³/s² | Same constant KSROP uses for Earth |
| `ngeo_deg` | 0–2 (ships with J2 only, see §8) | — |

### `input.dat` — Simulation parameters (no drag line — the Moon is airless)

```
nrev  istep  tole
n_force(1)=lunar geo   n_force(2)=Earth 3rd-body   n_force(3)=Sun 3rd-body
CR  AM(m2/kg)  IPSR(0/1)  ISHAD(0/1/2)
```

### `input.opm` — Initial state (CCSDS OPM v2.0, `CENTER_NAME = MOON`)

Same format as KSROP's Earth-centered OPM, Moon-centered. `REF_FRAME =
EME2000` denotes axes parallel to the Earth-equatorial J2000 frame,
translated (not rotated) to the Moon's center — **not** the Moon's own
body-fixed rotating frame.

### Output: same CCSDS OEM/OPM v2.0 formats as KSROP itself.

---

## 8. Known Issues / Limitations

- ~~Gravity field is zonal-only (J2), and that's a real ceiling, not just
  a missing data file~~ — **fixed 2026-08-08**: `input/GRAIL_tess22.dat`
  ships the real degree-2 order-2 term, C̄₂₂ = 3.4673798×10⁻⁵, S̄₂₂ ≈
  −2.5×10⁻¹⁰ (consistent with zero, as expected for the Moon's
  tidally-locked long axis) — GRAIL Primary Mission solution (Konopliv
  et al. 2013, *JGR: Planets* 118, Table 4), wired via KSROP's general
  `(n,m)` tesseral support ([KSROP#30](https://github.com/hari251086/KSROP/issues/30)).
  A new `moon_rotation_angle_deg` (`src/moon_rotation.F`) supplies the
  lunar prime-meridian angle this needs — see ALGORITHM.md and Version
  History for the source, and an important caveat: this uses only the
  *secular* part of the IAU rotation model, omitting the Moon's real
  (multi-degree-amplitude) physical libration, a larger relative gap
  than the analogous Earth/Mars simplifications. The Moon's real field
  remains **mascon-dominated** beyond (2,2) — a fuller GRGM-lineage
  table (e.g. GRGM900C) would need no further code change
  (`geo_coeff_tess_general` already supports arbitrary `(n,m)`), just
  more coefficient rows.
- **No atmospheric drag** — the Moon is airless; this is by design, not a
  gap.
- **SRP shadow radius correctly reuses the Moon's own radius** (the
  `R_Earth` common-block slot, populated from `const_moon.dat`, is consumed
  by KSROP's `shadfncyl`/`shadfncone` as-is — they were already
  body-agnostic, no change needed).
- **No GMAT (or other independent) cross-validation yet** — KSROP's own
  Earth-centered full-force model was validated against GMAT to ~1.9 km
  over 2 GTO revolutions (see KSROP's `project_ksrop_gmat_validation`
  history); an equivalent lunar-centered validation campaign has not been
  done here. The test suite validates internal consistency (conservation
  laws, bounded drift) but not absolute accuracy against an independent
  tool.
- Earth/Sun third-body Legendre degree ships at 2 (`const_moon.dat` line 2)
  — matches KSROP's own Sun-degree convention; adjustable per run.

---

## 9. Version History

| Date | Change |
|---|---|
| 2026-08-07 | Initial repo: `driver_KS_lunar.F` (Moon-centered driver, reusing KSROP v2.3.0's KS engine/force-model math unmodified via the third-body slot repurposing convention), lunar J2 coefficient file (JPL/GRAIL-sourced), unit tests (6 checks) and multi-regime integration test (110 checks across 10 lunar orbit regimes, all passing). Found and fixed a real gfortran stack-layout bug (reproduced on both Windows and Linux CI) by isolating loop-control scalars into a COMMON block — both `gfortran` and `ifx` now build and run clean. Companion to KSROP issue #26. |
| 2026-08-08 | **Tesseral (2,2) rollout from KSROP#30**: bumped the KSROP dependency to v2.8.0 (general `(n,m)` tesseral/mascon support). Added `moon_rotation_angle_deg` (`src/moon_rotation.F`) — the Moon's prime-meridian rotation angle, `W = 38.3213 + 13.17635815·d` (`d`=days from J2000), sourced directly from NASA NAIF's public `pck00011.tpc` SPICE kernel (IAU 2009 rotation model); the rate matches the Moon's known synchronous rotation (360°/27.321661-day sidereal month = 13.17640°/day) to 5 significant figures. **Deliberately secular-only**: the full IAU model has 13 periodic physical-libration correction terms (also in the same kernel, `BODY301_NUT_PREC_PM`) that are *not* implemented — a materially larger relative gap than the analogous Earth/Mars simplifications, since real lunar libration amplitude is degrees, not arcseconds. Added `input/GRAIL_tess22.dat` with the real GRAIL Primary Mission (2,2) coefficients (Konopliv et al. 2013, Table 4), and wired `geo_coeff_tess_general`/`rotate_tess_coeffs`/`tess_general_force` into `driver_KS_lunar.F`'s propagation loop, mirroring KSROP's own driver wiring pattern exactly. Verified genuinely active via a before/after comparison at 100 km altitude (a=1837.4 km) over 10 revolutions: real vs. zeroed (2,2) coefficients diverge by ~12 m — consistent in scale with KSROP's own Earth LEO validation for this same non-resonant regime. 116/116 existing checks still pass unchanged. |

---

## 10. Dependencies / References

- **KSROP** (git tag `v2.8.0`): the KS-regularized propagation engine,
  force-model math, and CCSDS I/O are all consumed as an fpm library
  dependency, unmodified — see `fpm.toml`.
- **Lunar physical constants**: JPL Horizons / GRAIL mission
  (GM = 4902.8001 km³/s², R = 1737.4 km,
  [Lunar Fact Sheet](https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html);
  JPL SSD [Planetary Physical Parameters](https://ssd.jpl.nasa.gov/planets/phys_par.html)).
- **Lunar J2**: 202.7×10⁻⁶, from GRAIL-era lunar gravity field solutions
  (e.g. Lemoine et al. 2014, GRGM900C, *Geophysical Research Letters*).
- **Frozen lunar orbit inclinations** (27°, 50°, 76°, 86°): Elipe & Lara
  (2003), "Frozen Orbits About the Moon"; Nie & Gurfil (2018) generalize
  and refine the same result.
- **Related**: [KSROP-Mars](https://github.com/hari251086/KSROP-Mars),
  the analogous Mars-centered extension, built the same day using the same
  slot-repurposing pattern.
