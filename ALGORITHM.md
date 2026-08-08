# ALGORITHM.md — KSROP-Lunar

## 1. Overview
KSROP-Lunar is a Moon-centered orbit propagator built on top of `KSROP`
(the Earth-orbiting KS-regularized propagator, consumed here as an fpm
library dependency, git tag `v2.3.0`). It exists to answer KSROP issue #26
(does KSROP's central-body-plus-third-body architecture generalize beyond
Earth?) and is a sibling to `KSROP-Mars`. It has no other dependency
relationship within `GitHub\`.

## 2. Problem Statement
Numerically integrate a spacecraft's trajectory forward in time around the
Moon, under lunar oblateness (zonal, currently J2-only — see §9), Earth and
Sun third-body gravity, and solar radiation pressure, given an initial
Moon-centered Cartesian state and epoch. "Correct" means the same thing it
means in KSROP itself: energy/angular-momentum/orbit-closure conservation
to machine precision for the unperturbed two-body case, and bounded
(not exactly conserved, but physically plausible) drift under full
dynamics — verified in `test_lunar_regimes.py` across 10 lunar orbit
regimes (low lunar orbit, eccentric, near-equatorial, four "frozen"
inclinations, higher circular, retrograde).

## 3. Inputs
Same field structure as KSROP's own `driver_KS.F`, with two differences:
`const_moon.dat`'s five constants are `mu_Moon R_Moon AU mu_Earth mu_Sun`
(not `mu_Earth R_Earth AU mu_Sun mu_Moon`), and `input.dat` has **no drag
line** (3 lines: `nrev/istep/tole`, force flags, SRP params — see
`README.md` §7). Initial state comes from `input/input.opm`
(`CENTER_NAME = MOON`).

## 4. Core Algorithm
1. **Initialization**: identical to KSROP's own driver — read
   `const_moon.dat`/`input.dat`/`input.opm`, `car2oe`/`oe2car` round-trip,
   `force_models` resolves on/off flags. The lunar gravity coefficient is
   loaded via `geo_coeff_body(ngeo_deg, c_j, 'input/GRAIL_lowdeg_zonal.dat')`
   — KSROP v2.3.0's file-parameterized entry point (added specifically to
   support this repo, see KSROP's own `README.md` revision history), rather
   than the Earth-hardcoded `geo_coeff` wrapper.
2. **Third-body ephemerides — the one genuinely new piece of logic**:
   KSROP's `solarnpv(dj,s)`/`lunarpv(dj,tm)` return the Sun's and Moon's
   position **as seen from Earth** (geocentric, EME2000). For a
   Moon-centered orbiter the two relevant third bodies are Earth and Sun,
   obtained by vector arithmetic on those same two calls, not a new
   ephemeris series:
   ```
   Earth position relative to Moon = -lunarpv(dj)
   Sun position relative to Moon   = solarnpv(dj) - lunarpv(dj)
   ```
   These feed into KSROP's existing two third-body "slots" (historically
   named `ts`/`amuS`/`nsun_deg` for the Sun and `tm`/`amuM`/`nmoon_deg` for
   the Moon) — repurposed here to mean Earth and Sun respectively. See
   `app/driver_KS_lunar.F`'s header comment (THIRD-BODY SLOT REPURPOSING)
   for the full mapping. `third_body_aux`/the `qsun`/`qmoon` KS-EOM force
   formulas are consumed **completely unmodified** from KSROP — they only
   ever operate on a generic position vector + GM, nothing Sun/Moon-specific
   internally.
3. **KS transform / integration loop / SRP**: identical to KSROP's own
   driver (`car2ks`/`ks2car`, RKG4 via `rkgil`, cannonball SRP with
   `shadfncyl`/`shadfncone`) — see KSROP's own `ALGORITHM.md` §4 for the
   full step-by-step, which applies here unchanged. The one behavioral
   difference: SRP's Sun-direction geometry uses the **`tm`** slot (Sun,
   in this repo's repurposing), not `ts` — the opposite of KSROP's
   Earth-centered driver, where the Sun literally is the `ts` slot.
4. **No drag**: this driver never opens an atmosphere table and never
   computes a drag force (the Moon is airless) — `qdrag`/`p_drag` and the
   entire per-revolution atmosphere setup block present in KSROP's
   `driver_KS.F` are simply absent here, not zeroed-out.
5. **Termination**: `h_alt = R(1) - R_Moon < 0` is treated as lunar surface
   impact (KSROP's Earth driver uses a `< 80 km` re-entry threshold, which
   has no lunar analog — the Moon has no atmosphere to define a re-entry
   boundary against, so impact is genuine surface contact, `h_alt < 0`).

## 5. Key Equations / Physics
Identical to KSROP's own `ALGORITHM.md` §5 (Sundman transform, KS regular
elements, angular-rate scaling Γ = w/w_Kep, oblateness potential) — none of
this changes for a different central body; only `amue`/`R_Earth` (holding
the Moon's GM/radius here, via the shared `/xy/` common block) and the two
third-body inputs change. See §4 above for the one real physics-adjacent
addition (the third-body vector construction).

## 6. Outputs
Same CCSDS OEM v2.0 trajectory / KS-elements debug dump / OPM initial
elements as KSROP, filenamed `KSROP_LUNAR_*` instead of `KSROP_*`.

## 7. Complexity & Performance
Same per-step cost structure as KSROP (dominated by the zonal-harmonic sum
degree, currently trivial at `ngeo_deg=2`). No parallelism, same as KSROP;
the `GitHub\CLAUDE.md` 4-core cap is not directly exercised (a single
propagation run is inherently sequential).

## 8. Validation & Accuracy
116 automated checks (6 Fortran unit tests + 110 Python integration
checks across 10 lunar orbit regimes, two-body and full-dynamics phases),
all passing as of 2026-08-07 (`gfortran` on Windows and Linux CI, and
`ifx`). Internal-consistency
validation only (conservation laws, bounded drift) — **no independent
cross-validation** (e.g. against GMAT with a real lunar force model) has
been done yet, unlike KSROP's own Earth-centered campaign (~1.9 km over 2
GTO revolutions vs GMAT, see KSROP's `project_ksrop_gmat_validation`
project memory). A GMAT (or other independent tool) lunar cross-check is
open follow-on work, not yet an issue.

## 9. Known Limitations
- **J2-only lunar gravity field, and it's a real ceiling** — see README
  §8; the Moon's real field is mascon-dominated (strong C22 and
  higher-degree/order terms) and KSROP's zonal-only force-model pipeline
  cannot consume those terms regardless of what coefficient file is
  supplied (KSROP#29).
- **No independent (GMAT) cross-validation** — see §8.
- **`EME2000`-labeled but Moon-centered frame** — axes parallel to
  Earth-equatorial J2000, translated (not rotated) to the Moon's center;
  not the Moon's own body-fixed rotating frame (e.g. not `IAU_MOON`/
  `MOON_ME`). Fine for inertial dynamics, but a consumer expecting a
  lunar-body-fixed frame would need an additional rotation this repo does
  not provide.

## 10. Dependencies
- **KSROP** (git tag `v2.3.0`): consumed as an fpm library dependency —
  the entire KS engine, integrator, force-model math, and CCSDS I/O are
  KSROP's own, unmodified. `geo_coeff_body` (KSROP v2.3.0+) was added
  upstream specifically to support this repo's lunar coefficient file.
- **No other `GitHub\` repo dependency.**
- **External data**: none required beyond the two small files this repo
  ships (`GRAIL_lowdeg_zonal.dat`, tiny; no large external download needed,
  unlike KSROP's own 231 MB EGM2008 file, since this repo intentionally
  ships only J2).
