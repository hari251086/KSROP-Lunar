#!/usr/bin/env python3
"""
test_lunar_regimes.py - Multi-case integration tests for KSROP-Lunar

Same two-phase pattern as KSROP's own test_initial_conditions.py, but
for a Moon-centered orbiter instead of Earth-centered:

Phase 1 - Two-body (all perturbations off):
  Verifies energy, angular momentum, orbit closure, and semi-major
  axis conservation to machine precision, at LUNAR mu/radius scale
  (mu ~ 4903 km3/s2, r ~ 1740-3000 km -- three orders of magnitude
  smaller than KSROP's usual Earth GTO/GEO scale).

Phase 2 - Full dynamics (lunar J2 + Earth/Sun third-body + SRP):
  Verifies the propagator completes without divergence, all states
  remain finite, altitude stays physical, and osculating energy /
  angular momentum / semi-major axis drift stay bounded.

Orbit regimes cover low lunar orbit, elliptical, polar, near-equatorial,
and the four inclinations (27, 50, 76, 86 deg) at which lunar orbits
are known to be long-term "frozen" (Elipe & Lara 2003; Nie & Gurfil
2018) -- the lunar-orbit-design analog of KSROP's Earth-side sweep
covering GTO/Molniya/SSO/critical-inclination.

Usage
-----
  python test_lunar_regimes.py [executable]
"""

import subprocess, os, shutil, math, sys, glob

MU      = 4902.8001       # km3/s2 (JPL/GRAIL)
R_MOON  = 1737.4          # km (volumetric mean radius)

# -------------------------------------------------------------------
# Vector/orbital helpers (identical to KSROP's own test harness)
# -------------------------------------------------------------------
def mag(v):
    return math.sqrt(sum(c*c for c in v))

def dot(a, b):
    return sum(x*y for x, y in zip(a, b))

def cross3(a, b):
    return [a[1]*b[2]-a[2]*b[1],
            a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0]]

def orbital_energy(x, xd):
    return dot(xd, xd)/2.0 - MU/mag(x)

def semimajor(x, xd):
    return -MU / (2.0*orbital_energy(x, xd))

def ang_momentum_mag(x, xd):
    return mag(cross3(x, xd))

def altitude(x):
    return mag(x) - R_MOON

def is_finite(v):
    return all(math.isfinite(c) for c in v)

def oe2cart(a, e, i_deg, raan_deg, aop_deg, nu_deg):
    """Keplerian elements to Cartesian state (km, km/s)."""
    i    = math.radians(i_deg)
    raan = math.radians(raan_deg)
    aop  = math.radians(aop_deg)
    nu   = math.radians(nu_deg)

    p = a * (1.0 - e*e)
    r = p / (1.0 + e*math.cos(nu))

    rx_pf = r * math.cos(nu)
    ry_pf = r * math.sin(nu)
    vx_pf = -math.sqrt(MU/p) * math.sin(nu)
    vy_pf =  math.sqrt(MU/p) * (e + math.cos(nu))

    cosO, sinO = math.cos(raan), math.sin(raan)
    cosw, sinw = math.cos(aop),  math.sin(aop)
    cosi, sini = math.cos(i),    math.sin(i)

    l1 = cosO*cosw - sinO*sinw*cosi
    l2 = sinO*cosw + cosO*sinw*cosi
    l3 = sinw*sini
    m1 = -cosO*sinw - sinO*cosw*cosi
    m2 = -sinO*sinw + cosO*cosw*cosi
    m3 = cosw*sini

    x  = [l1*rx_pf + m1*ry_pf, l2*rx_pf + m2*ry_pf, l3*rx_pf + m3*ry_pf]
    xd = [l1*vx_pf + m1*vy_pf, l2*vx_pf + m2*vy_pf, l3*vx_pf + m3*vy_pf]
    return x, xd

# -------------------------------------------------------------------
# OEM parser (identical to KSROP's own)
# -------------------------------------------------------------------
def read_oem(path):
    rows = []
    in_data = False
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s == 'DATA_START':
                in_data = True
                continue
            if s == 'DATA_STOP':
                in_data = False
                continue
            if not in_data:
                continue
            vals = s.split()
            if len(vals) >= 7:
                x  = [float(v) for v in vals[1:4]]
                xd = [float(v) for v in vals[4:7]]
                rows.append((x, xd))
    return rows

# -------------------------------------------------------------------
# Assertion helper
# -------------------------------------------------------------------
class Counter:
    def __init__(self):
        self.npass = 0
        self.nfail = 0
    @property
    def total(self):
        return self.npass + self.nfail

def check(ctr, name, ok, detail=''):
    tag = 'PASS' if ok else 'FAIL'
    msg = f'      [{tag}]  {name}'
    if not ok and detail:
        msg += f'  ({detail})'
    print(msg)
    if ok:
        ctr.npass += 1
    else:
        ctr.nfail += 1
    return ok

# -------------------------------------------------------------------
# Test cases: (name, a_km, e, i_deg, raan_deg, aop_deg, nu_deg, nrev)
# a is measured from the Moon's center; altitudes noted in comments.
# -------------------------------------------------------------------
ORBITS = [
    ("Low lunar orbit (100 km alt, equatorial)",
                                    1837.4,  0.001,   0.0,   0.0,   0.0,  45.0, 1),
    ("Low lunar orbit, polar (LRO-like)",
                                    1837.4,  0.001,  90.0,   0.0,   0.0,  45.0, 1),
    ("Elliptical, i=45",           2200.0,  0.20,   45.0,  60.0,  90.0,  30.0, 1),
    ("Near-equatorial, i=5",       2000.0,  0.05,    5.0, 180.0,  90.0,  60.0, 1),
    ("Frozen-candidate i=27",      2000.0,  0.05,   27.0,   0.0,   0.0,  90.0, 1),
    ("Frozen-candidate i=50",      2000.0,  0.05,   50.0,  90.0,   0.0,  90.0, 1),
    ("Frozen-candidate i=76",      2000.0,  0.05,   76.0,  90.0,   0.0,  90.0, 1),
    ("Frozen-candidate i=86",      2000.0,  0.05,   86.0,  90.0,   0.0,  90.0, 1),
    ("Higher circular (500 km alt)",
                                    2237.4,  0.001,  30.0,   0.0,   0.0, 120.0, 1),
    ("Retrograde, i=150",          2100.0,  0.10,  150.0,  45.0,  45.0,  30.0, 1),
]

# -------------------------------------------------------------------
# Input-file writers
# -------------------------------------------------------------------
def write_opm(x, xd):
    with open('input/input.opm', 'w') as f:
        f.write('CCSDS_OPM_VERS = 2.0\n')
        f.write('CREATION_DATE  = 2026-08-07T00:00:00.000\n')
        f.write('ORIGINATOR     = KSROP-Lunar\n\n')
        f.write('META_START\n')
        f.write('OBJECT_NAME    = LUNAR_ORBITER\n')
        f.write('CENTER_NAME    = MOON\n')
        f.write('REF_FRAME      = EME2000\n')
        f.write('TIME_SYSTEM    = UTC\n')
        f.write('META_STOP\n\n')
        f.write('STATE_VECTOR\n')
        f.write('EPOCH          = 2026-08-07T00:00:00.000\n')
        f.write(f'X              = {x[0]:20.9f} [km]\n')
        f.write(f'Y              = {x[1]:20.9f} [km]\n')
        f.write(f'Z              = {x[2]:20.9f} [km]\n')
        f.write(f'X_DOT          = {xd[0]:20.12f} [km/s]\n')
        f.write(f'Y_DOT          = {xd[1]:20.12f} [km/s]\n')
        f.write(f'Z_DOT          = {xd[2]:20.12f} [km/s]\n')

def write_two_body(nrev, istep=360):
    with open('input/const_moon.dat', 'w') as f:
        f.write('4902.8001D0 1737.4D0 1.495978707d08 0.0d0 0.0d0\n')
        f.write('0 0 0\n')
        f.write('4.56d-6\n')
    with open('input/input.dat', 'w') as f:
        f.write(f'{nrev} {istep} 1d-15\n')
        f.write('0 0 0\n')
        f.write('1.2 0.01 0 1\n')

def write_full_dynamics(nrev, istep=360):
    """Lunar J2 + Earth/Sun third-body (deg 2) + SRP -- see README."""
    with open('input/const_moon.dat', 'w') as f:
        f.write('4902.8001D0 1737.4D0 1.495978707d08'
                ' 3.986004415D5 1.32712440018d11\n')
        f.write('2 2 2\n')
        f.write('4.56d-6\n')
    with open('input/input.dat', 'w') as f:
        f.write(f'{nrev} {istep} 1d-15\n')
        f.write('1 1 1\n')
        f.write('1.2 0.01 1 1\n')

# -------------------------------------------------------------------
# Run the driver and return OEM rows (or None on failure)
# -------------------------------------------------------------------
def run_driver(exe, label):
    for f in glob.glob('output/KSROP_LUNAR_*.oem'):
        os.remove(f)
    try:
        result = subprocess.run(
            [exe], capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        print(f'      TIMEOUT ({label})')
        return None
    if result.returncode != 0:
        print(f'      RUNTIME ERROR ({label}): {result.stderr[:200]}')
        return None
    oem_files = sorted(glob.glob('output/KSROP_LUNAR_*.oem'))
    if not oem_files:
        print(f'      NO OEM OUTPUT ({label})')
        return None
    return read_oem(oem_files[-1])

# -------------------------------------------------------------------
# Phase 1: Two-body checks (exact conservation)
# -------------------------------------------------------------------
def run_two_body_checks(ctr, rows, expected_rows):
    x0, xd0 = rows[0]
    xf, xdf = rows[-1]
    ok = True

    ok &= check(ctr, f'Row count = {expected_rows}',
                len(rows) == expected_rows, f'got {len(rows)}')

    E0 = orbital_energy(x0, xd0)
    Ef = orbital_energy(xf, xdf)
    dE = abs(Ef - E0) / abs(E0)
    ok &= check(ctr, f'Energy conserved        (rel err {dE:.2e})',
                dE < 1e-8)

    h0 = ang_momentum_mag(x0, xd0)
    hf = ang_momentum_mag(xf, xdf)
    dh = abs(hf - h0) / abs(h0)
    ok &= check(ctr, f'Ang. momentum conserved (rel err {dh:.2e})',
                dh < 1e-8)

    dr = mag([xf[j]-x0[j] for j in range(3)])
    dv = mag([xdf[j]-xd0[j] for j in range(3)])
    ok &= check(ctr, f'Orbit closure            (dr={dr:.4e}, dv={dv:.4e})',
                dr < 0.1 and dv < 1e-5)

    a0 = semimajor(x0, xd0)
    af = semimajor(xf, xdf)
    da = abs(af - a0) / a0
    ok &= check(ctr, f'SMA conserved            (rel err {da:.2e})',
                da < 1e-6)
    return ok

# -------------------------------------------------------------------
# Phase 2: Full-dynamics checks (bounded drift)
# -------------------------------------------------------------------
def run_full_dynamics_checks(ctr, rows, expected_rows, a0_nom):
    x0, xd0 = rows[0]
    xf, xdf = rows[-1]
    ok = True

    ok &= check(ctr, f'Row count = {expected_rows}',
                len(rows) == expected_rows, f'got {len(rows)}')

    all_finite = all(is_finite(x) and is_finite(xd) for x, xd in rows)
    ok &= check(ctr, 'All states finite (no NaN/Inf)', all_finite)

    min_alt = min(altitude(x) for x, xd in rows)
    ok &= check(ctr, f'Min altitude > 0 km       (h={min_alt:.1f} km)',
                min_alt > 0.0)

    E0 = orbital_energy(x0, xd0)
    Ef = orbital_energy(xf, xdf)
    dE = abs(Ef - E0) / abs(E0)
    ok &= check(ctr, f'Energy drift bounded     (rel {dE:.4e})',
                dE < 1e-2, f'E0={E0:.6f}, Ef={Ef:.6f}')

    h0 = ang_momentum_mag(x0, xd0)
    hf = ang_momentum_mag(xf, xdf)
    dh = abs(hf - h0) / abs(h0)
    ok &= check(ctr, f'Ang. mom. drift bounded  (rel {dh:.4e})',
                dh < 1e-2, f'h0={h0:.4f}, hf={hf:.4f}')

    a0 = semimajor(x0, xd0)
    af = semimajor(xf, xdf)
    da = abs(af - a0_nom) / a0_nom
    ok &= check(ctr, f'SMA physically plausible (da/a={da:.4e})',
                da < 5e-2, f'a_nom={a0_nom:.1f}, af={af:.1f}')

    return ok

# -------------------------------------------------------------------
# Main
# -------------------------------------------------------------------
def main():
    exe = sys.argv[1] if len(sys.argv) > 1 else 'driver_KS_lunar.exe'
    exe = os.path.abspath(exe)
    if not os.path.isfile(exe):
        print(f'ERROR: executable not found: {exe}')
        sys.exit(1)

    backup = {}
    for fname in ['input/const_moon.dat', 'input/input.dat', 'input/input.opm']:
        bak = fname + '.bak_ic'
        if os.path.isfile(fname):
            shutil.copy(fname, bak)
            backup[fname] = bak

    ctr = Counter()
    failed_cases = []
    istep = 360

    # ==============================================================
    # PHASE 1: Two-body
    # ==============================================================
    print('=' * 64)
    print(' Phase 1 - Two-Body (exact conservation, lunar mu/radius)')
    print('=' * 64)

    for name, a, e, i, raan, aop, nu, nrev in ORBITS:
        x, xd = oe2cart(a, e, i, raan, aop, nu)
        write_opm(x, xd)
        write_two_body(nrev, istep)

        print(f'\n  --- {name} ---')
        print(f'      a={a:.1f} km  e={e:.4f}  i={i:.1f} deg'
              f'  alt={a - R_MOON:.0f} km')

        rows = run_driver(exe, name)
        if rows is None or len(rows) < 2:
            ctr.nfail += 5
            failed_cases.append(f'2B:{name}')
            continue

        if not run_two_body_checks(ctr, rows, nrev*istep+1):
            failed_cases.append(f'2B:{name}')

    # ==============================================================
    # PHASE 2: Full dynamics
    # ==============================================================
    print('\n' + '=' * 64)
    print(' Phase 2 - Full Dynamics (lunar J2 + Earth/Sun + SRP)')
    print('=' * 64)

    for name, a, e, i, raan, aop, nu, nrev in ORBITS:
        x, xd = oe2cart(a, e, i, raan, aop, nu)
        write_opm(x, xd)
        write_full_dynamics(nrev, istep)

        print(f'\n  --- {name} ---')
        print(f'      a={a:.1f} km  e={e:.4f}  i={i:.1f} deg')

        rows = run_driver(exe, name)
        if rows is None or len(rows) < 2:
            ctr.nfail += 5
            failed_cases.append(f'FD:{name}')
            continue

        if not run_full_dynamics_checks(ctr, rows, nrev*istep+1, a):
            failed_cases.append(f'FD:{name}')

    # ==============================================================
    # Summary
    # ==============================================================
    for fname, bak in backup.items():
        if os.path.isfile(bak):
            shutil.copy(bak, fname)
            os.remove(bak)

    print('\n' + '=' * 64)
    print(f' Total checks: {ctr.total}   Passed: {ctr.npass}'
          f'   Failed: {ctr.nfail}')
    if failed_cases:
        print(f' Failed: {", ".join(failed_cases)}')
    print('=' * 64)

    sys.exit(0 if ctr.nfail == 0 else 1)


if __name__ == '__main__':
    main()
