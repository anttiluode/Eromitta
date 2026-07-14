"""
EXP2b -- AUTOPSY OF THE FIRED KILL (post-hoc, labeled as such)
==============================================================
EXP2's registered P1 predicted gradient screening is sterile. The kill
condition fired the other way: nvort(gradient) = 2114 vs nvort(intensity)
= 402. Before that inversion is allowed to be a result, it must survive
four artifact checks. Steep fronts are exactly what gradient screening
makes, and steep fronts are exactly what aliases on a grid.

CHECKS (thresholds fixed before running this file):

A1  CORE AMPLITUDE. Real vortices sit at amplitude zeros. Mean |psi|^2 on
    winding plaquettes must be < 0.5x the field mean for both arms.
    If gradient-arm "vortices" sit on ordinary-amplitude pixels, they are
    phase noise, not defects.

A2  CHARGE BALANCE. |sum(wind)| / nvort < 0.05 for both arms (defects
    nucleate in +/- pairs; net topological charge ~ 0).

A3  TIMESTEP ROBUSTNESS. Rerun gradient arm at dt=0.03 (half). Spurious
    windings from integration error are dt-sensitive. PASS if nvort changes
    by < 40%.

A4  RESOLUTION ROBUSTNESS. Rerun gradient arm at N=192 (protocol scales
    envelope with N). If defect AREAL DENSITY is grid-pinned (artifact),
    count scales ~N^2 with fixed per-pixel density AND cores fail A1.
    Record density and re-run A1 at 192. PASS if A1 holds at 192.

A5  SPECTRAL PILE-UP (recorded, not pass/fail): fraction of |psi| power
    in the top half of k-space, both arms. Large pile-up = grid-scale
    dynamics; interpret counts with suspicion.

Do not hype. Do not lie. Just show.
"""
import numpy as np, json, os, sys
sys.path.insert(0, os.path.dirname(__file__))
from exp2_four_screenings import run

def core_stats(r):
    I = r['I']
    mask = r['wind'] != 0
    return float(I[mask].mean() / I.mean()), int(np.abs(r['wind'].sum())), int(mask.sum())

out = {}

rB = run('intensity')
rC = run('gradient')
for name, r in [('B_intensity', rB), ('C_gradient', rC)]:
    ratio, netQ, n = core_stats(r)
    out[name] = dict(nvort=n, core_I_over_mean_I=ratio, net_charge=netQ,
                     A1_pass=bool(ratio < 0.5),
                     A2_pass=bool(netQ < 0.05 * max(n, 1)))

rC_dt = run('gradient', dt=0.03)
out['A3_dt_half'] = dict(nvort=rC_dt['nvort'],
                         change=abs(rC_dt['nvort'] - rC['nvort']) / rC['nvort'],
                         A3_pass=bool(abs(rC_dt['nvort'] - rC['nvort']) / rC['nvort'] < 0.40))

rC_192 = run('gradient', N=192)
ratio192, netQ192, n192 = core_stats(rC_192)
out['A4_N192'] = dict(nvort=n192,
                      density_128=rC['nvort'] / 128**2,
                      density_192=n192 / 192**2,
                      core_I_over_mean_I=ratio192,
                      A4_pass=bool(ratio192 < 0.5))

def hi_k_frac(r):
    F = np.abs(np.fft.fftshift(np.fft.fft2(r['I'] - r['I'].mean())))**2
    N = F.shape[0]
    x = np.arange(N) - N // 2
    KX, KY = np.meshgrid(x, x, indexing='ij')
    K = np.sqrt(KX**2 + KY**2)
    return float(F[K > N / 4].sum() / F.sum())

out['A5_hi_k_power'] = dict(B_intensity=hi_k_frac(rB), C_gradient=hi_k_frac(rC))

base = os.path.join(os.path.dirname(__file__), '..')
with open(os.path.join(base, 'results', 'exp2b_autopsy.json'), 'w') as f:
    json.dump(out, f, indent=2)
print(json.dumps(out, indent=2))
