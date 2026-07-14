"""
EXP2d -- PERSISTENT DEFECTS, dt-CHECKED AT MATCHED PHYSICAL TIME (post-hoc)
===========================================================================
Two counters measure two different things:

  SNAPSHOT count (exp2c): every instantaneous 2pi winding, including
  short-lived turbulent froth.
  PERSISTENT count (ArrowField's counter): windings of the phase of the
  TIME-AVERAGED field over the window -- only defects that hold their
  position for ~T_window survive the averaging.

ArrowField's 402 is a persistent count. The exp2b dt-check that seemed to
show fragility (2114 -> 588) was confounded: halving dt halved the physical
duration. Here: persistent counts for all four arms at dt in {0.06, 0.03}
with matched physical time (T_burn=96, T_window=42), psi accumulated
incrementally (no giant buffers).

The question this settles: does the gradient > intensity inversion hold
for defects that PERSIST, and is the persistent count dt-robust once the
confound is removed?

Do not hype. Do not lie. Just show.
"""
import numpy as np, json, os, sys
from scipy.signal import convolve2d
sys.path.insert(0, os.path.dirname(__file__))
from exp2_four_screenings import lap, biharm, grad2

def run_persist(mode, kappa=5.0, amp=2.0, N=128, dt=0.06, damping=0.001,
                pcub=0.2, g=0.02, T_burn=96.0, T_window=42.0, seed=1):
    burn = int(round(T_burn / dt)); window = int(round(T_window / dt))
    rng = np.random.default_rng(seed)
    x = np.arange(N); X, Y = np.meshgrid(x, x, indexing='ij')
    c = N // 2; r = N / 15.0
    env = amp * np.exp(-((X - c)**2 + (Y - c)**2) / (2 * r**2))
    ph = 0.9 * rng.standard_normal((N, N))
    ph = convolve2d(ph, np.ones((5, 5)) / 25., mode='same', boundary='wrap')
    psi = env * np.exp(1j * 2.5 * ph)
    psi_old = psi.copy()
    psi_sum = np.zeros_like(psi)
    for t in range(burn + window):
        I_loc = np.abs(psi)**2
        Vp = -psi + pcub * I_loc * psi
        if mode == 'linear':
            a = lap(psi) - Vp - g * biharm(psi)
        elif mode == 'intensity':
            a = lap(psi) / (1.0 + kappa * I_loc + 1e-9) - Vp - g * biharm(psi)
        elif mode == 'gradient':
            a = lap(psi) / (1.0 + kappa * grad2(psi) + 1e-9) - Vp - g * biharm(psi)
        elif mode == 'mixed':
            a = lap(psi) / (1.0 + kappa * (I_loc + grad2(psi)) + 1e-9) - Vp - g * biharm(psi)
        v = psi - psi_old
        psi_old, psi = psi, psi + (1.0 - damping * dt) * v + (dt**2) * a
        if t >= burn:
            psi_sum += psi
    th = np.angle(psi_sum / window)
    def w(a_, b_):
        d = b_ - a_
        return (d + np.pi) % (2 * np.pi) - np.pi
    d1 = w(th, np.roll(th, -1, 0))
    d2 = w(np.roll(th, -1, 0), np.roll(np.roll(th, -1, 0), -1, 1))
    d3 = w(np.roll(np.roll(th, -1, 0), -1, 1), np.roll(th, -1, 1))
    d4 = w(np.roll(th, -1, 1), th)
    wind = np.round((d1 + d2 + d3 + d4) / (2 * np.pi)).astype(int)
    return int((wind != 0).sum())

if __name__ == '__main__':
    out = {}
    for mode in ['linear', 'intensity', 'gradient', 'mixed']:
        for dt in [0.06, 0.03]:
            n = run_persist(mode, dt=dt)
            out[f"{mode}_dt{dt}"] = n
            print(f"{mode:10s} dt={dt}: persistent nvort = {n}")
    ratio = out['gradient_dt0.03'] / max(out['intensity_dt0.03'], 1)
    out['gradient_over_intensity_dt0.03'] = round(ratio, 2)
    print(f"\ngradient / intensity (dt=0.03, persistent): {ratio:.2f}")
    base = os.path.join(os.path.dirname(__file__), '..')
    with open(os.path.join(base, 'results', 'exp2d_persistence.json'), 'w') as f:
        json.dump(out, f, indent=2)
