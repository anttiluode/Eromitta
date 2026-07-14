"""
EXP2c -- THE COMPARISON, RE-RUN ON REPAIRED INSTRUMENTS (post-hoc)
==================================================================
The autopsy (exp2b) found two instrument faults, one inherited, one mine:

  FAULT 1 (inherited from ArrowField's counter): vortices were counted on
  the phase of the TIME-AVERAGED field. Moving defects smear; core-amplitude
  checks against time-averaged intensity are meaningless. FIX: count on the
  FINAL SNAPSHOT; check core amplitude on the same snapshot.

  FAULT 2 (mine): the dt-robustness check halved dt without doubling steps,
  so it compared different physical durations. FIX: fix physical time
  T_burn = 96, T_window = 42; steps = T/dt.

Re-measured quantities, per arm {linear, intensity, gradient, mixed} x
dt {0.06, 0.03}:
    nvort_snap   snapshot vortex count
    core_ratio   |psi|^2 at cores / field mean (real defects: << 1)
    sigma        field structure

Post-hoc status: the EXP2 registration is spent. Whatever comes out here
is reported as the repaired measurement, not as a confirmed prediction.

Do not hype. Do not lie. Just show.
"""
import numpy as np, json, os, sys
from scipy.signal import convolve2d
sys.path.insert(0, os.path.dirname(__file__))
from exp2_four_screenings import lap, biharm, grads, grad2

def run_snap(mode, kappa=5.0, amp=2.0, N=128, dt=0.06, damping=0.001,
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
    sig_acc, nacc = 0.0, 0
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
            sig_acc += np.abs(psi).std(); nacc += 1
    th = np.angle(psi)
    def w(a_, b_):
        d = b_ - a_
        return (d + np.pi) % (2 * np.pi) - np.pi
    d1 = w(th, np.roll(th, -1, 0))
    d2 = w(np.roll(th, -1, 0), np.roll(np.roll(th, -1, 0), -1, 1))
    d3 = w(np.roll(np.roll(th, -1, 0), -1, 1), np.roll(th, -1, 1))
    d4 = w(np.roll(th, -1, 1), th)
    wind = np.round((d1 + d2 + d3 + d4) / (2 * np.pi)).astype(int)
    mask = wind != 0
    I = np.abs(psi)**2
    core = float(I[mask].mean() / I.mean()) if mask.any() else float('nan')
    return dict(nvort=int(mask.sum()), core_ratio=core,
                sigma=float(sig_acc / max(nacc, 1)),
                I=I, th=th, wind=wind)

if __name__ == '__main__':
    out, maps = {}, {}
    for mode in ['linear', 'intensity', 'gradient', 'mixed']:
        for dt in [0.06, 0.03]:
            r = run_snap(mode, dt=dt)
            key = f"{mode}_dt{dt}"
            out[key] = dict(nvort=r['nvort'], core_ratio=round(r['core_ratio'], 3)
                            if r['nvort'] else None, sigma=round(r['sigma'], 3))
            maps[key] = r
            print(f"{key:20s} nvort={r['nvort']:5d}  core|psi|^2/mean="
                  f"{out[key]['core_ratio']}  sigma={out[key]['sigma']}")
    base = os.path.join(os.path.dirname(__file__), '..')
    with open(os.path.join(base, 'results', 'exp2c_repaired.json'), 'w') as f:
        json.dump(out, f, indent=2)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    show = ['intensity_dt0.03', 'gradient_dt0.03', 'mixed_dt0.03']
    fig, axes = plt.subplots(2, 3, figsize=(12, 7.5))
    for j, name in enumerate(show):
        m = maps[name]
        axes[0, j].imshow(m['I'], cmap='inferno')
        axes[0, j].set_title(f"{name}\n|psi|^2 (snapshot)  nvort={m['nvort']}")
        axes[1, j].imshow(m['th'], cmap='twilight')
        yy, xx = np.nonzero(m['wind'])
        axes[1, j].plot(xx, yy, 'w.', ms=1.5)
        axes[1, j].set_title(f"phase + defects  core ratio={m['core_ratio']:.2f}")
        for ax in (axes[0, j], axes[1, j]):
            ax.set_xticks([]); ax.set_yticks([])
    plt.suptitle("EXP2c -- repaired instruments: snapshot counts, dt=0.03, matched physical time")
    plt.tight_layout()
    plt.savefig(os.path.join(base, 'figs', 'fig2_repaired.png'), dpi=110)
    print("figure written")
