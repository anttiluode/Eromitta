"""
EXP2 -- FOUR SCREENINGS: what must the medium SENSE to make matter?
===================================================================
ArrowField's headline: intensity self-slowing c^2 = c0^2/(1+beta|psi|^2)
generates topological defects (402 vs 0). Bianconi's WARM-UP entropic
action gives a different screening: h(w) = a/(1+a*w) with w = |GRAD psi|^2
-- the medium slows on gradients, not intensity. Her FULL theory adds the
term (m^2 + xi*R)|Phi><Phi| to the induced metric -- an INTENSITY term.

So Bianconi's own two-step structure makes a prediction about phiworld:
her warm-up (gradient-only) should be sterile, and her full induced metric
(gradient + mass*intensity) should make matter -- because phiworld's beta
term IS the mass term of the entropic metric.

Protocol: IDENTICAL to ArrowField complex_phiworld.py (same init, same
seed, same steps, same counters), only the screening argument varies.

    A  linear      c^2 = 1                                   (control)
    B  intensity   c^2 = 1/(1 + 5|psi|^2)                    (phiworld)
    C  gradient    c^2 = 1/(1 + 5|grad psi|^2)               (Bianconi warm-up,
                                                              c^2-form)
    Cf gradient    psi_tt = div( grad psi / (1+5|grad psi|^2) )  (true
                   flux/Euler-Lagrange form of her Eq. 21)
    D  mixed       c^2 = 1/(1 + 5(|psi|^2 + |grad psi|^2))   (full scalar
                                                              sector: mass +
                                                              gradient)
    C4 gradient    same as C at DOUBLE drive amplitude       (amplitude
                                                              control: can
                                                              energy substitute?)

=== REGISTERED PREDICTIONS (before any run) ===

P1  GRADIENT SCREENING IS STERILE. Kerr-type self-focusing needs the speed
    to fall with INTENSITY: bright regions become waveguides that trap their
    own light. A smooth bright condensate is invisible to gradient screening,
    so it has no trapping mechanism -- fronts steepen, nothing condenses.
    Predict: nvort(C) < 0.10 * nvort(B), and same for Cf.
    KILL: nvort(C or Cf) >= 0.5 * nvort(B) -- then gradient screening IS in
    the defect class and my mechanism story is wrong.

P2  THE MASS TERM RESCUES IT. Predict nvort(D) >= 0.33 * nvort(B): adding
    the intensity component to the entropic metric restores defect
    generation; the gradient component does not destroy it.
    KILL: nvort(D) < 0.10 * nvort(B) -- then "phiworld is the mass-dominated
    limit of Bianconi's induced metric" fails DYNAMICALLY even though it
    holds at the level of the formula, and we say so.

P3  AMPLITUDE CANNOT SUBSTITUTE (again). Predict nvort(C4, amp x2) stays
    < 0.10 * nvort(B), mirroring ArrowField's beta=0 amplitude control.

P4  STRUCTURE. Predict field sigma: sigma(C) < sigma(B). Gradient screening
    slows transport where fronts are steep but cannot trap amplitude, so it
    should fail to hold structure the way beta=0 fails (sigma ~ 0.36).
    Uncertain; registered anyway.

Interpretation stakes: if P1+P2 hold, Bianconi's OWN diagnosis of her
warm-up ("first limitation: no mass term") is exactly phiworld's diagnosis
of beta=0 ("no self-slowing, no things") -- the two limitations are the
same limitation, seen from entropy and from simulation.

Do not hype. Do not lie. Just show.
"""
import numpy as np, json, os, time
from scipy.signal import convolve2d

LAP = np.array([[0, 1, 0], [1, -4, 1], [0, 1, 0]], dtype=np.float64)

def lap(f):
    return (convolve2d(f.real, LAP, mode='same', boundary='wrap')
            + 1j * convolve2d(f.imag, LAP, mode='same', boundary='wrap'))

def biharm(f):
    return lap(lap(f))

def grads(f):
    gx = 0.5 * (np.roll(f, -1, 0) - np.roll(f, 1, 0))
    gy = 0.5 * (np.roll(f, -1, 1) - np.roll(f, 1, 1))
    return gx, gy

def grad2(f):
    gx, gy = grads(f)
    return np.abs(gx)**2 + np.abs(gy)**2

def div(hx, hy):
    return (0.5 * (np.roll(hx, -1, 0) - np.roll(hx, 1, 0))
            + 0.5 * (np.roll(hy, -1, 1) - np.roll(hy, 1, 1)))

def count_vortices(buf):
    th = np.angle(buf.mean(axis=0))
    def w(a, b):
        d = b - a
        return (d + np.pi) % (2 * np.pi) - np.pi
    d1 = w(th, np.roll(th, -1, 0))
    d2 = w(np.roll(th, -1, 0), np.roll(np.roll(th, -1, 0), -1, 1))
    d3 = w(np.roll(np.roll(th, -1, 0), -1, 1), np.roll(th, -1, 1))
    d4 = w(np.roll(th, -1, 1), th)
    wind = np.round((d1 + d2 + d3 + d4) / (2 * np.pi)).astype(int)
    return int((wind != 0).sum()), wind

def run(mode, kappa=5.0, amp=2.0, N=128, dt=0.06, damping=0.001, pcub=0.2,
        g=0.02, burn=1600, window=700, seed=1):
    """Identical protocol to ArrowField complex_phiworld.py run()."""
    rng = np.random.default_rng(seed)
    x = np.arange(N); X, Y = np.meshgrid(x, x, indexing='ij')
    c = N // 2; r = N / 15.0
    env = amp * np.exp(-((X - c)**2 + (Y - c)**2) / (2 * r**2))
    ph = 0.9 * rng.standard_normal((N, N))
    ph = convolve2d(ph, np.ones((5, 5)) / 25., mode='same', boundary='wrap')
    psi = env * np.exp(1j * 2.5 * ph)
    psi_old = psi.copy()
    buf = None
    for t in range(burn + window):
        I_loc = np.abs(psi)**2
        Vp = -psi + pcub * I_loc * psi
        if mode == 'linear':
            a = lap(psi) - Vp - g * biharm(psi)
        elif mode == 'intensity':
            c2 = 1.0 / (1.0 + kappa * I_loc + 1e-9)
            a = c2 * lap(psi) - Vp - g * biharm(psi)
        elif mode == 'gradient':
            c2 = 1.0 / (1.0 + kappa * grad2(psi) + 1e-9)
            a = c2 * lap(psi) - Vp - g * biharm(psi)
        elif mode == 'gradient_flux':
            gx, gy = grads(psi)
            h = 1.0 / (1.0 + kappa * (np.abs(gx)**2 + np.abs(gy)**2) + 1e-9)
            a = div(h * gx, h * gy) - Vp - g * biharm(psi)
        elif mode == 'mixed':
            c2 = 1.0 / (1.0 + kappa * (I_loc + grad2(psi)) + 1e-9)
            a = c2 * lap(psi) - Vp - g * biharm(psi)
        else:
            raise ValueError(mode)
        v = psi - psi_old
        new = psi + (1.0 - damping * dt) * v + (dt**2) * a
        psi_old, psi = psi, new
        if t == burn:
            buf = np.zeros((window, N, N), dtype=complex)
        if t >= burn:
            buf[t - burn] = psi
    nvort, wind = count_vortices(buf)
    return dict(nvort=nvort, sigma=float(np.abs(buf).std()),
                I=(np.abs(buf)**2).mean(axis=0),
                th=np.angle(buf.mean(axis=0)), wind=wind)

if __name__ == '__main__':
    arms = [('A_linear',        dict(mode='linear')),
            ('B_intensity',     dict(mode='intensity')),
            ('C_gradient',      dict(mode='gradient')),
            ('Cf_gradient_flux',dict(mode='gradient_flux')),
            ('D_mixed',         dict(mode='mixed')),
            ('C4_gradient_amp4',dict(mode='gradient', amp=4.0))]
    results, maps = {}, {}
    for name, kw in arms:
        t0 = time.time()
        r = run(**kw)
        results[name] = dict(nvort=r['nvort'], sigma=r['sigma'],
                             seconds=round(time.time() - t0, 1))
        maps[name] = r
        print(f"{name:20s}  nvort={r['nvort']:5d}  sigma={r['sigma']:.3f}"
              f"  ({results[name]['seconds']}s)")

    nB = max(results['B_intensity']['nvort'], 1)
    verdict = dict(
        P1_gradient_sterile=bool(results['C_gradient']['nvort'] < 0.10 * nB
                                 and results['Cf_gradient_flux']['nvort'] < 0.10 * nB),
        P1_KILL=bool(results['C_gradient']['nvort'] >= 0.5 * nB
                     or results['Cf_gradient_flux']['nvort'] >= 0.5 * nB),
        P2_mass_rescues=bool(results['D_mixed']['nvort'] >= 0.33 * nB),
        P2_KILL=bool(results['D_mixed']['nvort'] < 0.10 * nB),
        P3_amplitude_no_substitute=bool(results['C4_gradient_amp4']['nvort'] < 0.10 * nB),
        P4_sigma=bool(results['C_gradient']['sigma'] < results['B_intensity']['sigma']))
    results['verdicts'] = verdict
    print(json.dumps(verdict, indent=2))

    base = os.path.join(os.path.dirname(__file__), '..')
    with open(os.path.join(base, 'results', 'exp2_screenings.json'), 'w') as f:
        json.dump(results, f, indent=2)

    # figure
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    show = ['B_intensity', 'C_gradient', 'D_mixed']
    fig, axes = plt.subplots(2, 3, figsize=(12, 7.5))
    for j, name in enumerate(show):
        m = maps[name]
        axes[0, j].imshow(m['I'], cmap='inferno')
        axes[0, j].set_title(f"{name}  |psi|^2   nvort={m['nvort']}")
        axes[1, j].imshow(m['th'], cmap='twilight')
        yy, xx = np.nonzero(m['wind'])
        axes[1, j].plot(xx, yy, 'w.', ms=1.5)
        axes[1, j].set_title("phase + defects")
        for ax in (axes[0, j], axes[1, j]):
            ax.set_xticks([]); ax.set_yticks([])
    plt.suptitle("EXP2 -- what the medium senses decides whether matter forms")
    plt.tight_layout()
    plt.savefig(os.path.join(base, 'figs', 'fig1_four_screenings.png'), dpi=110)
    print("figure written")
