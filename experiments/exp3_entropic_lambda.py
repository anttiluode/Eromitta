"""
EXP3 -- IS THE COSMOLOGICAL CONSTANT THE ENTROPY COST OF THE FOSSIL?
====================================================================
Bianconi's emergent cosmological constant (her Eq. 64) is

    Lambda_G ~ Tr( Gtilde - I - ln Gtilde )        >= 0,  = 0 iff Gtilde = I

-- the Bregman residue x - 1 - ln x of the dressing field: the entropy cost
of the G-field's failure to relax to the identity. Rajapinta exp3 showed
the crystal forgets everything EXCEPT its topological charge. ArrowField
showed defect count obeys Kibble-Zurek: N ~ tau_Q^(-sigma). Chaining the
three: the part of Lambda_G that cannot relax should be the part pinned to
defects, so the residual vacuum energy should be SET BY THE QUENCH RATE
on the same power law as the defect census.

Protocol. Complex field, intensity self-slowing (the Clockfield arm,
kappa=5), quenched through its U(1) transition:

    psi_tt = lap(psi)/(1+5|psi|^2) - eps(t) psi - g|psi|^2 psi
             - gamma psi_t + noise
    eps(t): +eps0 -> -eps0 linearly over 2*tau_Q, then HELD at -eps0 for a
    fixed relax time T_relax = 40 (same absolute coarsening time for all
    runs; deviation from ArrowField's threshold-crossing protocol, chosen
    so final states are comparable -- noted, not hidden).

Local dressing (flat space, scalar sector, small-coupling regime a=0.1):
    u(x)      = |psi|^2 + |grad psi|^2      (matter density seen by the
                                             induced metric)
    s(x)      = 1/(1 + a u)                 (G-field eigenvalue, Eq. 60)
    lam(x)    = s - 1 - ln s                (pointwise Lambda density)
    lam_bg    = lam at the spatial MEDIAN of u   (uniform-condensate floor)
    Lam_topo  = sum |lam - lam_bg|          (topological excess)

=== REGISTERED PREDICTIONS (before any run; 4 tau_Q x 3 seeds) ===

R1  KZ REPLICATES. Final snapshot defect count N(tau_Q) fits a power law
    N ~ tau_Q^(-sigma_N) with sigma_N in [0.3, 1.0], r^2 >= 0.8.
    (Sanity anchor to ArrowField Addendum 3.)

R2  LAMBDA COUNTS DEFECTS. Lam_topo(tau_Q) fits a power law with exponent
    sigma_L matching sigma_N within +/- 0.15.
    KILL: |sigma_L - sigma_N| > 0.3 or r^2(Lam) < 0.8 -- then the entropic
    residue is NOT a defect census and the chain to Bianconi's Lambda breaks.

R3  FIXED CHARGE PER DEFECT. Lam_topo / N constant across tau_Q within a
    factor of 2 (each fossil carries a fixed entropy cost).

R4  CONCENTRATION. Fraction of Lam_topo mass within r<=3 px of defect
    cores exceeds the area fraction of those disks by > 2x
    (the residue lives at the thaw line, not in the bulk).

Do not hype. Do not lie. Just show.
"""
import numpy as np, json, os
from scipy.signal import convolve2d

LAP = np.array([[0, 1, 0], [1, -4, 1], [0, 1, 0]], dtype=np.float64)
def lap(f):
    return (convolve2d(f.real, LAP, mode='same', boundary='wrap')
            + 1j * convolve2d(f.imag, LAP, mode='same', boundary='wrap'))

def grad2(f):
    gx = 0.5 * (np.roll(f, -1, 0) - np.roll(f, 1, 0))
    gy = 0.5 * (np.roll(f, -1, 1) - np.roll(f, 1, 1))
    return np.abs(gx)**2 + np.abs(gy)**2

def snapshot_wind(psi):
    th = np.angle(psi)
    def w(a_, b_):
        d = b_ - a_
        return (d + np.pi) % (2 * np.pi) - np.pi
    d1 = w(th, np.roll(th, -1, 0))
    d2 = w(np.roll(th, -1, 0), np.roll(np.roll(th, -1, 0), -1, 1))
    d3 = w(np.roll(np.roll(th, -1, 0), -1, 1), np.roll(th, -1, 1))
    d4 = w(np.roll(th, -1, 1), th)
    return np.round((d1 + d2 + d3 + d4) / (2 * np.pi)).astype(int)

def quench(tau_Q, seed, N=128, dt=0.08, eps0=1.0, gc=0.25, gamma=0.3,
           kappa=5.0, a_dress=0.1, T_relax=40.0, noise=2e-3):
    rng = np.random.default_rng(seed)
    psi = 1e-2 * (rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N)))
    psi_old = psi.copy()
    n_ramp = int(round(2 * tau_Q / dt)); n_relax = int(round(T_relax / dt))
    for t in range(n_ramp + n_relax):
        eps = eps0 - 2 * eps0 * min(t * dt, 2 * tau_Q) / (2 * tau_Q)
        I = np.abs(psi)**2
        a = (lap(psi) / (1 + kappa * I) - eps * psi - gc * I * psi)
        v = (1.0 - gamma * dt) * (psi - psi_old)
        kick = noise * np.sqrt(dt) * (rng.standard_normal((N, N))
                                      + 1j * rng.standard_normal((N, N)))
        psi_old, psi = psi, psi + v + dt**2 * a + kick
    wind = snapshot_wind(psi)
    n_def = int((wind != 0).sum())
    u = np.abs(psi)**2 + grad2(psi)
    s = 1.0 / (1.0 + a_dress * u)
    lam = s - 1 - np.log(s)
    lam_bg = float(np.interp(0, [0], [0])) * 0  # placeholder clarity
    s_bg = 1.0 / (1.0 + a_dress * np.median(u))
    lam_bg = s_bg - 1 - np.log(s_bg)
    dlam = np.abs(lam - lam_bg)
    Lam_topo = float(dlam.sum())
    # concentration
    if n_def:
        mask = convolve2d((wind != 0).astype(float), np.ones((7, 7)),
                          mode='same', boundary='wrap') > 0
        conc = float(dlam[mask].sum() / Lam_topo)
        area = float(mask.mean())
    else:
        conc, area = float('nan'), float('nan')
    return dict(n_def=n_def, Lam_topo=Lam_topo, conc=conc, area=area,
                lam=lam, wind=wind, I=np.abs(psi)**2)

def fit_power(x, y):
    lx, ly = np.log(x), np.log(np.maximum(y, 1e-12))
    A = np.vstack([lx, np.ones_like(lx)]).T
    coef, res, *_ = np.linalg.lstsq(A, ly, rcond=None)
    yhat = A @ coef
    r2 = 1 - ((ly - yhat)**2).sum() / ((ly - ly.mean())**2).sum()
    return float(-coef[0]), float(r2)

if __name__ == '__main__':
    taus = [8, 16, 32, 64]
    seeds = [1, 2, 3]
    table, sample = {}, None
    for tq in taus:
        row = [quench(tq, s) for s in seeds]
        if tq == taus[0]:
            sample = row[0]
        table[tq] = dict(
            n_def=[r['n_def'] for r in row],
            Lam=[r['Lam_topo'] for r in row],
            conc=[r['conc'] for r in row],
            area=[r['area'] for r in row])
        print(f"tau_Q={tq:3d}  n_def={table[tq]['n_def']}  "
              f"Lam_topo={[round(v,1) for v in table[tq]['Lam']]}  "
              f"conc={[round(c,2) for c in table[tq]['conc']]} "
              f"(area={[round(a,2) for a in table[tq]['area']]})")

    x = np.array(taus, float)
    Nm = np.array([np.mean(table[t]['n_def']) for t in taus])
    Lm = np.array([np.mean(table[t]['Lam']) for t in taus])
    sN, r2N = fit_power(x, Nm)
    sL, r2L = fit_power(x, Lm)
    charge = Lm / np.maximum(Nm, 1)
    concs = np.array([np.mean(table[t]['conc']) for t in taus])
    areas = np.array([np.mean(table[t]['area']) for t in taus])

    verdict = dict(
        sigma_N=round(sN, 3), r2_N=round(r2N, 3),
        sigma_L=round(sL, 3), r2_L=round(r2L, 3),
        R1=bool(0.3 <= sN <= 1.0 and r2N >= 0.8),
        R2=bool(abs(sL - sN) <= 0.15 and r2L >= 0.8),
        R2_KILL=bool(abs(sL - sN) > 0.3 or r2L < 0.8),
        R3_charge_ratio=round(float(charge.max() / charge.min()), 2),
        R3=bool(charge.max() / charge.min() <= 2.0),
        R4_conc_over_area=[round(float(c / a), 2) for c, a in zip(concs, areas)],
        R4=bool(np.all(concs / areas > 2.0)))
    print(json.dumps(verdict, indent=2))

    base = os.path.join(os.path.dirname(__file__), '..')
    with open(os.path.join(base, 'results', 'exp3_lambda.json'), 'w') as f:
        json.dump(dict(table={str(k): v for k, v in table.items()},
                       verdict=verdict), f, indent=2)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
    axes[0].imshow(sample['lam'], cmap='magma')
    yy, xx = np.nonzero(sample['wind'])
    axes[0].plot(xx, yy, 'c.', ms=2)
    axes[0].set_title(f"lam(x) = s-1-ln s, tau_Q={taus[0]}\ncyan = defects")
    axes[0].set_xticks([]); axes[0].set_yticks([])
    axes[1].loglog(x, Nm, 'o-', label=f"N_def, sigma={sN:.2f} (r2={r2N:.2f})")
    axes[1].loglog(x, Lm / Lm[0] * Nm[0], 's--',
                   label=f"Lam_topo (scaled), sigma={sL:.2f} (r2={r2L:.2f})")
    axes[1].set_xlabel("tau_Q"); axes[1].legend()
    axes[1].set_title("defect census vs entropic residue")
    axes[2].semilogx(x, charge, 'd-')
    axes[2].set_xlabel("tau_Q"); axes[2].set_title("Lam_topo / N_def (charge per fossil)")
    plt.tight_layout()
    plt.savefig(os.path.join(base, 'figs', 'fig3_entropic_lambda.png'), dpi=110)
    print("figure written")
