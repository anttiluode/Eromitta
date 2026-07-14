"""
EXP1 -- THE DERIVATIVE TOWER
============================
Claim to verify (the candidate resolution of the ArrowField "open seam"):

    The three saturating forms in the ecosystem are the SAME LOG at three
    derivative levels:

        level 0   K(x)      =  ln(1 + kx)          entropy / Kahler potential
                                                    (Bianconi's L, Kahler repo's K)
        level 1   K'(x)     =  k/(1 + kx)          power -1: Bianconi's h(w),
                                                    phiworld's c^2(|psi|^2)
        level 2   dd~K      =  k/(1 + kb)^2        power -2: Gamma, the
                                                    Fubini-Study conformal factor

    If this holds, phiworld's exponent (-1) and the Clockfield Gamma's
    exponent (-2) are not competing laws. They are the stiffness and the
    metric of one entropic potential.

=== REGISTERED CHECKS (pass/fail thresholds fixed before running) ===

T1  BIANCONI EIGENVALUE STRUCTURE. For a random Lorentzian metric g and a
    random complex gradient v_mu, the mixed tensor (g + a*M) g^{-1} with
    M_{mu nu} = vbar_mu v_nu must have eigenvalues {1,1,1, 1 + a*|v|^2_g},
    so that  -Tr ln(G g^{-1}) = -ln(1 + a*|v|^2_g)  exactly.
    PASS: max |eig residual| < 1e-9 over 200 random draws (spacelike |v|^2_g>0).
    Also RECORD (not pass/fail): fraction of unrestricted draws where
    1 + a*|v|^2_g <= 0 (timelike gradients) -- the regime where Bianconi's
    positivity requirement bites. Honesty item, not a check.

T2  KAHLER HESSIAN IDENTITY. K(phi,phibar) = ln(1 + t*phi*phibar).
    Wirtinger finite differences of d^2K/dphi dphibar at 200 random complex
    points must equal t/(1 + t*b)^2 with b=|phi|^2.
    PASS: max relative error < 1e-4 (h=1e-3 central differences).

T3  FIRST-DERIVATIVE = SCREENING. dK/dw at w = |grad psi|^2 equals
    Bianconi's h(w) = a/(1+aw) (her Eq. 22) and, with w -> |psi|^2,
    phiworld's c^2/c0^2 = 1/(1+beta*|psi|^2) up to the constant a vs beta.
    PASS: max abs error < 1e-12 (this is symbolic-trivial; recorded so the
    tower is complete in one runnable file, not scattered across READMEs).

Do not hype. Do not lie. Just show.
"""
import numpy as np, json, os

rng = np.random.default_rng(7)
out = {}

# ---------------------------------------------------------------- T1
eta = np.diag([-1.0, 1.0, 1.0, 1.0])
alpha = 0.37
max_resid = 0.0
timelike = 0
total_unrestricted = 0
draws = 0
while draws < 200:
    Lm = rng.standard_normal((4, 4))
    if abs(np.linalg.det(Lm)) < 0.1:
        continue
    g = Lm.T @ eta @ Lm                      # random Lorentzian metric
    ginv = np.linalg.inv(g)
    v = rng.standard_normal(4) + 1j * rng.standard_normal(4)   # nabla psi (lower)
    w = np.real(np.conj(v) @ ginv @ v)       # |grad psi|^2_g  (real: g symmetric)
    total_unrestricted += 1
    if 1 + alpha * w <= 0:
        timelike += 1
        continue
    M = np.outer(np.conj(v), v)              # M_{mu nu} = vbar_mu v_nu
    G = g + alpha * M
    N = G @ ginv                             # mixed tensor N_mu^nu
    eig = np.sort_complex(np.linalg.eigvals(N))
    target = np.sort_complex(np.array([1, 1, 1, 1 + alpha * w], dtype=complex))
    max_resid = max(max_resid, float(np.abs(eig - target).max()))
    # and the Lagrangian itself:
    Ltrace = -np.log(np.linalg.det(N)).real  # -Tr ln = -ln det
    Ldirect = -np.log(1 + alpha * w)
    max_resid = max(max_resid, abs(Ltrace - Ldirect))
    draws += 1

out['T1'] = dict(max_residual=max_resid,
                 passed=bool(max_resid < 1e-9),
                 timelike_fraction_unrestricted=timelike / total_unrestricted)

# ---------------------------------------------------------------- T2
tau = 1.3
h = 1e-3   # balances truncation O(h^2) vs roundoff O(eps/h^2) for 2nd differences
max_rel = 0.0
for _ in range(200):
    phi = (rng.standard_normal() + 1j * rng.standard_normal())
    K = lambda p: np.log(1 + tau * (p * np.conj(p)).real)
    # Wirtinger: d2K/dphi dphibar = 1/4 (Kxx + Kyy)
    x, y = phi.real, phi.imag
    Kxx = (K(complex(x + h, y)) - 2 * K(phi) + K(complex(x - h, y))) / h**2
    Kyy = (K(complex(x, y + h)) - 2 * K(phi) + K(complex(x, y - h))) / h**2
    hess = 0.25 * (Kxx + Kyy)
    b = (phi * np.conj(phi)).real
    exact = tau / (1 + tau * b)**2           # Gamma * tau  (Fubini-Study)
    max_rel = max(max_rel, abs(hess - exact) / abs(exact))
out['T2'] = dict(max_relative_error=max_rel, passed=bool(max_rel < 1e-4))

# ---------------------------------------------------------------- T3
ws = np.abs(rng.standard_normal(200))
a = 0.37
dK = a / (1 + a * ws)                         # d/dw ln(1+aw)
h_bianconi = a / (1 + a * ws)                 # her Eq. 22
err = float(np.abs(dK - h_bianconi).max())
out['T3'] = dict(max_error=err, passed=bool(err < 1e-12))

os.makedirs(os.path.join(os.path.dirname(__file__), '..', 'results'), exist_ok=True)
with open(os.path.join(os.path.dirname(__file__), '..', 'results', 'exp1_tower.json'), 'w') as f:
    json.dump(out, f, indent=2)

print(json.dumps(out, indent=2))
print()
print("THE TOWER:  ln(1+kx)  --d/dx-->  k/(1+kx)  --dd~ on |phi|^2-->  k/(1+kb)^2")
print("            entropy/Kahler K     h(w), c^2(psi)                 Gamma, g_FS")
