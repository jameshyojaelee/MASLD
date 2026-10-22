"""Independent Decimal native BB dispersion/eta mixed derivatives.

Closed rising-product derivatives replace float mode ratios and complex
perturbations in this reference. No candidate derivative helper is used.
Gaussian normalization cancels in the posterior derivative moments.
"""
from decimal import Decimal
import math


def exact(value):
    return value if isinstance(value, Decimal) else Decimal.from_float(float(value))


def require(ok, message):
    if not ok:
        raise ValueError(message)


def prefixes(shape, n, derivatives=False):
    zero, one = Decimal(0), Decimal(1)
    logs, reciprocal, squared = [zero], [zero], [zero]
    for j in range(n):
        term = shape+Decimal(j)
        logs.append(logs[-1]+term.ln())
        if derivatives:
            inverse = one/term
            reciprocal.append(reciprocal[-1]+inverse)
            squared.append(squared[-1]+inverse*inverse)
    return logs, reciprocal, squared


def prepare(rows, p, G, Z, bb_vector):
    zero, one, two = Decimal(0), Decimal(1), Decimal(2)
    require(rows and len(p) == 2*G+Z+5, 'Native synthetic dimensions differ')
    q = [exact(x) for x in p]
    prepared = []
    for row in rows:
        n = int(row.n)
        lo = max(int(row.minor), (int(row.percent)*n+99)//100)
        support = tuple(range(lo, n-lo+1))
        require(0 < n <= 5000 and row.n == n and row.n >= row.depth and row.a in support
                and len(row.z) == Z and row.orientation in (-1, 1), 'Native row/support differs')
        stage = exact(row.stage)
        eta0 = q[row.gene]*(one+q[2*G]*stage+sum(
            (q[2*G+1+j]*exact(row.z[j]) for j in range(Z)), zero))
        eta0 -= Decimal(row.orientation)*(q[2*G+1+Z]*stage+exact(row.omega))
        rr = q[G+row.gene]+sum((q[-3+j]*exact(x)
            for j, x in enumerate((row.stage, row.background, row.duplication))), zero)
        pref, palt, px = (exact(x) for x in (row.p_ref, row.p_alt, row.p_x))
        prior = [two*pref*palt, pref*pref+two*pref*px, palt*palt+two*palt*px]
        if row.orientation == -1:
            prior = [prior[0], prior[2], prior[1]]
        require(all(x >= zero for x in prior) and sum(prior, zero) > zero, 'Invalid genotype prior')
        choose = [Decimal(math.comb(n, k)).ln() for k in range(n+1)]
        e, phi = exact(row.e), exact(row.phi)
        require(zero < e < one and zero <= phi < one, 'Invalid homozygote components')
        def homo(mean):
            if phi == zero:
                return [Decimal(math.comb(n, k))*mean**k*(one-mean)**(n-k) for k in support]
            t = (one-phi)/phi
            return [x.exp() for x in bb_vector(support, n, mean*t, (one-mean)*t, choose)]
        homref, homalt = homo(e), homo(one-e)
        hom = [prior[1]*a+prior[2]*b for a, b in zip(homref, homalt)]
        t, tr = shapes(eta0, rr)[2:4]
        t_prefixes = prefixes(t, n, True)
        prepared.append(dict(n=n, observed_index=int(row.a)-lo, support=support, eta0=eta0,
            stage=stage, omega_slope=-Decimal(row.orientation)*stage, r=rr,
            hetero_prior=prior[0], hom=hom, logchoose=choose, row=row,
            t=t, tr=tr, t_prefixes=t_prefixes))
    require(sum(d['n']+1 for d in prepared) <= 100000, 'Synthetic gene support cap exceeded')
    return prepared


def shapes(eta, r):
    one = Decimal(1)
    sigmoid = one/(one+(-r).exp())
    complement = one/(one+r.exp())
    rho = exact(1e-6)+exact(1-2e-6)*sigmoid
    rho_r = exact(1-2e-6)*sigmoid*complement
    t = (one-rho)/rho
    tr = -rho_r/(rho*rho)
    m, q = one/(one+(-eta).exp()), one/(one+eta.exp())
    return t*m, t*q, t, tr, m, q


def row_logmass(d, eta, r, bb_vector):
    """Probability-only finite-difference comparator; no derivative prefixes."""
    a, b = shapes(eta, r)[:2]
    logs = bb_vector(d['support'], d['n'], a, b, d['logchoose'])
    masses = [d['hetero_prior']*lp.exp()+hom for lp, hom in zip(logs, d['hom'])]
    return masses[d['observed_index']].ln()-sum(masses, Decimal(0)).ln()


def row_terms(d, eta, r):
    """Selected logmass, eta score, r score, eta/r mixed log curvature."""
    zero = Decimal(0)
    a, b, t, tr, m, q = shapes(eta, r)
    la, ia, ja = prefixes(a, d['n'], True)
    lb, ib, jb = prefixes(b, d['n'], True)
    lt, it, _ = d['t_prefixes'] if r == d['r'] else prefixes(t, d['n'], True)
    ae, be = t*m*q, -t*m*q
    ar, br = tr*m, tr*q
    aer, ber = tr*m*q, -tr*m*q
    entries = []
    for k, hom in zip(d['support'], d['hom']):
        nk = d['n']-k
        lp = d['logchoose'][k]+la[k]+lb[nk]-lt[-1]
        se = ae*ia[k]+be*ib[nk]
        sr = ar*ia[k]+br*ib[nk]-tr*it[-1]
        her = aer*ia[k]-ae*ar*ja[k]+ber*ib[nk]-be*br*jb[nk]
        het = d['hetero_prior']*lp.exp()
        entries.append((het+hom, het*se, het*sr, het*(her+se*sr)))
    den = [sum((entry[j] for entry in entries), zero) for j in range(4)]
    obs = entries[d['observed_index']]
    def log_derivatives(values):
        mass, de, dr, der = values
        require(mass.is_finite() and mass > 0, 'Invalid native selected mixture mass')
        return mass.ln(), de/mass, dr/mass, der/mass-de*dr/(mass*mass)
    numerator, denominator = log_derivatives(obs), log_derivatives(den)
    answer = tuple(a-b for a, b in zip(numerator, denominator))
    require(all(x.is_finite() for x in answer), 'Nonfinite Decimal mixed derivative')
    return answer


def conditional_moments(prepared, u):
    """Sum row log likelihood and r_g/omega_S scores and mixed curvature."""
    ell = score_r = score_omega = mixed = Decimal(0)
    for d in prepared:
        value, se, sr, her = row_terms(d, d['eta0']+u*d['stage'], d['r'])
        ell += value
        score_r += sr
        score_omega += d['omega_slope']*se
        mixed += d['omega_slope']*her
    return ell, score_r, score_omega, mixed
