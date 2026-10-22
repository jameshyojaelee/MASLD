"""Decimal Gauss-Legendre moments for the fixed synthetic mixed derivatives.

No candidate derivatives or numerical integration routines are used. Raw
moments must be summed across the complete domain before forming covariance.
Convergence checks do not certify quadrature error.
"""
from decimal import Decimal
import math


RULES = ((24, '0.25', 16), (26, '0.125', 16), (26, '0.125', 32))
CHUNKS = 16


def require(ok, message):
    if not ok:
        raise ValueError(message)


def legendre(order, x):
    previous, current = Decimal(1), x
    for k in range(2, order+1):
        previous, current = current, ((2*k-1)*x*current-(k-1)*previous)/k
    derivative = order*(x*current-previous)/(x*x-1)
    return current, derivative


def rule(order, digits):
    """Refine float root seeds in the active Decimal context, then check moments."""
    require(order in (16,32), 'Undeclared GL order')
    tolerance = Decimal(10)**(-digits+10)
    positive = []
    for j in range(1,order//2+1):
        x = Decimal.from_float(math.cos(math.pi*(j-0.25)/(order+0.5)))
        for iteration in range(30):
            value, derivative = legendre(order,x)
            change = value/derivative
            x -= change
            if abs(change) <= tolerance:
                break
        else:
            raise ValueError('Decimal GL root did not converge')
        value, derivative = legendre(order,x)
        require(abs(value) <= tolerance*100 and 0 < x < 1, 'Decimal GL residual failed')
        weight = 2/((1-x*x)*derivative*derivative)
        require(weight > 0, 'Nonpositive GL weight')
        positive.append((x,weight))
    nodes = sorted([(-x,w) for x,w in positive]+positive)
    require(len({x for x,w in nodes}) == order, 'Duplicate GL root')
    checks = []
    for degree in range(2*order):
        found = sum((w*x**degree for x,w in nodes),Decimal(0))
        expected = Decimal(0) if degree%2 else Decimal(2)/Decimal(degree+1)
        error = abs(found-expected)
        require(error <= tolerance*1000, 'GL polynomial moment failed')
        checks.append(str(error))
    return nodes,checks


def panels(rule_index, chunk):
    require(rule_index in range(len(RULES)) and chunk in range(CHUNKS),
            'Undeclared rule or chunk')
    bound,spacing,order = RULES[rule_index]
    step = Decimal(spacing)
    count = int(Decimal(2*bound)/step)
    require(Decimal(count)*step == 2*bound, 'Noninteger panel family')
    begin,end = count*chunk//CHUNKS,count*(chunk+1)//CHUNKS
    return bound,step,order,count,begin,end


def integrate(prepared, variance, shift, rule_index, chunk, digits, conditional_moments,
              progress=None):
    bound,spacing,order,count,begin,end = panels(rule_index,chunk)
    nodes,checks = rule(order,digits)
    rootv = variance.sqrt()
    half = spacing/2
    moments = [Decimal(0) for _ in range(6)]
    calls = 0
    for panel in range(begin,end):
        center = -Decimal(bound)+(Decimal(panel)+Decimal('0.5'))*spacing
        for x,w in nodes:
            z = center+half*x
            ell,sr,sw,hrw = conditional_moments(prepared,rootv*z)
            gv = (z*z-1)/(2*variance)
            weight = half*w*(ell-shift-z*z/2).exp()
            terms = (Decimal(1),sr,sw,gv,hrw+sr*sw,sr*gv)
            for j,value in enumerate(terms):
                moments[j] += weight*value
            calls += 1
        if progress is not None:
            progress(panel,calls)
    require(moments[0] > 0 and all(x.is_finite() for x in moments), 'Invalid raw moments')
    return dict(digits=digits,rule_index=rule_index,chunk=chunk,bound=bound,
        spacing=str(spacing),order=order,total_panels=count,first_panel=begin,
        last_panel_exclusive=end,node_calls=calls,moments=list(map(str,moments)),
        moment_order=['mass','r_score','omega_score','variance_score',
                      'r_omega_curvature_plus_score_product','r_score_times_variance_score'],
        GL_polynomial_absolute_errors=checks,no_chunk_covariance=True)


def whole_line_tail_bounds(rows, variance, shift, bound):
    """Analytic absolute raw-moment bounds, separate from quadrature error.

For this regularized logistic rho, |t_r/t| <= 1. A raw BB eta score is
bounded by n, r score by 3n, and eta/r curvature by 2n. Mixture log
curvature adds a covariance bounded by 6n^2; selection subtracts another
such log normalizer. These deliberately loose bounds hold for every u.
The selected row probability is at most one. Gaussian Mills bounds then
bound the discarded moments. Decimal evaluation is not interval arithmetic.
"""
    zero = Decimal(0)
    b = sum((Decimal(6*int(row.n)) for row in rows),zero)
    w = sum((2*Decimal(int(row.n))*abs(Decimal.from_float(float(row.stage)))
             for row in rows),zero)
    c = sum((Decimal(4*int(row.n)+12*int(row.n)**2)*
             abs(Decimal.from_float(float(row.stage))) for row in rows),zero)
    z = Decimal(bound)
    envelope = (-shift-z*z/2).exp()
    gaussian_mass = 2*envelope/z
    gaussian_square_plus_one = 2*envelope*(z+2/z)
    return list(map(str,(gaussian_mass,b*gaussian_mass,w*gaussian_mass,
        gaussian_square_plus_one/(2*variance),(c+b*w)*gaussian_mass,
        b*gaussian_square_plus_one/(2*variance))))
