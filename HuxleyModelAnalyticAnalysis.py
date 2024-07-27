#%%
from sympy import *
# sympy plotting backends
from spb import *


#%%
P1, P2 = symbols('P_1, P_2', real = True)
V1, V2 = symbols('V_1, V_2', real = True)
f = symbols('f', integer = True, positive = True)

# A)
Wa = 0
dUa = f * (P2-P1) * V1
Qa = dUa - Wa

# B)
Wb = P2 * (V1 - V2)
dUb = f * P2 * (V2 - V1)
Qb = dUb - Wb

# C)
Wc = 0
dUc = f * V2 * (P1 - P2)
Qc = dUc - Wc

# D)
Wd = P1 * (V2 - V1)
dUd = f * P1 * (V1 - V2)
Qd = dUd - Wd

#%%
# net work
Wnet = Wa + Wb + Wc + Wd
Wnet.factor()

#%%
# net heat
Qnet = Qa + Qb + Qc + Qd
Qnet.simplify().factor()

#%% net energy
dUnet = dUa + dUb + dUc + dUd
dUnet.simplify()

#%%

# %%
x = symbols('x', real = True)

# these are some muscle properties
f1 = symbols('f_1', positive = True)
g0, g1, g2 = symbols('g_0, g_1, g_2', positive = True)
h = symbols('h', positive = True)

# for calculating force, this is an important parameter
# Gamma = msk/(2 ell)
Gamma = symbols('Gamma')

# velocity (first assuming that v is positive)
vp = symbols('v', positive = True)
vn = symbols('v', negative = True)


# the factor varphi shows up everywhere, so we define it
varphi = f1/(f1 + g1)
psi_p = (f1+g1)*h/(2*vp)
psi_n = (f1+g1)*h/(2*vn)

def_vals = {f1: 15.0,
            g0: 170.0,
            g1: 8.0,
            g2: 25.0,
            h: 1.0}

def eval_default(expr):
    rv = expr
    for var in def_vals.keys():
        rv = rv.subs(var, def_vals[var])
    return rv

# %%
# this is for shortening
n_pos = Piecewise(
    (varphi*(1-exp(-psi_p))*exp(g0*x/vp), x <= 0),
    (varphi *(1-exp(((x/h)**2 - 1)*psi_p)), x <= h),
    (0, True)
)

term0 = varphi * (1-exp(-psi_p))*exp(g0*x/vp)  # -oo to 0
term1 = varphi * (1-exp(((x/h)**2 - 1)*psi_p)) # from 0 to h

force_pos = integrate(x * term0, (x, -oo, 0)) + integrate(x*term1, (x, 0, h))
force_pos

term1_n = varphi * (1 - exp((x/h)**2 * psi_n)) # from 0 to h
term2_n = varphi * (1 - exp(psi_n)) * exp((x-h)/(2*vn) * (2*g1 + g2 * (x/h - 1)))

force_neg = integrate(x*term1_n, (x, 0, h)) #+ integrate(x*term2_n, (x, h, oo))

# the second term we do by hand to get
u = vn/h
a = g2/(2*u)
b = (g2 - g1)/u
z1 = (1 - (b/(2*a))) / (1/sqrt(2*a))
term2_n_integrated = h**2/(2*a) * exp(-b**2/(2*a)) * exp(-z1**2 / 2) + (b/(2*a)) * exp(-b**2/(2*a)) * h**2 / sqrt(2*a) * sqrt(pi/2) * (1 - erf(z1/sqrt(2)))
#term2_n_integrated *= (f1/(f1+g1))*

#plot_piecewise(eval_default(n_pos.subs(v,10*h)), xlim = (-2*def_vals[h], 2*def_vals[h]))

#plot(eval_default(force), xlim = (0, 80.0))
#term2_n


# %%
