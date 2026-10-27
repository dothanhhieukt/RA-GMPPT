# Averaged model of the PV-battery buck charger (development copy of avg_plant.m)
import numpy as np
from scipy.optimize import brentq

q = 1.602176634e-19; kB = 1.380649e-23
# ASMS-230M, NREL SAM single-diode parameters (as in the Simulink PV Array block)
IL0 = 8.388; I00 = 4.26e-10; NID = 1.019; RSH = 148.4; RS = 0.378; NCELL = 60
A_ISC = 0.0996e-2; B_VOC = -0.370e-2; TREF = 298.15; EG = 1.22   # effective band gap fitted to the Simulink GMPP at 50 and 65 C
VF_BP = 0.8; RON_BP = 1e-3          # bypass diode (Simulink default diode)

def mod_params(G, T, egEff=None):
    Tk = T + 273.15
    Vt = kB*Tk/q
    a = NID*NCELL*Vt
    Iph = IL0*G/1000*(1 + A_ISC*(T-25))
    eg = EG if egEff is None else egEff
    I0 = I00*(Tk/TREF)**3*np.exp(q*eg/(NID*NCELL*kB)*NCELL*(1/TREF - 1/Tk)) if False else \
         I00*(Tk/TREF)**3*np.exp(q*eg/(NID*kB)*(1/TREF - 1/Tk))
    return Iph, I0, a, RSH*1000/max(G, 1)

def vmod_of_i(I, Iph, I0, a, rsh=RSH):
    """module voltage for current array I; bypass diode clamps the voltage at -(Vf + Ron*I)"""
    I = np.asarray(I, float)
    Vg = np.linspace(-0.79, 45.0, 4000)
    Ig = np.full_like(Vg, Iph)
    for _ in range(50):
        ex = np.exp(np.clip((Vg + Ig*RS)/a, -60, 60))
        f = Iph - I0*(ex-1) - (Vg+Ig*RS)/rsh - Ig
        df = -I0*ex*RS/a - RS/rsh - 1
        Ig = Ig - f/df
    m = Ig > -1e-6
    Vm = Vg[m][::-1]; Im = Ig[m][::-1]          # increasing current
    V = np.interp(I, Im, Vm, left=Vm[0], right=-10.0)
    Vbp = -(VF_BP + RON_BP*np.maximum(I, 0))
    return np.maximum(V, Vbp)

class PVString:
    def __init__(self, G, T, n=6000, egEff=None):
        G = np.atleast_1d(G).astype(float); T = np.broadcast_to(np.atleast_1d(T).astype(float), G.shape)
        self.G, self.T = G, T
        Imax = max(IL0*g/1000*(1+A_ISC*(t-25)) for g, t in zip(G, T))*1.02 + 0.05
        Igrid = np.linspace(0, Imax, n)
        Vs = np.zeros(n)
        for g, t in zip(G, T):
            if g <= 0: Vs += -(VF_BP + RON_BP*Igrid); continue
            Iph, I0, a, rsh = mod_params(g, t, egEff)
            Vs += vmod_of_i(Igrid, Iph, I0, a, rsh)
        o = np.argsort(Vs)
        self.V = Vs[o]; self.I = Igrid[o]
        self.Voc = self.V[self.I > 1e-9].max() if np.any(self.I > 1e-9) else 0
        # ensure monotone: I decreasing in V
        self.I = np.maximum.accumulate(self.I[::-1])[::-1]
    def i(self, v):
        return np.interp(v, self.V, self.I, left=self.I[0], right=0.0)
    def slope(self, v, dv=0.05):
        return (self.i(v+dv) - self.i(v-dv))/(2*dv)
    def pv(self, vmin=0):
        V = np.linspace(max(vmin, 0), self.V.max(), 4000); return V, V*self.i(V)

# battery Samsung INR18650-35E pack
SOCtab = np.array([0.0000,0.0413,0.0855,0.1150,0.1593,0.1888,0.2625,0.3363,0.4100,0.4838,0.5575,0.6313,0.7050,0.7788,0.8525,0.9263,1.0000])
OCVtab = np.array([2.807,3.007,3.207,3.307,3.427,3.447,3.517,3.587,3.657,3.727,3.807,3.867,3.917,3.997,4.067,4.107,4.167])
def battery(pack, soc0, Qscale=1.0):
    Ns = 13 if pack == 48 else 20
    return dict(Ns=Ns, Q=6*3.39*Qscale, R=Ns*0.045/6, soc=soc0, Vcv=54.6 if pack == 48 else 84.0)
def ocv(b): return b['Ns']*np.interp(b['soc'], SOCtab, OCVtab)

# converter
L = 1e-3; CIN = 100e-6; RD = 3.3; CD = 150e-6; VF_DB = 0.5; VF_DF = 0.5; RON = 0.02; CPV = 1e-6/4

class Plant:
    def __init__(self, pack, soc0, Qscale=1.0, v0=None):
        self.b = battery(pack, soc0, Qscale)
        self.v = 150.0 if v0 is None else v0; self.vd = self.v; self.iL = 0.0
    def vbat(self, i): return ocv(self.b) + self.b['R']*i
    def step(self, D, pv, dt, nsub=5):
        h = dt/nsub
        for _ in range(nsub):
            v = self.v
            ipv = pv.i(v); g = pv.slope(v)
            iin = D*self.iL
            id_ = (v - self.vd)/RD
            C = CIN + 4*CPV/4
            f = ipv - iin - id_
            dv = h*f/(C - h*(g - 1/RD))           # semi-implicit on PV slope and Rd
            self.v = v + dv
            self.vd += h*(self.v - self.vd)/(RD*CD)
            # inductor (diode-blocked: iL >= 0 -> DCM approximated by clamp)
            vb = self.vbat(self.iL)
            vL = D*(self.v - VF_DB - RON*self.iL) - (1-D)*VF_DF - vb
            self.iL = max(0.0, self.iL + h*vL/L)
            self.b['soc'] += h*self.iL/(3600*self.b['Q'])
        return self.v, pv.i(self.v), self.vbat(self.iL), self.iL

class FastPV:
    """uniform-grid lookup of I(V) for the time loop"""
    def __init__(self, pvs, dv=0.01):
        self.dv = dv; self.n = int(160/dv) + 2
        V = np.arange(self.n)*dv
        self.Itab = pvs.i(V).tolist(); self.pvs = pvs
    def i_g(self, v):
        if v <= 0: return self.Itab[0], 0.0
        x = v/self.dv; k = int(x)
        if k >= self.n - 1: return 0.0, 0.0
        i0 = self.Itab[k]; i1 = self.Itab[k+1]
        return i0 + (x-k)*(i1-i0), (i1-i0)/self.dv
