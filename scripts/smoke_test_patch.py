"""Smoke test for patched eval_pcf_forward."""
import sys, importlib.util
from mpmath import mp, mpf

spec = importlib.util.spec_from_file_location("relay", "siarc_t1t6_relay.py")
mod = importlib.util.module_from_spec(spec)
sys.modules["relay"] = mod
spec.loader.exec_module(mod)

# Golden ratio: a(n)=1, b(n)=1. K = phi ~ 1.618
a = lambda n: mpf(1)
b = lambda n: mpf(1)
K = mod.eval_pcf_forward(a, b, 200, 50)
mp.dps = 50
phi = (1 + mpf(5)**mpf("0.5")) / 2
print("Golden ratio test:")
print("  eval_pcf_forward =", K)
print("  phi              =", phi)
print("  match:", abs(K - phi) < mpf("1e-40"))

# a(1)=0 test: a(n) = n(n-1), b(n)=1. a(1)=0, b(0)=1 => K=1
a2 = lambda n: mpf(n**2 - n)
b2 = lambda n: mpf(1)
K2 = mod.eval_pcf_forward(a2, b2, 200, 50)
print("\na(1)=0 test:")
print("  K =", K2, ", b(0) = 1, match:", abs(K2 - 1) < mpf("1e-40"))

# d=2 fam=67 from the survey: a=[-4,3,1], b=[-2,2,-1]
# a(1)=0, b(0)=-2 => K=-2
a3 = lambda n: mpf(-4 + 3*n + n**2)
b3 = lambda n: mpf(-2 + 2*n - n**2)
K3 = mod.eval_pcf_forward(a3, b3, 200, 50)
print("\nd=2 fam=67 test:")
print("  K =", K3, ", b(0) = -2, match:", abs(K3 - (-2)) < mpf("1e-40"))
