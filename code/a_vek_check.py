# a_vek_check.py -- Talbot-Emulation: Kontrolle der vektorisierten Bausteine gegen die skalaren (Pre-Reg v4, A-7).
# Das ist eine Code-Kontrolle, KEIN A2-Ergebnis: es wird keine Entropie-Abweichung, kein Fehlalarm, kein sigma*/f* berechnet.
# Verglichen werden nur Wahrscheinlichkeitsvektoren und Schaetzwerte bei IDENTISCHEN Eingaben.
import numpy as np
from a_kern import lift, schaetzer, schaetzer_batch, bin_wahrscheinlichkeiten_batch, entropie_batch, H, WORTE


def skalar_mit_phasen(k, a, wort, phis_zeile):
    F, pdiag, U, L = lift(k)
    v = U[:, a].copy(); stufe = 0
    for c in reversed(wort):
        if c == "S":
            v = F @ v
        else:
            v = pdiag * np.exp(1j * phis_zeile[stufe]) * v; stufe += 1
    p = np.abs(v)**2
    return p / p.sum()


rng = np.random.default_rng(424242)
worst_p = 0.0
n = 25
for k in (1, 2, 3, 4):
    F, pdiag, U, L = lift(k)
    for a in range(k + 1):
        for w in WORTE:
            nT = w.count("T")
            phis = [rng.normal(0.0, 0.3, (n, L)) for _ in range(nT)]
            pb = bin_wahrscheinlichkeiten_batch(k, a, w, 0.0, rng, n, phis=phis if nT else None)
            for i in range(n):
                ps = skalar_mit_phasen(k, a, w, [ph[i] for ph in phis])
                worst_p = max(worst_p, float(np.max(np.abs(pb[i] - ps))))
print(f"Vektorisierung: max |p_batch - p_skalar| = {worst_p:.2e}  (k=1..4, alle Labels, sieben Woerter, 25 Zeilen)")

worst_h = 0.0
for _ in range(200):
    p = rng.dirichlet(np.ones(12) * rng.uniform(0.05, 2.0))
    worst_h = max(worst_h, abs(float(entropie_batch(p[None, :])[0]) - H(p)))
print(f"entropie_batch gegen H: max Abweichung {worst_h:.2e}")

worst_s = 0.0
for _ in range(200):
    L = int(rng.integers(6, 15))
    p = rng.dirichlet(np.ones(L) * rng.uniform(0.05, 2.0))
    c = rng.multinomial(10_000, p)
    h1, s1 = schaetzer(c)
    h2, s2 = schaetzer_batch(c[None, :])
    worst_s = max(worst_s, abs(h1 - h2[0]), abs(s1 - s2[0]))
print(f"schaetzer_batch gegen schaetzer: max Abweichung {worst_s:.2e}  (200 Zufallszaehlungen)")

# 2D-Multinomial: Zeilensummen = N, Mittelwerte nahe p (nur eine Funktionskontrolle)
p2 = bin_wahrscheinlichkeiten_batch(2, 0, "ST", 0.05, rng, 4000)
cnt = rng.multinomial(10_000, p2)
print("2D-Multinomial: Zeilensummen alle 10000:", bool((cnt.sum(axis=1) == 10_000).all()), " Form", cnt.shape)

# Negativkontrolle: die Stufenreihenfolge vertauschen (Wort "ST" statt "TS" in der Batch-Fassung) muss gegen die skalare Fassung auffallen
F, pdiag, U, L = lift(3)
phis = [rng.normal(0.0, 0.3, (5, L))]
pb_falsch = bin_wahrscheinlichkeiten_batch(3, 1, "TS", 0.0, rng, 5, phis=phis)     # Batch: erst S, dann T
ps_ref = np.array([skalar_mit_phasen(3, 1, "ST", [ph]) for ph in phis[0]])       # skalar: erst T, dann S
nk_abw = float(np.max(np.abs(pb_falsch - ps_ref)))
print(f"Negativkontrolle (vertauschte Stufenreihenfolge): Abweichung {nk_abw:.2e}, sichtbar: {nk_abw > 1e-3}")

ok = worst_p < 1e-12 and worst_h < 1e-12 and worst_s < 1e-12 and bool((cnt.sum(axis=1) == 10_000).all()) and nk_abw > 1e-3
print("KONTROLLE:", "BESTANDEN" if ok else "NICHT bestanden")
