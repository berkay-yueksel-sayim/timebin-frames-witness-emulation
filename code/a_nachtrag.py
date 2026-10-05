# a_nachtrag.py -- Talbot-Emulation: Nachtrag zu zwei Ergaenzungen der zweiten Route (1: exaktes Quantil, 3: Kontrollzeile).
#  (1) Lauf-1-MC gegen das EXAKTE 99,865-%-Quantil auf dem GANZEN sigma-Gitter (40 x 2 = 80 Werte). Exakt heisst hier:
#      Verteilung von d(Delta), Delta ~ N(0, 2 sigma^2), auf einem feinen Delta-Gitter mit Normalgewichten, Quantil aus der
#      gewichteten Verteilungsfunktion -- OHNE Monotonie-Annahme (Hinweis der zweiten Route). Nur Frame ST, nur Delta (in a_quantil_exakt.py
#      selbst nachgeprueft). Maximum ueber die blinden Labels.
#  (3) Kontrollzeile sigma = 0 hochstatistisch mit EIGENER Saat: 40000 Wiederholungen, k = 2 und k = 4.
#  Kein Eingriff in den gelockten Lauf; Lauf 1 wird nur gelesen (md5 geprueft). Ergebnis ist ein NACHTRAG, nicht Teil der Pre-Reg.
import json, math, hashlib, pathlib
import numpy as np
from a_kern import lift, H, WORTE
from a_lauf import blind_labels, fehlalarm_naiv, PREREG

HIER = pathlib.Path(__file__).parent
P = HIER / "a_tests_results_LAUF1_2026-10-03.json"
assert hashlib.md5(P.read_bytes()).hexdigest() == "aba6a738ddc6bd4c4196b1f1fab7146d"
r = json.loads(P.read_text(encoding="utf-8"))["A2"]

def p_st(k, a, delta):
    F, pdiag, U, L = lift(k)
    m = a + 1
    phi = np.zeros(L); phi[m] = delta
    v = F @ (pdiag * np.exp(1j * phi) * U[:, a].astype(complex))
    p = np.abs(v)**2
    return p / p.sum()

def quantil_exakt(k, sigma, q=0.99865, n=40001):
    blind, h = blind_labels(k); wi = WORTE.index("ST")
    best = 0.0
    s = math.sqrt(2) * sigma
    x = np.linspace(-8 * s, 8 * s, n)
    w = np.exp(-0.5 * (x / s)**2); w /= w.sum()
    for a in blind:
        d = np.array([abs(H(p_st(k, a, xi)) - 1 - h[a, wi]) for xi in x])
        o = np.argsort(d); cdf = np.cumsum(w[o])
        best = max(best, float(d[o][np.searchsorted(cdf, q)]))
    return best

print("(1) Lauf-1-MC gegen exaktes Quantil, ganzes Gitter")
alle = []
for k in (2, 4):
    rel = []
    for key, v in r[str(k)]["b_kurve"].items():
        s = float(key)
        ex = quantil_exakt(k, s)
        rel.append((s, (v["delta_max"] - ex) / ex))
    arr = np.array([x[1] for x in rel])
    ueber2 = [f"{s:.2f} ({100*x:+.2f} %)" for s, x in rel if abs(x) > 0.02]
    print(f"  k={k}: {len(rel)} Werte | Mittel {100*arr.mean():+.2f} % | Streuung {100*arr.std(ddof=1):.2f} % | "
          f"max |rel| {100*np.abs(arr).max():.2f} % | ueber 2 %: {len(ueber2)} -> {ueber2}")
    alle += list(arr)
alle = np.array(alle)
print(f"  gesamt: {len(alle)} Werte, max |rel| {100*np.abs(alle).max():.2f} %, Mittel {100*alle.mean():+.2f} %")

print("(3) Kontrollzeile sigma = 0, 40000 Wiederholungen, eigene Saat 31337+")
for k in (2, 4):
    blind, _ = blind_labels(k)
    for lab in blind:
        rate = fehlalarm_naiv(k, lab, 0.0, np.random.default_rng(31337 + 10 * k + lab), PREREG["NPH"], 40000)
        se = math.sqrt(rate * (1 - rate) / 40000)
        lauf1 = r[str(k)]["a_rate_je_blindes_label"][str(lab)]["kontrolle_sigma_0"]
        print(f"  k={k} Label {lab}: {rate:.4f} +- {se:.4f} | Lauf 1 (1000 Wdh.) {lauf1:.3f} = {(lauf1 - rate)/math.sqrt(rate*(1-rate)/1000):+.1f} SE(1000)")
