# a_quantil_exakt.py -- Talbot-Emulation: Aufloesung einer Abweichung zwischen Monte Carlo und zweiter Route (k=2, sigma=0,02: MC 0,0208 gegen 0,0213, rel 2,37 %).
# Eigene halbanalytische Route, nach einem Hinweis der zweiten Route gebaut, mit eigenem Code:
#  (1) PRUEFEN statt annehmen: wirkt Phasenrauschen fuer die blinden Labels nur im Frame ST, und dort nur ueber
#      Delta = phi_(+m) - phi_(-m)? -> Test: zufaellige Phasen auf ALLEN Bins gegen nur Delta auf einem Bin.
#  (2) Ist d(Delta) = |H_bin - 1 - h| monoton in |Delta|? Dann ist das 99,865-%-Quantil von d = d(q), mit
#      q = sqrt(2) * sigma * z,  P(|Z| <= z) = 0,99865  (exakt, ohne Ziehungsfehler).
#  (3) Erwarteter Ziehungsfehler meines MC-Quantils bei 1e5 Ziehungen (Ordnungsstatistik), zum Einordnen der 2,37 %.
import json, math, hashlib, pathlib
import numpy as np
from a_kern import lift, H, WORTE, bin_wahrscheinlichkeiten_batch, entropie_batch, label_entropien
from a_lauf import blind_labels

HIER = pathlib.Path(__file__).parent
r = json.loads((HIER / "a_tests_results_LAUF1_2026-10-03.json").read_text(encoding="utf-8"))["A2"]

def p_mit_phasen(k, a, wort, phi_je_stufe):
    F, pdiag, U, L = lift(k)
    v = U[:, a].astype(complex); st = 0
    for c in reversed(wort):
        if c == "S":
            v = F @ v
        else:
            v = pdiag * np.exp(1j * phi_je_stufe[st]) * v; st += 1
    p = np.abs(v)**2
    return p / p.sum()

# z mit P(|Z| <= z) = 0,99865 per Bisektion (ohne scipy)
lo, hi = 0.0, 10.0
for _ in range(200):
    mid = (lo + hi) / 2
    if math.erf(mid / math.sqrt(2)) < 0.99865: lo = mid
    else: hi = mid
z = (lo + hi) / 2
print(f"z (zweiseitig 99,865 %) = {z:.6f}")

rng = np.random.default_rng(99)
for k in (2, 4):
    blind, h = blind_labels(k)
    F, pdiag, U, L = lift(k)
    print(f"\n== k={k}, blinde Labels {blind}")
    # (1) welche Frames reagieren ueberhaupt auf Phasen?
    for a in blind:
        for wi, w in enumerate(WORTE):
            nT = w.count("T")
            if nT == 0: continue
            dev = max(abs(H(p_mit_phasen(k, a, w, [rng.normal(0, 0.3, L) for _ in range(nT)])) - 1 - h[a, wi]) for _ in range(200))
            print(f"   Label {a} Frame {w or 'id':3s}: max |H-1-h| bei sigma=0,3 ueber 200 Ziehungen = {dev:.2e}")
        # (1b) im Frame ST: haengt p nur von Delta = phi_(+m) - phi_(-m) ab?  m = a+1, -m = L-m
        m = a + 1
        phi = rng.normal(0, 0.3, L)
        nur_delta = np.zeros(L); nur_delta[m] = phi[m] - phi[L - m]
        diff = np.max(np.abs(p_mit_phasen(k, a, "ST", [phi]) - p_mit_phasen(k, a, "ST", [nur_delta])))
        print(f"   Label {a} Frame ST: alle Bins zufaellig gegen nur Delta auf Bin +m: max |dp| = {diff:.2e}")
    # (2) Monotonie und exakter Quantilwert, nur Frame ST
    wi = WORTE.index("ST")
    for s in (0.02, 0.05, 0.10, 0.14, 0.15):
        werte = []
        for a in blind:
            m = a + 1
            def d_von(delta):
                ph = np.zeros(L); ph[m] = delta
                return abs(H(p_mit_phasen(k, a, "ST", [ph])) - 1 - h[a, wi])
            q = math.sqrt(2) * s * z
            gitter = np.linspace(0, q, 400)
            dg = np.array([d_von(x) for x in gitter])
            monoton = bool(np.all(np.diff(dg) >= -1e-15))
            werte.append((d_von(q), monoton))
        exakt = max(w for w, _ in werte)
        mc = r[str(k)]["b_kurve"][f"{s:.3f}"]["delta_max"]
        print(f"   sigma {s:.2f}: exakt {exakt:.5f} (monoton bis q: {[m for _, m in werte]})  |  MC {mc:.5f}  rel {100*(mc-exakt)/exakt:+.2f} %")

# (3) Ziehungsfehler des MC-Quantils (k=2, sigma=0,02): 300 Wiederholungen mit je 1e5 Ziehungen, nur Frame ST
k, s = 2, 0.02
blind, h = blind_labels(k); a = blind[0]; wi = WORTE.index("ST")
qs = []
rng2 = np.random.default_rng(7)
for _ in range(300):
    p = bin_wahrscheinlichkeiten_batch(k, a, "ST", s, rng2, 100_000)
    qs.append(float(np.quantile(np.abs(entropie_batch(p) - 1 - h[a, wi]), 0.99865)))
qs = np.array(qs)
print(f"\nMC-Quantil k=2 sigma=0,02, 300 x 1e5 Ziehungen: Mittel {qs.mean():.5f}, Streuung {qs.std(ddof=1):.5f} "
      f"(rel {100*qs.std(ddof=1)/qs.mean():.2f} %); Anteil <= 0,02079 (mein Lauf): {np.mean(qs <= 0.020794):.3f}")
