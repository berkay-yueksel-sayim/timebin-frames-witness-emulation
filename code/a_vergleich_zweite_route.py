# a_vergleich_zweite_route.py -- Talbot-Emulation: Vergleich Lauf 1 (Monte Carlo) gegen eine zweite, unabhaengig
# geschriebene Route (halbanalytisch + eigenes MC), nach Kriterien, die VOR dem Sehen der Lauf-1-Zahlen festgelegt wurden:
#   (c) deterministisch: f* gleich, Dmax(f) auf 1e-4 gleich
#   (b) Dmax(sigma) innerhalb ~2 % (relativ); sigma* = 0,14 ODER 0,15 bei k=2 ist Uebereinstimmung; k=4 zensiert
#   (a) Abstand <= 0,045 je Wert
# Die Werte der zweiten Route (R2) stehen hier von Hand abgeschrieben; ihr Code liegt als a2_second_route.py bei, die
# Eintraege "k" ihrer Ergebnisdatei a2_second_route_result.json sind dieselben Zahlen.
# (English: the values of the second route (R2) are transcribed here by hand; its code is deposited as a2_second_route.py,
# and the "k" entries of its result file a2_second_route_result.json hold the same numbers.)
import json, hashlib, pathlib

HIER = pathlib.Path(__file__).parent
P = HIER / "a_tests_results_LAUF1_2026-10-03.json"
assert hashlib.md5(P.read_bytes()).hexdigest() == "aba6a738ddc6bd4c4196b1f1fab7146d", "nicht Lauf 1"
r = json.loads(P.read_text(encoding="utf-8"))["A2"]

R2 = {
 "a": {"2": {"1": (0.714, 0.462)}, "4": {"1": (0.618, 0.454), "3": (0.623, 0.438)}},
 "b_dmax": {"2": {0.02: 0.0213, 0.05: 0.0988, 0.10: 0.2885, 0.13: 0.4173, 0.14: 0.4606, 0.15: 0.5036, 0.20: 0.7052, 0.30: 0.9674},
            "4": {0.05: 0.0500, 0.10: 0.1305, 0.15: 0.2016, 0.20: 0.2438, 0.30: 0.2516, 0.40: 0.2516}},
 "b_stern": {"2": 0.14, "4": None},
 "c_dmax": {"2": {0.01: 0.0756, 0.05: 0.2901, 0.07: 0.3801, 0.10: 0.5032},
            "4": {0.01: 0.0889, 0.035: 0.2579, 0.05: 0.3466, 0.10: 0.6073}},
 "c_stern": {"2": 0.070, "4": 0.035},
 "delta_best": {"2": 0.5000, "4": 0.4150},
}

befunde = []
print("== (a)  Kriterium |MC - R2| <= 0,045 ==")
for k, labs in R2["a"].items():
    for lab, (r2_rate, r2_kontr) in labs.items():
        mine = r[k]["a_rate_je_blindes_label"][lab]
        for name, a, b in (("sigma=0,05", mine["sigma_0.05"], r2_rate), ("Kontrolle sigma=0", mine["kontrolle_sigma_0"], r2_kontr)):
            d = abs(a - b); ok = d <= 0.045
            print(f"  k={k} Label {lab} {name:18s}: MC {a:.3f}  R2 {b:.3f}  |d| {d:.3f}  {'OK' if ok else 'BEFUND'}")
            if not ok: befunde.append(("a", k, lab, name, a, b))

print("== (b)  Kriterium relativ <= 2 %; sigma* 0,14 oder 0,15 (k=2); k=4 zensiert ==")
for k, tab in R2["b_dmax"].items():
    for s, r2v in tab.items():
        mine = r[k]["b_kurve"][f"{s:.3f}"]["delta_max"]
        rel = abs(mine - r2v) / r2v
        ok = rel <= 0.02
        print(f"  k={k} sigma {s:.2f}: MC {mine:.4f}  R2 {r2v:.4f}  rel {100*rel:5.2f} %  {'OK' if ok else 'BEFUND'}")
        if not ok: befunde.append(("b", k, s, mine, r2v))
    print(f"  k={k} delta_best: MC {r[k]['bester_abstand_rauschfrei']:.4f}  R2 {R2['delta_best'][k]:.4f}")
s2 = r["2"]["b_sigma_stern"]; z4 = r["4"]["b_zensiert_kein_bruch_im_gitter"]
print(f"  sigma* k=2: MC {s2}  (Uebereinstimmung, wenn 0,14 oder 0,15) -> {'OK' if s2 in (0.14, 0.15) else 'BEFUND'}")
print(f"  k=4 zensiert: MC {z4}  R2 zensiert -> {'OK' if z4 else 'BEFUND'}")

print("== (c)  Kriterium f* gleich, Dmax(f) |d| <= 1e-4 ==")
for k, tab in R2["c_dmax"].items():
    for f, r2v in tab.items():
        mine = r[k]["c_kurve"][f"{f:.3f}"]["delta_max"]
        d = abs(mine - r2v); ok = d <= 1e-4 + 5e-5      # R2 druckt auf 4 Stellen: Rundung bis 5e-5 zulassen
        print(f"  k={k} f {f:.3f}: MC {mine:.5f}  R2 {r2v:.4f}  |d| {d:.1e}  {'OK' if ok else 'BEFUND'}")
        if not ok: befunde.append(("c", k, f, mine, r2v))
    fs = r[k]["c_f_stern"]
    print(f"  f* k={k}: MC {fs}  R2 {R2['c_stern'][k]} -> {'OK' if abs(fs - R2['c_stern'][k]) < 1e-9 else 'BEFUND'}")

print("\nBEFUNDE ausserhalb der vorab festgelegten Grenzen:", len(befunde))
for b in befunde:
    print("  ", b)
