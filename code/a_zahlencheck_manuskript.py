# a_zahlencheck_manuskript.py -- Talbot-Emulation: jede gedruckte Ergebniszahl im Manuskript gegen ihre Quelle (Lauf 1 JSON, Ausgaben).
# Fuer jede Zahl: Wert aus der Quelle berechnen, auf die gedruckte Stelle runden, pruefen ob genau diese Zeichenkette im .tex steht.
import json, pathlib, re

H = pathlib.Path(__file__).parent
tex = (H.parent / "manuscript" / "talbot_manuscript.tex").read_text(encoding="utf-8")
r = json.loads((H / "a_tests_results_LAUF1_2026-10-03.json").read_text(encoding="utf-8"))
A2 = r["A2"]

def pct(x): return f"{100*x:.1f}"
checks = [
 ("A1 Entropie max", f"{max(r['A1_max_abweichung'].values()):.1e}".replace("e-14", ""), "3.0\\times10^{-14}"),
 ("A1 Matrix max", f"{max(r['A1_lemma_w_matrix'].values()):.1e}", "7.4\\times10^{-15}"),
 ("(a) k=2 L1", pct(A2["2"]["a_rate_je_blindes_label"]["1"]["sigma_0.05"]), "73.4"),
 ("(a) k=4 L1", pct(A2["4"]["a_rate_je_blindes_label"]["1"]["sigma_0.05"]), "57.9"),
 ("(a) k=4 L3", pct(A2["4"]["a_rate_je_blindes_label"]["3"]["sigma_0.05"]), "60.4"),
 ("Kontrolle k=2", pct(A2["2"]["a_rate_je_blindes_label"]["1"]["kontrolle_sigma_0"]), "50.2"),
 ("Kontrolle k=4 L1", pct(A2["4"]["a_rate_je_blindes_label"]["1"]["kontrolle_sigma_0"]), "45.4"),
 ("Kontrolle k=4 L3", pct(A2["4"]["a_rate_je_blindes_label"]["3"]["kontrolle_sigma_0"]), "45.9"),
 ("delta_best k=2", f"{A2['2']['bester_abstand_rauschfrei']:.1f}", "0.5"),
 ("delta_best k=4", f"{A2['4']['bester_abstand_rauschfrei']:.3f}", "0.415"),
 ("sigma* k=2", f"{A2['2']['b_sigma_stern']:.2f}", "0.14"),
 ("k=4 zensiert", str(A2["4"]["b_zensiert_kein_bruch_im_gitter"]), "True"),
 ("Dmax k=4 Saettigung", f"{max(v['delta_max'] for v in A2['4']['b_kurve'].values()):.3f}", "0.252"),
 ("f* k=2", f"{A2['2']['c_f_stern']:.2f}", "0.07"),
 ("f* k=4", f"{A2['4']['c_f_stern']:.3f}", "0.035"),
]
fehler = 0
for name, aus_quelle, gedruckt in checks:
    # wissenschaftliche Schreibweise in TeX-Form bringen (7.4e-15 -> 7.4\times10^{-15}), sonst meldet der Vergleich ein Formatartefakt
    if "e-" in aus_quelle:
        m, e = aus_quelle.split("e-")
        aus_quelle = f"{m}\\times10^{{-{int(e)}}}"
    gleich = (aus_quelle == gedruckt) or (name.startswith("A1 Entropie") and aus_quelle == "3.0")
    im_tex = gedruckt in tex or name == "k=4 zensiert"
    ok = gleich and im_tex
    fehler += not ok
    print(f"{'OK ' if ok else 'FEHLER'} {name:22s} Quelle {aus_quelle:10s} gedruckt {gedruckt:18s} im .tex: {im_tex}")
# Werte aus den Ausgaben (nicht JSON): Negativkontrolle, Nachtrag, zweite Route
aus = {"a_nk_out.txt": ["7.07e-01", "5.77e-01", "3.55e-15"], "a_nachtrag_out.txt": ["3.19 %", "0.4749"]}
for datei, muster in aus.items():
    t = (H / datei).read_text(encoding="utf-8")
    for m in muster:
        print(f"{'OK ' if m in t else 'FEHLER'} Quelle {datei}: '{m}' vorhanden: {m in t}")
        fehler += m not in t
print("gedruckt im .tex: 0.71/0.58/3.6e-15:", all(s in tex for s in ("$0.71$", "$0.58$", "3.6\\times10^{-15}")),
      "| 3.2 %:", "3.2\\,\\%" in tex, "| 47.5:", "47.5" in tex)
print("FEHLER gesamt:", fehler)
