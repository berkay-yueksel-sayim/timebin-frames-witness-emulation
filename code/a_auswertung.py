# a_auswertung.py -- Talbot-Emulation: liest NUR a_tests_results_LAUF1_2026-10-03.json (Lauf 1) und wertet nach Pre-Reg v4 §2 aus.
# Zweiter Weg fuer die ENTSCHEIDUNGSLOGIK: sigma*/f* und die Urteile werden hier aus den Kurven NEU bestimmt,
# nicht aus den Booleans des Laufs uebernommen; danach wird beides verglichen.
import json, hashlib, pathlib

HIER = pathlib.Path(__file__).parent
P = HIER / "a_tests_results_LAUF1_2026-10-03.json"
SOLL_MD5 = "aba6a738ddc6bd4c4196b1f1fab7146d"
ist = hashlib.md5(P.read_bytes()).hexdigest()
print("JSON md5", ist, "== Lauf 1:", ist == SOLL_MD5)
r = json.loads(P.read_text(encoding="utf-8"))


def bruch(kurve):
    keys = sorted(kurve, key=float)
    for i, k in enumerate(keys):
        if not kurve[k]["zertifizierbar"]:
            return (float(keys[i - 1]) if i > 0 else 0.0), False, keys
    return None, True, keys


a1_ok = max(r["A1_max_abweichung"].values()) < 1e-12 and max(r["A1_lemma_w_matrix"].values()) < 1e-12
print(f"A1: max Entropie {max(r['A1_max_abweichung'].values()):.2e}, max Matrix {max(r['A1_lemma_w_matrix'].values()):.2e} -> "
      f"{'BESTAETIGT' if a1_ok else 'WIDERLEGT'} (Lauf sagt {r['A1_bestaetigt']})")

urteil = {"a": [], "b": [], "c": []}
for k, e in r["A2"].items():
    raten = e["a_rate_je_blindes_label"]
    a_ok = all(v["sigma_0.05"] >= 0.5 for v in raten.values())
    kontrolle_hoch = any(v["kontrolle_sigma_0"] >= 0.5 for v in raten.values())
    s_stern, zens_b, kb = bruch(e["b_kurve"])
    b_ok = (not zens_b) and 0.02 <= s_stern <= 0.30
    f_stern, zens_c, kc = bruch(e["c_kurve"])
    c_ok = (not zens_c) and 0.01 <= f_stern <= 0.30
    urteil["a"].append(a_ok); urteil["b"].append(b_ok); urteil["c"].append(c_ok)
    print(f"\nk={k}: blind {e['blind']}, bester rauschfreier Abstand {e['bester_abstand_rauschfrei']:.4f}")
    for lab, v in raten.items():
        print(f"  (a) Label {lab}: sigma=0,05 -> {v['sigma_0.05']:.3f} | Kontrolle sigma=0 -> {v['kontrolle_sigma_0']:.3f}")
    print(f"  (a) -> {'bestaetigt' if a_ok else 'widerlegt'} (Lauf: {e['a_bestaetigt']}); "
          f"Kontrollzeile >= 50 %: {kontrolle_hoch} -> Zuschreibung an Phasenrauschen {'NICHT erlaubt' if kontrolle_hoch else 'erlaubt'}")
    print(f"  (b) sigma* = {s_stern} zensiert {zens_b} -> {'bestaetigt' if b_ok else 'widerlegt'} (Lauf: {e['b_sigma_stern']}, {e['b_bestaetigt']})")
    for s in (kb[0], kb[len(kb) // 2], kb[-1]):
        print(f"      Delta_max({s}) = {e['b_kurve'][s]['delta_max']:.4f}")
    if not zens_b:
        i = kb.index(f"{s_stern:.3f}") if f"{s_stern:.3f}" in kb else None
        if i is not None:
            print(f"      am Bruch: Delta_max({kb[i]}) = {e['b_kurve'][kb[i]]['delta_max']:.4f}, Delta_max({kb[i+1]}) = {e['b_kurve'][kb[i+1]]['delta_max']:.4f}")
    print(f"  (c) f* = {f_stern} zensiert {zens_c} -> {'bestaetigt' if c_ok else 'widerlegt'} (Lauf: {e['c_f_stern']}, {e['c_bestaetigt']})")
    if not zens_c:
        i = kc.index(f"{f_stern:.3f}")
        for kk in kc[i:i + 2]:
            print(f"      f={kk}: Delta_max {e['c_kurve'][kk]['delta_max']:.4f}, bester verrauschter Abstand {e['c_kurve'][kk]['bester_abstand_verrauscht']:.4f}")


def gesamt(liste):
    return "bestaetigt" if all(liste) else ("teilweise" if any(liste) else "widerlegt")


print("\nGESAMT nach Pre-Reg §2 (beide k = bestaetigt, eines = teilweise):")
for h in ("a", "b", "c"):
    print(f"  A2({h}): {gesamt(urteil[h])}  {urteil[h]}")
