# a_reproduzieren.py -- Talbot-Emulation: rechnet die pre-registrierten Tests A1 und A2 OHNE die Lock-Pruefung nach
# und vergleicht das Ergebnis mit Lauf 1 (a_tests_results_LAUF1_2026-10-03.json).
# Warum: a_tests.py prueft vor dem Rechnen die gelockte Pre-Registrierung; deren Datei liegt diesem Record (v1.0) noch
# nicht bei, nur ihr SHA-256 (A_PREREG_LOCK.sha256, Manuskript). Die Rechen-Aufrufe hier sind dieselben wie in
# a_tests.py (gleiche Funktionen, gleiche Festlegungen a_lauf.PREREG, gleiche feste Saaten).
# Laufzeit: Minuten bis etwa eine Stunde, je nach Rechner.
import json, pathlib
from a_lauf import PREREG, run_a1, run_a2_k

HIER = pathlib.Path(__file__).parent
neu = {}
a1, a1_matrix = run_a1([5, 7, 8, 9, 11, 12])
neu["A1_max_abweichung"] = a1
neu["A1_lemma_w_matrix"] = a1_matrix
neu["A1_bestaetigt"] = bool(max(a1.values()) < 1e-12 and max(a1_matrix.values()) < 1e-12)
neu["A2"] = {str(k): run_a2_k(k, PREREG) for k in (2, 4)}
neu = json.loads(json.dumps(neu))   # gleiche Darstellung wie in der gespeicherten Datei (Schluessel als Text)
(HIER / "a_tests_results_reproduziert.json").write_text(json.dumps(neu, indent=1), encoding="utf-8")

lauf1 = json.loads((HIER / "a_tests_results_LAUF1_2026-10-03.json").read_text(encoding="utf-8"))
print("Vergleich mit Lauf 1 (ohne Umgebung und Pre-Reg-Kopf):")
gleich_alle = True
for schluessel in ("A1_bestaetigt", "A2"):
    gleich = neu[schluessel] == lauf1[schluessel]
    gleich_alle &= gleich
    print(f"  {schluessel:20s}: {'identisch' if gleich else 'ABWEICHUNG'}")
for schluessel in ("A1_max_abweichung", "A1_lemma_w_matrix"):
    # Rundungsfehler-Groessen: die Bestaetigung haengt nur an < 1e-12
    groesste = max(abs(neu[schluessel][k] - lauf1[schluessel][k]) for k in lauf1[schluessel])
    print(f"  {schluessel:20s}: groesste Differenz zu Lauf 1 {groesste:.1e} (beide < 1e-12: "
          f"{max(neu[schluessel].values()) < 1e-12 and max(lauf1[schluessel].values()) < 1e-12})")
print("ERGEBNIS:", "A2 und A1-Urteil identisch mit Lauf 1" if gleich_alle else "ABWEICHUNG -- siehe oben")
