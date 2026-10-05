# a_tests.py -- Talbot-Emulation: die pre-registrierten Tests A1 (neue k) und A2 (Fehlerbudget), Fassung zu Pre-Reg v4.
# LAEUFT NUR NACH DEM LOCK: bricht ab, wenn A_PREREG_LOCK.sha256 fehlt oder nicht zur gelockten Datei passt.
# Die Rechenfunktionen stehen in a_lauf.py; die Groessen in a_lauf.PREREG sind die der Pre-Reg v4.
import hashlib, json, pathlib, platform, sys
import numpy as np
from a_lauf import PREREG, run_a1, run_a2_k

HIER = pathlib.Path(__file__).parent
LOCK = HIER / "A_PREREG_LOCK.sha256"
if not LOCK.exists():
    print("ABBRUCH: Pre-Reg nicht gelockt (A_PREREG_LOCK.sha256 fehlt)."); sys.exit(2)
soll, datei = LOCK.read_text(encoding="utf-8").split()[:2]
ist = hashlib.sha256((HIER / datei).read_bytes().replace(b"\r\n", b"\n")).hexdigest()
if ist != soll:
    print("ABBRUCH: gelockte Pre-Reg wurde veraendert:", datei); sys.exit(2)
print("Pre-Reg gelockt und unveraendert:", datei, soll[:16])

out = {"prereg": datei, "prereg_sha256": soll,
       "umgebung": {"python": platform.python_version(), "numpy": np.__version__, "system": platform.platform()},
       "festlegungen": {k: v for k, v in PREREG.items() if not k.endswith("GITTER")} |
                       {"sigma_gitter": [PREREG["SIGMA_GITTER"][0], PREREG["SIGMA_GITTER"][-1], len(PREREG["SIGMA_GITTER"])],
                        "f_gitter": [PREREG["F_GITTER"][0], PREREG["F_GITTER"][-1], len(PREREG["F_GITTER"])]}}

a1, a1_matrix = run_a1([5, 7, 8, 9, 11, 12])
out["A1_max_abweichung"] = a1
out["A1_lemma_w_matrix"] = a1_matrix
out["A1_bestaetigt"] = bool(max(a1.values()) < 1e-12 and max(a1_matrix.values()) < 1e-12)
print("A1:", a1, "Matrix:", a1_matrix, "bestaetigt:", out["A1_bestaetigt"])

out["A2"] = {}
for k in (2, 4):
    e = run_a2_k(k, PREREG)
    out["A2"][str(k)] = e
    print(f"A2 k={k}: blind {e['blind']}; (a) {e['a_rate_je_blindes_label']}; (b) sigma* = {e['b_sigma_stern']} "
          f"(zensiert: {e['b_zensiert_kein_bruch_im_gitter']}); (c) f* = {e['c_f_stern']} (zensiert: {e['c_zensiert_kein_bruch_im_gitter']})")
(HIER / "a_tests_results.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
print("geschrieben: a_tests_results.json")
