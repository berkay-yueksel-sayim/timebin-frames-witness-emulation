# a_lauf.py -- Talbot-Emulation: die Rechen-Funktionen fuer A1 und A2 (Pre-Reg v4). KEIN Top-Level-Lauf.
# Aufgerufen von a_tests.py (nur nach dem Lock, Groessen aus der Pre-Reg) und von a_probelauf.py (Code-Test auf einem
# k, das NICHT in A1/A2 steht, mit reduzierten Groessen).
import itertools
import numpy as np
from a_kern import (su2k, lift, H, WORTE, label_entropien, bin_wahrscheinlichkeiten_batch,
                    schaetzer_batch, entropie_batch)

# Groessen und Bereiche der Pre-Reg v4 (A-1 ... A-7)
PREREG = {"NPH": 10_000, "WDH": 1000, "SIGMA_A": 0.05, "ZIEH": 100_000,
          "SIGMA_GITTER": [round(0.01 * n, 2) for n in range(1, 41)],
          "F_GITTER": [round(0.005 * n, 3) for n in range(1, 81)],
          "SIGMA_BEREICH": (0.02, 0.30), "F_BEREICH": (0.01, 0.30)}


def run_a1(ks, laenge=6):
    a1 = {}; a1_matrix = {}
    for k in ks:
        S, T = su2k(k); F, pdiag, U, L = lift(k); P = np.diag(pdiag)
        Ud = U.conj().T
        a1_matrix[str(k)] = max(float(np.max(np.abs(1j * Ud @ F @ U - S))),
                                float(np.max(np.abs(np.exp(-1j * np.pi / L) * Ud @ P @ U - T))))
        dev = 0.0
        for n in range(laenge + 1):
            for w in itertools.product("FP", repeat=n):
                Mo = np.eye(L, dtype=complex); Ml = np.eye(k + 1, dtype=complex)
                for c in w:
                    Mo = (F if c == "F" else P) @ Mo
                    Ml = (S if c == "F" else T) @ Ml
                for a in range(k + 1):
                    dev = max(dev, abs(H(np.abs(Mo @ U[:, a])**2) - H(np.abs(Ml[:, a])**2) - 1))
        a1[str(k)] = dev
    return a1, a1_matrix


def blind_labels(k):
    h = label_entropien(k)
    return [a for a in range(k + 1) if np.max(np.abs(h[a] - np.round(h[a]))) < 1e-9], h


def fehlalarm_naiv(k, label, sigma, rng, nph, wdh):
    """A2(a): Rate, mit der die NAIVE Regel dieses blinde Label zertifiziert. Je Wiederholung UND je Frame frische
    Phasen (A-4); in einem Frame teilen sich die nph Photonen EINE Phasen-Realisierung; zertifiziert, wenn
    irgendeiner der sieben Frames ausloest."""
    ausgeloest = np.zeros(wdh, bool)
    for w in WORTE:
        p = bin_wahrscheinlichkeiten_batch(k, label, w, sigma, rng, wdh)
        Hh, sd = schaetzer_batch(rng.multinomial(nph, p))
        ausgeloest |= (np.abs(Hh - np.round(Hh)) - 3 * sd > 0)
    return float(ausgeloest.mean())


def delta_max_phase(k, sigma, blind, h, rng, zieh):
    """99,865-%-Quantil von |H_bin - 1 - h| bei unendlich vielen Photonen, Maximum ueber blinde Labels und die vier
    Frames mit Phasenstufe (ohne Phasenstufe kein Phasenrauschen)."""
    m = 0.0
    for a in blind:
        for wi, w in enumerate(WORTE):
            if "T" not in w:
                continue
            p = bin_wahrscheinlichkeiten_batch(k, a, w, sigma, rng, zieh)
            d = np.abs(entropie_batch(p) - 1 - h[a, wi])
            m = max(m, float(np.quantile(d, 0.99865)))
    return m


def stern_aus_gitter(gitter, zert):
    """A-6: letzter Gitterwert VOR dem ersten Bruch; 0.0, wenn schon der erste Gitterwert bricht;
    (None, True), wenn die Bedingung im ganzen Gitter nie bricht (zensiert)."""
    bruch = next((i for i, z in enumerate(zert) if not z), None)
    if bruch is None:
        return None, True
    return (gitter[bruch - 1] if bruch > 0 else 0.0), False


def run_a2_k(k, cfg, seedbasis=0):
    nph, wdh, zieh = cfg["NPH"], cfg["WDH"], cfg["ZIEH"]
    blind, h = blind_labels(k)
    nicht_blind = [a for a in range(k + 1) if a not in blind]
    bester_abstand = max(float(np.max(np.abs(h[a] - np.round(h[a])))) for a in nicht_blind)   # rauschfrei, alle sieben Frames
    eintrag = {"blind": blind, "nicht_blind": nicht_blind, "bester_abstand_rauschfrei": bester_abstand}

    # (a) naive Regel, JEDES blinde Label einzeln (A-3); Kontrollzeile sigma = 0
    rate = {}
    for lab in blind:
        rate[str(lab)] = {
            "sigma_0.05": fehlalarm_naiv(k, lab, cfg["SIGMA_A"], np.random.default_rng(seedbasis + 2000 + 100 * k + lab), nph, wdh),
            "kontrolle_sigma_0": fehlalarm_naiv(k, lab, 0.0, np.random.default_rng(seedbasis + 2500 + 100 * k + lab), nph, wdh)}
    eintrag["a_rate_je_blindes_label"] = rate
    eintrag["a_bestaetigt"] = bool(all(r["sigma_0.05"] >= 0.5 for r in rate.values()))
    eintrag["a_kontrolle_sigma_0_unter_50"] = bool(all(r["kontrolle_sigma_0"] < 0.5 for r in rate.values()))

    # (b) geschuetzte Regel, Phasenrauschen, viele Photonen: rauschfreier bester Abstand gegen Delta_max(sigma)
    rng = np.random.default_rng(seedbasis + 3000 + k)
    kurve_b = {}; zert_b = []
    for s in cfg["SIGMA_GITTER"]:
        dm = delta_max_phase(k, s, blind, h, rng, zieh)
        z = bool(bester_abstand > dm); zert_b.append(z)
        kurve_b[f"{s:.3f}"] = {"delta_max": dm, "zertifizierbar": z}
    s_stern, zens_b = stern_aus_gitter(cfg["SIGMA_GITTER"], zert_b)
    bruch_b = next((i for i, z in enumerate(zert_b) if not z), None)
    lo, hi = cfg["SIGMA_BEREICH"]
    eintrag.update({"b_sigma_stern": s_stern, "b_zensiert_kein_bruch_im_gitter": zens_b,
                    "b_spaeter_wieder_zertifizierbar": bool(bruch_b is not None and any(zert_b[bruch_b:])),
                    "b_kurve": kurve_b, "b_bestaetigt": bool((not zens_b) and lo <= s_stern <= hi)})

    # (c) geschuetzte Regel, Dunkelzaehlungen, viele Photonen: VERRAUSCHTER bester Abstand (alle Frames, alle nicht-blinden
    #     Labels) gegen Delta_max(f) (alle sieben Frames, alle blinden Labels); deterministisch
    L = 2 * (k + 2)
    ideal = {(a, w): bin_wahrscheinlichkeiten_batch(k, a, w, 0.0, None, 1)[0] for a in range(k + 1) for w in WORTE}
    def h_dunkel(a, wi, f):
        return float(H((1 - f) * ideal[(a, WORTE[wi])] + f / L)) - 1
    kurve_c = {}; zert_c = []
    for f in cfg["F_GITTER"]:
        dmf = max(abs(h_dunkel(a, wi, f) - h[a, wi]) for a in blind for wi in range(len(WORTE)))
        vals = [h_dunkel(a, wi, f) for a in nicht_blind for wi in range(len(WORTE))]
        abst = max(abs(x - round(x)) for x in vals)
        z = bool(abst > dmf); zert_c.append(z)
        kurve_c[f"{f:.3f}"] = {"delta_max": dmf, "bester_abstand_verrauscht": abst, "zertifizierbar": z}
    f_stern, zens_c = stern_aus_gitter(cfg["F_GITTER"], zert_c)
    bruch_c = next((i for i, z in enumerate(zert_c) if not z), None)
    lo, hi = cfg["F_BEREICH"]
    eintrag.update({"c_f_stern": f_stern, "c_zensiert_kein_bruch_im_gitter": zens_c,
                    "c_spaeter_wieder_zertifizierbar": bool(bruch_c is not None and any(zert_c[bruch_c:])),
                    "c_kurve": kurve_c, "c_bestaetigt": bool((not zens_c) and lo <= f_stern <= hi)})
    return eintrag
