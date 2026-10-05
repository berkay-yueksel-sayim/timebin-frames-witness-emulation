# a_tabelle_woerter.py -- Talbot-Emulation: Werte fuer die Bauteil-Tabelle (Pre-Reg §3, nur Tabelle, keine Hypothese).
# Je Frame-Wort: Zahl der Fourier- und Phasenstufen, ob eine Phasenstufe VOR einer Fourierstufe liegt (nur dann wirkt
# Phasenrauschen auf die Intensitaeten), und je Label h_a(M) mit Abstand zur naechsten ganzen Zahl. k = 2 und 4.
from a_kern import WORTE, label_entropien
import numpy as np

for k in (2, 4):
    h = label_entropien(k)
    L = 2 * (k + 2)
    print(f"k={k}  (L={L} Bins)")
    for wi, w in enumerate(WORTE):
        nF, nP = w.count("S"), w.count("T")
        # Anwendungsreihenfolge: "XY" = erst Y, dann X -> von rechts nach links lesen
        folge = w[::-1]
        p_vor_f = any(folge[i] == "T" and "S" in folge[i + 1:] for i in range(len(folge)))
        werte = " ".join(f"{h[a, wi]:.4f}" for a in range(k + 1))
        abst = max(abs(h[a, wi] - round(h[a, wi])) for a in range(k + 1))
        print(f"  {w or 'id':3s} F={nF} P={nP} Phase-vor-Fourier={str(p_vor_f):5s} h_a: {werte}   max Abstand {abst:.4f}")
    nb = [a for a in range(k + 1) if np.max(np.abs(h[a] - np.round(h[a]))) > 1e-9]
    erstes = {a: next(w or 'id' for wi, w in enumerate(WORTE) if abs(h[a, wi] - round(h[a, wi])) > 1e-9) for a in nb}
    print(f"  nicht-blind: {nb}; erstes Wort, das sie zeigt: {erstes}")
