# a_paare_wortlaenge.py -- Talbot-Emulation: kuerzestes Frame-Wort, das ein NICHT-blindes Paar (i, jbar) als nicht-ganzzahlig zeigt.
# Breitensuche ueber Woerter in S, T nach Laenge; je Paar das erste Wort mit Abstand > 1e-6 und dieser Abstand.
import itertools
import numpy as np
from a_kern import su2k, H

for k in (2, 4):
    S, T = su2k(k); S = S.astype(complex); n = k + 1
    offen = {(i, j) for i in range(n) for j in range(n)}
    blind_bekannt = {2: {(1, 1)}, 4: {(1, 1), (1, 3), (3, 1), (3, 3)}}[k]   # aus a_paare.py (Gruppe)
    offen -= blind_bekannt
    gefunden = {}
    for L in range(0, 9):
        for w in itertools.product("ST", repeat=L):
            M = np.eye(n, dtype=complex)
            for c in reversed(w):
                M = (S if c == "S" else T) @ M
            h = [H(np.abs(M[:, a])**2) for a in range(n)]
            for (i, j) in list(offen):
                d = abs(h[i] + h[j] - round(h[i] + h[j]))
                if d > 1e-6:
                    gefunden[(i, j)] = ("".join(w) or "id", L, round(d, 4))
                    offen.discard((i, j))
        if not offen:
            break
    laengste = max(v[1] for v in gefunden.values())
    print(f"k={k}: laengstes noetiges Wort = {laengste}; offen geblieben: {sorted(offen)}")
    for p in sorted(gefunden, key=lambda q: -gefunden[q][1])[:6]:
        print(f"   Paar {p}: erstes Wort {gefunden[p][0]!r} (Laenge {gefunden[p][1]}), Abstand {gefunden[p][2]}")
