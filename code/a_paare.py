# a_paare.py -- Talbot-Emulation: blinde Labels und blinde Produkt-Paare (i, jbar) fuer k = 2 und k = 4, EIGENE direkte Rechnung.
# Methode: Rahmengruppe = Abschluss von <S, T> bis auf globale Phase (BFS, Schluessel nach Phasen-Normierung auf die erste
# Komponente mit |u| > 1e-6, gerundet). Fuer die doppelte Theorie gilt M_D = M (x) Mbar (Lemma F, Mueger 2003), also
# Paarwert = h_i(M) + h_j(M) (|Mbar e_j|^2 = |M e_j|^2). Blind = in JEDEM Frame ganzzahlig (Toleranz 1e-9);
# zur Begruendung der Toleranz wird der kleinste Abstand zur naechsten ganzen Zahl unter den NICHT-blinden Faellen ausgegeben.
import itertools
import numpy as np
from a_kern import su2k, H, WORTE

TOL = 1e-9

def schluessel(M):
    v = M.flatten()
    i = int(np.argmax(np.abs(v) > 1e-6))
    w = v * np.conj(v[i]) / abs(v[i])
    return tuple(np.round(np.concatenate([w.real, w.imag]), 8))

def rahmengruppe(S, T):
    E = np.eye(S.shape[0], dtype=complex)
    gesehen = {schluessel(E): E}; rand = [E]
    while rand:
        neu = []
        for M in rand:
            for G in (S, T):
                X = G @ M
                s = schluessel(X)
                if s not in gesehen:
                    gesehen[s] = X; neu.append(X)
        rand = neu
        if len(gesehen) > 100000: raise RuntimeError("keine Abgeschlossenheit")
    return list(gesehen.values())

def entropien(Ms, n):
    return np.array([[H(np.abs(M[:, a])**2) for a in range(n)] for M in Ms])

def wortmatrizen(S, T):
    out = []
    for w in WORTE:
        M = np.eye(S.shape[0], dtype=complex)
        for c in reversed(w):
            M = (S if c == "S" else T) @ M
        out.append(M)
    return out

def abstand(x):
    return np.abs(x - np.round(x))


for k in (2, 4):
    S, T = su2k(k); n = k + 1
    G = rahmengruppe(S.astype(complex), T)
    hG = entropien(G, n)                          # |G| x n
    hW = entropien(wortmatrizen(S.astype(complex), T), n)
    einzeln_G = [a for a in range(n) if abstand(hG[:, a]).max() < TOL]
    einzeln_W = [a for a in range(n) if abstand(hW[:, a]).max() < TOL]
    paare_G, paare_W, luecke_G = [], [], []
    for i, j in itertools.product(range(n), repeat=2):
        dG = abstand(hG[:, i] + hG[:, j]).max()
        dW = abstand(hW[:, i] + hW[:, j]).max()
        (paare_G if dG < TOL else luecke_G).append((i, j) if dG < TOL else dG)
        if dW < TOL: paare_W.append((i, j))
    print(f"k={k}: |Rahmengruppe| bis auf Phase = {len(G)}")
    print(f"   blinde Einzel-Labels: Gruppe {einzeln_G} | sieben Woerter {einzeln_W}")
    print(f"   blinde Paare (i, jbar): Gruppe {paare_G}")
    print(f"                            sieben Woerter {paare_W}")
    print(f"   kleinster Abstand zur ganzen Zahl unter NICHT-blinden Paaren (Gruppe): {min(luecke_G):.4f}  (Toleranz {TOL})")
