# a_nk.py -- Talbot-Emulation: Negativ- und Positivkontrolle des WERKZEUGS (keine pre-registrierte Hypothese).
# Frage: erkennt die Validierung einen Fehler? Mit der GERADEN Einbettung (cos statt sin) darf weder
# Lemma W noch die Auslese-Identitaet H_bin = h + 1 gelten. Nur k = 2, 4 (in Pre-Reg §0 schon bekannt);
# die Test-k {5,7,8,9,11,12} werden hier NICHT angefasst.
import itertools
import numpy as np
from a_kern import su2k, lift, H


def einbettung(k, gerade):
    N = k + 2; L = 2 * N
    U = np.zeros((L, k + 1), complex)
    for m in range(1, N):
        U[m, m - 1] += 1 / np.sqrt(2)
        U[(-m) % L, m - 1] += (1 if gerade else -1) / np.sqrt(2)
    return U


def pruefe(k, gerade):
    S, T = su2k(k); F, pdiag, _, L = lift(k); P = np.diag(pdiag)
    U = einbettung(k, gerade)
    lemma_w = float(np.max(np.abs(1j * U.conj().T @ F @ U - S)))
    a1 = 0.0
    for n in range(1, 5):
        for w in itertools.product("FP", repeat=n):
            Mo = np.eye(L, dtype=complex); Ml = np.eye(k + 1, dtype=complex)
            for c in w:
                Mo = (F if c == "F" else P) @ Mo
                Ml = (S if c == "F" else T) @ Ml
            for a in range(k + 1):
                a1 = max(a1, abs(H(np.abs(Mo @ U[:, a])**2) - H(np.abs(Ml[:, a])**2) - 1))
    return lemma_w, a1


ok = True
for k in (2, 4):
    pk = pruefe(k, gerade=False)
    nk = pruefe(k, gerade=True)
    print(f"k={k}  PK ungerade: Lemma W {pk[0]:.2e}, A1 {pk[1]:.2e}   |   NK gerade: Lemma W {nk[0]:.2e}, A1 {nk[1]:.2e}")
    ok &= pk[0] < 1e-12 and pk[1] < 1e-12 and nk[0] > 1e-2 and nk[1] > 1e-2
print("Kontrollen:", "BESTANDEN (PK haelt, NK scheitert sichtbar)" if ok else "NICHT bestanden")
