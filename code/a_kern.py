# a_kern.py -- Talbot-Emulation: gemeinsame Bausteine (eigene Implementierung nach Formeln)
import numpy as np

WORTE = ["", "S", "T", "SS", "ST", "TS", "TT"]   # Woerter der Laenge <= 2; "XY" = erst Y, dann X

def su2k(k):
    N = k + 2
    a = np.arange(k + 1)
    S = np.sqrt(2.0 / N) * np.sin(np.pi * np.outer(a + 1, a + 1) / N)
    T = np.diag(np.exp(2j * np.pi * a * (a + 2) / (4.0 * N)))
    return S, T

def lift(k):
    N = k + 2; L = 2 * N
    x = np.arange(L)
    F = np.exp(-2j * np.pi * np.outer(x, x) / L) / np.sqrt(L)
    pdiag = np.exp(1j * np.pi * x**2 / L)
    U = np.zeros((L, k + 1), complex)
    for m in range(1, N):
        U[m, m - 1] += 1 / np.sqrt(2); U[(-m) % L, m - 1] -= 1 / np.sqrt(2)
    return F, pdiag, U, L

def H(p):
    p = np.asarray(p, float); p = p[p > 1e-15]; return float(-(p * np.log2(p)).sum())

def label_entropien(k):
    """rauschfrei: h_a(M) fuer alle Labels a und die sieben Woerter"""
    S, T = su2k(k)
    out = np.zeros((k + 1, len(WORTE)))
    for wi, w in enumerate(WORTE):
        M = np.eye(k + 1, dtype=complex)
        for c in reversed(w):
            M = (S if c == "S" else T) @ M
        for a in range(k + 1):
            out[a, wi] = H(np.abs(M[:, a])**2)
    return out

def bin_wahrscheinlichkeiten(k, a, wort, sigma, rng):
    """optisch, mit Phasenfehlern: jede P-Stufe bekommt diag(exp(i phi_x)), phi_x ~ N(0, sigma), je Bin und Stufe neu"""
    F, pdiag, U, L = lift(k)
    v = U[:, a].copy()
    for c in reversed(wort):
        if c == "S":
            v = F @ v
        else:
            phi = rng.normal(0.0, sigma, L) if sigma > 0 else np.zeros(L)
            v = pdiag * np.exp(1j * phi) * v
    p = np.abs(v)**2
    return p / p.sum()

def schaetzer(counts):
    """Plug-in-Entropie der Bins mit Miller-Madow-Korrektur, minus 1 Bit; sigma aus der Delta-Methode"""
    n = counts.sum(); q = counts[counts > 0] / n; K = len(q)
    Hb = float(-(q * np.log2(q)).sum()) + (K - 1) / (2 * n * np.log(2))
    var = (float((q * np.log2(q)**2).sum()) - float((q * np.log2(q)).sum())**2) / n
    return Hb - 1.0, float(np.sqrt(max(var, 0.0)))

# ---- vektorisierte Fassungen (Pre-Reg v4, A-7: 10^5 Ziehungen) ----------------------------------------
def entropie_batch(p):
    """H in Bit je Zeile; wie H(), nur ueber die letzte Achse"""
    p = np.asarray(p, float)
    ok = p > 1e-15
    return -(np.where(ok, p, 0.0) * np.log2(np.where(ok, p, 1.0))).sum(axis=-1)

def bin_wahrscheinlichkeiten_batch(k, a, wort, sigma, rng, n, phis=None):
    """n unabhaengige Rausch-Realisierungen auf einmal -> (n, L). Je P-Stufe und je Zeile frische Phasen phi_x ~ N(0, sigma);
    phis (Liste je P-Stufe von (n, L)-Feldern) nur fuer die Kontrolle gegen die skalare Fassung."""
    F, pdiag, U, L = lift(k)
    v = np.tile(U[:, a], (n, 1)).astype(complex)
    stufe = 0
    for c in reversed(wort):
        if c == "S":
            v = v @ F.T                      # F ist symmetrisch
        else:
            if phis is not None:
                phi = phis[stufe]
            elif sigma > 0:
                phi = rng.normal(0.0, sigma, (n, L))
            else:
                phi = np.zeros((n, L))
            stufe += 1
            v = pdiag * np.exp(1j * phi) * v
    p = np.abs(v)**2
    return p / p.sum(axis=1, keepdims=True)

def schaetzer_batch(counts):
    """wie schaetzer(), zeilenweise: K = Zahl der Bins mit Zaehlung > 0 (Pre-Reg A-2)"""
    counts = np.asarray(counts)
    n = counts.sum(axis=1)
    q = counts / n[:, None]
    pos = counts > 0
    K = pos.sum(axis=1)
    lq = np.log2(np.where(pos, q, 1.0))
    Hb = -(q * lq).sum(axis=1) + (K - 1) / (2 * n * np.log(2))
    var = ((q * lq**2).sum(axis=1) - ((q * lq).sum(axis=1))**2) / n
    return Hb - 1.0, np.sqrt(np.maximum(var, 0.0))
