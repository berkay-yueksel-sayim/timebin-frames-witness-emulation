# a2_second_route.py -- independent second route for hypothesis A2 (parts a/b/c), model orders k = 2 and k = 4.
# Usage:  python a2_second_route.py [output_directory]      (default: current directory; writes a2_second_route_result.json)
# Requires: numpy, scipy.  No input files, no network, no absolute paths.
#
# Basis: only the definitions of the pre-registration (frames, model, certification rules, noise models).
# Route: (b) and (c) semi-analytic -- phase noise acts only where a Fourier stage follows a phase stage, and there only through the
# difference of two bin phases (a single random variable); the dark-count background is deterministic. The leading-order
# approximation is printed alongside. (a) is a separate Monte-Carlo implementation (finite photon numbers) plus an asymptotic
# formula for the degenerate frames.
# Where the definition leaves something open, the assumptions A-1 ... A-6 below are stated and, where a number changes, run as variants.
# The keys of the result file are kept exactly as in the reference run so the two result files can be compared entry by entry.
import hashlib, json, math, os, sys, time
from statistics import NormalDist
import numpy as np
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")
OUTDIR = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else os.path.abspath(".")
FRAMES = ["id", "S", "T", "SS", "ST", "TS", "TT"]          # 'XY' = first Y, then X
Q = 0.99865
TOL_BLIND = 1e-9


def ent(p):
    p = np.asarray(p, float)
    q = p[p > 0]
    return float(-(q * np.log2(q)).sum())


def modell(k):
    N = k + 2; L = 2 * N; n = k + 1
    a = np.arange(n)
    S = math.sqrt(2 / N) * np.sin(np.pi * np.outer(a + 1, a + 1) / N)
    T = np.diag(np.exp(2j * np.pi * a * (a + 2) / (4 * N)))
    x = np.arange(L)
    F = np.exp(-2j * np.pi * np.outer(x, x) / L) / math.sqrt(L)
    P = np.exp(1j * np.pi * x.astype(float) ** 2 / L)       # diagonal of the phase stage
    U = np.zeros((L, n), complex)
    for i in range(n):
        m = i + 1
        U[m % L, i] += 1 / math.sqrt(2)
        U[(-m) % L, i] -= 1 / math.sqrt(2)
    return dict(k=k, N=N, L=L, n=n, S=S, T=T, F=F, P=P, U=U)


def label_vert(M, w, a):
    v = np.zeros(M["n"], complex); v[a] = 1
    for c in reversed(w if w != "id" else ""):
        v = (M["S"] if c == "S" else M["T"]) @ v
    return np.abs(v) ** 2


def bin_vert(M, w, a, rng=None, sigma=0.0, phi_fest=None):
    """Bin distribution in frame w. Phase noise: every application of P becomes P*diag(e^{i phi_x}), phi_x ~ N(0, sigma^2) per bin
    and stage. phi_fest: a fixed phase vector for the (single) P application -- used by the one-dimensional route."""
    v = M["U"][:, a].copy()
    for c in reversed(w if w != "id" else ""):
        if c == "S":
            v = M["F"] @ v
        else:
            v = M["P"] * v
            if phi_fest is not None:
                v = v * np.exp(1j * phi_fest)
            elif sigma > 0:
                v = v * np.exp(1j * rng.normal(0, sigma, M["L"]))
    p = np.abs(v) ** 2
    return p / p.sum()


def main():
    t0 = time.time()
    os.makedirs(OUTDIR, exist_ok=True)
    out = {"k": {}}
    print("INDEPENDENT SECOND ROUTE, hypothesis A2 -- basis: definitions of the pre-registration")
    print("ASSUMPTIONS where the definition is open (each run as a variant where it changes the number):")
    print("  A-1  N_ph = 10^4 applies PER FRAME (variant: 10^4 in total, 10^4/7 per frame)")
    print("  A-2  Miller-Madow K = number of bins with count > 0 (variant: K = L)")
    print("  A-3  (a) 'one blind label': every blind label reported separately; for k=4 additionally 'any of them'")
    print("  A-4  (a) per repetition all seven frames are measured anew, phase noise redrawn per repetition and frame")
    print("  A-5  (c) 'best frame' = frame with the largest NOISY distance (variant: best frame without noise)")
    print("  A-6  (b),(c) 'up to sigma*/f*' = largest grid value at which the condition holds; whether it holds for all smaller ones is printed")
    zq = NormalDist().inv_cdf(0.5 + Q / 2)

    for k in (2, 4):
        M = modell(k)
        L, n = M["L"], M["n"]
        print("\n" + "=" * 100 + "\nk = %d  (N = %d, L = %d, labels 0..%d)" % (k, M["N"], L, k))
        # positive control: Lemma W and the A1 core
        eS = np.abs(M["S"] - 1j * M["U"].conj().T @ M["F"] @ M["U"]).max()
        eT = np.abs(M["T"] - np.exp(-1j * np.pi / L) * M["U"].conj().T @ np.diag(M["P"]) @ M["U"]).max()
        a1 = max(abs(ent(bin_vert(M, w, a)) - ent(label_vert(M, w, a)) - 1) for w in FRAMES for a in range(n))
        print("PK+  Lemma W: |S - iU'FU| = %.1e, |T - e^{-i pi/L} U'PU| = %.1e ; A1 core max|H_bin - H_label - 1| = %.1e" % (eS, eT, a1))
        if max(eS, eT, a1) > 1e-12:
            print("ABORT: the model does not satisfy Lemma W / A1"); sys.exit(2)
        h = {(w, a): ent(label_vert(M, w, a)) for w in FRAMES for a in range(n)}
        abst = {a: max(abs(h[w, a] - round(h[w, a])) for w in FRAMES) for a in range(n)}
        bester = {a: max(FRAMES, key=lambda w: abs(h[w, a] - round(h[w, a]))) for a in range(n)}
        blind = [a for a in range(n) if abst[a] < TOL_BLIND]
        nicht = [a for a in range(n) if a not in blind]
        print("Blind labels:", blind, "(definition: k=2 {1}, k=4 {1,3})")
        for a in range(n):
            print("   label %d: h per frame %s | delta_a = %.6f (best frame %s)" % (a, ["%.4f" % h[w, a] for w in FRAMES], abst[a], bester[a]))
        soll = {2: [1], 4: [1, 3]}[k]
        if blind != soll:
            print("ABORT: blind labels differ from the definition"); sys.exit(2)
        rk = {"blind": blind, "delta": {str(a): abst[a] for a in nicht}, "bester_frame": {str(a): bester[a] for a in nicht}}

        # ---------------- (b) phase noise, protected rule, infinitely many photons ----------------
        # Structure: P precedes F only in the word 'ST' (first T, then S); in T, TT, TS the phase stage acts as a diagonal BEFORE the
        # measurement and has no effect. In the word ST, U e_a carries only the bins +-m (m = a+1), so only Delta = phi_m - phi_{-m}
        # ~ N(0, 2 sigma^2) matters.
        def abw_st(a, d):
            phi = np.zeros(L); phi[(a + 1) % L] = d
            return abs(ent(bin_vert(M, "ST", a, phi_fest=phi)) - 1 - h["ST", a])
        # structural control: in all other frames the deviation is exactly 0
        rng = np.random.default_rng(20261003)
        rest = max(abs(ent(bin_vert(M, w, a, rng, 0.3)) - 1 - h[w, a]) for w in FRAMES if w != "ST" for a in blind for _ in range(50))
        print("\n(b) structural control: max|h~ - h| outside ST at sigma = 0.3 (50 draws):", "%.1e" % rest)
        # check monotonicity in |Delta|
        dd = np.linspace(0, 3.2, 3201)
        mono = {a: bool(np.all(np.diff([abw_st(a, d) for d in dd]) >= -1e-15)) for a in blind}
        print("    deviation in the ST frame monotone in |Delta| on [0; 3.2] rad:", mono)
        gitter_s = [round(0.01 * i, 2) for i in range(1, 41)]
        dmax_s, dmax_mc, dmax_fo = [], [], []
        # quantile WITHOUT the monotonicity assumption: |Delta| half-normal with spread sqrt(2)*sigma; weighted grid on [0; 9*sqrt(2)*sigma]
        def quantil_st(a, s, npkt=6001):
            sd = math.sqrt(2) * s
            d = np.linspace(0, 9 * sd, npkt)
            w = np.exp(-0.5 * (d / sd) ** 2); w[0] *= 0.5; w[-1] *= 0.5; w = w / w.sum()
            g = np.array([abw_st(a, x) for x in d])
            o = np.argsort(g)
            cdf = np.cumsum(w[o])
            return float(g[o][np.searchsorted(cdf, Q)])
        for s in gitter_s:
            dq = zq * math.sqrt(2) * s
            dmax_s.append(max(quantil_st(a, s) for a in blind))
            # leading order: the zero bins n0 each receive Delta^2/(2L) -> n0*(D2/2L)*log2(2L/D2) + B*D2
            fo = 0.0
            for a in blind:
                p0 = bin_vert(M, "ST", a)
                n0 = int((p0 < 1e-15).sum())
                D2 = dq ** 2
                y = np.arange(L); c = 2 * np.pi * (a + 1) * y / M["N"]
                nz = p0 > 1e-15
                B = -(np.log2(p0[nz]) * np.cos(c[nz]) / (2 * L)).sum() - ((1 + np.cos(c[nz])) / (2 * L * math.log(2))).sum()
                fo = max(fo, n0 * D2 / (2 * L) * math.log2(2 * L / D2) + B * D2)
            dmax_fo.append(fo)
        # counter-route following the LITERAL definition: >= 2000 draws, independent phi per bin and stage, all frames
        for s in (0.02, 0.05, 0.10, 0.20, 0.30):
            werte = []
            for a in blind:
                for w in FRAMES:
                    z = [abs(ent(bin_vert(M, w, a, rng, s)) - 1 - h[w, a]) for _ in range(4000)]
                    werte.append(float(np.quantile(z, Q)))
            dmax_mc.append((s, max(werte)))
        print("    Delta_max(sigma): semi-analytic | leading order | literal MC (4000 draws)")
        for s, dm, fo in zip(gitter_s, dmax_s, dmax_fo):
            mc = [v for t, v in dmax_mc if abs(t - s) < 1e-9]
            if s in (0.01, 0.02, 0.03, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40) or mc:
                print("     sigma=%.2f  %.6f | %.6f | %s" % (s, dm, fo, ("%.6f" % mc[0]) if mc else "-"))
        best_delta = max(abst[a] for a in nicht)
        gilt = [best_delta > dm for dm in dmax_s]
        s_stern = max([s for s, g in zip(gitter_s, gilt) if g], default=None)
        luecke = s_stern is not None and not all(gilt[:gitter_s.index(s_stern) + 1])
        print("    best non-blind distance max delta_a = %.6f (label %s)" % (best_delta, [a for a in nicht if abst[a] == best_delta]))
        zensiert = s_stern == gitter_s[-1]
        print("    sigma* = %s rad  | condition for all smaller grid values: %s | prediction A2(b): sigma* in [0.02 ; 0.30]  -> %s"
              % (s_stern, not luecke, "WITHIN RANGE" if s_stern is not None and 0.02 <= s_stern <= 0.30 else "OUTSIDE"))
        if zensiert:
            print("    CENSORED: the condition never breaks anywhere on the grid -- max Delta_max(sigma) = %.6f at sigma = %.2f lies below delta = %.6f"
                  % (max(dmax_s), gitter_s[int(np.argmax(dmax_s))], best_delta))
        rk["b"] = {"sigma_stern": s_stern, "zensiert": zensiert, "luecke": luecke, "delta_best": best_delta,
                   "dmax": dict(zip(map(str, gitter_s), dmax_s)), "dmax_fuehrend": dict(zip(map(str, gitter_s), dmax_fo)),
                   "dmax_mc_woertlich": {str(s): v for s, v in dmax_mc}, "monoton": {str(a): v for a, v in mono.items()}}

        # ---------------- (c) dark counts, protected rule, deterministic ----------------
        gitter_f = [round(0.005 * i, 3) for i in range(1, 81)]
        pb = {(w, a): bin_vert(M, w, a) for w in FRAMES for a in range(n)}
        def htil(w, a, f):
            return ent((1 - f) * pb[w, a] + f / L) - 1
        dmax_f, gilt1, gilt2, fo_f = [], [], [], []
        for f in gitter_f:
            dm = max(abs(htil(w, a, f) - h[w, a]) for w in FRAMES for a in blind)
            dmax_f.append(dm)
            v1 = max(abs(htil(w, a, f) - round(htil(w, a, f))) for w in FRAMES for a in nicht)
            v2 = max(abs(htil(bester[a], a, f) - round(htil(bester[a], a, f))) for a in nicht)
            gilt1.append(v1 > dm); gilt2.append(v2 > dm)
            # leading order in the frame with the most zero bins: n0*(f/L)*log2(L/f) + f*[...]
            fo = 0.0
            for w in FRAMES:
                for a in blind:
                    p = pb[w, a]; nz = p > 1e-15; n0 = int((~nz).sum())
                    lin = ((p[nz] - 1 / L) * np.log2(p[nz])).sum() + (n0 / L) / math.log(2)
                    fo = max(fo, abs(n0 * f / L * math.log2(L / f) + f * lin))
            fo_f.append(fo)
        def stern(g):
            s = max([f for f, x in zip(gitter_f, g) if x], default=None)
            return s, (s is not None and all(g[:gitter_f.index(s) + 1]))
        f1, ok1 = stern(gilt1); f2, ok2 = stern(gilt2)
        print("\n(c) Delta_max(f) exact | leading order (zero bins + linear):")
        for f, dm, fo in zip(gitter_f, dmax_f, fo_f):
            if f in (0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4):
                print("     f=%.3f  %.6f | %.6f" % (f, dm, fo))
        print("    f* (A-5: noisy best frame) = %s | throughout: %s" % (f1, ok1))
        print("    f* (variant: best frame without noise) = %s | throughout: %s" % (f2, ok2))
        print("    prediction A2(c): f* in [0.01 ; 0.30] -> %s" % ("WITHIN RANGE" if f1 is not None and 0.01 <= f1 <= 0.30 else "OUTSIDE"))
        rk["c"] = {"f_stern_A5": f1, "durchgehend_A5": ok1, "f_stern_variante": f2, "durchgehend_variante": ok2,
                   "dmax": dict(zip(map(str, gitter_f), dmax_f))}

        # ---------------- (a) naive rule, phase noise sigma = 0.05, N_ph = 10^4, 1000 repetitions ----------------
        def lauf(sigma, nph_frame, k_modus, seed):
            r = np.random.default_rng(seed)
            treffer = {a: 0 for a in blind}; irgend = 0; je_frame = {(a, w): 0 for a in blind for w in FRAMES}
            for _ in range(1000):
                zert_irgend = False
                for a in blind:
                    zert = False
                    for w in FRAMES:
                        p = bin_vert(M, w, a, r, sigma)
                        c = r.multinomial(nph_frame, p)
                        ph = c / nph_frame
                        K = int((c > 0).sum()) if k_modus == "beob" else L
                        Hp = ent(ph)
                        Hh = Hp + (K - 1) / (2 * nph_frame * math.log(2)) - 1
                        q = ph[ph > 0]
                        var = ((q * np.log2(q) ** 2).sum() - Hp ** 2) / nph_frame
                        sh = math.sqrt(max(var, 0.0))
                        if abs(Hh - round(Hh)) - 3 * sh > 0:
                            zert = True; je_frame[a, w] += 1
                    treffer[a] += zert; zert_irgend |= zert
                irgend += zert_irgend
            return treffer, irgend, je_frame
        print("\n(a) naive rule -- share of the 1000 repetitions with a false certification (prediction: >= 50 %)")
        varianten = [("A-1/A-2 wie angenommen", 0.05, 10000, "beob"), ("ohne Phasenrauschen (sigma=0)", 0.0, 10000, "beob"),
                     ("Variante K = L", 0.05, 10000, "L"), ("Variante N_ph insgesamt (10^4/7 je Frame)", 0.05, 10000 // 7, "beob")]
        anzeige = ["A-1/A-2 as assumed", "without phase noise (sigma=0)", "variant K = L", "variant N_ph in total (10^4/7 per frame)"]
        rk["a"] = {}
        for i, (name, s, nph, km) in enumerate(varianten):
            tr, irg, jf = lauf(s, nph, km, 7000 + 100 * k + i)
            je = {a: tr[a] / 1000 for a in blind}
            print("   %-44s per label %s | any blind label %.3f | per frame %s" % (
                anzeige[i], {a: "%.3f" % v for a, v in je.items()}, irg / 1000,
                {"%d:%s" % (a, w): jf[a, w] for a in blind for w in FRAMES if jf[a, w]}))
            rk["a"][name] = {"je_label": {str(a): v for a, v in je.items()}, "irgendein": irg / 1000,
                             "je_frame": {"%d:%s" % (a, w): jf[a, w] for a in blind for w in FRAMES}}
        # asymptotic formula for the degenerate frames (bin distribution uniform over K bins, no noise):
        # Hhat - h = (K-1-X)/(2N ln2), 3*sigmahat = 3*sqrt(X)/(N ln2), X ~ chi2_{K-1}  -> certified if |K-1-X| > 6 sqrt(X)
        def p_entartet(Kb):
            # |K-1-X| > 6 sqrt(X): for X < (sqrt(K+8)-3)^2 or X > (3+sqrt(K+8))^2
            lo = (math.sqrt(Kb - 1 + 9) - 3) ** 2; hi = (3 + math.sqrt(Kb - 1 + 9)) ** 2
            return stats.chi2.cdf(lo, Kb - 1) + stats.chi2.sf(hi, Kb - 1)
        for a in blind:
            fr = {w: int((pb[w, a] > 1e-15).sum()) for w in FRAMES}
            gleich = {w: bool(np.allclose(pb[w, a][pb[w, a] > 1e-15], 1 / fr[w])) for w in FRAMES}
            pf = {w: p_entartet(fr[w]) for w in FRAMES if gleich[w]}
            ges = 1 - np.prod([1 - v for v in pf.values()])
            print("   asymptotic, label %d (degenerate frames only, no phase noise): per frame %s -> at least one %.3f"
                  % (a, {w: "K=%d:%.3f" % (fr[w], v) for w, v in pf.items()}, ges))
            rk["a"]["asymptotisch_label_%d" % a] = {"je_frame": {w: v for w, v in pf.items()}, "gesamt": float(ges)}
        out["k"][str(k)] = rk

    out["sekunden"] = round(time.time() - t0, 1)
    pfad = os.path.join(OUTDIR, "a2_second_route_result.json")
    json.dump(out, open(pfad, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    print("\nwritten: a2_second_route_result.json | md5 of the result entries %s | %.0f s"
          % (hashlib.md5(json.dumps(out["k"], sort_keys=True).encode("utf-8")).hexdigest()[:8], out["sekunden"]))
    print("SCOPE: definitions of the pre-registration as stated for k = 2 and k = 4; if a definition changes, the numbers must be recomputed.")


main()
