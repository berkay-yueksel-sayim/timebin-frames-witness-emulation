# a_fig1_hero.py -- Talbot-Emulation: Abbildung 1 des Manuskripts.
# (a) Abbildung Label -> Zeit-Bins am Beispiel k = 2, Eingangs-Label a = 0, Frame S: jede Label-Wahrscheinlichkeit wird halbiert
#     auf die Bins +-(b+1) verteilt, also H_bin = h + 1 (rechnet mit a_kern, prueft die Identitaet). Direkte Rechnung, kein Lauf.
# (b) Phasenfehler: Rand Delta_max(sigma) gegen den rauschfreien besten Abstand, k = 2 und 4 (aus Lauf 1, b_kurve).
# (c) Untergrund: Rand Delta_max(f) gegen den um den Untergrund verschobenen besten Abstand, k = 2 und 4 (aus Lauf 1, c_kurve).
#     Die gestrichelten Kurven enden am ersten Gitterpunkt, an dem sie den Rand erreichen: dahinter kommt der groesste Abstand aus
#     Frames, in denen jedes Label ganzzahlige Entropie hat, und liegt deshalb genau auf dem Rand (kein Signal mehr).
# Schreibt a_fig1_hero.pdf und a_fig1_hero.png; prueft sigma*, f*, Saettigung und Abstaende gegen Lauf 1.
import hashlib, json, pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from a_kern import su2k, bin_wahrscheinlichkeiten, H

HIER = pathlib.Path(__file__).parent
P = HIER / "a_tests_results_LAUF1_2026-10-03.json"
assert hashlib.md5(P.read_bytes()).hexdigest() == "aba6a738ddc6bd4c4196b1f1fab7146d", "nicht Lauf 1"
A2 = json.loads(P.read_text(encoding="utf-8"))["A2"]

# ---------- (a) Beispiel ----------
k, a = 2, 0
S, T = su2k(k)
p_label = np.abs(S[:, a]) ** 2
p_bin = bin_wahrscheinlichkeiten(k, a, "S", 0.0, np.random.default_rng(0))
L = len(p_bin)
h_label, h_bin = H(p_label), H(p_bin)
erwartet = np.zeros(L)
for b in range(k + 1):
    erwartet[b + 1] += p_label[b] / 2; erwartet[(-(b + 1)) % L] += p_label[b] / 2
print(f"(a) k={k}, Label {a}, Frame S: p_label {np.round(p_label, 4)} -> h = {h_label:.4f} bit; "
      f"p_bin {np.round(p_bin, 4)} -> H_bin = {h_bin:.4f} bit; H_bin - h - 1 = {h_bin - h_label - 1:.1e} "
      f"| Bins wie +-(b+1)-Regel: {np.allclose(p_bin, erwartet, atol=1e-12)}")

# ---------- (b), (c) aus Lauf 1 ----------
kurven = {}
for kk in ("2", "4"):
    e = A2[kk]
    sb = np.array(sorted(float(s) for s in e["b_kurve"]))
    db = np.array([e["b_kurve"][f"{s:.3f}"]["delta_max"] for s in sb])
    fc = np.array(sorted(float(s) for s in e["c_kurve"]))
    dc = np.array([e["c_kurve"][f"{s:.3f}"]["delta_max"] for s in fc])
    bc = np.array([e["c_kurve"][f"{s:.3f}"]["bester_abstand_verrauscht"] for s in fc])
    schnitt = int(np.argmax(dc >= bc))                        # erster Gitterpunkt, an dem der Rand das Signal erreicht
    kurven[kk] = dict(sb=sb, db=db, best=e["bester_abstand_rauschfrei"], fc=fc, dc=dc, bc=bc, schnitt=schnitt,
                      sstern=e["b_sigma_stern"], fstern=e["c_f_stern"])
    zert_b = [e["b_kurve"][f"{s:.3f}"]["zertifizierbar"] for s in sb]
    print(f"(b/c) k={kk}: bester Abstand {e['bester_abstand_rauschfrei']:.4f}; sigma* {e['b_sigma_stern']}; "
          f"Delta_max(sigma) max {db.max():.4f}; letztes zertifizierbares sigma {sb[np.array(zert_b)].max() if any(zert_b) else None}; "
          f"f* {e['c_f_stern']}; erstes f mit Rand >= Signal {fc[schnitt]} (gestrichelt bis dort); "
          f"spaeter wieder zertifizierbar: {e['c_spaeter_wieder_zertifizierbar']}")

# ---------- Abbildung ----------
plt.rcParams.update({"font.family": "STIXGeneral", "mathtext.fontset": "stix", "font.size": 8,
                     "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "xtick.major.size": 2.5, "ytick.major.size": 2.5, "axes.edgecolor": "#898781",
                     "xtick.color": "#52514e", "ytick.color": "#52514e", "axes.labelcolor": "#0b0b0b"})
BLAU, ORANGE, RAND, TEXT1, TEXT2 = "#2a78d6", "#eb6834", "#898781", "#0b0b0b", "#52514e"
FARBE = {"2": BLAU, "4": ORANGE}
TEXTFARBE = {"2": "#1c5cab", "4": "#b9430f"}        # dunklere Stufen derselben Farbfamilie fuer Beschriftungen (Kontrast)
fig = plt.figure(figsize=(7.0, 2.75))
TITEL_Y = 0.965

# (a) zwei kleine Achsen mit GLEICHER y-Skala: Labels oben, Bins unten; Label b und seine Bins +-(b+1) in derselben Stufe
STUFEN = ["#86b6ef", "#2a78d6", "#104281"]                      # Ordinal-Rampe, validiert gegen Weiss
farbe_bin = ["white"] * L
for b in range(k + 1):
    farbe_bin[b + 1] = farbe_bin[(-(b + 1)) % L] = STUFEN[b]
xa0 = 0.045
aa = fig.add_axes([xa0, 0.56, 0.10, 0.25])
ab = fig.add_axes([xa0, 0.15, 0.19, 0.25])
aa.bar(np.arange(k + 1), p_label, width=0.62, color=STUFEN[:k + 1])
aa.set_xticks(range(k + 1)); aa.set_xlim(-0.6, k + 0.6)
aa.set_xlabel("label $b$", labelpad=1)
aa.text(1.12, 0.66, f"$h = {h_label:.1f}$ bit", transform=aa.transAxes, va="center", ha="left", color=TEXT1)
aa.text(1.12, 0.32, "frame $S$, $k=2$,\ninput label $a=0$", transform=aa.transAxes, va="center", ha="left",
        color=TEXT2, fontsize=7)
ab.bar(np.arange(L), p_bin, width=0.62, color=farbe_bin)
ab.set_xticks(range(L)); ab.set_xlim(-0.6, L - 0.4)
ab.set_xlabel(f"time bin $x$ ($L={L}$)", labelpad=1)
ab.text(1.04, 0.66, "$H_{\\rm bin} = h + 1$", transform=ab.transAxes, va="center", ha="left", color=TEXT1)
ab.text(1.04, 0.36, f"$= {h_bin:.1f}$ bit", transform=ab.transAxes, va="center", ha="left", color=TEXT1)
ab.text(0.5, 1.04, "shade = label $b$, split onto bins $\\pm(b+1)$", transform=ab.transAxes, va="bottom", ha="center",
        color=TEXT2, fontsize=7)
for ax in (aa, ab):
    ax.set_ylim(0, 0.6); ax.set_yticks([0, 0.25, 0.5]); ax.set_yticklabels(["0", "1/4", "1/2"])
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
fig.text(0.012, TITEL_Y, "(a)  each label $\\to$ two half-weight bins", fontsize=8.5, ha="left", va="top")

# gemeinsamer Linienschluessel fuer (b) und (c), neutral grau (die Farben = k sind direkt beschriftet)
fig.legend(handles=[Line2D([0], [0], color=TEXT2, lw=1.3), Line2D([0], [0], color=TEXT2, lw=0.9, ls=(0, (4, 2.5)))],
           labels=["margin $\\Delta_{\\max}$ (blind labels)", "signal (best non-blind label)"], loc="center",
           bbox_to_anchor=(0.665, 0.875), ncol=2, frameon=False, fontsize=7, handlelength=2.2, columnspacing=1.6,
           borderaxespad=0, borderpad=0)

# (b) Phasenfehler
gb = fig.add_axes([0.375, 0.15, 0.27, 0.64])
for kk in ("2", "4"):
    c = kurven[kk]
    gb.plot(c["sb"], c["db"], color=FARBE[kk], lw=1.3)
    gb.axhline(c["best"], color=FARBE[kk], lw=0.9, ls=(0, (4, 2.5)))
    gb.text(0.405, c["best"] + 0.012, f"$k={kk}$ signal", color=TEXTFARBE[kk], fontsize=7, ha="right", va="bottom")
c2 = kurven["2"]
y_s = float(np.interp(c2["sstern"], c2["sb"], c2["db"]))
gb.plot([c2["sstern"]], [y_s], "o", ms=3.2, color=BLAU, mec="white", mew=0.6, zorder=5)
gb.annotate(f"$\\sigma^*={c2['sstern']:.2f}$ rad\n(last certified grid point)", xy=(c2["sstern"], y_s),
            xytext=(0.19, 0.345), color=TEXT2, fontsize=7, va="center",
            arrowprops=dict(arrowstyle="-", color=RAND, lw=0.5, relpos=(0.0, 0.5)))
gb.text(0.405, kurven["4"]["db"].max() + 0.012, f"$k=4$ margin saturates at {kurven['4']['db'].max():.3f}", color=TEXT2,
        fontsize=7, ha="right", va="bottom")
gb.set_xlim(0, 0.41); gb.set_ylim(0, 0.58); gb.set_yticks([0, 0.25, 0.5])
gb.set_xlabel("phase error $\\sigma$ (rad)"); gb.set_ylabel("bit")
for s in ("top", "right"):
    gb.spines[s].set_visible(False)
fig.text(0.335, TITEL_Y, "(b)  phase errors", fontsize=8.5, ha="left", va="top")

# (c) Untergrund: gestrichelt nur bis zum ersten Schnitt mit dem Rand
gc = fig.add_axes([0.715, 0.15, 0.255, 0.64])
ziel = {"2": (0.058, 0.17), "4": (0.058, 0.07)}           # Beschriftungen getauscht: Leitlinien kreuzen kein Label
for kk in ("2", "4"):
    c = kurven[kk]
    gc.plot(c["fc"], c["dc"], color=FARBE[kk], lw=1.3)
    m = c["schnitt"] + 1
    gc.plot(c["fc"][:m], c["bc"][:m], color=FARBE[kk], lw=0.9, ls=(0, (4, 2.5)))
    y = float(np.interp(c["fstern"], c["fc"], c["dc"]))
    gc.plot([c["fstern"]], [y], "o", ms=3.2, color=FARBE[kk], mec="white", mew=0.6, zorder=5)
    gc.annotate(f"$k={kk}$: $f^*={c['fstern']:g}$", xy=(c["fstern"], y), xytext=ziel[kk], color=TEXTFARBE[kk], fontsize=7,
                va="center", arrowprops=dict(arrowstyle="-", color=RAND, lw=0.5, relpos=(0.0, 0.5)))
gc.set_xlim(0, 0.1); gc.set_ylim(0, 0.58); gc.set_yticks([0, 0.25, 0.5])
gc.set_xticks([0, 0.025, 0.05, 0.075, 0.1]); gc.set_xticklabels(["0", "0.025", "0.05", "0.075", "0.1"])
gc.set_xlabel("background fraction $f$")
for s in ("top", "right"):
    gc.spines[s].set_visible(False)
fig.text(0.675, TITEL_Y, "(c)  background counts", fontsize=8.5, ha="left", va="top")

fig.savefig("a_fig1_hero.pdf", metadata={"CreationDate": None})   # ohne Datum: byte-gleich reproduzierbar
fig.savefig("a_fig1_hero.png", dpi=300)
print("geschrieben: a_fig1_hero.pdf, a_fig1_hero.png")
