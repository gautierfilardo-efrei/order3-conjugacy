"""Render the five figures of the manuscript from the recorded outputs.
Inputs (JSON, from the order3-conjugacy / CPWI repositories): algo_results.json, solver_variants.json,
round_algebra_results.json, wreath_summary.json, round_solver_profile.json.
Usage: python make_figures.py [input_dir] [output_dir]   (defaults: data/ and . ; writes fig_*.pdf and fig_*.png)
Figure width = 5.2 in, the text width of the elsarticle 12pt preprint layout, so fonts are not scaled."""
import json, math, os, sys
import numpy as np, matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch

HERE = os.path.dirname(os.path.abspath(__file__))
IN = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "data"); OUT = sys.argv[2] if len(sys.argv) > 2 else HERE; os.makedirs(OUT, exist_ok=True)
W = 5.2; FS = 7.0; GREY = "#7f7f7f"
mpl.rcParams.update({"font.size": FS, "axes.labelsize": FS, "axes.titlesize": FS+0.5, "xtick.labelsize": FS-0.5,
    "ytick.labelsize": FS-0.5, "legend.fontsize": FS-0.5, "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "savefig.dpi": 300, "axes.titlelocation": "left", "axes.titleweight": "normal"})
load = lambda n: json.load(open(os.path.join(IN, n)))
A = load("algo_results.json"); SV = load("solver_variants.json"); RA = load("round_algebra_results.json")
S = load("wreath_summary.json"); P = load("round_solver_profile.json")
def letter(ax, l): ax.text(-0.12, 1.06, l, transform=ax.transAxes, fontsize=FS+2, fontweight="bold", va="bottom")
def save(fig, name):
    fig.savefig(os.path.join(OUT, name+".pdf"), bbox_inches="tight"); fig.savefig(os.path.join(OUT, name+".png"), bbox_inches="tight"); plt.close(fig)

# ---------- Figure: statistics of the solution count
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(W, 2.3))
ks = np.arange(0, 7); pois = np.exp(-1)/np.array([math.factorial(k) for k in ks])
cols = ["#c7dcef", "#8fbbdb", "#4a90c4", "#1f5f99"]
for i, N in enumerate((6, 7, 8, 9)):
    h = A["partA_exhaustive"][str(N)]["nsol_hist_all_pairs"] if N < 9 else A["partA_N9"]["nsol_hist_all_pairs"]
    tot = sum(h.values()); p = [h.get(str(k), 0)/tot for k in ks]
    ax1.plot(ks+(i-1.5)*0.08, p, "o", ms=2.8, color=cols[i], label=f"all pairs in Sym({N})")
ax1.plot(ks, pois, "_", ms=9, mew=1.3, color="k", label="Poisson(1)")
for ax in (ax1, ax2):
    ax.set_yscale("log"); ax.set_ylim(2e-4, 1); ax.set_yticks([1e-3, 1e-2, 1e-1, 1]); ax.set_yticklabels(["0.001", "0.01", "0.1", "1"]); ax.margins(x=0.06)
ax1.set_xlabel("number of solutions $x$ of $x\\,a\\,x^2=c$"); ax1.set_ylabel("fraction of pairs $(a,c)$")
ax1.set_title("Near Poisson(1), heavier tail"); ax1.legend(frameon=False, loc="lower left", borderaxespad=0.2)
C = RA["C"]["m5_all_3600_AB"]
pr = [C["pmf_random_target"].get(str(k), 0) for k in ks]; pp = [C["pmf_planted_target"].get(str(k), 0) for k in ks]
sb = [k*np.exp(-1)/math.factorial(k) for k in ks]
ax2.plot(ks-0.06, pr, "o", ms=3, color="#4a90c4", label="random target, $m=5$")
ax2.plot(ks+0.06, pp, "s", ms=3, color="#d62728", label="planted target $c=x_0ax_0^2$")
ax2.plot(ks, pois, "_", ms=9, mew=1.3, color="#4a90c4", alpha=.8, label="Poisson(1)"); ax2.plot(ks, sb, "_", ms=9, mew=1.3, color="#d62728", alpha=.8, label="size-biased Poisson(1)")
ax2.set_xlabel("number of solutions"); ax2.set_ylabel("fraction of targets"); ax2.set_title("Planted targets: mean 2 solutions"); ax2.legend(frameon=False, loc="lower left", handlelength=1.2, borderaxespad=0.2)
letter(ax1, "a"); letter(ax2, "b"); fig.tight_layout(); save(fig, "fig_statistics")

# ---------- Figure: cost law and peeling depth
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(W, 2.4), gridspec_kw=dict(width_ratios=[1.35, 1]))
col = {"first": "#1f77b4", "mcf": "#d62728"}; lab = {"first": "first free point", "mcf": "most constrained first"}; mm = np.arange(6, 37)
for h in ("first", "mcf"):
    med = S["fits"][h]["medians"]; ms = sorted(int(k) for k in med); ys = [np.log2(med[str(m)]) for m in ms]
    ax1.plot(ms, ys, "o", ms=3, color=col[h], label=lab[h]); L = S["fits"][h]["linear"]; M = S["fits"][h]["mlogm"]
    ax1.plot(mm, M["alpha"]*mm*np.log2(mm)+M["beta"], "-", lw=1, color=col[h], alpha=.85); ax1.plot(mm, L["slope"]*mm+L["intercept"], "--", lw=0.9, color=col[h], alpha=.6)
ax1.plot([], [], "-", color="k", lw=1, label=r"fit $\alpha\, m\log_2 m+\beta$"); ax1.plot([], [], "--", color="k", lw=0.9, label=r"fit $s\,m+i$")
ax1.axvline(36, color=GREY, lw=0.8, ls=":"); ax1.text(35.3, 2, "m = 36", ha="right", va="bottom", fontsize=FS-1, color=GREY)
ax1.set_xlabel(r"degree $m$ of $\mathrm{Alt}(m)$"); ax1.set_ylabel(r"$\log_2$ median search nodes")
ax1.set_title("Fits agree on data, diverge at m = 36"); ax1.legend(frameon=False, loc="upper left"); ax1.set_ylim(0, 75); ax1.margins(x=0.04)
for name, c, lb in (("first", "#1f77b4", "first free point"), ("chains", "#9467bd", "most partially known chains"), ("mcf", "#d62728", "most constrained first")):
    pm = P[name]["per_m"]; ms = sorted(int(k) for k in pm)
    for m in ms: ax2.plot([m]*len(pm[str(m)]["planted_depth"]), pm[str(m)]["planted_depth"], ".", color=c, alpha=.22, ms=3.5)
    ax2.plot(ms, [pm[str(m)]["k_mean"] for m in ms], "-o", ms=2.8, color=c, lw=1, label=lb)
mm2 = np.arange(5, 17); ax2.plot(mm2, 0.38*mm2+0.6, "--", color=GREY, lw=0.8)
ax2.set_xlabel(r"degree $m$ of $\mathrm{Alt}(m)$"); ax2.set_ylabel("free decisions $k(m)$")
ax2.set_title(r"$k(m)\approx 0.38\,m+0.6$"); ax2.legend(frameon=False, loc="upper left"); ax2.margins(0.06); ax2.set_xticks([6,8,10,12,14,16])
letter(ax1, "a"); letter(ax2, "b"); fig.tight_layout(); save(fig, "fig_cost_law")

# ---------- Figure: solver variants
lab = {"a_baseline_first": "first free point", "b1_mcf_lookahead": "most constrained first", "b2_chains": "most partially known chains",
       "c_first_ctype": "first free point + cycle-type test", "f_reference_solve_round": "coupled $(U,V)$ solver", "e_z3_onehot": "SMT (one-hot SAT), seconds"}
col = {"a_baseline_first": "#1f77b4", "b1_mcf_lookahead": "#d62728", "b2_chains": "#9467bd", "c_first_ctype": "#7f7f7f", "f_reference_solve_round": "#2c2c2c", "e_z3_onehot": "#e377c2"}
mk = {"a_baseline_first": "o", "b1_mcf_lookahead": "o", "b2_chains": "o", "c_first_ctype": "x", "f_reference_solve_round": "^"}
fig, ax = plt.subplots(figsize=(W, 3.0)); axb = ax.twinx()
for v in ("f_reference_solve_round", "a_baseline_first", "c_first_ctype", "b2_chains", "b1_mcf_lookahead"):
    rows = [s for s in SV["summary"] if s["variant"] == v and not s["exhausted_or_timeout"]]
    ax.plot([s["m"] for s in rows], [s["log2_median"] for s in rows], mk[v], ms=3.2 if mk[v] != "x" else 4, color=col[v], label=lab[v], mfc="none" if v == "f_reference_solve_round" else col[v])
rows = [s for s in SV["summary"] if s["variant"] == "e_z3_onehot" and not s["exhausted_or_timeout"]]
axb.plot([s["m"] for s in rows], [np.log2(s["median_seconds"]) for s in rows], "s", ms=3, color=col["e_z3_onehot"], alpha=.8, label=lab["e_z3_onehot"])
ax.set_xlabel(r"degree $m$ of $\mathrm{Alt}(m)$"); ax.set_ylabel(r"$\log_2$ median search nodes to first solution")
axb.set_ylabel(r"$\log_2$ median seconds (SMT, one machine)", color=col["e_z3_onehot"]); axb.tick_params(axis="y", colors=col["e_z3_onehot"])
axb.spines["right"].set_visible(True); axb.spines["right"].set_color(col["e_z3_onehot"]); axb.spines["top"].set_visible(False)
ax.set_title("The branching rule changes the slope; the cycle-type test changes nothing")
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = axb.get_legend_handles_labels(); ax.legend(h1+h2, l1+l2, frameon=False, loc="upper left"); ax.margins(0.05); ax.set_ylim(0, 26)
fig.tight_layout(); save(fig, "fig_variants")

# ---------- Figure: the reduction (schematic)
fig, ax = plt.subplots(figsize=(W, 2.3)); ax.set_xlim(0, 10); ax.set_ylim(-0.75, 3.65); ax.axis("off"); L = 5; R = 0.45; FS = 6.3
blocks = [(1.6, 1.75, "w", "$e_w = s(w)$", "#4a90c4"), (5.0, 1.75, "x", "$e_x = s(x)+3B$", "#d62728"), (8.4, 1.75, "y", "$e_y = s(y)+T-4B$", "#9467bd")]
for cx, cy, name, lb, c in blocks:
    ax.add_patch(Circle((cx, cy), R, fc="white", ec=c, lw=1.3))
    for k in range(L):
        th = 2*np.pi*k/L+np.pi/2; ax.plot(cx+0.32*np.cos(th), cy+0.32*np.sin(th), ".", color=c, ms=3.5)
    ax.annotate("", xy=(cx+0.32*np.cos(np.pi/2+2*np.pi/L), cy+0.32*np.sin(np.pi/2+2*np.pi/L)), xytext=(cx+0.32*np.cos(np.pi/2), cy+0.32*np.sin(np.pi/2)), arrowprops=dict(arrowstyle="->", color=c, lw=0.7, connectionstyle="arc3,rad=-0.4"))
    ax.text(cx, cy-R-0.1, f"block $B_{{{name}}}$ ($L$-cycle of $a$)", ha="center", va="top", fontsize=FS-0.5)
    ax.text(cx, cy+R+0.08, lb, ha="center", va="bottom", fontsize=FS-0.5, color=c)
for (x1, y1), (x2, y2), rad in (((1.6+R, 1.75), (5.0-R, 1.75), -0.3), ((5.0+R, 1.75), (8.4-R, 1.75), -0.3), ((8.4, 0.85), (1.6, 0.85), -0.3)):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), connectionstyle=f"arc3,rad={rad}", arrowstyle="->", mutation_scale=8, lw=1.0, color="#2c2c2c"))
ax.text(3.3, 2.35, r"$x:\ z\mapsto m_x z+t_w$", ha="center", va="bottom", fontsize=FS-0.5); ax.text(6.7, 2.35, r"$z\mapsto m_y z+t_x$", ha="center", va="bottom", fontsize=FS-0.5)
ax.text(5.0, 0.15, r"$z\mapsto m_w z+t_y$", ha="center", va="top", fontsize=FS-0.5)
ax.text(5.0, 3.6, r"$a$: $z\mapsto z+1$ on every block;  $c$: $z\mapsto z+m_u$ on $B_u$, $m_u=g^{e_u}$;  $x^3=1$ forces $x$ to permute blocks by $\sigma$, $\sigma^3=1$", ha="center", va="top", fontsize=FS-0.5)
ax.text(5.0, -0.55, r"$x^3=\mathrm{id}$ on the triple iff $m_w m_x m_y\equiv1$ iff $e_w+e_x+e_y\equiv 0\ (\mathrm{mod}\ T)$ iff $s(w)+s(x)+s(y)=B$", ha="center", va="top", fontsize=FS-0.5)
save(fig, "fig_reduction")

# ---------- Figure: the wreath normal form (schematic)
fig, ax = plt.subplots(figsize=(W, 1.8)); ax.set_xlim(0, 10); ax.set_ylim(0, 3.2); ax.axis("off"); f = 5.9
def blk(x0, y0, w, h, label, c):
    ax.add_patch(FancyBboxPatch((x0, y0), w, h, boxstyle="round,pad=0.02,rounding_size=0.08", fc=c, ec="#2c2c2c", lw=0.7)); ax.text(x0+w/2, y0+h/2, label, ha="center", va="center", fontsize=f)
ax.text(0.1, 3.05, "round of the layered walk", fontsize=f, va="top")
blk(0.1, 1.6, 1.1, 0.6, "$U$", "#dbe9f6"); blk(0.1, 0.7, 1.1, 0.6, "$V$", "#f6dbdb")
ax.annotate("", xy=(2.6, 1.9), xytext=(1.25, 1.9), arrowprops=dict(arrowstyle="->", lw=0.8)); ax.annotate("", xy=(2.6, 1.0), xytext=(1.25, 1.0), arrowprops=dict(arrowstyle="->", lw=0.8))
ax.text(1.92, 2.0, "$UAVU$", ha="center", va="bottom", fontsize=f); ax.text(1.92, 1.1, "$VBUV$", ha="center", va="bottom", fontsize=f)
blk(2.65, 1.6, 1.1, 0.6, "$U'$", "#dbe9f6"); blk(2.65, 0.7, 1.1, 0.6, "$V'$", "#f6dbdb")
ax.text(1.9, 0.15, r"unknowns $(U,V)\in\mathrm{Alt}(m)^2$, constants $A,B$", ha="center", fontsize=f)
ax.annotate("", xy=(5.3, 1.45), xytext=(4.3, 1.45), arrowprops=dict(arrowstyle="-|>", lw=1.1, color="#2c2c2c")); ax.text(4.8, 1.6, "Theorem", ha="center", va="bottom", fontsize=f)
ax.text(5.4, 3.05, r"one-variable equation in $\mathrm{Alt}(m)\wr C_2\subseteq\mathrm{Sym}(2m)$", fontsize=f, va="top")
ax.text(6.4, 2.3, r"block 1 $=\{0,\dots,m-1\}$, block 2 $=\{m,\dots,2m-1\}$", ha="center", fontsize=f)
blk(5.4, 1.6, 2.0, 0.6, "block 1: $U$ acts", "#dbe9f6"); blk(5.4, 0.7, 2.0, 0.6, "block 2: $V$ acts", "#f6dbdb")
ax.add_patch(FancyArrowPatch((7.5, 1.9), (7.5, 1.0), connectionstyle="arc3,rad=-0.9", arrowstyle="<->", mutation_scale=8, lw=0.9, color="#2c2c2c")); ax.text(7.95, 1.45, r"$\sigma$", fontsize=f, va="center")
ax.text(8.45, 2.05, r"$x=(U,V)\,\sigma$", fontsize=f, va="center"); ax.text(8.45, 1.45, r"$a=(B,A)$", fontsize=f, va="center"); ax.text(8.45, 0.85, r"$c=(U',V')\,\sigma$", fontsize=f, va="center")
ax.text(7.4, 0.15, r"$(U,V)$ solves the round  $\Longleftrightarrow$  $x\,a\,x^{2}=c$ in $\mathrm{Alt}(m)\wr C_2$", ha="center", fontsize=f)
save(fig, "fig_wreath")
print("wrote", sorted(os.listdir(OUT)))
