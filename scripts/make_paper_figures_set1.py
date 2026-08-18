import os
import csv
import glob
import math
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import TwoSlopeNorm

RESULTS = "results"
OUTDIR = os.path.join("figures", "paper_BP_set1")
os.makedirs(OUTDIR, exist_ok=True)

DEFECT = (4, 3)

def locate(candidates):
    for p in candidates:
        if p and os.path.exists(p):
            return p
    return None

def g(pattern):
    xs = glob.glob(pattern)
    return xs[0] if xs else None

def read_csv(path):
    if path is None or not os.path.exists(path):
        return []
    with open(path, newline="") as f:
        return list(csv.DictReader(f))

def fval(x, default=np.nan):
    try:
        return float(x)
    except Exception:
        return default

def savefig(fig, name):
    png = os.path.join(OUTDIR, name + ".png")
    pdf = os.path.join(OUTDIR, name + ".pdf")
    fig.tight_layout()
    fig.savefig(png, dpi=220, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return png, pdf

def load_global(path):
    rows = read_csv(path)
    times = [fval(r["time"]) for r in rows]
    delta = [fval(r["delta_flux"]) for r in rows]
    return rows, np.array(times), np.array(delta)

def common_xy(path_ref, path_other):
    rows1 = read_csv(path_ref)
    rows2 = read_csv(path_other)
    m1 = {round(fval(r["time"]), 8): fval(r["delta_flux"]) for r in rows1}
    m2 = {round(fval(r["time"]), 8): fval(r["delta_flux"]) for r in rows2}
    ts = sorted(set(m1.keys()) & set(m2.keys()))
    return np.array(ts), np.array([m1[t] for t in ts]), np.array([m2[t] for t in ts])

baseline_dir = locate([
    os.path.join(RESULTS, "paper_K2p2_set1_light", "K2p2_D16_dt005_Lie")
])
d8_dir = locate([
    os.path.join(RESULTS, "paper_K2p2_set1_light", "K2p2_D8_dt005_Lie")
])
d24_dir = locate([
    os.path.join(RESULTS, "paper_K2p2_set1_light", "K2p2_D24_dt005_Lie")
])
dt010_dir = locate([
    os.path.join(RESULTS, "paper_K2p2_set1_light", "K2p2_D16_dt010_Lie")
])
dt0025_dir = locate([
    os.path.join(RESULTS, "paper_K2p2_set1_light", "K2p2_D16_dt0025_Lie")
])
strang005_dir = locate([
    os.path.join(RESULTS, "paper_K2p2_set1_light", "K2p2_D16_dt005_Strang")
])
strang0025_dir = locate([
    os.path.join(RESULTS, "paper_K2p2_set1_light", "K2p2_D16_dt0025_Strang")
])

baseline_global = os.path.join(baseline_dir, "global.csv") if baseline_dir else None
baseline_links = os.path.join(baseline_dir, "links.csv") if baseline_dir else None
baseline_profiles = os.path.join(baseline_dir, "profiles.csv") if baseline_dir else None

d8_global = os.path.join(d8_dir, "global.csv") if d8_dir else None
d24_global = os.path.join(d24_dir, "global.csv") if d24_dir else None
dt010_global = os.path.join(dt010_dir, "global.csv") if dt010_dir else None
dt0025_global = os.path.join(dt0025_dir, "global.csv") if dt0025_dir else None
strang005_global = os.path.join(strang005_dir, "global.csv") if strang005_dir else None
strang0025_global = os.path.join(strang0025_dir, "global.csv") if strang0025_dir else None

localX_csv = locate([
    os.path.join(RESULTS, "light_qpu_X_ZZ_K2p2_8x6", "local_X.csv")
])
zz_csv = locate([
    os.path.join(RESULTS, "light_qpu_X_ZZ_K2p2_8x6", "electric_strings_ZZ.csv")
])
xx_csv = locate([
    os.path.join(RESULTS, "light_qpu_XX_K2p2_8x6", "connected_XX.csv")
])

wilson_base_csv = locate([
    os.path.join(RESULTS, "light_wilson_K2p2_8x6", "wilson_averages.csv")
])
wilson_d24_csv = locate([
    os.path.join(RESULTS, "light_wilson_K2p2_8x6_D24_Lie_dt005", "wilson_averages.csv")
])
wilson_strang_csv = locate([
    os.path.join(RESULTS, "light_wilson_K2p2_8x6_D16_Strang_dt0025", "wilson_averages.csv")
])

flux_valid_csv = locate([
    os.path.join(RESULTS, "fast_exact_vs_bp_4x4_K2p2_dt005", "global.csv"),
    g(os.path.join(RESULTS, "fast_exact_vs_bp_4x4_K2p2_dt005*", "global.csv"))
])
qpu_valid_csv = locate([
    os.path.join(RESULTS, "exact_vs_bp_qpuobs_4x4.csv")
])
wilson_valid_csv = locate([
    os.path.join(RESULTS, "exact_vs_bp_wilson_4x4", "wilson_exact_vs_bp.csv")
])

made = []
summary_lines = []

# ------------------------------------------------------------------
# Figure 1: Delta flux baseline + controls
# ------------------------------------------------------------------
if baseline_global and os.path.exists(baseline_global):
    _, t0, d0 = load_global(baseline_global)
    fig = plt.figure(figsize=(7.2, 4.8))
    plt.plot(t0, d0, marker="o", linewidth=2.2, label="Lie dt=0.05 D=16 (baseline)")
    if d8_global and os.path.exists(d8_global):
        _, t, d = load_global(d8_global)
        plt.plot(t, d, marker=".", linestyle="--", label="Lie dt=0.05 D=8")
    if d24_global and os.path.exists(d24_global):
        _, t, d = load_global(d24_global)
        plt.plot(t, d, marker=".", linestyle="--", label="Lie dt=0.05 D=24")
    if dt010_global and os.path.exists(dt010_global):
        _, t, d = load_global(dt010_global)
        plt.plot(t, d, marker="s", linestyle=":", label="Lie dt=0.10 D=16")
    if dt0025_global and os.path.exists(dt0025_global):
        _, t, d = load_global(dt0025_global)
        plt.plot(t, d, marker="^", linestyle="-.", label="Lie dt=0.025 D=16")
    if strang005_global and os.path.exists(strang005_global):
        _, t, d = load_global(strang005_global)
        plt.plot(t, d, marker="d", linestyle=":", label="Strang dt=0.05 D=16")
    if strang0025_global and os.path.exists(strang0025_global):
        _, t, d = load_global(strang0025_global)
        plt.plot(t, d, marker="v", linestyle="-.", label="Strang dt=0.025 D=16")
    plt.axvline(0.4, linestyle="--", linewidth=1.0)
    plt.xlabel("time")
    plt.ylabel(r"$\Delta \Phi(t)$")
    plt.title(r"Differential electric flux, $K/\Gamma=2.2$")
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=8, ncol=2)
    made += list(savefig(fig, "01_delta_flux_controls"))

# ------------------------------------------------------------------
# Figure 2: control differences vs baseline
# ------------------------------------------------------------------
if baseline_global and os.path.exists(baseline_global):
    fig = plt.figure(figsize=(7.2, 4.8))
    anyline = False
    for label, path in [
        ("D=8", d8_global),
        ("D=24", d24_global),
        ("Lie dt=0.10", dt010_global),
        ("Lie dt=0.025", dt0025_global),
        ("Strang dt=0.05", strang005_global),
        ("Strang dt=0.025", strang0025_global),
    ]:
        if path and os.path.exists(path):
            t, yref, yoth = common_xy(baseline_global, path)
            if len(t) > 0:
                plt.plot(t, np.abs(yoth - yref), marker="o", label=label)
                anyline = True
    if anyline:
        plt.axvline(0.4, linestyle="--", linewidth=1.0)
        plt.xlabel("time")
        plt.ylabel(r"$|\Delta\Phi_{\rm control}-\Delta\Phi_{\rm baseline}|$")
        plt.title("Control-channel deviation from baseline")
        plt.grid(True, alpha=0.3)
        plt.legend(fontsize=8)
        made += list(savefig(fig, "02_control_differences"))
    else:
        plt.close(fig)

# ------------------------------------------------------------------
# Figure 3: spatial heatmaps
# ------------------------------------------------------------------
if baseline_links and os.path.exists(baseline_links):
    rows = read_csv(baseline_links)
    bytime = defaultdict(list)
    allvals = []
    xmin = ymin = 1e9
    xmax = ymax = -1e9
    for r in rows:
        t = round(fval(r["time"]), 8)
        x1 = fval(r["v1_x"]); y1 = fval(r["v1_y"])
        x2 = fval(r["v2_x"]); y2 = fval(r["v2_y"])
        midx = fval(r["mid_x"]); midy = fval(r["mid_y"])
        dv = fval(r["delta_flux"])
        bytime[t].append((x1, y1, x2, y2, midx, midy, dv))
        allvals.append(dv)
        xmin = min(xmin, x1, x2, midx)
        xmax = max(xmax, x1, x2, midx)
        ymin = min(ymin, y1, y2, midy)
        ymax = max(ymax, y1, y2, midy)
    avail = sorted(bytime.keys())
    targets = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
    chosen = [min(avail, key=lambda x: abs(x - t)) for t in targets]
    amax = max(abs(min(allvals)), abs(max(allvals)))
    norm = TwoSlopeNorm(vmin=-amax, vcenter=0.0, vmax=amax)

    fig, axes = plt.subplots(2, 3, figsize=(11.5, 6.5))
    axes = axes.flatten()
    for ax, t in zip(axes, chosen):
        segs = [[(a, b), (c, d)] for (a, b, c, d, mx, my, dv) in bytime[t]]
        vals = [dv for (_, _, _, _, _, _, dv) in bytime[t]]
        lc = LineCollection(segs, cmap="coolwarm", norm=norm, linewidths=4)
        lc.set_array(np.array(vals))
        ax.add_collection(lc)
        ax.scatter([DEFECT[0]], [DEFECT[1]], marker="*", s=80)
        ax.set_xlim(xmin - 0.5, xmax + 0.5)
        ax.set_ylim(ymin - 0.5, ymax + 0.5)
        ax.set_aspect("equal")
        ax.set_title(f"t = {t:.2f}")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    cbar = fig.colorbar(lc, ax=axes.tolist(), shrink=0.9)
    cbar.set_label(r"$\Delta \phi_\ell$")
    fig.suptitle("Spatial differential electric flux")
    made += list(savefig(fig, "03_spatial_heatmaps"))

# ------------------------------------------------------------------
# Figure 4: radial profiles
# ------------------------------------------------------------------
if baseline_profiles and os.path.exists(baseline_profiles):
    rows = read_csv(baseline_profiles)
    bytime = defaultdict(list)
    for r in rows:
        t = round(fval(r["time"]), 8)
        rb = fval(r["radius_bin"])
        mad = fval(r["mean_abs_delta_flux"])
        md = fval(r["mean_delta_flux"])
        bytime[t].append((rb, md, mad))
    fig = plt.figure(figsize=(7.0, 4.8))
    for t in [0.1, 0.2, 0.3, 0.4, 0.5]:
        tt = min(bytime.keys(), key=lambda x: abs(x - t))
        pts = sorted(bytime[tt], key=lambda z: z[0])
        r = [p[0] for p in pts]
        y = [p[2] for p in pts]
        plt.plot(r, y, marker="o", label=f"t={tt:.2f}")
    plt.xlabel("radius bin")
    plt.ylabel(r"mean $|\Delta \phi_\ell|$")
    plt.title("Radial profile of differential flux")
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=8)
    made += list(savefig(fig, "04_radial_profiles"))

# ------------------------------------------------------------------
# Figure 5: local X shell averages
# ------------------------------------------------------------------
if localX_csv and os.path.exists(localX_csv):
    rows = read_csv(localX_csv)
    shell = defaultdict(lambda: defaultdict(list))
    for r in rows:
        t = round(fval(r["time"]), 8)
        x = int(fval(r["x"]))
        y = int(fval(r["y"]))
        dx = x - DEFECT[0]
        dy = y - DEFECT[1]
        rr = round(math.hypot(dx, dy), 3)
        shell[rr][t].append(fval(r["delta_x"]))
    fig = plt.figure(figsize=(7.0, 4.8))
    for rr in sorted(shell.keys())[:6]:
        ts = sorted(shell[rr].keys())
        vals = [np.mean(shell[rr][t]) for t in ts]
        plt.plot(ts, vals, marker="o", label=f"r={rr:g}")
    plt.xlabel("time")
    plt.ylabel(r"shell-averaged $\Delta\langle X\rangle$")
    plt.title("Local magnetic response")
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=8, ncol=2)
    made += list(savefig(fig, "05_localX_shells"))

# ------------------------------------------------------------------
# Figure 6: electric strings ZZ
# ------------------------------------------------------------------
if zz_csv and os.path.exists(zz_csv):
    rows = read_csv(zz_csv)
    bydist = defaultdict(lambda: defaultdict(list))
    for r in rows:
        d = int(fval(r["distance"]))
        t = round(fval(r["time"]), 8)
        bydist[d][t].append(fval(r["delta_zz"]))
    fig = plt.figure(figsize=(7.0, 4.8))
    for d in sorted(bydist.keys()):
        ts = sorted(bydist[d].keys())
        vals = [np.mean(bydist[d][t]) for t in ts]
        plt.plot(ts, vals, marker="o", label=f"d={d}")
    plt.xlabel("time")
    plt.ylabel(r"$\Delta\langle Z_{r_0}Z_r\rangle$")
    plt.title("Short electric strings")
    plt.grid(True, alpha=0.3)
    plt.legend()
    made += list(savefig(fig, "06_strings_ZZ"))

# ------------------------------------------------------------------
# Figure 7: connected XX
# ------------------------------------------------------------------
if xx_csv and os.path.exists(xx_csv):
    rows = read_csv(xx_csv)
    bydist = defaultdict(lambda: defaultdict(list))
    for r in rows:
        d = int(fval(r["distance"]))
        t = round(fval(r["time"]), 8)
        bydist[d][t].append(fval(r["delta_connected"]))
    fig = plt.figure(figsize=(7.0, 4.8))
    for d in sorted(bydist.keys()):
        ts = sorted(bydist[d].keys())
        vals = [np.mean(bydist[d][t]) for t in ts]
        plt.plot(ts, vals, marker="o", label=f"d={d}")
    plt.xlabel("time")
    plt.ylabel(r"$\Delta C^{XX}_{0r}$")
    plt.title("Connected magnetic correlators")
    plt.grid(True, alpha=0.3)
    plt.legend()
    made += list(savefig(fig, "07_connected_XX"))

# ------------------------------------------------------------------
# Figure 8: Wilson loops baseline
# ------------------------------------------------------------------
if wilson_base_csv and os.path.exists(wilson_base_csv):
    rows = read_csv(wilson_base_csv)
    byshape = defaultdict(list)
    for r in rows:
        byshape[r["shape"]].append((fval(r["time"]), fval(r["delta_mean"])))
    fig = plt.figure(figsize=(7.0, 4.8))
    for shape in ["1x1", "2x1", "1x2", "2x2"]:
        pts = sorted(byshape[shape], key=lambda z: z[0])
        plt.plot([p[0] for p in pts], [p[1] for p in pts], marker="o", label=shape)
    plt.xlabel("time")
    plt.ylabel(r"$\Delta W$")
    plt.title("Wilson loops (dual-area X products)")
    plt.grid(True, alpha=0.3)
    plt.legend(title="shape")
    made += list(savefig(fig, "08_wilson_loops"))

# ------------------------------------------------------------------
# Figure 9: Wilson control differences
# ------------------------------------------------------------------
if wilson_base_csv and os.path.exists(wilson_base_csv):
    base = read_csv(wilson_base_csv)
    def wmap(rows):
        m = {}
        for r in rows:
            m[(r["shape"], round(fval(r["time"]), 8))] = fval(r["delta_mean"])
        return m
    m0 = wmap(base)
    fig = plt.figure(figsize=(7.0, 4.8))
    drawn = False
    for label, path in [("D24 Lie dt=0.05", wilson_d24_csv), ("D16 Strang dt=0.025", wilson_strang_csv)]:
        if path and os.path.exists(path):
            m = wmap(read_csv(path))
            times = sorted(set(t for (_, t) in m0.keys()) & set(t for (_, t) in m.keys()))
            vals = []
            keep_times = []
            for t in times:
                diffs = []
                for shape in ["1x1", "2x1", "1x2", "2x2"]:
                    if (shape, t) in m0 and (shape, t) in m:
                        diffs.append(abs(m[(shape, t)] - m0[(shape, t)]))
                if diffs:
                    keep_times.append(t)
                    vals.append(np.mean(diffs))
            if vals:
                plt.plot(keep_times, vals, marker="o", label=label)
                drawn = True
    if drawn:
        plt.xlabel("time")
        plt.ylabel("mean absolute Wilson deviation")
        plt.title("Wilson-loop control-channel deviation")
        plt.grid(True, alpha=0.3)
        plt.legend(fontsize=8)
        made += list(savefig(fig, "09_wilson_control_deviation"))
    else:
        plt.close(fig)

# ------------------------------------------------------------------
# Figure 10: exact-vs-BP validation summary
# ------------------------------------------------------------------
fig = plt.figure(figsize=(7.2, 4.8))
anything = False

if flux_valid_csv and os.path.exists(flux_valid_csv):
    rows = read_csv(flux_valid_csv)
    t = [fval(r["time"]) for r in rows]
    e = [fval(r["abs_error"]) for r in rows]
    plt.plot(t, e, marker="o", label=r"Flux $\Delta\Phi$ error")
    anything = True

if qpu_valid_csv and os.path.exists(qpu_valid_csv):
    rows = read_csv(qpu_valid_csv)
    fam = defaultdict(lambda: defaultdict(list))
    for r in rows:
        t = round(fval(r["time"]), 8)
        fam[r["observable"]][t].append(fval(r["abs_error"]))
    label_map = {"X": "Local X error", "ZZ": "String ZZ error", "CXX": "Connected XX error"}
    for obs in ["X", "ZZ", "CXX"]:
        if obs in fam:
            ts = sorted(fam[obs].keys())
            vals = [max(fam[obs][tt]) for tt in ts]
            plt.plot(ts, vals, marker="o", label=label_map.get(obs, obs))
            anything = True

if wilson_valid_csv and os.path.exists(wilson_valid_csv):
    rows = read_csv(wilson_valid_csv)
    bytime = defaultdict(list)
    for r in rows:
        t = round(fval(r["time"]), 8)
        bytime[t].append(fval(r["delta_abs_error"]))
    ts = sorted(bytime.keys())
    vals = [max(bytime[t]) for t in ts]
    plt.plot(ts, vals, marker="o", label="Wilson-loop max error")
    anything = True

if anything:
    plt.axvline(0.4, linestyle="--", linewidth=1.0)
    plt.xlabel("time")
    plt.ylabel("absolute exact-vs-BP error")
    plt.title("Exact-vs-BP validation summary on 4x4")
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=8)
    made += list(savefig(fig, "10_validation_summary"))
else:
    plt.close(fig)

# ------------------------------------------------------------------
# Figure 11: baseline diagnostics summary text plot
# ------------------------------------------------------------------
summary_lines.append("Frozen BP baseline:")
summary_lines.append("  K/Gamma = 2.2, D = 16, dt = 0.05, Lie-Trotter, 8x6 dual lattice")
summary_lines.append("Controlled window:")
summary_lines.append("  t <= 0.4 strictly controlled; 0.4 < t <= 0.5 shown with explicit BP uncertainty")
summary_lines.append("")
summary_lines.append("Files detected:")
for name, path in [
    ("baseline_global", baseline_global),
    ("baseline_links", baseline_links),
    ("baseline_profiles", baseline_profiles),
    ("local_X", localX_csv),
    ("strings_ZZ", zz_csv),
    ("connected_XX", xx_csv),
    ("wilson_base", wilson_base_csv),
    ("flux_validation", flux_valid_csv),
    ("qpuobs_validation", qpu_valid_csv),
    ("wilson_validation", wilson_valid_csv),
]:
    summary_lines.append(f"  {name}: {'OK' if path and os.path.exists(path) else 'MISSING'}")

summary_lines.append("")
summary_lines.append("Generated figure files:")
for p in made:
    summary_lines.append("  " + os.path.basename(p))

with open(os.path.join(OUTDIR, "figure_manifest.txt"), "w") as f:
    f.write("\n".join(summary_lines) + "\n")

print("Created figure set in:", OUTDIR)
for p in made:
    print(" -", p)
print(" -", os.path.join(OUTDIR, "figure_manifest.txt"))
