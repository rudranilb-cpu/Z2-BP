import csv
import math
import os
from collections import defaultdict

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import TwoSlopeNorm

BASE = "results/fast_spatial_8x6"
GLOBAL_CSV = os.path.join(BASE, "global.csv")
LINKS_CSV = os.path.join(BASE, "links.csv")
PROFILES_CSV = os.path.join(BASE, "profiles.csv")
OUTDIR = os.path.join("figures", "fast_spatial_8x6")
os.makedirs(OUTDIR, exist_ok=True)

TARGET_TIMES = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]

def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))

global_rows = read_csv(GLOBAL_CSV)
link_rows = read_csv(LINKS_CSV)
profile_rows = read_csv(PROFILES_CSV)

# ---------- helpers ----------
def nearest_time(target, available):
    return min(available, key=lambda t: abs(t - target))

# ---------- parse global ----------
g_time = []
g_delta = []
g_front_rms = []
g_front_q80 = []
for r in global_rows:
    g_time.append(float(r["time"]))
    g_delta.append(float(r["delta_flux"]))
    g_front_rms.append(float(r["front_rms"]))
    g_front_q80.append(float(r["front_q80"]))

# ---------- parse links ----------
links_by_time = defaultdict(list)
all_times = set()
all_delta = []
xmin = 1e9
xmax = -1e9
ymin = 1e9
ymax = -1e9

for r in link_rows:
    t = round(float(r["time"]), 8)
    all_times.add(t)
    x1 = float(r["v1_x"])
    y1 = float(r["v1_y"])
    x2 = float(r["v2_x"])
    y2 = float(r["v2_y"])
    mx = float(r["mid_x"])
    my = float(r["mid_y"])
    delta = float(r["delta_flux"])
    links_by_time[t].append((x1, y1, x2, y2, mx, my, delta))
    all_delta.append(delta)
    xmin = min(xmin, x1, x2, mx)
    xmax = max(xmax, x1, x2, mx)
    ymin = min(ymin, y1, y2, my)
    ymax = max(ymax, y1, y2, my)

available_times = sorted(all_times)
amax = max(abs(min(all_delta)), abs(max(all_delta)))
norm = TwoSlopeNorm(vmin=-amax, vcenter=0.0, vmax=amax)

# ---------- parse profiles ----------
profiles_by_time = defaultdict(list)
profile_times = set()
for r in profile_rows:
    t = round(float(r["time"]), 8)
    rb = float(r["radius_bin"])
    md = float(r["mean_delta_flux"])
    mad = float(r["mean_abs_delta_flux"])
    n = int(r["nlinks"])
    profiles_by_time[t].append((rb, md, mad, n))
    profile_times.add(t)

profile_times = sorted(profile_times)

# ---------- plot 1: total differential flux ----------
plt.figure(figsize=(6.5, 4.5))
plt.plot(g_time, g_delta, marker="o")
plt.xlabel("time")
plt.ylabel(r"$\Delta \Phi(t)$")
plt.title("8x6 BP run: total differential flux")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "global_delta_flux.png"), dpi=180)
plt.close()

# ---------- plot 2: front metrics ----------
plt.figure(figsize=(6.5, 4.5))
plt.plot(g_time, g_front_rms, marker="o", label="front_rms")
plt.plot(g_time, g_front_q80, marker="s", label="front_q80")
plt.xlabel("time")
plt.ylabel("front measure")
plt.title("8x6 BP run: front measures")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "front_measures.png"), dpi=180)
plt.close()

# ---------- plot 3: heatmaps at requested times ----------
chosen_times = [nearest_time(t, available_times) for t in TARGET_TIMES]

for t in chosen_times:
    rows = links_by_time[t]
    segs = [[(x1, y1), (x2, y2)] for (x1, y1, x2, y2, mx, my, delta) in rows]
    vals = [delta for (_, _, _, _, _, _, delta) in rows]

    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    lc = LineCollection(segs, cmap="coolwarm", norm=norm, linewidths=4)
    lc.set_array(vals)
    ax.add_collection(lc)

    # also mark link midpoints lightly
    mxs = [mx for (_, _, _, _, mx, my, delta) in rows]
    mys = [my for (_, _, _, _, mx, my, delta) in rows]
    ax.scatter(mxs, mys, s=8)

    ax.set_xlim(xmin - 0.5, xmax + 0.5)
    ax.set_ylim(ymin - 0.5, ymax + 0.5)
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"Differential link flux map, t = {t:.2f}")
    cbar = fig.colorbar(lc, ax=ax)
    cbar.set_label(r"$\Delta \phi_\ell$")
    plt.tight_layout()
    fname = os.path.join(OUTDIR, f"heatmap_t{t:.2f}.png")
    plt.savefig(fname, dpi=180)
    plt.close()

# ---------- plot 4: signed radial profiles ----------
plt.figure(figsize=(7.0, 5.0))
for t in chosen_times:
    pts = sorted(profiles_by_time[t], key=lambda z: z[0])
    r = [p[0] for p in pts]
    md = [p[1] for p in pts]
    plt.plot(r, md, marker="o", label=f"t={t:.2f}")
plt.xlabel("radius bin")
plt.ylabel("mean differential flux")
plt.title("Radial profile: signed differential flux")
plt.grid(True, alpha=0.3)
plt.legend(ncol=2, fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "radial_profile_signed.png"), dpi=180)
plt.close()

# ---------- plot 5: absolute radial profiles ----------
plt.figure(figsize=(7.0, 5.0))
for t in chosen_times:
    pts = sorted(profiles_by_time[t], key=lambda z: z[0])
    r = [p[0] for p in pts]
    mad = [p[2] for p in pts]
    plt.plot(r, mad, marker="o", label=f"t={t:.2f}")
plt.xlabel("radius bin")
plt.ylabel("mean absolute differential flux")
plt.title("Radial profile: absolute differential flux")
plt.grid(True, alpha=0.3)
plt.legend(ncol=2, fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "radial_profile_abs.png"), dpi=180)
plt.close()

print("Created figures in:", OUTDIR)
for fn in sorted(os.listdir(OUTDIR)):
    print(" -", fn)
