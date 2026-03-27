"""
visualize.py  --  Beam FEA result visualiser
Reads results/ CSV files produced by the C++ solver and renders:
  - Animated heatmap of the beam (stress / deflection / axial)
  - Time history plot of max deflection
  - Saves animation to results/beam_animation.gif
"""

import csv, os, sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.animation as animation
from matplotlib.gridspec import GridSpec

RESULTS_DIR = "results"

# ── helpers ──────────────────────────────────────────────────────────────────

def read_csv(fname):
    path = os.path.join(RESULTS_DIR, fname)
    if not os.path.exists(path):
        sys.exit(f"[error] {path} not found. Run the solver first.")
    with open(path) as f:
        reader = csv.reader(f)
        header = next(reader)
        rows   = [[float(v) for v in row] for row in reader]
    return header, np.array(rows)

def read_metadata():
    path = os.path.join(RESULTS_DIR, "metadata.csv")
    meta = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:    meta[row["key"]] = float(row["value"])
            except: meta[row["key"]] = row["value"]
    return meta

# ── load data ─────────────────────────────────────────────────────────────────

print("[vis] Loading results...")
meta    = read_metadata()
L       = meta["L"]
n_nodes = int(meta["n_nodes"])
n_elem  = int(meta["n_elem"])
n_frames= int(meta["n_frames"])

_, defl_data  = read_csv("deflection.csv")  # [n_frames, 1+n_nodes]
_, axial_data = read_csv("axial.csv")       # [n_frames, 1+n_nodes]
_, stress_data= read_csv("stress.csv")      # [n_frames, 1+n_elem]

times    = defl_data[:, 0]
defl     = defl_data[:, 1:]          # [n_frames, n_nodes]
axial    = axial_data[:, 1:]
stress   = stress_data[:, 1:]        # [n_frames, n_elem]

node_x  = np.linspace(0, L, n_nodes)
elem_x  = 0.5 * (node_x[:-1] + node_x[1:])   # element midpoints

print(f"[vis] {n_frames} frames  |  {n_nodes} nodes  |  {n_elem} elements")
print(f"[vis] t = [{times[0]:.3f}, {times[-1]:.3f}] s")

# ── colour maps ───────────────────────────────────────────────────────────────
CMAPS = {
    "Transverse deflection": (defl,  node_x,  "m",   "RdYlBu_r"),
    "Axial displacement":    (axial, node_x,  "m",   "PuOr"),
    "Bending stress":        (stress,elem_x,  "Pa",  "coolwarm"),
}

# ── figure layout ─────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(13, 9), facecolor="#0f0f0f")
fig.suptitle("Beam FEA — Transient Analysis", color="white",
             fontsize=14, fontweight="bold", y=0.98)

gs = GridSpec(4, 3, figure=fig, hspace=0.55, wspace=0.35,
              top=0.93, bottom=0.07, left=0.07, right=0.97)

# top row: 3 beam heatmaps
axes_beam = [fig.add_subplot(gs[0:2, i]) for i in range(3)]
# bottom row: time-history for each quantity
axes_hist = [fig.add_subplot(gs[2:4, i]) for i in range(3)]

for ax in axes_beam + axes_hist:
    ax.set_facecolor("#1a1a1a")
    for spine in ax.spines.values(): spine.set_edgecolor("#444")
    ax.tick_params(colors="#aaa", labelsize=8)
    ax.xaxis.label.set_color("#aaa")
    ax.yaxis.label.set_color("#aaa")
    ax.title.set_color("white")

# ── beam drawing ──────────────────────────────────────────────────────────────
BEAM_HALF = 0.08    # half-height of beam rectangle in data (normalised) coords

beam_artists = []   # (patchcollection per column)
time_lines   = []
hist_vals    = []   # track max-abs per frame for time histories

labels = list(CMAPS.keys())
for col, (label, (data, xpos, unit, cmap_name)) in enumerate(CMAPS.items()):
    ax  = axes_beam[col]
    axh = axes_hist[col]

    vmax = np.max(np.abs(data)) or 1.0
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    cmap = plt.get_cmap(cmap_name)

    n_pts = len(xpos)
    rects = []
    from matplotlib.patches import Rectangle
    from matplotlib.collections import PatchCollection

    dx = xpos[1] - xpos[0] if len(xpos) > 1 else 1.0
    patches = []
    for i in range(n_pts):
        x0 = xpos[i] - dx/2 if i == 0 else (xpos[i-1]+xpos[i])/2
        x1 = xpos[i] + dx/2 if i == n_pts-1 else (xpos[i]+xpos[i+1])/2
        patches.append(Rectangle((x0, -BEAM_HALF), x1-x0, 2*BEAM_HALF))

    pc = PatchCollection(patches, cmap=cmap, norm=norm, alpha=0.95)
    pc.set_array(data[0])
    ax.add_collection(pc)
    beam_artists.append((pc, data, cmap, norm))

    # deflection overlay line (only for deflection column)
    if col == 0:
        defl_scale = BEAM_HALF / (vmax + 1e-30) * 0.7
        (defl_line,) = ax.plot(node_x, data[0] * defl_scale, color="#fff",
                                lw=1.0, alpha=0.7, zorder=5)
    else:
        defl_line = None

    ax.set_xlim(-0.05*L, 1.05*L)
    ax.set_ylim(-BEAM_HALF*2, BEAM_HALF*2)
    ax.set_xlabel("x [m]", fontsize=8)
    ax.set_yticks([])
    ax.set_title(label, fontsize=9, pad=4)
    cb = plt.colorbar(pc, ax=ax, orientation="horizontal", pad=0.25, shrink=0.9)
    cb.set_label(unit, color="#aaa", fontsize=7)
    cb.ax.tick_params(labelsize=6, colors="#aaa")
    cb.outline.set_edgecolor("#444")

    # time history
    max_abs = np.max(np.abs(data), axis=1)
    hist_vals.append(max_abs)
    axh.plot(times, max_abs, color="#555", lw=0.8, alpha=0.6)
    (tl,) = axh.plot([], [], color="#f08030", lw=1.5)
    (dot,) = axh.plot([], [], "o", color="#f08030", ms=4)
    time_lines.append((tl, dot, max_abs))
    axh.set_xlim(times[0], times[-1])
    axh.set_ylim(0, max_abs.max()*1.15 or 1.0)
    axh.set_xlabel("t [s]", fontsize=8)
    axh.set_ylabel(f"|max| [{unit}]", fontsize=8)
    axh.set_title(f"max |{label.split()[0].lower()}|", fontsize=9, pad=4)

    beam_artists[-1] = (pc, data, cmap, norm, defl_line)

# timestamp
ts_text = fig.text(0.5, 0.005, "t = 0.000 s", ha="center", va="bottom",
                   color="#aaa", fontsize=10)

# ── animation ─────────────────────────────────────────────────────────────────

def update(frame):
    artists_out = []
    for col, (pc, data, cmap, norm, defl_line) in enumerate(beam_artists):
        pc.set_array(data[frame])
        artists_out.append(pc)
        if defl_line is not None:
            vmax = norm.vmax or 1.0
            defl_scale = BEAM_HALF / (vmax + 1e-30) * 0.7
            defl_line.set_ydata(data[frame] * defl_scale)
            artists_out.append(defl_line)

        tl, dot, mvals = time_lines[col]
        tl.set_data(times[:frame+1], mvals[:frame+1])
        dot.set_data([times[frame]], [mvals[frame]])
        artists_out += [tl, dot]

    ts_text.set_text(f"t = {times[frame]:.3f} s")
    artists_out.append(ts_text)
    return artists_out

n_animate = min(n_frames, 200)   # cap at 200 frames for gif speed
frame_idx = np.round(np.linspace(0, n_frames-1, n_animate)).astype(int)

ani = animation.FuncAnimation(fig, lambda i: update(frame_idx[i]),
                               frames=n_animate, interval=40, blit=True)

out_gif = os.path.join(RESULTS_DIR, "beam_animation.gif")
print(f"[vis] Saving animation to {out_gif}  (this may take ~30s)...")
ani.save(out_gif, writer="pillow", fps=25, dpi=110)
print(f"[vis] Saved: {out_gif}")

plt.tight_layout(rect=[0,0.02,1,0.97])
plt.savefig(os.path.join(RESULTS_DIR, "beam_final_frame.png"), dpi=130,
            facecolor=fig.get_facecolor())
print(f"[vis] Saved: {RESULTS_DIR}/beam_final_frame.png")
plt.show()
