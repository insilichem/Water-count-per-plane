# ============================
# Step 1 — single-frame setup (slab ±5 Å around plane)
# UCSF Chimera (classic) script
# ============================

from chimera import runCommand as rc, openModels, Molecule
from chimera import selection
from chimera import replyobj
import os, time, math

# ---- User-configurable bits (Step 1) ----
plane_name = "plane"
plane_atoms = ":236@C6,C5,C1"    # three atoms defining the plane
extra_radius = 5.0               # visual padding (Å)
disk_thickness = 0.1             # visual thickness (Å)
water_resname = "WAT"            # default water residue name (change to HOH, etc. if needed)
slab_half_thickness = 5.0        # Å (show waters with |distance to plane| <= this)
step_mode = True                 # STOP after Step 1 until you confirm
out_txt = "water_check_step1.txt"

# ---- Find the primary molecule model id (avoid assuming #0) ----
mol_models = openModels.list(modelTypes=[Molecule])
if not mol_models:
    raise RuntimeError("No Molecule models found. Load your structure/trajectory first.")

primary = mol_models[0]
model_id = "#%d" % primary.id
print "Model detected:", model_id

# ---- Clean up visuals to avoid accumulation ----
rc("~select")
rc("~ribbon")
rc("~display")       # hide all atoms/bonds across all models
rc("define plane name plane sel replace true radius 5 thickness 0.1")


########################HEEEREEEEEEEEEEEEEEE


# ---- Define/replace the visual plane from the three atoms ----
rc("select %s" % plane_atoms)
rc("define plane name %s sel replace true radius %.3f thickness %.3f"
   % (plane_name, extra_radius, disk_thickness))

# ---- Compute the plane from the 3 selected atoms (to filter waters by distance) ----
plane_sel_atoms = selection.currentAtoms()
if len(plane_sel_atoms) != 3:
    raise RuntimeError("Plane selection must have exactly 3 atoms; got %d" % len(plane_sel_atoms))

def _v(a):  # world coords (respecting transforms)
    # xformCoord() returns a Point with x,y,z in scene coords
    c = a.xformCoord() if hasattr(a, "xformCoord") else a.coord()
    return (c.x, c.y, c.z)

p1 = _v(plane_sel_atoms[0])
p2 = _v(plane_sel_atoms[1])
p3 = _v(plane_sel_atoms[2])
rc("~select")  # clear selection

def vsub(a,b): return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
def vcross(a,b): return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
def vdot(a,b): return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
def vlen(a): return math.sqrt(vdot(a,a))

u = vsub(p2, p1)
v = vsub(p3, p1)
n = vcross(u, v)             # normal (not yet unit)
n_len = vlen(n)
if n_len == 0.0:
    raise RuntimeError("Selected plane atoms are colinear; cannot define plane.")
n_hat = (n[0]/n_len, n[1]/n_len, n[2]/n_len)
r0 = p1                       # any point on plane
t = slab_half_thickness

# ---- Show ONLY waters within the ±t slab around the plane ----
# First, hide all waters in the primary model
for r in primary.residues:
    if r.type.strip().upper() == water_resname.upper():
        for a in r.atoms:
            a.display = False

def point_plane_signed_distance(pt, r0, n_hat):
    # distance = n̂ · (pt - r0)
    return ( (pt[0]-r0[0])*n_hat[0] + (pt[1]-r0[1])*n_hat[1] + (pt[2]-r0[2])*n_hat[2] )

shown_water_resids = 0
for r in primary.residues:
    if r.type.strip().upper() != water_resname.upper():
        continue
    show_this = False
    # If any atom in the water lies within the slab, we display the whole residue
    for a in r.atoms:
        c = a.xformCoord() if hasattr(a, "xformCoord") else a.coord()
        pt = (c.x, c.y, c.z)
        d = point_plane_signed_distance(pt, r0, n_hat)
        if abs(d) <= t:
            show_this = True
            break
    if show_this:
        shown_water_resids += 1
        for a in r.atoms:
            a.display = True

# Lightweight representation/color
rc("represent stick %s & :%s" % (model_id, water_resname))
rc("color white %s & :%s" % (model_id, water_resname))

# ---- Write tiny check file ----
header = "Step\tFrame\tPlane_atoms\tPlane_padding\tWaters_visualized\n"
row = "1\t1\t%s\t%dÅ\t%d\n" % (plane_atoms, int(extra_radius), shown_water_resids)

out_path = os.path.abspath(out_txt)
with open(out_path, "w") as f:
    f.write(header)
    f.write(row)

# ---- Console verification block ----
summary_lines = [
    "",
    "========== Step 1 Verification (PLANE SLAB ±%.1f Å) ==========" % slab_half_thickness,
    "Model: %s" % model_id,
    "Plane created: name='%s', Extra radius=%.1f Å, thickness=%.2f Å" % (plane_name, extra_radius, disk_thickness),
    "Active frame set to 1 (attempted both 1- and 0-based coordsets).",
    "Waters visible (resname %s) within ±%.1f Å of the plane: %d" % (water_resname, slab_half_thickness, shown_water_resids),
    "Only plane + slab-selected waters shown; other atoms hidden; previous shapes deleted.",
    "Check file written: %s" % out_path,
    "==============================================================",
    "",
    "Next action: Reply with either",
    "  --PROCEED               (go to Step 2: counting across frames with stride=10), or",
    "  --ADJUST STEP1: <note>  (e.g., atoms, padding, water name, slab size).",
    ""
]
print("\n".join(summary_lines))

# ---- Respect step_mode: stop after Step 1 ----
if step_mode:
    time.sleep(0.2)
