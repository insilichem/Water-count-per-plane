# =============================================================================
# Water Counting Script v2.3
# Counts water molecules on either side of a plane defined by 3 SELECTED atoms.
# Only counts waters within a specified radius of the plane center.
# Iterates through the entire trajectory.
# Determines which side is most populated overall.
#
# Usage:
#   1. Open your trajectory in UCSF Chimera.
#   2. Select exactly 3 atoms to define the plane (Ctrl+Click).
#   3. Run this script via File -> Open, or command: open script-plane_water_v2.py
# =============================================================================

from chimera import runCommand as rc, openModels, Molecule
from chimera import selection
import os
import math

# ---- Configuration ----
water_resname = "WAT"               # Residue name for water
output_filename = "water_count_v2.txt"
search_radius = 5.0                 # Angstroms (radius from plane centroid)

# ---- Initialization ----

# 1. Validate Plane Selection
plane_sel_atoms = selection.currentAtoms()
if len(plane_sel_atoms) != 3:
    raise RuntimeError("Error: Select 3 atoms to define the plane. (Found %d atoms selected)" % len(plane_sel_atoms))

# Store references to the atoms
p_atom1 = plane_sel_atoms[0]
p_atom2 = plane_sel_atoms[1]
p_atom3 = plane_sel_atoms[2]

# 2. Find the primary molecule
mol_models = openModels.list(modelTypes=[Molecule])
if not mol_models:
    raise RuntimeError("No Molecule models found. Load your structure/trajectory first.")

# Use the molecule the first plane atom belongs to
mol = p_atom1.molecule
print "Processing molecule: %s (ID: #%d)" % (mol.name, mol.id)

# 3. Identify water residues
water_residues = []
for r in mol.residues:
    if r.type.strip().upper() == water_resname.upper():
        water_residues.append(r)

print "Found %d water residues (%s)." % (len(water_residues), water_resname)

# ---- Analysis Loop ----
frames = sorted(mol.coordSets.keys())
print "Analyzing %d frames with search radius %.1f A..." % (len(frames), search_radius)

results_data = [] # List to store (frame, side_a, side_b)
grand_total_a = 0
grand_total_b = 0

for frame in frames:
    mol.activeCoordSet = mol.coordSets[frame]

    # Helper to get world coordinates
    def get_coord(atom):
        c = atom.xformCoord()
        return (c.x, c.y, c.z)

    # 1. Get Plane Coordinates
    p1 = get_coord(p_atom1)
    p2 = get_coord(p_atom2)
    p3 = get_coord(p_atom3)

    # 2. Compute Centroid of the 3 atoms
    cx = (p1[0] + p2[0] + p3[0]) / 3.0
    cy = (p1[1] + p2[1] + p3[1]) / 3.0
    cz = (p1[2] + p2[2] + p3[2]) / 3.0
    centroid = (cx, cy, cz)

    # 3. Compute Plane Normal
    u = (p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2])
    v = (p3[0]-p1[0], p3[1]-p1[1], p3[2]-p1[2])

    nx = u[1]*v[2] - u[2]*v[1]
    ny = u[2]*v[0] - u[0]*v[2]
    nz = u[0]*v[1] - u[1]*v[0]

    n_len = math.sqrt(nx*nx + ny*ny + nz*nz)

    if n_len == 0.0:
        print "Warning: Frame %d has collinear plane atoms. Skipping." % frame
        results_data.append((frame, -1, -1)) # Mark as invalid
        continue

    n_hat = (nx/n_len, ny/n_len, nz/n_len)
    r0 = p1 # Point on plane (for signed distance calculation)

    # 4. Count Waters within Radius
    side_a = 0
    side_b = 0

    for r in water_residues:
        # Calculate geometric center of the water residue
        sx, sy, sz = 0.0, 0.0, 0.0
        n_atoms = 0
        for a in r.atoms:
            c = a.xformCoord()
            sx += c.x
            sy += c.y
            sz += c.z
            n_atoms += 1

        if n_atoms == 0:
            continue

        wx = sx / n_atoms
        wy = sy / n_atoms
        wz = sz / n_atoms

        # Check Distance to Centroid
        dist_sq = (wx - cx)**2 + (wy - cy)**2 + (wz - cz)**2
        if dist_sq > search_radius**2:
            continue

        # Signed distance to plane
        d = (wx - r0[0])*n_hat[0] + (wy - r0[1])*n_hat[1] + (wz - r0[2])*n_hat[2]

        if d >= 0:
            side_a += 1
        else:
            side_b += 1

    # Accumulate
    grand_total_a += side_a
    grand_total_b += side_b

    # Store
    results_data.append((frame, side_a, side_b))

    # Progress
    if frame % 10 == 0 or frame == frames[0] or frame == frames[-1]:
        print "Frame %d: Side A = %d, Side B = %d" % (frame, side_a, side_b)

# ---- Result Writing ----

# Determine winner
winner_text = ""
if grand_total_a > grand_total_b:
    winner_text = "The side-a is the most water populated"
elif grand_total_b > grand_total_a:
    winner_text = "The side-b is the most water populated"
else:
    winner_text = "Both sides are equally populated"

summary_line = "# %s (Total A: %d, Total B: %d)" % (winner_text, grand_total_a, grand_total_b)
frame_count_line = "# Amount of frames calculated: %d" % len(frames)

print "\nAnalysis Complete."
print summary_line
print frame_count_line

# Write to file
out_path = os.path.abspath(output_filename)
try:
    f_out = open(out_path, "w")
    f_out.write(summary_line + "\n")
    f_out.write(frame_count_line + "\n")
    f_out.write("Frame\tSideA\tSideB\n")

    for row in results_data:
        frame, a, b = row
        if a == -1:
            f_out.write("%d\tNaN\tNaN\n" % frame)
        else:
            f_out.write("%d\t%d\t%d\n" % (frame, a, b))

    f_out.close()
    print "Results saved to %s" % out_path
except IOError:
    raise RuntimeError("Could not open output file: %s" % out_path)
