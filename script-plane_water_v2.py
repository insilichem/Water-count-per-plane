# =============================================================================
# Water Counting Script v2
# Counts water molecules on either side of a plane defined by 3 atoms.
# Iterates through the entire trajectory.
#
# Usage:
#   Open your trajectory in UCSF Chimera.
#   Run this script via File -> Open, or command: open script-plane_water_v2.py
# =============================================================================

from chimera import runCommand as rc, openModels, Molecule
from chimera import selection
import os
import math

# ---- Configuration ----
plane_atoms_str = ":236@C6,C5,C1"   # Atoms defining the plane
water_resname = "WAT"               # Residue name for water
output_filename = "water_count_v2.txt"

# ---- Initialization ----
# Find the primary molecule
mol_models = openModels.list(modelTypes=[Molecule])
if not mol_models:
    raise RuntimeError("No Molecule models found. Load your structure/trajectory first.")

mol = mol_models[0]
print "Processing molecule: %s (ID: #%d)" % (mol.name, mol.id)

# Select plane atoms
rc("select %s" % plane_atoms_str)
plane_sel_atoms = selection.currentAtoms()
if len(plane_sel_atoms) != 3:
    rc("~select")
    raise RuntimeError("Plane selection must have exactly 3 atoms; found %d. Check 'plane_atoms_str'." % len(plane_sel_atoms))

# Store references to the atoms (these persist across frames)
p_atom1 = plane_sel_atoms[0]
p_atom2 = plane_sel_atoms[1]
p_atom3 = plane_sel_atoms[2]

rc("~select") # Clear selection

# Identify water residues to avoid searching every frame
water_residues = []
for r in mol.residues:
    if r.type.strip().upper() == water_resname.upper():
        water_residues.append(r)

print "Found %d water residues (%s)." % (len(water_residues), water_resname)

# Prepare output file
out_path = os.path.abspath(output_filename)
try:
    f_out = open(out_path, "w")
    f_out.write("Frame\tSideA\tSideB\n")
    print "Writing results to: %s" % out_path
except IOError:
    raise RuntimeError("Could not open output file: %s" % out_path)

# ---- Analysis Loop ----
# Get all available frames
# coordSets is a map: key=frame_number (int), value=CoordSet object
frames = sorted(mol.coordSets.keys())
print "Analyzing %d frames..." % len(frames)

for frame in frames:
    # Set the active coordinate set for the molecule
    mol.activeCoordSet = mol.coordSets[frame]

    # 1. Get coordinates of the plane defining atoms for this frame
    # We use xformCoord() to get coordinates including any transformations,
    # though for relative calculations coord() would work if consistent.
    # xformCoord() is safer if the user has manipulated the view but we want world space relative.
    # Note: If the molecule tumbles in the trajectory, the plane moves with it.

    def get_coord(atom):
        c = atom.xformCoord()
        return (c.x, c.y, c.z)

    p1 = get_coord(p_atom1)
    p2 = get_coord(p_atom2)
    p3 = get_coord(p_atom3)

    # 2. Compute Plane Normal
    # u = p2 - p1
    u = (p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2])
    # v = p3 - p1
    v = (p3[0]-p1[0], p3[1]-p1[1], p3[2]-p1[2])

    # n = u x v
    nx = u[1]*v[2] - u[2]*v[1]
    ny = u[2]*v[0] - u[0]*v[2]
    nz = u[0]*v[1] - u[1]*v[0]

    n_len = math.sqrt(nx*nx + ny*ny + nz*nz)

    if n_len == 0.0:
        print "Warning: Frame %d has collinear plane atoms. Skipping." % frame
        f_out.write("%d\tNaN\tNaN\n" % frame)
        continue

    n_hat = (nx/n_len, ny/n_len, nz/n_len)
    r0 = p1 # Point on plane

    # 3. Count Waters
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

        # Center of water
        wx = sx / n_atoms
        wy = sy / n_atoms
        wz = sz / n_atoms

        # Signed distance to plane: d = n_hat . (water_pos - r0)
        d = (wx - r0[0])*n_hat[0] + (wy - r0[1])*n_hat[1] + (wz - r0[2])*n_hat[2]

        # Classify
        if d >= 0:
            side_a += 1
        else:
            side_b += 1

    # Output result for this frame
    # Use modulo to avoid spamming the console too much, or just print all
    if frame % 10 == 0 or frame == frames[0] or frame == frames[-1]:
        print "Frame %d: Side A = %d, Side B = %d" % (frame, side_a, side_b)

    f_out.write("%d\t%d\t%d\n" % (frame, side_a, side_b))

    # Flush occasionally to ensure data is written
    if frame % 50 == 0:
        f_out.flush()

f_out.close()
print "Analysis complete. Results saved to %s" % output_filename
