# =============================================================================
# Water Counting Script v3.1 (with Live Visualization)
#
# 1. Analyzes the entire trajectory and saves statistics to a file.
# 2. Activates a REAL-TIME VISUALIZATION mode that draws lines from the
#    plane centroid to water molecules within the search radius during playback.
#
# Usage:
#   1. Open your trajectory in UCSF Chimera.
#   2. Select exactly 3 atoms to define the plane (Ctrl+Click).
#   3. Run this script via File -> Open, or command: open script-plane_water_v2.py
# =============================================================================

import chimera
from chimera import runCommand as rc, openModels, Molecule, Element
from chimera import selection, triggers
import os
import math
import sys

# ---- Configuration ----
water_resname = "WAT"               # Residue name for water
output_filename = "water_count_v2.txt"
search_radius = 5.0                 # Angstroms (radius from plane centroid)

# Visualization Colors
color_centroid = "yellow"
color_side_a = "cyan"
color_side_b = "orange"

# =============================================================================
# Part 1: Initialization & Validation
# =============================================================================

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

mol = p_atom1.molecule
print "Processing molecule: %s (ID: #%d)" % (mol.name, mol.id)

# 3. Identify water residues
water_residues = []
for r in mol.residues:
    if r.type.strip().upper() == water_resname.upper():
        water_residues.append(r)

print "Found %d water residues (%s)." % (len(water_residues), water_resname)

# =============================================================================
# Part 2: Batch Analysis (Statistics File)
# =============================================================================

print "Starting batch analysis..."
frames = sorted(mol.coordSets.keys())

results_data = [] # List to store (frame, side_a, side_b)
grand_total_a = 0
grand_total_b = 0

# Store current frame to restore later
initial_frame = mol.activeCoordSet.id

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
        results_data.append((frame, -1, -1)) # Mark as invalid
        continue

    n_hat = (nx/n_len, ny/n_len, nz/n_len)
    r0 = p1

    # 4. Count Waters within Radius
    side_a = 0
    side_b = 0

    for r in water_residues:
        sx, sy, sz = 0.0, 0.0, 0.0
        n_atoms = 0
        for a in r.atoms:
            c = a.xformCoord()
            sx += c.x
            sy += c.y
            sz += c.z
            n_atoms += 1

        if n_atoms == 0: continue

        wx, wy, wz = sx/n_atoms, sy/n_atoms, sz/n_atoms

        # Check Distance to Centroid
        dist_sq = (wx - cx)**2 + (wy - cy)**2 + (wz - cz)**2
        if dist_sq > search_radius**2:
            continue

        # Signed distance
        d = (wx - r0[0])*n_hat[0] + (wy - r0[1])*n_hat[1] + (wz - r0[2])*n_hat[2]

        if d >= 0: side_a += 1
        else:      side_b += 1

    grand_total_a += side_a
    grand_total_b += side_b
    results_data.append((frame, side_a, side_b))

    if frame % 50 == 0:
        sys.stdout.write(".")
        sys.stdout.flush()

print "\nAnalysis complete."

# Restore frame
if initial_frame in mol.coordSets:
    mol.activeCoordSet = mol.coordSets[initial_frame]

# ---- Write File ----
winner_text = ""
if grand_total_a > grand_total_b:
    winner_text = "The side-a is the most water populated"
elif grand_total_b > grand_total_a:
    winner_text = "The side-b is the most water populated"
else:
    winner_text = "Both sides are equally populated"

summary_line = "# %s (Total A: %d, Total B: %d)" % (winner_text, grand_total_a, grand_total_b)
frame_count_line = "# Amount of frames calculated: %d" % len(frames)

print summary_line
print frame_count_line

out_path = os.path.abspath(output_filename)
try:
    f_out = open(out_path, "w")
    f_out.write(summary_line + "\n")
    f_out.write(frame_count_line + "\n")
    f_out.write("Frame\tSideA\tSideB\n")
    for row in results_data:
        frame, a, b = row
        if a == -1: f_out.write("%d\tNaN\tNaN\n" % frame)
        else:       f_out.write("%d\t%d\t%d\n" % (frame, a, b))
    f_out.close()
    print "Results saved to %s" % out_path
except IOError:
    print "Error writing to file: %s" % out_path


# =============================================================================
# Part 3: Live Visualization (Trigger)
# =============================================================================

class WaterVisualizer:
    def __init__(self, mol, p_atoms, water_residues, radius):
        self.mol = mol
        self.p_atoms = p_atoms
        self.water_residues = water_residues
        self.radius_sq = radius * radius
        self.handler_name = "WaterVisUpdate"
        self.pbg_name = "WaterPlaneLinks"

        # Colors
        self.c_side_a = chimera.Color.lookup(color_side_a)
        self.c_side_b = chimera.Color.lookup(color_side_b)
        self.c_centroid = chimera.Color.lookup(color_centroid)

        # Create Dummy Molecule for Centroid (needed for PseudoBonds)
        self.dummy_mol = Molecule()
        self.dummy_mol.name = "Plane Centroid"
        self.dummy_res = self.dummy_mol.newResidue("CTR", "1", 1)
        self.dummy_atom = self.dummy_mol.newAtom("C", Element("C"))
        self.dummy_res.addAtom(self.dummy_atom)
        openModels.add([self.dummy_mol])

        # Style the dummy atom
        self.dummy_atom.drawMode = chimera.Atom.Sphere
        self.dummy_atom.color = self.c_centroid
        self.dummy_atom.radius = 0.5

        # Initialize PseudoBondGroup
        self.pbg = chimera.misc.getPseudoBondGroup(self.pbg_name)

        # Register trigger
        triggers.addHandler('NewFrame', self.update, None)
        print "Live visualization enabled. Play trajectory to see updates."

        # Force one update now
        self.update('NewFrame', None, None)

    def cleanup(self):
        # Remove trigger
        try:
            triggers.deleteHandler('NewFrame', self.update)
        except:
            pass

        # Remove dummy molecule
        if self.dummy_mol:
            openModels.close([self.dummy_mol])
            self.dummy_mol = None

        # Clear PBG
        if self.pbg:
            self.pbg.deleteAll()
            # self.pbg = None # Keep reference if we want to reuse, but usually fine

    def update(self, trigger, myData, triggersData):
        if not self.mol or not self.mol.activeCoordSet:
            return

        # 1. Calculate Centroid
        def get_c(atom):
            c = atom.xformCoord()
            return (c.x, c.y, c.z)

        p1 = get_c(self.p_atoms[0])
        p2 = get_c(self.p_atoms[1])
        p3 = get_c(self.p_atoms[2])

        cx = (p1[0] + p2[0] + p3[0]) / 3.0
        cy = (p1[1] + p2[1] + p3[1]) / 3.0
        cz = (p1[2] + p2[2] + p3[2]) / 3.0

        # Update dummy atom position
        self.dummy_atom.setCoord(chimera.Point(cx, cy, cz))

        # 2. Calculate Normal (for side determination)
        u = (p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2])
        v = (p3[0]-p1[0], p3[1]-p1[1], p3[2]-p1[2])
        nx = u[1]*v[2] - u[2]*v[1]
        ny = u[2]*v[0] - u[0]*v[2]
        nz = u[0]*v[1] - u[1]*v[0]

        n_len = math.sqrt(nx*nx + ny*ny + nz*nz)
        if n_len == 0: return

        n_hat = (nx/n_len, ny/n_len, nz/n_len)
        r0 = p1

        # 3. Update PseudoBonds
        self.pbg.deleteAll()

        for r in self.water_residues:
            # Connect to the first atom (usually Oxygen)
            if not r.atoms: continue
            water_atom = r.atoms[0]

            # Use Center of Mass for filtering logic (to match batch analysis)
            sx, sy, sz = 0.0, 0.0, 0.0
            n = 0
            for a in r.atoms:
                c = a.xformCoord()
                sx += c.x; sy += c.y; sz += c.z
                n += 1
            wx, wy, wz = sx/n, sy/n, sz/n

            # Distance Check
            if (wx-cx)**2 + (wy-cy)**2 + (wz-cz)**2 > self.radius_sq:
                continue

            # Create Bond
            pb = self.pbg.newPseudoBond(self.dummy_atom, water_atom)

            # Color Check
            d = (wx-r0[0])*n_hat[0] + (wy-r0[1])*n_hat[1] + (wz-r0[2])*n_hat[2]
            if d >= 0:
                pb.color = self.c_side_a
            else:
                pb.color = self.c_side_b

# Clean up existing instance if script re-run
if hasattr(chimera, 'water_visualizer_instance'):
    chimera.water_visualizer_instance.cleanup()

# Instantiate and attach
chimera.water_visualizer_instance = WaterVisualizer(mol, [p_atom1, p_atom2, p_atom3], water_residues, search_radius)
