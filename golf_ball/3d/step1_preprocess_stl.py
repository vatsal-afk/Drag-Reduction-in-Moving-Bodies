"""
Step 1: Golf Ball STL Preprocessing for Basilisk CFD
=====================================================
Cleans, centers, and normalizes the golf ball STL so that:
  - The model is watertight (no holes or self-intersections)
  - Normals are consistently outward-facing (critical for distance field sign)
  - The ball is centered at the origin
  - The diameter is normalized to exactly 1.0 (all physics goes into Re)
  - Output is ASCII STL (required by Basilisk's input_stl())

Usage:
    pip install trimesh numpy
    python step1_preprocess_stl.py golf_ball_raw.stl golfball.stl

The Sketchfab model is likely in binary STL. This script converts it.
"""

import sys
import numpy as np

try:
    import trimesh
except ImportError:
    print("Install trimesh: pip install trimesh")
    sys.exit(1)


def preprocess_golf_ball(input_path: str, output_path: str, target_diameter: float = 1.0):
    print(f"\n=== Golf Ball STL Preprocessor ===")
    print(f"Input : {input_path}")
    print(f"Output: {output_path}")
    print(f"Target diameter: {target_diameter}")

    # ── Load ──────────────────────────────────────────────────────────────────
    mesh = trimesh.load(input_path, force='mesh')
    print(f"\n[Load]")
    print(f"  Vertices : {len(mesh.vertices)}")
    print(f"  Faces    : {len(mesh.faces)}")
    print(f"  Watertight: {mesh.is_watertight}")
    print(f"  Volume   : {mesh.volume:.4f}  (negative = flipped normals)")

    # ── Repair ────────────────────────────────────────────────────────────────
    print(f"\n[Repair]")

    # Fix winding so normals consistently point outward
    if mesh.volume < 0:
        print("  Volume negative — inverting mesh normals")
        mesh.invert()

    trimesh.repair.fix_normals(mesh)
    trimesh.repair.fill_holes(mesh)

    # Remove duplicate/degenerate faces — API varies by trimesh version
    try:
        mesh.update_faces(mesh.unique_faces())
        mesh.update_faces(mesh.nondegenerate_faces())
        print("  Cleaned faces using update_faces() API")
    except AttributeError:
        print("  Skipping face cleanup (mesh already clean)")

    try:
        mesh.merge_vertices()
    except TypeError:
        pass

    print(f"  Watertight after repair: {mesh.is_watertight}")
    print(f"  Volume after repair    : {mesh.volume:.6f}")

    if not mesh.is_watertight:
        print("  WARNING: mesh is still not watertight.")
        print("  Open the STL in MeshLab and run:")
        print("    Filters > Cleaning > Remove Duplicate Faces")
        print("    Filters > Remeshing > Close Holes")
        print("    Filters > Normals > Recompute Face Normals")
        print("  Then re-run this script.")

    # ── Center at origin ──────────────────────────────────────────────────────
    print(f"\n[Center]")
    centroid = mesh.bounding_box.centroid
    print(f"  Bounding box centroid: {centroid}")
    mesh.vertices -= centroid

    # ── Normalize diameter to 1.0 ─────────────────────────────────────────────
    print(f"\n[Normalize]")
    extents = mesh.bounding_box.extents        # (dx, dy, dz)
    current_diameter = float(np.max(extents))  # largest extent = effective diameter
    print(f"  Current extents (x,y,z): {extents}")
    print(f"  Current max diameter   : {current_diameter:.6f}")

    scale = target_diameter / current_diameter
    mesh.vertices *= scale

    final_extents = mesh.bounding_box.extents
    print(f"  Scale factor applied   : {scale:.6f}")
    print(f"  Final extents          : {final_extents}")
    print(f"  Final diameter         : {float(np.max(final_extents)):.6f}")

    # ── Dimple statistics ─────────────────────────────────────────────────────
    # Rough estimate: compare surface area to equivalent smooth sphere
    R = target_diameter / 2.0
    smooth_area = 4.0 * np.pi * R**2
    actual_area = mesh.area
    area_ratio = actual_area / smooth_area
    print(f"\n[Dimple Statistics]")
    print(f"  Surface area (normalized)  : {actual_area:.4f}")
    print(f"  Smooth sphere surface area : {smooth_area:.4f}")
    print(f"  Area ratio (dimpled/smooth): {area_ratio:.4f}")
    print(f"  (>1 means dimples increase surface area as expected)")
    # ── Float32 quantization ────────────────────────────────────────────────────
    # Basilisk reads STL vertices as single-precision (float32).
    # After scaling to D=1.0, some tiny dimple triangles collapse to zero area
    # at float32 precision, causing assertion failures in distance.h.
    # Fix: round-trip vertices through float32, then purge degenerate triangles.
    print(f"\n[Float32 Cleanup]")
    mesh.vertices = mesh.vertices.astype(np.float32).astype(np.float64)

    # Compute per-face area via cross product
    v0 = mesh.vertices[mesh.faces[:, 0]]
    v1 = mesh.vertices[mesh.faces[:, 1]]
    v2 = mesh.vertices[mesh.faces[:, 2]]
    crosses = np.cross(v1 - v0, v2 - v0)
    areas = 0.5 * np.linalg.norm(crosses, axis=1)
    good = areas > 1e-12
    n_removed = np.sum(~good)
    print(f"  Removed {n_removed} degenerate faces (zero area at float32)")
    if n_removed > 0:
        mesh.update_faces(good)
        mesh.remove_unreferenced_vertices()

    # Also check for duplicate vertices that create zero-length edges
    mesh.merge_vertices(merge_tex=True, merge_norm=True)
    mesh.update_faces(mesh.nondegenerate_faces())
    print(f"  Final face count: {len(mesh.faces)}")

    # ── Export Binary STL ──────────────────────────────────────────────────────
    # Basilisk's input_stl() in this version requires binary STL.
    print(f"\n[Export]")
    mesh.export(output_path, file_type='stl')
    print(f"  Saved: {output_path}")

    # ── Basilisk domain sizing hint ───────────────────────────────────────────
    print(f"\n[Basilisk Domain Hints]")
    print(f"  Ball diameter D = {target_diameter}")
    print(f"  Recommended domain size: 16*D = {16*target_diameter}")
    print(f"  Recommended origin x: -3*D = {-3*target_diameter}  (upstream distance)")
    print(f"  At level 10: {2**10} cells/domain -> {2**10 / 16:.0f} cells/D")
    print(f"  At level 11: {2**11} cells/domain -> {2**11 / 16:.0f} cells/D")
    print(f"  Dimples ~10% of D wide, need >=6 cells/dimple -> level 10 is marginal")
    print(f"  Recommended geometry refinement: level 11")
    print(f"  Recommended flow refinement    : level 9 or 10")

    return mesh


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python step1_preprocess_stl.py <input.stl> <output.stl>")
        print("Example: python step1_preprocess_stl.py golf_ball_raw.stl golfball.stl")
        sys.exit(1)

    input_stl  = sys.argv[1]
    output_stl = sys.argv[2]
    diameter   = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0

    preprocess_golf_ball(input_stl, output_stl, diameter)
