/**
 * VTU (VTK Unstructured Grid) output for 3D Basilisk simulations.
 * Writes an XML-format .vtu file readable by ParaView.
 *
 * Each leaf cell is output as a VTK_HEXAHEDRON (type 12) in 3D
 * or VTK_QUAD (type 9) in 2D.
 *
 * Usage:
 *   output_vtu(scalars, vectors, "snapshot", t);
 *
 * scalars: list of scalar fields (e.g., {p, cs})
 * vectors: list of vector fields (e.g., {u})
 * fname:   base filename (will produce fname.vtu)
 * tstep:   simulation time (stored as metadata)
 */

#ifndef OUTPUT_VTU_H
#define OUTPUT_VTU_H

static void output_vtu (scalar * scalars, vector * vectors,
                        const char * fname, double tstep)
{
  char name[256];
  snprintf(name, sizeof(name), "%s.vtu", fname);
  FILE * fp = fopen(name, "w");
  if (!fp) {
    fprintf(stderr, "output_vtu: cannot open %s\n", name);
    return;
  }

  /* ── Count leaf cells ───────────────────────────────────────────── */
  long nc = 0;
  foreach(reduction(+:nc))
    nc++;

  /* ── Header ─────────────────────────────────────────────────────── */
  fputs("<?xml version=\"1.0\"?>\n", fp);
  fputs("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
        "byte_order=\"LittleEndian\">\n", fp);
  fputs("  <UnstructuredGrid>\n", fp);

#if dimension == 3
  int verts_per_cell = 8;
  int cell_type = 12; /* VTK_HEXAHEDRON */
#else
  int verts_per_cell = 4;
  int cell_type = 9; /* VTK_QUAD */
#endif

  fprintf(fp, "    <Piece NumberOfPoints=\"%ld\" NumberOfCells=\"%ld\">\n",
          nc * verts_per_cell, nc);

  /* ── Points ─────────────────────────────────────────────────────── */
  fputs("      <Points>\n", fp);
  fputs("        <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
        "format=\"ascii\">\n", fp);

  foreach() {
#if dimension == 3
    /* 8 vertices of the hexahedron */
    for (int dz = 0; dz <= 1; dz++)
      for (int dy = 0; dy <= 1; dy++)
        for (int dx = 0; dx <= 1; dx++)
          fprintf(fp, "          %g %g %g\n",
                  x + (dx - 0.5)*Delta,
                  y + (dy - 0.5)*Delta,
                  z + (dz - 0.5)*Delta);
#else
    for (int dy = 0; dy <= 1; dy++)
      for (int dx = 0; dx <= 1; dx++)
        fprintf(fp, "          %g %g 0\n",
                x + (dx - 0.5)*Delta,
                y + (dy - 0.5)*Delta);
#endif
  }

  fputs("        </DataArray>\n", fp);
  fputs("      </Points>\n", fp);

  /* ── Cells ──────────────────────────────────────────────────────── */
  fputs("      <Cells>\n", fp);

  /* connectivity */
  fputs("        <DataArray type=\"Int64\" Name=\"connectivity\" "
        "format=\"ascii\">\n", fp);
  {
    long idx = 0;
    foreach() {
#if dimension == 3
      /* Hexahedron vertex ordering per VTK convention */
      fprintf(fp, "          %ld %ld %ld %ld %ld %ld %ld %ld\n",
              idx, idx+1, idx+3, idx+2,
              idx+4, idx+5, idx+7, idx+6);
#else
      fprintf(fp, "          %ld %ld %ld %ld\n",
              idx, idx+1, idx+3, idx+2);
#endif
      idx += verts_per_cell;
    }
  }
  fputs("        </DataArray>\n", fp);

  /* offsets */
  fputs("        <DataArray type=\"Int64\" Name=\"offsets\" "
        "format=\"ascii\">\n", fp);
  {
    long off = 0;
    foreach() {
      off += verts_per_cell;
      fprintf(fp, "          %ld\n", off);
    }
  }
  fputs("        </DataArray>\n", fp);

  /* types */
  fputs("        <DataArray type=\"UInt8\" Name=\"types\" "
        "format=\"ascii\">\n", fp);
  foreach()
    fprintf(fp, "          %d\n", cell_type);
  fputs("        </DataArray>\n", fp);

  fputs("      </Cells>\n", fp);

  /* ── Cell data ──────────────────────────────────────────────────── */
  fputs("      <CellData>\n", fp);

  /* Level (always useful for AMR visualization) */
  fputs("        <DataArray type=\"Int32\" Name=\"level\" "
        "format=\"ascii\">\n", fp);
  foreach()
    fprintf(fp, "          %d\n", level);
  fputs("        </DataArray>\n", fp);

  /* Scalar fields */
  for (scalar s in scalars) {
    fprintf(fp, "        <DataArray type=\"Float64\" Name=\"%s\" "
            "format=\"ascii\">\n", s.name);
    foreach()
      fprintf(fp, "          %g\n", val(s));
    fputs("        </DataArray>\n", fp);
  }

  /* Vector fields */
  for (vector v in vectors) {
    fprintf(fp, "        <DataArray type=\"Float64\" Name=\"%s\" "
            "NumberOfComponents=\"3\" format=\"ascii\">\n", v.x.name);
    foreach()
#if dimension == 3
      fprintf(fp, "          %g %g %g\n", val(v.x), val(v.y), val(v.z));
#else
      fprintf(fp, "          %g %g 0\n", val(v.x), val(v.y));
#endif
    fputs("        </DataArray>\n", fp);
  }

  fputs("      </CellData>\n", fp);

  /* ── Footer ─────────────────────────────────────────────────────── */
  fputs("    </Piece>\n", fp);
  fputs("  </UnstructuredGrid>\n", fp);
  fputs("</VTKFile>\n", fp);
  fclose(fp);

  fprintf(stderr, "# VTU: wrote %s (%ld cells, t=%g)\n", name, nc, tstep);
}

#endif /* OUTPUT_VTU_H */
