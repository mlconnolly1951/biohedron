#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

/* msdensity.c Copyright 1995 by Michael L. Connolly */

/* February 8, 2000 */


struct surface *add_densities (struct surface *den1, struct surface *den2)
{
	struct surface *den;
	den = combine_densities (den1, den2, (double) ( 1.0));
	return (den);
}


struct surface *subtract_densities (struct surface *den1, struct surface *den2)
{
	struct surface *den;
	den = combine_densities (den1, den2, (double) (-1.0));
	return (den);
}


struct surface *combine_densities (struct surface *den1, struct surface *den2, double sgn)
{
	int k, in1, in2, inx1, iny1, inz1, inx2, iny2, inz2;
	long n1, n2, nboth, nneither;
	long dx, dx1, dx2, dy, dy1, dy2, dz, dz1, dz2;
	long x, y, z, idx, idx1, idx2;
	long ncube1, ncube2, n_cube;
	long origin1[3], width1[3], limit1[3];
	long origin2[3], width2[3], limit2[3];
	long origin[3], width[3], limit[3];
	float d, d1, d2;
	float *densities1, *densities2, *densities;
	double cube_width, average;
	char message[MAXLINE];
	struct surface *den;

	if (den1 == NULL) return (NULL);
	else if (den2 == NULL) return (NULL);
	else if (den1 -> cube_width != den2 -> cube_width) {
		sprintf (message, "combine_densities: cube widths not equal: %6.3f %6.3f",
			den1 -> cube_width, den2 -> cube_width);
		set_error1 (message);
		return (NULL);
	}
	cube_width = den1 -> cube_width;
	
	/* copy to local variables */
	for (k = 0; k < 3; k++) {
		origin1[k] = den1 -> origin[k];
		width1[k] = den1 -> width[k];
		limit1[k] = den1 -> limit[k];
		origin2[k] = den2 -> origin[k];
		width2[k] = den2 -> width[k];
		limit2[k] = den2 -> limit[k];
	}
	densities1 = den1 -> densities;
	densities2 = den2 -> densities;
	if (densities1 == NULL || densities2 == NULL) {
		set_error1 ("combine_densities: null input density");
		return (NULL);
	}
	ncube1 = den1 -> n_cube;
	ncube2 = den2 -> n_cube;
	if (ncube1 <= 0 || ncube2 <= 0) {
		set_error1 ("combine_densities: 0 n_cube input density");
		return (NULL);
	}

	/* define a bounds large enough to encompass both */
	for (k = 0; k < 3; k++) {
		if (origin1[k] < origin2[k]) origin[k] = origin1[k];
		else origin[k] = origin2[k];
		if (limit1[k] > limit2[k]) limit[k] = limit1[k];
		else limit[k] = limit2[k];
	}
	for (k = 0; k < 3; k++) {
		width[k] = limit[k] - origin[k];
		if (width[k] <= 0) {
			set_error1 ("combine_densities: zero-width density");
			return (NULL);
		}
	}
	n_cube = width[0] * width[1] * width[2];
	sprintf (message, "%8ld %8ld %8ld cubes in arguments and result, respectively",
		 ncube1, ncube2, n_cube);
	inform (message);
	
	/* create new density */
	den = create_density (cube_width, (int *) origin, (int *) width);
	if (den == NULL) {
		set_error1 ("combine_densities: not enough memory for new density");
		return (NULL);
	}
	densities = den -> densities;

	/* compute data */
	n1 = n2 = nboth = nneither = 0; average = 0.0;
	for (z = origin[2]; z < limit[2]; z++) {
		inz1 = (z >= origin1[2] && z < limit1[2]);
		inz2 = (z >= origin2[2] && z < limit2[2]);
		for (y = origin[1]; y < limit[1]; y++) {
			iny1 = (y >= origin1[1] && y < limit1[1]);
			iny2 = (y >= origin2[1] && y < limit2[1]);
			for (x = origin[0]; x < limit[0]; x++) {
				inx1 = (x >= origin1[0] && x < limit1[0]);
				inx2 = (x >= origin2[0] && x < limit2[0]);
				in1 = inx1 && iny1 && inz1;
				in2 = inx2 && iny2 && inz2;
				if (in1) {
					dx1 = x - origin1[0];
					dy1 = y - origin1[1];
					dz1 = z - origin1[2];
					idx1 = dz1 * (width1[0] * width1[1]) + dy1 * width1[0] + dx1;
					d1 = *(densities1 + idx1);
				}
				else d1 = 0.0;
				if (in2) {
					dx2 = x - origin2[0];
					dy2 = y - origin2[1];
					dz2 = z - origin2[2];
					idx2 = dz2 * (width2[0] * width2[1]) + dy2 * width2[0] + dx2;
					d2 = *(densities2 + idx2);
				}
				else d2 = 0.0;
				if (in1 && in2) nboth++;
				else if (in1) n1++;
				else if (in2) n2++;
				else nneither++;
				d = d1 + sgn * d2;
				dx = x - origin[0];
				dy = y - origin[1];
				dz = z - origin[2];
				idx = dz * (width[0] * width[1]) + dy * width[0] + dx;
				*(densities+idx) = d;
				average += d;
			}
		}
	}
	average /= n_cube;
	sprintf (message, "%8ld cubes in both given densities", nboth);
	inform (message);
	sprintf (message, "%8ld cubes in first given density", n1);
	inform (message);
	sprintf (message, "%8ld cubes in second given density", n2);
	inform (message);
	sprintf (message, "%8ld cubes in neither given density", nneither);
	inform (message);
	sprintf (message, "%8.4f average combined density", average);
	inform (message);
	return (den);
}

long collect_den_vtx0 (struct surface *den, double ctrlev, struct denvtx **denvtx_handle)
{
	long ix, ix0, ix1, iy, iy0, iy1, iz, iz0, iz1;
	long idx000, idx001, idx010, idx011;
	long idx100, idx101, idx110, idx111;
	double den000, den001, den010, den011;
	double den100, den101, den110, den111;
	long size1, size2, size3;
	long n_vertex;
	struct denvtx *denvtxs, *dv, *dv0, *dv1, *dv2, *dv3;

	size1 = den -> width[0];
	size2 = size1 * den -> width[1];
	size3 = size2 * den -> width[2];
	denvtxs = (struct denvtx *) allocate_objects (DENVTX, size3);
	if (denvtxs == NULL) {
		set_error1 ("no memory for denvtxs");
		return (0);
	}
	/* find surface vertices of density lattice */
	n_vertex = 0;
	for (iz = den -> origin[2]+1; iz < den -> limit[2]; iz++) {
		for (iy = den -> origin[1]+1; iy < den -> limit[1]; iy++) {
			for (ix = den -> origin[0]+1; ix < den -> limit[0]; ix++) {
				ix0 = ix - den -> origin[0] - 1; ix1 = ix0 + 1;
				iy0 = iy - den -> origin[1] - 1; iy1 = iy0 + 1;
				iz0 = iz - den -> origin[2] - 1; iz1 = iz0 + 1;
				idx000 = iz0 * size2 + iy0 * size1 + ix0;
				idx001 = iz0 * size2 + iy0 * size1 + ix1;
				idx010 = iz0 * size2 + iy1 * size1 + ix0;
				idx011 = iz0 * size2 + iy1 * size1 + ix1;
				idx100 = iz1 * size2 + iy0 * size1 + ix0;
				idx101 = iz1 * size2 + iy0 * size1 + ix1;
				idx110 = iz1 * size2 + iy1 * size1 + ix0;
				idx111 = iz1 * size2 + iy1 * size1 + ix1;
				den000 = *(den -> densities + idx000);
				den001 = *(den -> densities + idx001);
				den010 = *(den -> densities + idx010);
				den011 = *(den -> densities + idx011);
				den100 = *(den -> densities + idx100);
				den101 = *(den -> densities + idx101);
				den110 = *(den -> densities + idx110);
				den111 = *(den -> densities + idx111);
				dv = denvtxs + idx000;
				dv -> normal[0] -= den000;
				dv -> normal[1] -= den000;
				dv -> normal[2] -= den000;
				dv -> normal[0] += den001;
				dv -> normal[1] -= den001;
				dv -> normal[2] -= den001;
				dv -> normal[0] -= den010;
				dv -> normal[1] += den010;
				dv -> normal[2] -= den010;
				dv -> normal[0] += den011;
				dv -> normal[1] += den011;
				dv -> normal[2] -= den011;
				dv -> normal[0] -= den100;
				dv -> normal[1] -= den100;
				dv -> normal[2] += den100;
				dv -> normal[0] += den101;
				dv -> normal[1] -= den101;
				dv -> normal[2] += den101;
				dv -> normal[0] -= den110;
				dv -> normal[1] += den110;
				dv -> normal[2] += den110;
				dv -> normal[0] += den111;
				dv -> normal[1] += den111;
				dv -> normal[2] += den111;
				dv -> average = (den000+den001+den010+den011+den100+den101+den110+den111) / 8;
				if      (den110 < ctrlev && ctrlev < den111) {
					dv0 = denvtxs + idx000; if (0 == dv0 -> number) dv0 -> number = ++n_vertex;
					dv1 = denvtxs + idx010; if (0 == dv1 -> number) dv1 -> number = ++n_vertex;
					dv2 = denvtxs + idx110; if (0 == dv2 -> number) dv2 -> number = ++n_vertex;
					dv3 = denvtxs + idx100; if (0 == dv3 -> number) dv3 -> number = ++n_vertex;
				}
				else if (den110 > ctrlev && ctrlev > den111) {
					dv0 = denvtxs + idx000; if (0 == dv0 -> number) dv0 -> number = ++n_vertex;
					dv1 = denvtxs + idx100; if (0 == dv1 -> number) dv1 -> number = ++n_vertex;
					dv2 = denvtxs + idx110; if (0 == dv2 -> number) dv2 -> number = ++n_vertex;
					dv3 = denvtxs + idx010; if (0 == dv3 -> number) dv3 -> number = ++n_vertex;
				}
				if (den101 < ctrlev && ctrlev < den111) {
					dv0 = denvtxs + idx000; if (0 == dv0 -> number) dv0 -> number = ++n_vertex;
					dv1 = denvtxs + idx001; if (0 == dv1 -> number) dv1 -> number = ++n_vertex;
					dv2 = denvtxs + idx101; if (0 == dv2 -> number) dv2 -> number = ++n_vertex;
					dv3 = denvtxs + idx100; if (0 == dv3 -> number) dv3 -> number = ++n_vertex;
				}
				else if (den101 > ctrlev && ctrlev > den111) {
					dv0 = denvtxs + idx000; if (0 == dv0 -> number) dv0 -> number = ++n_vertex;
					dv1 = denvtxs + idx100; if (0 == dv1 -> number) dv1 -> number = ++n_vertex;
					dv2 = denvtxs + idx101; if (0 == dv2 -> number) dv2 -> number = ++n_vertex;
					dv3 = denvtxs + idx001; if (0 == dv3 -> number) dv3 -> number = ++n_vertex;
				}
				if (den011 < ctrlev && ctrlev < den111) {
					dv0 = denvtxs + idx000; if (0 == dv0 -> number) dv0 -> number = ++n_vertex;
					dv1 = denvtxs + idx001; if (0 == dv1 -> number) dv1 -> number = ++n_vertex;
					dv2 = denvtxs + idx011; if (0 == dv2 -> number) dv2 -> number = ++n_vertex;
					dv3 = denvtxs + idx010; if (0 == dv3 -> number) dv3 -> number = ++n_vertex;
				}
				else if (den011 > ctrlev && ctrlev > den111) {
					dv0 = denvtxs + idx000; if (0 == dv0 -> number) dv0 -> number = ++n_vertex;
					dv1 = denvtxs + idx010; if (0 == dv1 -> number) dv1 -> number = ++n_vertex;
					dv2 = denvtxs + idx011; if (0 == dv2 -> number) dv2 -> number = ++n_vertex;
					dv3 = denvtxs + idx001; if (0 == dv3 -> number) dv3 -> number = ++n_vertex;
				}
			}
		}
	}
	if (n_vertex <= 0) {
		set_error1 ("no vertices at surface of density");
		return (0);
	}
	*denvtx_handle = denvtxs;
	return(n_vertex);
}



int polygonize_density (struct surface *den, double ctrlev)
{
	long n_phnvtx, n_polygon, n_edge;
	struct denvtx *denvtxs;
	char message[MAXLINE];

	n_phnvtx = collect_den_vtx0 (den, ctrlev, &denvtxs);
	sprintf (message, "%8ld density surface vertices", n_phnvtx);
	inform (message);
	if (n_phnvtx == 0) return (0);
	if (error()) return (0);
	n_phnvtx = collect_den_vtx1 (den, denvtxs, n_phnvtx);
	if (n_phnvtx == 0) return (0);
	if (error()) return (0);
	n_polygon = collect_den_pgn (den, denvtxs, ctrlev);
	sprintf (message, "%8ld density surface polygons", n_polygon);
	inform (message);
	if (n_polygon == 0) return (0);
	if (error()) return (0);
	free_objects (DENVTX, (short *) denvtxs);
	denvtxs = NULL;
	n_edge = collect_edges (den);
	sprintf (message, "%8ld density surface edges", n_edge);
	inform (message);
	if (n_edge == 0) return (0);
	do_polyhedron_bounds (den);
	if (error ()) return (0);
	return (1);
}

long collect_den_vtx1 (struct surface *den, struct denvtx *denvtxs, long n_phnvtx)
{
	long ix, ix0, iy, iy0, iz, iz0;
	long idx000, i;
	int k;
	long size1, size2;
	struct denvtx *dv;
	struct phnvtx **phnvtx_handles, **iv;
	struct phnvtx *vtx, *pv;

	size1 = den -> width[0];
	size2 = size1 * den -> width[1];
	phnvtx_handles = (struct phnvtx **)
		allocate_pointers (PHNVTX, n_phnvtx);
	if (phnvtx_handles == NULL) {
		set_error1 ("no memory for phnvtx_handles");
		return (0);
	}
	for (i = 0; i < n_phnvtx; i++) {
		pv = allocate_phnvtx ();
		if (pv == NULL) return (0);
		*(phnvtx_handles + i) = pv;
	}
	for (iz = den -> origin[2]+1; iz <= den -> limit[2]; iz++) {
		for (iy = den -> origin[1]+1; iy <= den -> limit[1]; iy++) {
			for (ix = den -> origin[0]+1; ix <= den -> limit[0]; ix++) {
				ix0 = ix - den -> origin[0] - 1;
				iy0 = iy - den -> origin[1] - 1;
				iz0 = iz - den -> origin[2] - 1;
				idx000 = iz0 * size2 + iy0 * size1 + ix0;
				dv = denvtxs + idx000;
				if (dv -> number <= 0) continue;
				iv = phnvtx_handles + dv -> number - 1;
				vtx = *iv;
				for (k = 0; k < 3; k++)
					vtx -> outward[k] = -dv -> normal[k];
				if (!normalize(vtx -> outward)) return (0);
				vtx -> center[0] = den -> cube_width * ix;
				vtx -> center[1] = den -> cube_width * iy;
				vtx -> center[2] = den -> cube_width * iz;
			}
		}
	}
	den -> phnvtx_handles = phnvtx_handles;
	den -> n_phnvtx = n_phnvtx;
	return (n_phnvtx);
}

long collect_den_pgn (struct surface *den, struct denvtx *denvtxs, double ctrlev)
{
	long ix, ix0, ix1, iy, iy0, iy1, iz, iz0, iz1;
	long idx000, idx001, idx010, idx011;
	long idx100, idx101, idx110, idx111;
	long ndv0, ndv1, ndv2, ndv3;
	long size1, size2;
	long n_polygon, i;
	double den011;
	double den101;
	double den110;
	double den111;
	char message[MAXLINE];
	struct denvtx *dv0, *dv1, *dv2, *dv3;
	struct polygon **polygon_handles, **ip;
	struct polygon *poly;

	size1 = den -> width[0];
	size2 = size1 * den -> width[1];
	/* set up faces */
	n_polygon = 0;
	for (iz = den -> origin[2]+1; iz < den -> limit[2]; iz++) {
		for (iy = den -> origin[1]+1; iy < den -> limit[1]; iy++) {
			for (ix = den -> origin[0]+1; ix < den -> limit[0]; ix++) {
				ix0 = ix - den -> origin[0] - 1; ix1 = ix0 + 1;
				iy0 = iy - den -> origin[1] - 1; iy1 = iy0 + 1;
				iz0 = iz - den -> origin[2] - 1; iz1 = iz0 + 1;
				idx000 = iz0 * size2 + iy0 * size1 + ix0;
				idx001 = iz0 * size2 + iy0 * size1 + ix1;
				idx010 = iz0 * size2 + iy1 * size1 + ix0;
				idx011 = iz0 * size2 + iy1 * size1 + ix1;
				idx100 = iz1 * size2 + iy0 * size1 + ix0;
				idx101 = iz1 * size2 + iy0 * size1 + ix1;
				idx110 = iz1 * size2 + iy1 * size1 + ix0;
				idx111 = iz1 * size2 + iy1 * size1 + ix1;
				den011 = *(den -> densities + idx011);
				den101 = *(den -> densities + idx101);
				den110 = *(den -> densities + idx110);
				den111 = *(den -> densities + idx111);
				if      (den110 < ctrlev && ctrlev < den111) n_polygon++;
				else if (den110 > ctrlev && ctrlev > den111) n_polygon++;
				if (den101 < ctrlev && ctrlev < den111) n_polygon++;
				else if (den101 > ctrlev && ctrlev > den111) n_polygon++;
				if (den011 < ctrlev && ctrlev < den111) n_polygon++;
				else if (den011 > ctrlev && ctrlev > den111) n_polygon++;
			}
		}
	}
	if (n_polygon <= 0) {
		set_error1 ("collect_den_pgn: no polygons at surface of density");
		return (0);
	}
	polygon_handles = (struct polygon **)
		allocate_pointers (POLYGON, n_polygon);
	if (polygon_handles == NULL) {
		set_error1 ("collect_den_pgn: not enough memory for density polygon pointers");
		return (0);
	}
	for (i = 0; i < n_polygon; i++) {
		poly = allocate_polygon ();
		if (poly == NULL) {
			set_error1 ("collect_den_pgn: allocate polygon fails");
			sprintf (message, "n_polygon = %8ld, this = %8ld", n_polygon, i);
			set_error2 (message);
			return (0);
		}
		*(polygon_handles + i) = poly;
	}
	ip = polygon_handles;
	for (iz = den -> origin[2]+1; iz < den -> limit[2]; iz++) {
		for (iy = den -> origin[1]+1; iy < den -> limit[1]; iy++) {
			for (ix = den -> origin[0]+1; ix < den -> limit[0]; ix++) {
				ix0 = ix - den -> origin[0] - 1; ix1 = ix0 + 1;
				iy0 = iy - den -> origin[1] - 1; iy1 = iy0 + 1;
				iz0 = iz - den -> origin[2] - 1; iz1 = iz0 + 1;
				idx000 = iz0 * size2 + iy0 * size1 + ix0;
				idx001 = iz0 * size2 + iy0 * size1 + ix1;
				idx010 = iz0 * size2 + iy1 * size1 + ix0;
				idx011 = iz0 * size2 + iy1 * size1 + ix1;
				idx100 = iz1 * size2 + iy0 * size1 + ix0;
				idx101 = iz1 * size2 + iy0 * size1 + ix1;
				idx110 = iz1 * size2 + iy1 * size1 + ix0;
				idx111 = iz1 * size2 + iy1 * size1 + ix1;
				den011 = *(den -> densities + idx011);
				den101 = *(den -> densities + idx101);
				den110 = *(den -> densities + idx110);
				den111 = *(den -> densities + idx111);
				if (ip - polygon_handles >= n_polygon) break;
				poly = *ip;
				if (poly == NULL) {
					set_error1 ("collect_den_pgn: null polygon pointer");
					return (0);
				}
				if      (den110 < ctrlev && ctrlev < den111) {
					dv0 = denvtxs + idx000; ndv0 = dv0 -> number;
					dv1 = denvtxs + idx010; ndv1 = dv1 -> number;
					dv2 = denvtxs + idx110; ndv2 = dv2 -> number;
					dv3 = denvtxs + idx100; ndv3 = dv3 -> number;
					if (ip - polygon_handles >= n_polygon) {
						set_error1 ("polygon array overflow");
						return (0);
					}
					poly -> n_side = 4;
					poly -> material_index = 0;
					if (ndv0 <= 0 || ndv1 <= 0 || ndv2 <= 0 || ndv3 <= 0) {
						sprintf (message, "collect_den_pgn: ix = %4ld iy = %4ld iz = %4ld null vertex number", ix0, iy0, iz0);
						set_error1 (message);
						return (0);
					}
					poly -> vertex_index[0] = ndv0 - 1;
					poly -> vertex_index[1] = ndv1 - 1;
					poly -> vertex_index[2] = ndv2 - 1;
					poly -> vertex_index[3] = ndv3 - 1;
					ip++;
				}
				else if (den110 > ctrlev && ctrlev > den111) {
					dv0 = denvtxs + idx000; ndv0 = dv0 -> number;
					dv1 = denvtxs + idx100; ndv1 = dv1 -> number;
					dv2 = denvtxs + idx110; ndv2 = dv2 -> number;
					dv3 = denvtxs + idx010; ndv3 = dv3 -> number;
					if (ip - polygon_handles >= n_polygon) {
						set_error1 ("polygon array overflow");
						return (0);
					}
					poly -> n_side = 4;
					poly -> material_index = 0;
					if (ndv0 <= 0 || ndv1 <= 0 || ndv2 <= 0 || ndv3 <= 0) {
						sprintf (message, "collect_den_pgn: ix = %4ld iy = %4ld iz = %4ld null vertex number", ix0, iy0, iz0);
						set_error1 (message);
						return (0);
					}
					poly -> vertex_index[0] = ndv0 - 1;
					poly -> vertex_index[1] = ndv1 - 1;
					poly -> vertex_index[2] = ndv2 - 1;
					poly -> vertex_index[3] = ndv3 - 1;
					ip++;
				}
				if (ip - polygon_handles >= n_polygon) break;
				poly = *ip;
				if (poly == NULL) return (0);
				if (den101 < ctrlev && ctrlev < den111) {
					dv0 = denvtxs + idx000; ndv0 = dv0 -> number;
					dv1 = denvtxs + idx001; ndv1 = dv1 -> number;
					dv2 = denvtxs + idx101; ndv2 = dv2 -> number;
					dv3 = denvtxs + idx100; ndv3 = dv3 -> number;
					if (ip - polygon_handles >= n_polygon) {
						set_error1 ("polygon array overflow");
						return (0);
					}
					poly -> n_side = 4;
					poly -> material_index = 0;
					if (ndv0 <= 0 || ndv1 <= 0 || ndv2 <= 0 || ndv3 <= 0) {
						sprintf (message, "collect_den_pgn: ix = %4ld iy = %4ld iz = %4ld null vertex number", ix0, iy0, iz0);
						set_error1 (message);
						return (0);
					}
					poly -> vertex_index[0] = ndv0 - 1;
					poly -> vertex_index[1] = ndv1 - 1;
					poly -> vertex_index[2] = ndv2 - 1;
					poly -> vertex_index[3] = ndv3 - 1;
					ip++;
				}
				else if (den101 > ctrlev && ctrlev > den111) {
					dv0 = denvtxs + idx000; ndv0 = dv0 -> number;
					dv1 = denvtxs + idx100; ndv1 = dv1 -> number;
					dv2 = denvtxs + idx101; ndv2 = dv2 -> number;
					dv3 = denvtxs + idx001; ndv3 = dv3 -> number;
					if (ip - polygon_handles >= n_polygon) {
						set_error1 ("polygon array overflow");
						return (0);
					}
					poly -> n_side = 4;
					poly -> material_index = 0;
					if (ndv0 <= 0 || ndv1 <= 0 || ndv2 <= 0 || ndv3 <= 0) {
						sprintf (message, "collect_den_pgn: ix = %4ld iy = %4ld iz = %4ld null vertex number", ix0, iy0, iz0);
						set_error1 (message);
						return (0);
					}
					poly -> vertex_index[0] = ndv0 - 1;
					poly -> vertex_index[1] = ndv1 - 1;
					poly -> vertex_index[2] = ndv2 - 1;
					poly -> vertex_index[3] = ndv3 - 1;
					ip++;
				}
				if (ip - polygon_handles >= n_polygon) break;
				poly = *ip;
				if (poly == NULL) return (0);
				if (den011 < ctrlev && ctrlev < den111) {
					dv0 = denvtxs + idx000; ndv0 = dv0 -> number;
					dv1 = denvtxs + idx001; ndv1 = dv1 -> number;
					dv2 = denvtxs + idx011; ndv2 = dv2 -> number;
					dv3 = denvtxs + idx010; ndv3 = dv3 -> number;
					if (ip - polygon_handles >= n_polygon) {
						set_error1 ("polygon array overflow");
						return (0);
					}
					poly -> n_side = 4;
					poly -> material_index = 0;
					if (ndv0 <= 0 || ndv1 <= 0 || ndv2 <= 0 || ndv3 <= 0) {
						sprintf (message, "collect_den_pgn: ix = %4ld iy = %4ld iz = %4ld null vertex number", ix0, iy0, iz0);
						set_error1 (message);
						return (0);
					}
					poly -> vertex_index[0] = ndv0 - 1;
					poly -> vertex_index[1] = ndv1 - 1;
					poly -> vertex_index[2] = ndv2 - 1;
					poly -> vertex_index[3] = ndv3 - 1;
					ip++;
				}
				else if (den011 > ctrlev && ctrlev > den111) {
					dv0 = denvtxs + idx000; ndv0 = dv0 -> number;
					dv1 = denvtxs + idx010; ndv1 = dv1 -> number;
					dv2 = denvtxs + idx011; ndv2 = dv2 -> number;
					dv3 = denvtxs + idx001; ndv3 = dv3 -> number;
					if (ip - polygon_handles >= n_polygon) {
						set_error1 ("polygon array overflow");
						return (0);
					}
					poly -> n_side = 4;
					poly -> material_index = 0;
					if (ndv0 <= 0 || ndv1 <= 0 || ndv2 <= 0 || ndv3 <= 0) {
						sprintf (message, "collect_den_pgn: ix = %4ld iy = %4ld iz = %4ld null vertex number", ix0, iy0, iz0);
						set_error1 (message);
						return (0);
					}
					poly -> vertex_index[0] = ndv0 - 1;
					poly -> vertex_index[1] = ndv1 - 1;
					poly -> vertex_index[2] = ndv2 - 1;
					poly -> vertex_index[3] = ndv3 - 1;
					ip++;
				}
			}
		}
	}
	if (ip - polygon_handles != n_polygon) {
		set_error1 ("inconsistent polygon count");
		return (0);
	}
	den -> polygon_handles = polygon_handles;
	den -> n_polygon = n_polygon;
	return (n_polygon);
}

