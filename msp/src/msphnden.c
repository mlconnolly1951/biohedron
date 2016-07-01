/* Molecular Surface Package */
/* Copyright 1995 by Michael L. Connolly */
/* March 15, 2000 */

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"


/* CALLING   *   *   *   *   *   CALLING   */


struct surface *polyhedron_to_density (struct surface *poly, double cube_width)
{
	int k, result;
	int ibounds[2][3], width[3];
	long dimension;
	long x, y, z, nx, ny, nz, idx;
	long n_vertex, n_edge, n_triangle;
	long dimensions[3];
	long *triangles, *edges;
	float *densities, *f;
	double *vertices;
	double bounds[2][3];
	struct msdata *msd;
	struct surface *den;
	
	result = do_poly_bounds (poly);
	if (!result || error ()) return (NULL);
	for (k = 0; k < 3; k++) {
		ibounds[0][k] = floor (poly -> bounds[0][k] / cube_width) - 1;
		ibounds[1][k] = ceil (poly -> bounds[1][k] / cube_width) + 1;
		width[k] = ibounds[1][k] - ibounds[0][k];
	}

	den = create_density (cube_width, ibounds[0], width);
	if (error ()) return (NULL);
	make_phn_simple (poly);
	if (error ()) return (NULL);
	vertices = poly -> vertex_centers;
	if (vertices == NULL) {
		set_error1 ("null polyhedron vertex_centers");
		return (NULL);
	}
	edges = poly -> edgvtx;
	if (edges == NULL) {
		set_error1 ("null polyhedron edges");
		return (NULL);
	}
	triangles = poly -> triedgvtx;
	if (triangles == NULL) {
		set_error1 ("null polyhedron triangles");
		return (NULL);
	}
	dimension = 3;
	n_vertex = poly -> n_phnvtx;
	if (n_vertex <= 0) {
		set_error1 ("no polyhedron vertices");
		return (NULL);
	}
	n_edge = poly -> n_phnedg;
	if (n_edge <= 0) {
		set_error1 ("no polyhedron edges");
		return (NULL);
	}
	n_triangle = poly -> n_phntri;
	if (n_triangle <= 0) {
		set_error1 ("no polyhedron triangles");
		return (NULL);
	}
	for (k = 0; k < 3; k++) {
		dimensions[k] = den -> width[k];
		bounds[0][k] = den -> origin[k] * cube_width;
		bounds[1][k] = den -> limit[k] * cube_width;
	}
	msd = doPlex (dimension, dimensions, bounds, n_vertex, n_edge, n_triangle, vertices, edges, triangles);
	if (msd == NULL) {
		inform ("doPlex returns null");
		return (NULL);
	}
	if (msd -> narray != 1) {
		set_error1 ("wrong number of arrays returned by doPlex");
		return (NULL);
	}
	/* transfer results to local arrays */
	nx = dimensions[0]; ny = dimensions[1]; nz = dimensions[2];
	densities = den -> densities;
	for (z = 0; z < nz; z++) {
		for (y = 0; y < ny; y++) {
			for (x = 0; x < nx; x++) {
				idx = z * (ny * nx) + y * nx + x;
				f = msd -> floats[0] + idx;	/* density first now */
				*(densities + idx) = *f;
			}
		}
	}
	if (!freemsdata (msd)) return (NULL);
	if (!free_phn_simple (poly)) return (NULL);
	return (den);
}

int free_phn_simple (struct surface *phn)
{
	if (phn == NULL) return (0);
	if (phn -> n_phnvtx == 0) return (0);
	if (phn -> n_phnedg == 0) return (0);
	if (phn -> n_phntri == 0) return (0);
	if (phn -> vertex_centers == NULL) return (0);
	if (phn -> vertex_normals == NULL) return (0);
	if (phn -> edgvtx == NULL) return (0);
	if (phn -> triedgvtx == NULL) return (0);
	free_doubles (phn -> vertex_centers, 0, VERTEX_CENTERS);
	free_doubles (phn -> vertex_normals, 0, VERTEX_NORMALS);
	free_longs (phn -> edgvtx, 0, EDGVTX);
	free_longs (phn -> triedgvtx, 0, TRIEDGVTX);
	return (1);
}


/* Convert a polyhedron to a density */

struct msdata *doPlex (long dimension, long dimensions[3], double bounds[2][3], long n_vertex, long n_edge, long n_triangle, double *vertices, long *edges, long *triangles)
{
	int j, k;
	long nx, ny, nz;
	double xwidth, ywidth, zwidth;
	struct Plex *plex;
	struct msdata *msd;
	
	plex = (struct Plex *) allocate_object (PLEX);
	if (plex == NULL) {
		set_error1 ("doPlex: not enough memory to allocate Plex");
		return (NULL);
	}
	plex -> dimension = (integer) dimension;
	for (k = 0; k < 3; k++)
		plex -> dimensions[k] = (integer) dimensions[k];
	if (dimension == 2) plex -> dimensions[2] = (integer) 1;
	nx = (plex -> dimensions[0]);
	ny = (plex -> dimensions[1]);
	nz = (plex -> dimensions[2]);
	plex -> n_cube = (nx * ny * nz);
	for (j = 0; j < 2; j++)
		for (k = 0; k < 3; k++)
			plex -> bounds[j][k] = bounds[j][k];
	plex -> maxhash = (integer) 256L;
	plex -> maxvertex = (256L + plex -> n_cube + 4 * n_vertex);
	plex -> maxedge = (4 * n_edge + plex -> maxvertex + 4 * n_triangle);
	plex -> maxtriangle = (n_triangle + plex -> maxedge);
	if (sizeof (integer) == 2) {
		if (plex -> n_cube > 32767) {
			set_error1 ("doPlex: n_cube > max integer");
			return (NULL);
		}
		if (plex -> maxvertex > 32767) {
			set_error1 ("doPlex: maxvertex > max integer");
			return (NULL);
		}
		if (plex -> maxedge > 32767) {
			set_error1 ("doPlex: maxedge > max integer");
			return (NULL);
		}
		if (plex -> maxtriangle > 32767) {
			set_error1 ("doPlex: maxtriangle > max integer");
			return (NULL);
		}
	}
	plex -> glassblock = 1024L;
	plex -> mvertex = n_vertex;
	plex -> medge = n_edge;
	plex -> mtriangle = n_triangle;
	plex -> n_vertex = n_vertex;
	plex -> n_edge = n_edge;
	plex -> n_triangle = n_triangle;
	plex -> vertices = vertices;
	plex -> edges = edges;
	plex -> triangles = triangles;
	if (plex -> dimension == 2 && plex -> n_triangle > 0) return (NULL);
	if (plex -> dimension == 2 && plex -> triangles != NULL) return (NULL);
	if (plex -> vertices == NULL) return (NULL);
	if (plex -> edges == NULL) return (NULL);
	if (plex -> n_vertex == 0) return (NULL);
	if (plex -> n_edge == 0) return (NULL);
	xwidth = (plex -> bounds[1][0] - plex -> bounds[0][0]);
	ywidth = (plex -> bounds[1][1] - plex -> bounds[0][1]);
	zwidth = (plex -> bounds[1][2] - plex -> bounds[0][2]);
	if (xwidth <= 0.0 || ywidth <= 0.0) return (NULL);
	if (plex -> dimension == 2)
		plex -> cubeVolume = xwidth * ywidth;
	else plex -> cubeVolume = xwidth * ywidth * zwidth;
	plex -> cubeVolume /= plex -> n_cube;
	if (plex -> cubeVolume <= 0.0) return (NULL);
	informd ("initializePlex");
	if (!initializePlex (plex)) return (NULL);
	informd ("computeVolume");
	if (!computeVolume (plex)) return (NULL);
	informd ("storePoly");
	if (!storePoly (plex)) return (NULL);
	informd ("classifyCubes");
	if (!classifyCubes (plex)) return (NULL);
	informd ("createProvinces");
	if (!createProvinces(plex)) return (NULL);
	informd ("polyhedronProvince");
	if (!polyhedronProvince(plex)) return (NULL);
	informd ("densityProvince");
	if (!densityProvince(plex)) return (NULL);
	informd ("createBorderVertices");
	if (!createBorderVertices (plex)) return (NULL);
	informd ("createBorderEdges");
	if (!createBorderEdges (plex)) return (NULL);
	informd ("createBorderTriangles");
	if (!createBorderTriangles (plex)) return (NULL);
	informd ("allApices");
	if (!allApices (plex)) return (NULL);
	informd ("provinceCenter");
	if (!provinceCenter (plex)) return (NULL);
	informd ("computeJoins");
	if (!computeJoins(plex)) return (NULL);
	informd ("computeSurfaceDensity");
	if (!computeSurfaceDensity(plex)) return (NULL);
	informd ("transferdata");
	msd = transferdata (plex);
	if (msd == NULL) return (NULL);
	informd ("freePlex");
	if (!freePlex (plex)) return (NULL);
	return (msd);
}



/* INITIALIZATION   *   *   *   *   *   INITIALIZATION   */

/* Initialize Plex: store parameters and allocate memory */

int initializePlex(struct Plex *plex)
{
	long nitem;
	char message[MAXLINE];
	
	nitem = MaxPlexType * (long) plex -> maxhash;
	plex -> plexis = (struct Plexi *) allocate_objects (PLEXI, nitem);
	if (plex -> plexis == NULL) {
		set_error1 ("initializePlex: not enough memory for plexis");
		return (0);
	}
	nitem = plex -> maxvertex;
	plex -> plexvertices = (struct PlexVertex *)
		allocate_objects (PLEXVERTEX, nitem);
	if (plex -> plexvertices == NULL) {
		set_error1 ("initializePlex: not enough memory for plexvertices");
		return (0);
	}
	nitem = plex -> maxedge;
	plex -> plexedges = (struct PlexEdge *) allocate_objects (PLEXEDGE, nitem);
	if (plex -> plexedges == NULL) {
		set_error1 ("initializePlex: not enough memory for plexedges");
		return (0);
	}
	nitem = plex -> maxtriangle;
	if (plex -> maxtriangle > 0) {
		plex -> plextriangles = (struct PlexTriangle *)
			allocate_objects (PLEXTRIANGLE, nitem);
		if (plex -> plextriangles == NULL) {
			set_error1 ("initializePlex: not enough memory for plextriangles");
			return (0);
		}
	}
	else plex -> plextriangles = NULL;
	nitem = plex -> n_cube;
	plex -> plexcubes = (struct PlexCube *) allocate_objects (PLEXCUBE, nitem);
	if (plex -> plexcubes == NULL) {
		set_error1 ("initializePlex: not enough memory for plexcubes");
		return (0);
	}
	sprintf (message, "%8ld plex cubes allocated", nitem);
	inform (message);
	return(1);
}

int initializeGlass (struct Plex *plex)
{
	int result;
	result = moreGlass (plex);
	return (result);
}

struct msdata *transferdata (struct Plex *plex)
{
	long nx, ny, nz, x, y, z, idx, n;
	long ndensity, nplexline, count;
	double postot, negtot, value;
	float *densities, *plexlines;
	struct msdata *msd;
	char message[MAXLINE];

	nx = plex -> dimensions[0];
	ny = plex -> dimensions[1];
	nz = plex -> dimensions[2];

	postot = 0.0; negtot = 0.0;
	ndensity = nx * ny * nz;
	densities = allocate_floats (ndensity);
	if (densities == NULL) {
		set_error1 ("transferData: memory failure for densities");
		return (NULL);
	}
	for (z = 0; z < nz; z++)
		for (y = 0; y < ny; y++)
			for (x = 0; x < nx; x++) {
				idx = z * (nx * ny) + y * nx + x;
				value = (plex -> plexcubes + idx) -> occupancy;
				*(densities + idx) = value;
				if (value > 0.0) postot += value * plex -> cubeVolume;
				else if (value < 0.0) negtot += value * plex -> cubeVolume;
			}

	if (plex -> dimension == 3) {
		sprintf (message, "%8.3f positive cube volume", postot);
		inform (message);
		sprintf (message, "%8.3f negative cube volume", negtot);
		inform (message);
		sprintf (message, "%8.3f total    cube volume", postot + negtot);
		inform (message);
		sprintf (message, "%8.3f polyhedron    volume", plex -> polyvolume);
		inform (message);
	}

	/* transfer data to msdata */
	msd = newmsdata ();
	if (msd == NULL) {
		set_error1 ("transferData: no memory for msd");
		return (NULL);
	}
	/* just add densities, since not enough RAM on SGI Indy for lines */
	if (!addfloats (msd, 1L, ndensity, densities)) {
		set_error1 ("transferData: addfloats fails for densities");
		return (NULL);
	}
	/* return */
	if (msd == msd) return (msd);

	count = 0;
	nplexline = plex -> nglass[PolyEdge];
	if (nplexline == 0) {
		set_error1 ("transferData: no PolyEdge");
		return (NULL);
	}
	sprintf (message, "%8ld polygon edges", nplexline);
	inform (message);
	plexlines = allocate_floats (nplexline * 6);
	if (plexlines == NULL){
		set_error1 ("transferData: memory failure for PolyEdge plexlines");
		return (NULL);
	}
	n = transfertype (plex, PolyEdge, plexlines, count);
	if (n != plex -> nglass[PolyEdge]) {
		set_error1 ("transferData: inconsistent PolyEdge count(1)");
		return (NULL);
	}
	count += n;
	if (count != nplexline) {
		set_error1 ("transferData: inconsistent PolyEdge count(2)");
		return (NULL);
	}
	if (!addfloats (msd, 6L, nplexline, plexlines)) {
		set_error1 ("transferData: addfloats fails for PolyEdge plexlines");
		return (NULL);
	}

	count = 0;
	nplexline = plex -> nglass[DenEdge] + plex -> nglass[BorderEdge];
	sprintf (message, "%8ld density and border edges", nplexline);
	inform (message);
	if (nplexline == 0) {
		plexlines = NULL;
	}
	else {
		plexlines = allocate_floats (nplexline * 6);
		if (plexlines == NULL) {
			set_error1 ("transferData: memory failure for DenEdge and BorderEdge plexlines");
			return (NULL);
		}
		if (plex -> nglass[DenEdge] > 0) {
			n = transfertype (plex, DenEdge, plexlines, count);
			if (n != plex -> nglass[DenEdge]) {
				set_error1 ("transferData: inconsistent DenEdge count");
				return (NULL);
			}
			count += n;
		}
		if (plex -> nglass[BorderEdge] > 0) {
			n = transfertype (plex, BorderEdge, plexlines, count);
			if (n != plex -> nglass[BorderEdge]) {
				set_error1 ("transferData: inconsistent BorderEdge count");
				return (NULL);
			}
			count += n;
		}
	}
	if (count != nplexline) {
		sprintf (message, "transferData: count (%8ld) != nplexline (%8ld)",
			count, nplexline);
		set_error1 (message);
		return (NULL);
	}
	if (!addfloats (msd, 6L, nplexline, plexlines)) {
		set_error1 ("transferData: addfloats fails for DenEdge & BorderEdge");
		return (NULL);
	}

	return (msd);
}

long transfertype (struct Plex *plex, enum PlexType type, float *plexlines, long start)
{
	struct Plexi *plexi;
	struct Glass *glass;
	struct PlexVertex *plexvtxs[2];
	long h, maxhash, idx;
	long vns[4];
	int k;

	maxhash = plex -> maxhash;
	idx = start;
	for (h = 0; h < maxhash; h++) {
		plexi = plex -> plexis + ((int) type) * maxhash + h;
		for (glass = plexi -> head; glass != NULL; glass = glass -> next) {
			vns[0] = glass -> vertexnumbers[0];
			vns[1] = glass -> vertexnumbers[1];
			for (k = 0; k < 2; k++)
				plexvtxs[k] = plex -> plexvertices + vns[k] - 1;
			for (k = 0; k < 3; k++) {
				*(plexlines + 6 * idx + k) = plexvtxs[0] -> center[k];
				*(plexlines + 6 * idx + 3 + k) = plexvtxs[1] -> center[k];
			}
			if (plex -> dimension == 2) {
				*(plexlines + 6 * idx + 2) = 0.0;
				*(plexlines + 6 * idx + 5) = 0.0;
			}
			idx++;
		}
	}
	return (idx-start);
}

int freePlex (struct Plex *plex)
{
	struct Glass *glasses;
	struct GlassBlock *gb, *gbnext;

	free_objects (PLEXI, (short *) (plex -> plexis));
	plex -> plexis = NULL;
	free_objects (GLASS, (short *) (plex -> glasses));
	plex -> glasses = NULL; plex -> freeglass = NULL;
	free_objects (PLEXVERTEX, (short *) (plex -> plexvertices));
	plex -> plexvertices = NULL; plex -> n_vertex = 0;
	free_objects (PLEXJOIN, (short *) (plex -> joins));
	plex -> joins = NULL; plex -> n_join = 0;
	free_objects (PLEXCUBE, (short *) (plex -> plexcubes));
	plex -> plexcubes = NULL; plex -> n_cube = 0;
	free_objects (PROVINCE, (short *) (plex -> provinces));
	plex -> provinces = NULL;
	for (gb = plex -> head_glassblock; gb != NULL; gb = gbnext) {
		gbnext = gb -> next;
		glasses = gb -> glasses;
		free_objects (GLASS, (short *) glasses);
		free_object (GLASS_BLOCK, (short *) gb);
	}
	free_object (PLEX, (short *) plex);
	free_cache (GLASS_BLOCK);
	free_cache (PLEX);
	return (1);
}



int typeDim (enum PlexType type)
{
	switch (type) {
	case PolyVertex:
	case DenVertex:
	case ProvVertex:
	case BorderVertex:
		return (1);
	case PolyEdge:
	case DenEdge:
	case BorderEdge:
	case ProvEdge:
		return (2);
	case PolyTriangle:
	case DenSquare:
	case BorderTriangle:
	case ProvTriangle:
		return (3);
	case ProvTetrahedron:
		return (4);
	}
	return (0);
}



/*     DENSITY     DENSITY     DENSITY     DENSITY     DENSITY     */

int computeSurfaceDensity (struct Plex *plex)
{
	int jj, j, k, dim, result;
	long n_join, t, idx, v;
	long nfull, npartial, nempty, npositive, nnegative;
	long ico[3];
	double simplex_vol, vpositive, vnegative;
	double polyvolume, denvolume, bordervolume, bordernegative;
	double borderborder, borderpoly, borderden;
	double fg[3], center[3];
	double vect[3][3], simplex[4][3];
	struct PlexCube *cube;
	struct PlexJoin *join;
	struct PlexVertex *plexvtx;
	enum PlexType type, subtype;
	char message[MAXLINE];
	
	dim = plex -> dimension;
	n_join = plex -> n_join;
	nfull = npartial = nempty = 0;
	npositive = nnegative = 0;
	vpositive = vnegative = 0.0;
	polyvolume = denvolume = bordervolume = bordernegative = 0.0;
	borderborder = borderden = borderpoly = 0.0;
	for (t = 0; t < n_join; t++) {
		join = plex -> joins + t;
		type = join -> type;
		subtype = join -> subtype;
		if (join -> orn == 0) {
			set_error1 ("computeSurfaceDensity: join -> orn == 0");
			return (0);
		}
		if (dim == 2) {
			for (k = 0; k < 3; k++)
				simplex[0][k] = join -> apex[k];
			for (j = 0; j < 2; j++) {
				v = join -> vertices[j];
				if (v == 0) break;
				plexvtx = plex -> plexvertices + v - 1;
				for (k = 0; k < 3; k++)
					simplex[j+1][k] = plexvtx -> center[k];
			}
			for (j = 0; j < 2; j++)
				for (k = 0; k < 3; k++)
					vect[j][k] = simplex[j+1][k] - simplex[0][k];
			simplex_vol = join -> orn * (vect[0][0] * vect[1][1] - vect[0][1] * vect[1][0]) / 2.0;
			if (simplex_vol > 0.0) npositive++;
			else if (simplex_vol < 0.0) nnegative++;
			if (simplex_vol > 0.0) vpositive += simplex_vol;
			else if (simplex_vol < 0.0) vnegative += simplex_vol;
			for (k = 0; k < 3; k++)
				center[k] = join -> apex[k];
			result = center2indices(plex, center, ico);
			if (!result) return (0);
			idx = ico[2] * plex -> dimensions[0] * plex -> dimensions[1] + ico[1] * plex -> dimensions[0] + ico[0];
			cube = plex -> plexcubes + idx;
			cube -> occupancy += simplex_vol / plex -> cubeVolume;
		}
		else if (dim == 3 && join -> vertices[3] == 0) {
			for (k = 0; k < 3; k++)
				simplex[0][k] = join -> apex[k];
			for (j = 0; j < 3; j++) {
				v = join -> vertices[j];
				if (v == 0) break;
				plexvtx = plex -> plexvertices + v - 1;
				for (k = 0; k < 3; k++)
					simplex[j+1][k] = plexvtx -> center[k];
			}
			for (j = 0; j < 3; j++)
				for (k = 0; k < 3; k++)
					vect[j][k] = simplex[j+1][k] - simplex[0][k];
			fg[0] = vect[0][1] * vect[1][2] - vect[0][2] * vect[1][1];
			fg[1] = vect[0][2] * vect[1][0] - vect[0][0] * vect[1][2];
			fg[2] = vect[0][0] * vect[1][1] - vect[0][1] * vect[1][0];
			simplex_vol = join -> orn * (vect[2][0] * fg[0] + vect[2][1] * fg[1] + vect[2][2] * fg[2]) / 6.0;
			if (simplex_vol > 0.0) npositive++;
			else if (simplex_vol < 0.0) nnegative++;
			if (simplex_vol > 0.0) vpositive += simplex_vol;
			else if (simplex_vol < 0.0) vnegative += simplex_vol;
			if (type == PolyTriangle) polyvolume += simplex_vol;
			else if (type == BorderTriangle) bordervolume += simplex_vol;
			if (simplex_vol < 0.0 && type == BorderTriangle) bordernegative += simplex_vol;
			if (simplex_vol < 0.0 && type == BorderTriangle && subtype == PolyEdge) borderpoly += simplex_vol;
			if (simplex_vol < 0.0 && type == BorderTriangle && subtype == BorderEdge) borderborder += simplex_vol;
			if (simplex_vol < 0.0 && type == BorderTriangle && subtype == DenEdge) borderden += simplex_vol;
			for (k = 0; k < 3; k++)
				center[k] = join -> apex[k];
			result = center2indices(plex, center, ico);
			if (!result) return (0);
			idx = ico[2] * plex -> dimensions[0] * plex -> dimensions[1] + ico[1] * plex -> dimensions[0] + ico[0];
			cube = plex -> plexcubes + idx;
			cube -> occupancy += simplex_vol / plex -> cubeVolume;
		}
		else if (dim == 3 && join -> vertices[3] > 0) {
			for (k = 0; k < 3; k++)
				simplex[0][k] = join -> apex[k];
			for (j = 0; j < 3; j++) {
				v = join -> vertices[j];
				if (v == 0) break;
				plexvtx = plex -> plexvertices + v - 1;
				for (k = 0; k < 3; k++)
					simplex[j+1][k] = plexvtx -> center[k];
			}
			for (j = 0; j < 3; j++)
				for (k = 0; k < 3; k++)
					vect[j][k] = simplex[j+1][k] - simplex[0][k];
			fg[0] = vect[0][1] * vect[1][2] - vect[0][2] * vect[1][1];
			fg[1] = vect[0][2] * vect[1][0] - vect[0][0] * vect[1][2];
			fg[2] = vect[0][0] * vect[1][1] - vect[0][1] * vect[1][0];
			simplex_vol = join -> orn *
				(vect[2][0] * fg[0] + vect[2][1] * fg[1] + vect[2][2] * fg[2]) / 6.0;
			if (simplex_vol > 0.0) npositive++;
			else if (simplex_vol < 0.0) nnegative++;
			if (simplex_vol > 0.0) vpositive += simplex_vol;
			else if (simplex_vol < 0.0) vnegative += simplex_vol;
			if (type == DenSquare) denvolume += simplex_vol;
			for (k = 0; k < 3; k++)
				center[k] = join -> apex[k];
			result = center2indices(plex, center, ico);
			if (!result) return (0);
			idx = ico[2] * plex -> dimensions[0] * plex -> dimensions[1] + ico[1] * plex -> dimensions[0] + ico[0];
			cube = plex -> plexcubes + idx;
			cube -> occupancy += simplex_vol / plex -> cubeVolume;
			for (k = 0; k < 3; k++)
				simplex[0][k] = join -> apex[k];
			for (j = 0; j < 3; j++) {
				jj = (j < 2) ? j + 2 : 0;
				v = join -> vertices[jj];
				if (v == 0) break;
				plexvtx = plex -> plexvertices + v - 1;
				for (k = 0; k < 3; k++)
					simplex[j+1][k] = plexvtx -> center[k];
			}
			for (j = 0; j < 3; j++)
				for (k = 0; k < 3; k++)
					vect[j][k] = simplex[j+1][k] - simplex[0][k];
			fg[0] = vect[0][1] * vect[1][2] - vect[0][2] * vect[1][1];
			fg[1] = vect[0][2] * vect[1][0] - vect[0][0] * vect[1][2];
			fg[2] = vect[0][0] * vect[1][1] - vect[0][1] * vect[1][0];
			simplex_vol = join -> orn *
				(vect[2][0] * fg[0] + vect[2][1] * fg[1] + vect[2][2] * fg[2]) / 6.0;
			if (simplex_vol > 0.0) npositive++;
			else if (simplex_vol < 0.0) nnegative++;
			if (simplex_vol > 0.0) vpositive += simplex_vol;
			else if (simplex_vol < 0.0) vnegative += simplex_vol;
			if (type == DenSquare) denvolume += simplex_vol;
			for (k = 0; k < 3; k++)
				center[k] = join -> apex[k];
			result = center2indices(plex, center, ico);
			if (!result) return (0);
			idx = ico[2] * plex -> dimensions[0] * plex -> dimensions[1] + ico[1] * plex -> dimensions[0] + ico[0];
			cube = plex -> plexcubes + idx;
			cube -> occupancy += simplex_vol / plex -> cubeVolume;
		}
	}
	/* force empty to 0, full to 1, and partial to between 0 and 1 */
	for (idx = 0; idx < plex -> n_cube; idx++) {
		cube = plex -> plexcubes + idx;
		if (cube -> type == FullCube) {
			cube -> occupancy = 1.0;
			nfull++;
		}
		else if (cube -> type == EmptyCube) {
			cube -> occupancy = 0.0;
			nempty++;
		}
		else if (cube -> type == PartialCube) {
			if (cube -> occupancy < 0.0) cube -> occupancy = 0.0;
			else if (cube -> occupancy > 1.0) cube -> occupancy = 1.0;
			npartial++;
		}
	}
	if (dim == 3) {
		sprintf (message, "%8ld full cubes", nfull);
		inform (message);
		sprintf (message, "%8ld partial cubes", npartial);
		inform (message);
		sprintf (message, "%8ld empty cubes", nempty);
		inform (message);
		sprintf (message, "%8.3f positive simplex volume (%6ld simplices)", vpositive, npositive);
		inform (message);
		sprintf (message, "%8.3f negative simplex volume (%6ld simplices)", vnegative, nnegative);
		inform (message);
		sprintf (message, "%8.3f PolyTriangle join volume", polyvolume);
		inform (message);
		sprintf (message, "%8.3f DenSquare join volume", denvolume);
		inform (message);
		sprintf (message, "%8.3f BorderTriangle join volume", bordervolume);
		inform (message);
		sprintf (message, "%8.3f BorderTriangle negative join volume", bordernegative);
		inform (message);
		sprintf (message, "%8.3f BorderTriangle PolyEdge negative join volume", borderpoly);
		inform (message);
		sprintf (message, "%8.3f BorderTriangle BorderEdge negative join volume", borderborder);
		inform (message);
		sprintf (message, "%8.3f BorderTriangle DenEdge negative join volume", borderden);
		inform (message);
	}
	return (1);
}

/* Molecular Surface Package */
/* Copyright 1995 by Michael L. Connolly */
