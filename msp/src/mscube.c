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


/* HASH TABLE    * * * * *    HASH TABLE   */

/* hash */

long hashGlass (long idnumbers[4], long maxhash)
{
	long s, h;
	s = idnumbers[0] + idnumbers[1] + idnumbers[2] + idnumbers[3];
	h = s % maxhash;
	return (h);
}

/* check whether id numbers match (up to acceptable permutation) */

int idmatch (enum PlexType type, long ida[4], long idb[4])
{
	switch (type) {
	case PolyVertex:
	case ProvVertex:
		return (ida[0] == idb[0]);
	case DenVertex:			/* x, y and z should not be permuted */
	case BorderVertex:
		if      (ida[0] == idb[0] && ida[1] == idb[1] && ida[2] == idb[2]) return (1);
		else return (0);
	case PolyEdge:
	case DenEdge:
	case BorderEdge:
	case ProvEdge:
		if (ida[0] == idb[0] && ida[1] == idb[1])
			return (1);
		else if (ida[0] == idb[1] && ida[1] == idb[0])
			return (-1);
		else return (0);
	case PolyTriangle:
	case BorderTriangle:
	case ProvTriangle:
		if      (ida[0] == idb[0] && ida[1] == idb[1] && ida[2] == idb[2]) return (1);
		else if (ida[0] == idb[1] && ida[1] == idb[2] && ida[2] == idb[0]) return (1);
		else if (ida[0] == idb[2] && ida[1] == idb[0] && ida[2] == idb[1]) return (1);
		else if (ida[0] == idb[2] && ida[1] == idb[1] && ida[2] == idb[0]) return (-1);
		else if (ida[0] == idb[1] && ida[1] == idb[0] && ida[2] == idb[2]) return (-1);
		else if (ida[0] == idb[0] && ida[1] == idb[2] && ida[2] == idb[1]) return (-1);
		else return (0);
	case DenSquare:
		if      (ida[0] == idb[0] && ida[1] == idb[1] && ida[2] == idb[2] && ida[3] == idb[3]) return (1);
		else if (ida[0] == idb[1] && ida[1] == idb[2] && ida[2] == idb[3] && ida[3] == idb[0]) return (1);
		else if (ida[0] == idb[2] && ida[1] == idb[3] && ida[2] == idb[0] && ida[3] == idb[1]) return (1);
		else if (ida[0] == idb[3] && ida[1] == idb[0] && ida[2] == idb[1] && ida[3] == idb[2]) return (1);
		else if (ida[0] == idb[3] && ida[1] == idb[2] && ida[2] == idb[1] && ida[3] == idb[0]) return (-1);
		else if (ida[0] == idb[2] && ida[1] == idb[1] && ida[2] == idb[0] && ida[3] == idb[3]) return (-1);
		else if (ida[0] == idb[1] && ida[1] == idb[0] && ida[2] == idb[3] && ida[3] == idb[2]) return (-1);
		else if (ida[0] == idb[0] && ida[1] == idb[3] && ida[2] == idb[2] && ida[3] == idb[1]) return (-1);
		else return (0);
	case ProvTetrahedron: /* fix later */
		if      (ida[0] == idb[0] && ida[1] == idb[1] && ida[2] == idb[2] && ida[3] == idb[3]) return (1);
		else if (ida[0] == idb[1] && ida[1] == idb[2] && ida[2] == idb[3] && ida[3] == idb[0]) return (1);
		else if (ida[0] == idb[2] && ida[1] == idb[3] && ida[2] == idb[0] && ida[3] == idb[1]) return (1);
		else if (ida[0] == idb[3] && ida[1] == idb[0] && ida[2] == idb[1] && ida[3] == idb[2]) return (1);
		else if (ida[0] == idb[3] && ida[1] == idb[2] && ida[2] == idb[1] && ida[3] == idb[0]) return (-1);
		else if (ida[0] == idb[2] && ida[1] == idb[1] && ida[2] == idb[0] && ida[3] == idb[3]) return (-1);
		else if (ida[0] == idb[1] && ida[1] == idb[0] && ida[2] == idb[3] && ida[3] == idb[2]) return (-1);
		else if (ida[0] == idb[0] && ida[1] == idb[3] && ida[2] == idb[2] && ida[3] == idb[1]) return (-1);
		else return (0);
	default:
		break;
	}
	return (0);
}

/* store entry (glass) into hashing table */

struct Glass *storeGlass (struct Plex *plex, enum PlexType type, long idnumbers[4], double center[3], double measure)
{
	int i;
	long idns[4], h, maxhash;
	struct Glass *glass;
	struct Plexi *plexi;
	
	maxhash = plex -> maxhash;
	for (i = 0; i < 4; i++) idns[i] = idnumbers[i];
	/* kludge */
	if (type == DenVertex || type == BorderVertex) idns[3] = 0;
	/* hash */
	h = hashGlass (idns, maxhash);
	/* index into head (plexi) array depends on type & hash value */
	plexi = plex -> plexis + ((long) type) * maxhash + h;
	/* grab an unused slot */
	glass = allocateGlass (plex);
	
	if (glass == NULL) {
		return (NULL);
	}
	glass -> next = plexi -> head;
	plexi -> head = glass;
	/* store info */
	for (i = 0; i < 4; i++)
		glass -> idnumbers[i] = idns[i];
	glass -> type = type;
	for (i = 0; i < 3; i++) glass -> center[i] = center[i];
	glass -> measure = measure;
	switch (type) {
	case ProvVertex:
	case PolyVertex:
		if (idns[0] == 0) {
			set_error1 ("storeGlass: zero vertex number (1)");
			return (NULL);
		}
		glass -> vertexnumbers[0] = idns[0];
		break;
	case PolyEdge:
	case DenEdge:
	case BorderEdge:
	case ProvEdge:
		if (idns[0] == 0 || idns[1] == 0) {
			set_error1 ("storeGlass: zero vertex number (2)");
			return (NULL);
		}
		glass -> vertexnumbers[0] = idns[0];
		glass -> vertexnumbers[1] = idns[1];
		break;
	case PolyTriangle:
	case BorderTriangle:
	case ProvTriangle:
		if (idns[0] == 0 || idns[1] == 0 || idns[2] == 0) {
			set_error1 ("storeGlass: zero vertex number (3)");
			return (NULL);
		}
		glass -> vertexnumbers[0] = idns[0];
		glass -> vertexnumbers[1] = idns[1];
		glass -> vertexnumbers[2] = idns[2];
		break;
	case DenSquare:
	case ProvTetrahedron:
		if (idns[0] == 0 || idns[1] == 0 || idns[2] == 0 || idns[3] == 0) {
			set_error1 ("storeGlass: zero vertex number (4)");
			return (NULL);
		}
		glass -> vertexnumbers[0] = idns[0];
		glass -> vertexnumbers[1] = idns[1];
		glass -> vertexnumbers[2] = idns[2];
		glass -> vertexnumbers[3] = idns[3];
		break;
	case DenVertex:			/* kludge */
	case BorderVertex:
		glass -> vertexnumbers[0] = idnumbers[3];
		break;
	default:
		return (NULL);
	}
	plex -> nglass[type]++;
	glass -> number = plex -> nglass[type];
	return (glass);
}

int moreGlass(struct Plex *plex)
{
	struct GlassBlock *gb;
	struct Glass *glasses, *first, *last, *g;

	if (plex -> glassblock == 0) return (0);
	glasses = (struct Glass *) allocate_objects (GLASS, plex -> glassblock);
	if (glasses == NULL) {
		set_error1 ("moreGlass: not enough memory for glasses");
		return (0);
	}
	gb = (struct GlassBlock *) allocate_object (GLASS_BLOCK);
	if (gb == NULL) return (0);
	gb -> glasses = glasses;
	if (plex -> tail_glassblock == NULL)
		plex -> head_glassblock = gb;
	else plex -> tail_glassblock -> next = gb;
	plex -> tail_glassblock = gb;
	first = glasses;
	last = first + plex -> glassblock - 1;
	for (g = first; g < last; g++)
		g -> next = g + 1;
	last -> next = NULL;
	plex -> glasses = glasses;
	plex -> freeglass = glasses;
	return (1);
}

struct Glass *allocateGlass(struct Plex *plex)
{
	struct Glass *glass;
	/* grab an unused slot */
	if (plex -> freeglass == NULL) {
		if (!moreGlass (plex)) {
			return (NULL);
		}
	}
	glass = plex -> freeglass;
	plex -> freeglass = glass -> next;
	return (glass);
}

/* given identifying integers, find hash entry */

struct Glass *lookupGlass (struct Plex *plex, enum PlexType type, long idnumbers[4])
{
	int j, result;
	long h, maxhash;
	long gdnumbers[4];
	struct Glass *glass;
	struct Plexi *plexi;
	
	maxhash = plex -> maxhash;
	h = hashGlass (idnumbers, maxhash);
	plexi = plex -> plexis + ((long) type) * maxhash + h;
	for (glass = plexi -> head; glass != NULL; glass = glass -> next) {
		for (j = 0; j < 4; j++)
			gdnumbers[j] = glass -> idnumbers[j];
		result = idmatch (type, gdnumbers, idnumbers);
		if (result != 0) return (glass);
	}
	return (NULL);
}

int add_edge (struct Plex *plex, long v0, long v1)
{
	int sgn, j;
	long idns[4], edns[4];
	struct Glass *eglass;
	idns[0] = v0; idns[1] = v1; idns[2] = 0; idns[3] = 0;
	eglass = lookupGlass (plex, PolyEdge, idns);
	if (eglass == NULL) return (0);
	for (j = 0; j < 4; j++)
		edns[j] = eglass -> idnumbers[j];
	sgn = idmatch (PolyEdge, idns, edns);
	if (sgn == 0) return (0);
	eglass -> coefficient += sgn;
	return (1);
}



/* create vertex in special array; also store in hashing table */

long makeVertex (struct Plex *plex, enum PlexType type, double center[3])
{
	int k;
	long idnumbers[4];
	long vnumber;
	struct PlexVertex *plexvtx;
	struct Glass *glass;
	
	if (plex -> n_vertex >= plex -> maxvertex) {
		set_error1 ("too many vertices");
		return (0L);
	}
	/* step 1 */
	plexvtx = plex -> plexvertices + plex -> n_vertex;
	plex -> n_vertex++;
	for (k = 0; k < 3; k++)
		plexvtx -> center[k] = center[k];
	plexvtx -> type = type;
	/* step 2 */
	for (k = 0; k < 4; k++) idnumbers[k] = 0;
	idnumbers[0] = plex -> n_vertex;
	glass = storeGlass (plex, type, idnumbers, center, 0.0);
	if (glass == NULL) return (0L);
	plexvtx -> glass = glass;
	vnumber = plex -> n_vertex;
	return (vnumber);
}


/* convert density cube indices to floating-point center coordinates */

int indices2center (struct Plex *plex, long indices[3], double center[3], int addhalf)
{
	int k;
	double cube_width;
	
	for (k = 0; k < 3; k++) {
		if (k >= plex -> dimension) continue;
		cube_width = (plex -> bounds[1][k] - plex -> bounds[0][k]) / plex -> dimensions[k];
		center[k] = plex -> bounds[0][k] + cube_width * (indices[k] + addhalf * 0.5);
	}
	if (plex -> dimension == 2)
		center[2] = 0.0;
	return(1);
}

/* convert floating-point center coordinates to density cube indices */

int center2indices (struct Plex *plex, double center[3], long indices[3])
{
	int k;
	double fraction;
	for (k = 0; k < 3; k++) {
		if (k >= plex -> dimension) continue;
		fraction = (center[k] - plex -> bounds[0][k]) / (plex -> bounds[1][k] - plex -> bounds[0][k]);
		indices[k] = floor (fraction * plex -> dimensions[k]);
		if (indices[k] < 0 || indices[k] >= plex -> dimensions[k]) return (0);
	}
	if (plex -> dimension == 2)
		indices[2] = 0;
	return (1);
}


/* CUBE DENSITY    * * *   CUBE DENSITY   * * *   CUBE DENSITY*/


/* classify each cube as interior (full), surface (partial) or exterior (empty) */

int classifyCubes (struct Plex *plex)
{
	if (!markSurfaceCubes(plex)) return (0);		/* mark cubes intersecting polyhedron */
	if (!connectCubes(plex)) return (0);			/* find connected components of non-surface cubes */
	if (!numberRoots(plex)) return (0);				/* number component representatives (roots) */
	if (!numberCubes(plex)) return (0);				/* assign component number to all non-surface cubes */
	if (!inoroutCubes(plex)) return (0);			/* distinguish inner and outer cubes */
	return(1);
}

/* flag each cube with a vertex, line or triangle intersecting it */

int markSurfaceCubes(struct Plex *plex)
{
	int result;
	long maxhash;
	struct Plexi *plexis;
	
	maxhash = plex -> maxhash;
	plexis = plex -> plexis;
	/* markOneType(plex, plexis + PolyVertex * maxhash); */
	if (plex -> dimension == 2)
		result = markOneType(plex, plexis + PolyEdge * maxhash);
	if (plex -> dimension == 3)
		result = markOneType(plex, plexis + PolyTriangle * maxhash);
	if (!result) return (0);
	return (1);
}

int markOneType (struct Plex *plex, struct Plexi *plexis)
{
	int result;
	long h, maxhash;
	struct Plexi *plexi;
	struct Glass *glass;
	
	/* for each element of the given type in the hash table, mark it surface */
	maxhash = plex -> maxhash;
	for (h = 0; h < maxhash; h++) {
		plexi = plexis + h;
		for (glass = plexi -> head; glass != NULL; glass = glass -> next) {
			result = markOne (plex, glass);
			if (!result) return (0);
		}
	}
	return(1);
}

/* rough check: center of element */
/* this method fails sometimes for edges and triangles; no longer matters */

int markOne (struct Plex *plex, struct Glass *glass)
{
	int result, k;
	enum CubeType type;
	long indices[3];
	double center[3];
	
	for (k = 0; k < 3; k++)
		center[k] = glass -> center[k];
	result = center2indices(plex, center, indices);
	if (!result) return(0);
	type = PartialCube;
	result = markCubeType(plex, indices, type);
	return(1);

}

/* mark cube with given value (empty, partial or full) */

int markCubeType (struct Plex *plex, long indices[3], enum CubeType type)
{
	long idx, xdim, ydim;
	struct PlexCube *cube;

	xdim = plex -> dimensions[0];
	ydim = plex -> dimensions[1];
	idx = indices[2] * (xdim * ydim) + indices[1] * xdim + indices[0];
	cube = plex -> plexcubes + idx;
	cube -> type = type;
	return (1);
}

/* identify connected components of non-surface density */

int connectCubes(struct Plex *plex)
{
	int result;
	long x, y, z, idx, nx, ny, nz;
	struct PlexCube *plexcubes;
	struct PlexCube *cube0, *cube1, *cube2, *cube3;
	
	plexcubes = plex -> plexcubes;
	nx = plex -> dimensions[0];
	ny = plex -> dimensions[1];
	nz = plex -> dimensions[2];
	for (z = 0; z < nz; z++)
		for (y = 0; y < ny; y++)
			for (x = 0; x < nx; x++) {
				idx = z * (nx * ny) + y * nx + x;
				if (idx >= plex -> n_cube)
					return (0);
				cube0 = plexcubes + idx;
				if (cube0 -> type == PartialCube) continue;
				if (x < nx-1) {
					cube1 = cube0 + 1;
					if (cube1 -> type == cube0 -> type) {
						result = joinCubes (cube0, cube1);
						if (!result) return (0);
					}
				}
				if (y < ny-1) {
					cube2 = cube0 + nx;
					if (cube2 -> type == cube0 -> type) {
						result = joinCubes (cube0, cube2);
						if (!result) return (0);
					}
				}
				if (plex -> dimension < 3) continue;
				if (z < nz-1) {
					cube3 = cube0 + (nx * ny);
					if (cube3 -> type == cube0 -> type) {
						result = joinCubes (cube0, cube3);
						if (!result) return (0);
					}
				}
			}
	return (1);
}

/* mark both cubes as belonging to the same connected component */

int joinCubes (struct PlexCube *cube1, struct PlexCube *cube2)
{
	struct PlexCube *rcube1, *rcube2;
	
	if (cube1 == NULL || cube2 == NULL) return (0);
	/* find root of each growing component tree */
	rcube1 = cubeRoot (cube1);
	rcube2 = cubeRoot (cube2);
	if (rcube1 == NULL || rcube2 == NULL) return (0);
	if (rcube1 == rcube2) return (1);	/* already joined */
	/* connect roots */
	rcube1 -> next = rcube2;
	return (1);
}

/* find root (component-representative) cube of given cube */

struct PlexCube *cubeRoot (struct PlexCube *start)
{
	struct PlexCube *cube;
	
	for (cube = start; cube != NULL; cube = cube -> next)
		if (cube -> next == NULL) return (cube);
	return (NULL);
}

/* assign component numbers to one representative cube of each component */

int numberRoots (struct Plex *plex)
{
	long idx, comp, n_cube;
	char message[MAXLINE];
	struct PlexCube *plexcubes, *cube;
	
	n_cube = plex -> n_cube;
	plexcubes = plex -> plexcubes;
	comp = 0;
	for (idx = 0; idx < n_cube; idx++) {
		cube = plexcubes + idx;
		if (cube -> type == PartialCube) continue;
		if (cube -> next != NULL) continue;
		/* we have a root cube */
		comp++;
		cube -> comp = comp;
	}
	plex -> n_component = comp;
	if (plex -> dimension == 3) {
		sprintf (message, "%8ld non-surface cube components", comp);
	}
	return (1);
}

/* assign component number of root cube to each other cube pointing at it */

int numberCubes(struct Plex *plex)
{
	long idx, comp, n_cube;
	struct PlexCube *plexcubes, *cube, *rcube;
	long ncomps[512];
	char message[MAXLINE];
	
	n_cube = plex -> n_cube;
	plexcubes = plex -> plexcubes;
	for (comp = 0; comp < 512; comp++)
		ncomps[comp] = 0;
	for (idx = 0; idx < n_cube; idx++) {
		cube = plexcubes + idx;
		if (cube -> type == PartialCube) continue;	/* surface cubes out of it */
		if (cube -> next == NULL) {
			comp = cube -> comp;
			ncomps[comp]++;
			continue;			/* root cubes already handled */
		}
		rcube = cubeRoot (cube);					/* find representative */
		if (rcube == NULL) continue;
		cube -> comp = rcube -> comp;
		comp = cube -> comp;
		ncomps[comp]++;
	}
	if (plex -> dimension == 3)
	for (comp = 0; comp < 512; comp++)
		if (ncomps[comp] > 0) {
			sprintf (message, "%8ld cubes in component %4ld",
				ncomps[comp], comp);
			inform (message);
		}
	return (1);
}

/* determine which components are inner, and which outer */

int inoroutCubes (struct Plex *plex)
{
	int inside, dim, onedge, result;
	long idx, n_cube, nin, nout, nsurf;
	long nx, ny, nz, x, y, z;
	long indices[3];
	double center[3];
	struct PlexCube *plexcubes, *cube, *rcube;
	char message[MAXLINE];
	
	n_cube = plex -> n_cube;
	dim = plex -> dimension;
	nin = 0; nout = 0; nsurf = 0;
	plexcubes = plex -> plexcubes;
	nx = plex -> dimensions[0];
	ny = plex -> dimensions[1];
	nz = plex -> dimensions[2];
	/* first check the root cubes, to save time */
	for (z = 0; z < nz; z++)
		for (y = 0; y < ny; y++)
			for (x = 0; x < nx; x++) { 
				idx = z * (nx * ny) + y * nx + x;
				rcube = plexcubes + idx;
				/* skip surface cubes */
				if (rcube -> type == PartialCube) continue;
				if (rcube -> next != NULL) continue;
				/* we have a root cube (one per connected component) */
				/* special case: cubes at edges of bounding box */
				onedge = (x == 0 || x == nx-1 || y == 0 || y == ny-1 ||
						(dim == 3 && (z == 0 || z == nz-1)));
				if (onedge) inside = 0;
				else {
					/* use winding number or solid angle as check */
					indices[0] = x; indices[1] = y; indices[2] = z;
					result = indices2center (plex, indices, center, 1);
					if (result == 0) return (0);
					inside = insidePoly (plex, center);
				}
				if (inside) rcube -> type = FullCube;
				else rcube -> type = EmptyCube;
			}
	for (idx = 0; idx < n_cube; idx++) {
		cube = plexcubes + idx;
		if (cube -> type == PartialCube) {
			nsurf++;
			continue;
		}
		/* copy info from root cube */
		rcube = cubeRoot (cube);
		if (rcube -> type == FullCube) {
			cube -> type = FullCube;
			nin++;
		}
		else if (rcube -> type == EmptyCube) {
			cube -> type = EmptyCube;
			nout++;
		}
	}
	sprintf (message, "%8ld cubes inside", nin);
	inform (message);
	sprintf (message, "%8ld cubes outside", nout);
	inform (message);
	sprintf (message, "%8ld cubes surface", nsurf);
	inform (message);
	return (1);
}



/* add one square to hash table, return pointer to its entry */

struct Glass *addSquare (struct Plex *plex, long x0, long y0, long z0, long x1, long y1, long z1, long x2, long y2, long z2, long x3, long y3, long z3)
{
	long vns[4];
	struct Glass *glass, *eglass;
	
	vns[0] = loraDenVertex(plex, x0, y0, z0);
	if (vns[0] == 0)
		return (NULL);
	vns[1] = loraDenVertex(plex, x1, y1, z1);
	if (vns[1] == 0)
		return (NULL);
	vns[2] = loraDenVertex(plex, x2, y2, z2);
	if (vns[2] == 0)
		return (NULL);
	vns[3] = loraDenVertex(plex, x3, y3, z3);
	if (vns[3] == 0)
		return (NULL);
	eglass = loraDenEdge (plex, vns[0], vns[1]);
	if (eglass == NULL)
		return (NULL);
	eglass = loraDenEdge (plex, vns[1], vns[2]);
	if (eglass == NULL)
		return (NULL);
	eglass = loraDenEdge (plex, vns[2], vns[3]);
	if (eglass == NULL)
		return (NULL);
	eglass = loraDenEdge (plex, vns[3], vns[0]);
	if (eglass == NULL)
		return (NULL);
	glass = loraDenSquare (plex, vns[0], vns[1], vns[2], vns[3]);
	return (glass);
}

/* lookup or add density (cube) vertex */

long loraDenVertex (struct Plex *plex, long x, long y, long z)
{
	long v;
	long idnumbers[4];
	struct Glass *glass;

	idnumbers[0] = x;
	idnumbers[1] = y;
	idnumbers[2] = z;
	idnumbers[3] = 0;
	glass = lookupGlass(plex, DenVertex, idnumbers);
	if (glass == NULL) {
		v = makeDenVertex (plex, x, y, z);
		if (v == 0)
			return (0);
		return (v);
	}
	/* glass points at valid entry */
	v = glass -> vertexnumbers[0];
	return(v);
}

/* lookup or add density (cube) edge */

struct Glass *loraDenEdge (struct Plex *plex, long v0, long v1)
{
	int k;
	long idnumbers[4];
	double measure;
	double center[3], center0[3], center1[3];
	struct PlexVertex *plexvtx0, *plexvtx1;
	struct Glass *glass;

	idnumbers[0] = v0;
	idnumbers[1] = v1;
	idnumbers[2] = 0;
	idnumbers[3] = 0;
	glass = lookupGlass (plex, DenEdge, idnumbers);
	if (glass == NULL) {
		plexvtx0 = plex -> plexvertices + (v0 - 1);
		plexvtx1 = plex -> plexvertices + (v1 - 1);
		for (k = 0; k < 3; k++)
			center[k] = (plexvtx0 -> center[k] + plexvtx1 -> center[k]) / 2;
		for (k = 0; k < 3; k++) {
			center0[k] = plexvtx0 -> center[k];
			center1[k] = plexvtx1 -> center[k];
		}
		measure = distance (center0, center1);
		glass = storeGlass (plex, DenEdge, idnumbers, center, measure);
		if (glass == NULL)
			return (NULL);
	}
	return (glass);
}

/* lookup or add density (cube) square */

struct Glass *loraDenSquare (struct Plex *plex, long v0, long v1, long v2, long v3)
{
	int k;
	long idnumbers[4];
	double measure;
	double center[3], center0[3], center1[3], center2[3];
	struct PlexVertex *plexvtx0, *plexvtx1, *plexvtx2, *plexvtx3;
	struct Glass *glass;

	if (v0 == 0 || v1 == 0 || v2 == 0 || v3 == 0) {
		set_error1 ("loraDenSquare: zero vertex arg");
		return (NULL);
	}
	idnumbers[0] = v0;
	idnumbers[1] = v1;
	idnumbers[2] = v2;
	idnumbers[3] = v3;
	glass = lookupGlass (plex, DenSquare, idnumbers);
	if (glass == NULL) {
		plexvtx0 = plex -> plexvertices + (v0 - 1);
		plexvtx1 = plex -> plexvertices + (v1 - 1);
		plexvtx2 = plex -> plexvertices + (v2 - 1);
		plexvtx3 = plex -> plexvertices + (v3 - 1);
		for (k = 0; k < 3; k++) center[k] = 0.0;
		for (k = 0; k < 3; k++) center[k] += plexvtx0 -> center[k];
		for (k = 0; k < 3; k++) center[k] += plexvtx1 -> center[k];
		for (k = 0; k < 3; k++) center[k] += plexvtx2 -> center[k];
		for (k = 0; k < 3; k++) center[k] += plexvtx3 -> center[k];
		for (k = 0; k < 3; k++) center[k] /= 4;
		for (k = 0; k < 3; k++) {
			center0[k] = plexvtx0 -> center[k];
			center1[k] = plexvtx1 -> center[k];
			center2[k] = plexvtx2 -> center[k];
		}
		measure = distance (center0, center1) * distance (center1, center2);
		glass = storeGlass (plex, DenSquare, idnumbers, center, measure);
		if (glass == NULL) return (NULL);
	}
	return (glass);
}

long makeDenVertex (struct Plex *plex, long x, long y, long z)
{
	int k;
	long idnumbers[4];
	long vnumber;
	enum PlexType type;
	double center[3];
	struct PlexVertex *plexvtx;
	struct Glass *glass;

	type = DenVertex;
	if (plex -> n_vertex >= plex -> maxvertex) {
		set_error1 ("makeDenVertex: overflow");
		return (0);
	}
	plexvtx = plex -> plexvertices + plex -> n_vertex;
	plex -> n_vertex++;
	idnumbers[0] = x;
	idnumbers[1] = y;
	idnumbers[2] = z;
	indices2center (plex, idnumbers, center, 0);
	for (k = 0; k < 3; k++)
		plexvtx -> center[k] = center[k];
	plexvtx -> type = type;
	idnumbers[3] = plex -> n_vertex;
	glass = storeGlass (plex, type, idnumbers, center, 0.0);
	if (glass == NULL) return (0L);
	plexvtx -> glass = glass;
	vnumber = plex -> n_vertex;
	return (vnumber);
}


/* POLYHEDRON   * * * *  POLYHEDRON  * * * *  POLYHEDRON  */

/* store polygon/polyhedron into Plex */

int storePoly (struct Plex *plex)
{
	char message[MAXLINE];

	if (!storePolyBefore(plex)) return (0);
	if (plex -> dimension == 2) {
		if (!cutPoly2(plex)) return (0);
	}
	else if (plex -> dimension == 3) {
		if (!cutPoly3(plex)) return (0);
	}
	if (!storePolyAfter(plex)) return (0);
	free_objects (PLEXEDGE, (short *) (plex -> plexedges));
	plex -> plexedges = NULL;
	if (plex -> plextriangles != NULL) {
		free_objects (PLEXTRIANGLE, (short *) (plex -> plextriangles));
		plex -> plextriangles = NULL;
	}
	if (plex -> dimension == 3) {
		sprintf (message, "%8ld edges after dicing", plex -> n_edge);
		inform (message);
		sprintf (message, "%8ld triangles after dicing", plex -> n_triangle);
		inform (message);
	}
	return (1);
}

int storePolyBefore (struct Plex *plex)
{
	int j, k;
	long t, e, v;
	long *edges;
	long *triangles;
	double *vertices;
	struct PlexEdge *plexedg;
	struct PlexTriangle *plextri;
	struct PlexVertex *plexvtx;
	
	vertices = plex -> vertices;
	edges = plex -> edges;
	triangles = plex -> triangles;
	/* store vertices (for both dim = 2 & 3) */
	if (plex -> mvertex > plex -> maxvertex) return (0);
	/* special vertex array with coordinate information */
	for (v = 0; v < plex -> mvertex; v++) {
		plexvtx = plex -> plexvertices + v;
		for (k = 0; k < 3; k++)
			plexvtx -> center[k] = *(vertices + 3 * v + k);
		plexvtx -> type = PolyVertex;
	};
	/* for polyhedra, we store edges and triangles only after subdivision */
	if (plex -> dimension == 3) return (1);
	/* store edges */
	if (plex -> medge > plex -> maxedge) return (0);
	for (e = 0; e < plex -> medge; e++) {
		plexedg = plex -> plexedges + e;
		for (j = 0; j < 2; j++)
			plexedg -> vns[j] = *(edges + 2 * e + j);
	}
	/* store triangles */
	if (plex -> n_triangle > plex -> maxtriangle) return (0);
	for (t = 0; t < plex -> n_triangle; t++) {
		plextri = plex -> plextriangles + t;
		for (j = 0; j < 3; j++)
			plextri -> ens[j] = *(triangles + 6 * t + j);
		for (j = 0; j < 3; j++)
			plextri -> vns[j] = *(triangles + 6 * t + 3 + j);
	}
	return (1);
}

int cutPoly2 (struct Plex *plex)
{
	int result, k;
	long nx, ny;
	long i, e, v0, v1;
	double p1[3], p2[3], p3[3], p4[3];
	struct PlexEdge *plexedg;
	struct PlexVertex *pv0, *pv1;
	
	nx = plex -> dimensions[0]; ny = plex -> dimensions[1];
	if (!defineGrid (plex)) return (0);
	/* go through grid lines perpendicular to x-axis */
	for (i = 0; i < nx; i++) {
		for (k = 0; k < 3; k++) {
			p1[k] = *(plex -> xgrid + 6 * i + k);
			p2[k] = *(plex -> xgrid + 6 * i + 3 + k);
		}
		for (e = 0; e < plex -> n_edge; e++) {
			plexedg = plex -> plexedges + e;
			v0 = plexedg -> vns[0];
			v1 = plexedg -> vns[1];
			pv0 = plex -> plexvertices + (v0-1);
			pv1 = plex -> plexvertices + (v1-1);
			for (k = 0; k < 3; k++) {
				p3[k] = pv0 -> center[k];
				p4[k] = pv1 -> center[k];
			}
			result = lines_intersect (p1, p2, p3, p4);
			if (!result) continue;
			result = do_fraction (plex, plexedg, p1, p2, p3, p4);
			if (!result) return (0);
		}
	}
	/* go through grid lines perpendicular to y-axis */
	for (i = 0; i < ny; i++) {
		for (k = 0; k < 3; k++) {
			p1[k] = *(plex -> ygrid + 6 * i + k);
			p2[k] = *(plex -> ygrid + 6 * i + 3 + k);
		}
		for (e = 0; e < plex -> n_edge; e++) {
			plexedg = plex -> plexedges + e;
			v0 = plexedg -> vns[0];
			v1 = plexedg -> vns[1];
			pv0 = plex -> plexvertices + (v0-1);
			pv1 = plex -> plexvertices + (v1-1);
			for (k = 0; k < 3; k++) {
				p3[k] = pv0 -> center[k];
				p4[k] = pv1 -> center[k];
			}
			result = lines_intersect (p1, p2, p3, p4);
			if (!result) continue;
			result = do_fraction (plex, plexedg, p1, p2, p3, p4);
			if (!result) return (0);
		}
	}
	return (1);
}

int do_fraction (struct Plex *plex, struct PlexEdge *plexedg, double *p1, double *p2, double *p3, double *p4)
{
	int k;
	double fraction;
	struct PlexVertex *newvtx;
	struct PlexEdge *newedg;
	char message[MAXLINE];

	if (plex -> n_vertex >= plex -> maxvertex) {
		set_error1 ("too many vertices");
		return (0);
	}
	if (plex -> n_edge >= plex -> maxedge) {
		sprintf (message, "do_fraction: too many edges (%6ld)",
			plex -> maxedge);
		set_error1 (message);
		return (0);
	}
	fraction =  line_intersection (p1, p2, p3, p4);
	/* new vertex */
	newvtx = plex -> plexvertices + plex -> n_vertex;
	for (k = 0; k < 3; k++)
		newvtx -> center[k] = (1.0 - fraction) * *(p1 + k) + fraction * *(p2 + k);
	newvtx -> type = PolyVertex;
	plex -> n_vertex++;
	/* new edge */
	newedg = plex -> plexedges + plex -> n_edge;
	newedg -> vns[0] = plex -> n_vertex;
	newedg -> vns[1] = plexedg -> vns[1];
	plexedg -> vns[1] = plex -> n_vertex;
	plex -> n_edge++;
	plex -> n_cut++;
	return (1);
}

int defineGrid (struct Plex *plex)
{
	int k;
	long x, y, z, nx, ny, nz;
	long indices[2][3];
	double center[3];

	nx = plex -> dimensions[0]; ny = plex -> dimensions[1]; nz = plex -> dimensions[2];

	plex -> xgrid = allocate_doubles (6 * nx, 0, XGRID);
	if (plex -> xgrid == NULL)
		return (0);
	for (k = 0; k < 3; k++) {
		indices[0][k] = 0;
		indices[1][k] = plex -> dimensions[k];
	}
	for (x = 0; x < nx; x++) {
		indices[0][0] = x; indices[1][0] = x;
		if (!indices2center (plex, indices[0], center, 0)) return (0);
		for (k = 0; k < 3; k++)
			*(plex -> xgrid + x * 6 + k) = center[k];
		if (!indices2center (plex, indices[1], center, 0)) return (0);
		for (k = 0; k < 3; k++)
			*(plex -> xgrid + x * 6 + 3 + k) = center[k];
	}

	plex -> ygrid = allocate_doubles (6 * ny, 0, YGRID);
	if (plex -> ygrid == NULL)
		return (0);
	for (k = 0; k < 3; k++) {
		indices[0][k] = 0;
		indices[1][k] = plex -> dimensions[k];
	}
	for (y = 0; y < ny; y++) {
		indices[0][1] = y; indices[1][1] = y;
		if (!indices2center (plex, indices[0], center, 0)) return (0);
		for (k = 0; k < 3; k++)
			*(plex -> ygrid + y * 6 + k) = center[k];
		if (!indices2center (plex, indices[1], center, 0)) return (0);
		for (k = 0; k < 3; k++)
			*(plex -> ygrid + y * 6 + 3 + k) = center[k];
	}

	if (plex -> dimension == 2) {
		plex -> zgrid = NULL;
		return (1);
	}
	plex -> zgrid = allocate_doubles (6 * nz, 0, ZGRID);
	if (plex -> zgrid == NULL)
		return (0);
	for (k = 0; k < 3; k++) {
		indices[0][k] = 0;
		indices[1][k] = plex -> dimensions[k];
	}
	for (z = 0; z < nz; z++) {
		indices[0][1] = z; indices[1][1] = z;
		if (!indices2center (plex, indices[0], center, 0)) return (0);
		for (k = 0; k < 3; k++)
			*(plex -> zgrid + z * 6 + k) = center[k];
		if (!indices2center (plex, indices[1], center, 0)) return (0);
		for (k = 0; k < 3; k++)
			*(plex -> zgrid + z * 6 + 3 + k) = center[k];
	}
	return (1);
}

int cutPoly3 (struct Plex *plex)
{
	if (!initEST (plex)) return (0);
	if (!initPST (plex)) return (0);
	if (!cutST (plex)) return (0);
	if (!finishST (plex)) return (0);
	if (!fetchST (plex)) return (0);
	if (!freeST (plex)) return (0);
	return (1);
}

int initEST (struct Plex *plex)
{
	int j;
	long nitem;
	long e, vn;
	long *le;
	struct subedge *se;
	char message[MAXLINE];
	
	nitem = plex -> medge;
	plex -> subedges = (struct subedge **) allocate_pointers (SUBEDGE, nitem);
	if (plex -> subedges == NULL) return (0);
	for (e = 0; e < plex -> medge; e++) {
		le = plex -> edges + 2 * e;
		se = allocate_subedge ();
		if (se == NULL) return (0);
		*(plex -> subedges + e) = se;
		for (j = 0; j < 2; j++) {
			vn = *(le + j);
			if (vn <= 0) {
				sprintf (message, "initEST: zero vertex number (%d)", j);
				set_error1 (message);
				sprintf (message, "edge number %6ld", e + 1);
				set_error2 (message);
				return (0);
			}
			se -> vns[j] = vn;
		}
	}
	return (1);
}

initPST (struct Plex *plex)
{
	int j, orn;
	long nitem;
	long e, t;
	long *lt;
	struct subpolygon *sp;
	char message[MAXLINE];
	
	nitem = plex -> mtriangle;
	plex -> subpolygons = (struct subpolygon **) allocate_pointers (SUBPOLYGON, nitem);
	if (plex -> subpolygons == NULL) return (0);
	for (t = 0; t < plex -> mtriangle; t++) {
		lt = plex -> triangles + 6 * t;
		sp = allocate_subpolygon ();
		if (sp == NULL) return (0);
		*(plex -> subpolygons + t) = sp;
		sp -> n_edge = 3;
		for (j = 0; j < 3; j++) {
			e = *(lt + j);
			if (e == 0) {
				sprintf (message, "initPST: zero edge number (%d)", j);
				set_error1 (message);
				sprintf (message, "triangle number %6ld", t + 1);
				set_error2 (message);
				return (0);
			}
			orn = (e > 0) ? 0 : 1;
			sp -> orn[j] = (short) orn;
			sp -> edges[j] = *(plex -> subedges + (abs (e) - 1));
		}
	}
	return (1);
}


/* cut edges and polygons */

int cutST (struct Plex *plex)
{
	int result, which, k;
	double center[3], axis[3];
	long i, n;
	long ico[3];
	
	for (which = 0; which < 3; which++) {
		n = plex -> dimensions[which];
		ico[0] = 0; ico[1] = 0; ico[2] = 0;
		for (k = 0; k < 3; k++)
			axis[k] = (k == which);
		for (i = 1; i < n; i++) {
			ico[which] = i;
			indices2center (plex, ico, center, 0);
			result = onePlane (plex, center, axis);
			if (!result) {
				set_error1 ("onePlane fails");
				return (0);
			}
		}
	}
	return (1);
}

int onePlane (struct Plex *plex, double center[3], double axis[3])
{
	long e, t;
	struct subedge *se;
	struct subpolygon *sp;
	
	for (e = 0; e < plex -> medge; e++) {
		se = *(plex -> subedges + e);
		if (!cutEST (plex, se, center, axis)) return (0);
	}
	for (t = 0; t < plex -> mtriangle; t++) {
		sp = *(plex -> subpolygons + t);
		if (!cutPE (plex, sp, center, axis)) return (0);
	}
	for (t = 0; t < plex -> mtriangle; t++) {
		sp = *(plex -> subpolygons + t);
		if (!cutPST (plex, sp, center, axis)) return (0);
	}
	return (1);
}

int cutPE (struct Plex *plex, struct subpolygon *sp, double center[3], double axis[3])
{
	struct subpolygon *sp0, *sp1;
	struct subedge *se;
	
	if (sp == NULL) return (1);
	se = sp -> se;
	sp0 = sp -> child[0];
	sp1 = sp -> child[1];
	if (se != NULL) {
		if (!cutEST(plex, se,  center, axis)) return (0);
	}
	if (sp0 != NULL && sp1 != NULL) {
		if (!cutPE (plex, sp0, center, axis)) return (0);
		if (!cutPE (plex, sp1, center, axis)) return (0);
	}
	return (1);
}

int cutEST (struct Plex *plex, struct subedge *se, double center[3], double axis[3])
{
	struct subedge *se0, *se1;
	
	if (se == NULL) return (1);
	se0 = se -> child[0];
	se1 = se -> child[1];
	if (se0 != NULL && se1 != NULL) {
		if (!cutEST (plex, se0, center, axis)) return (0);
		if (!cutEST (plex, se1, center, axis)) return (0);
		return (1);
	}
	if (!cutsubedge (plex, se, center, axis)) return (0);
	return (1);
}

int cutPST (struct Plex *plex, struct subpolygon *sp, double center[3], double axis[3])
{
	struct subpolygon *sp0, *sp1;
	
	if (sp == NULL) return (1);
	sp0 = sp -> child[0];
	sp1 = sp -> child[1];
	if (sp0 != NULL && sp1 != NULL) {
		if (!cutPST (plex, sp0, center, axis)) {
			set_error1 ("cutPST returns null (1)");
			return (0);
		}
		if (!cutPST (plex, sp1, center, axis)) {
			set_error1 ("cutPST returns null (2)");
			return (0);
		}
		return (1);
	}
	if (!cutsubpolygon (sp)) {
		set_error1 ("cutPST returns null (3)");
		return (0);
	}
	return (1);
}

int cutsubedge (struct Plex *plex, struct subedge *se, double center[3], double axis[3])
{
	int k;
	double vect0[3], vect1[3], xpnt[3];
	double dot0, dot1, f;
	struct PlexVertex *pvtx0, *pvtx1;
	
	if (se == NULL) return (1);
	pvtx0 = plex -> plexvertices + (se -> vns[0] - 1);
	pvtx1 = plex -> plexvertices + (se -> vns[1] - 1);
	for (k = 0; k < 3; k++) {
		if (pvtx0 -> center[k] == center[k]) return (1);
		if (pvtx1 -> center[k] == center[k]) return (1);
	}
	for (k = 0; k < 3; k++) {
		vect0[k] = pvtx0 -> center[k] - center[k];
		vect1[k] = pvtx1 -> center[k] - center[k];
	}
	dot0 = dot_product (vect0, axis);
	dot1 = dot_product (vect1, axis);
	if (dot0 * dot1 >= 0.0) return (1);
	/* vertices lie on opposite sides */
	f = dot0 / (dot0 - dot1);
	for (k = 0; k < 3; k++) {
		xpnt[k] = (1.0 - f) * pvtx0 -> center[k] + f * pvtx1 -> center[k];
		if (axis[k] != 0.0) xpnt[k] = center[k];
	}
	if (!bisectEdge (plex, se, xpnt)) return (0);
	return (1);
}

int bisectEdge (struct Plex *plex, struct subedge *se, double xpnt[3])
{
	int k;
	struct subedge *se0, *se1;
	struct PlexVertex *pvtx;
	
	if (se == NULL) return (1);
	se0 = allocate_subedge ();
	if (se0 == NULL) {
		set_error1 ("bisectEdge: not enough memory");
		return (0);
	}
	se1 = allocate_subedge ();
	if (se1 == NULL) {
		set_error1 ("bisectEdge: not enough memory");
		return (0);
	}
	if (plex -> n_vertex >= plex -> maxvertex) {
		set_error1 ("too many vertices");
		return (0);
	}
	pvtx = plex -> plexvertices + plex -> n_vertex;
	plex -> n_vertex++;
	for (k = 0; k < 3; k++)
		pvtx -> center[k] = xpnt[k];
	pvtx -> type = PolyVertex;
	se0 -> vns[0] = se -> vns[0];
	se0 -> vns[1] = plex -> n_vertex;
	se1 -> vns[0] = se0 -> vns[1];
	se1 -> vns[1] = se -> vns[1];
	se -> child[0] = se0;
	se -> child[1] = se1;
	se0 -> parent = se;
	se1 -> parent = se;
	return (1);
}

/* if non-leaf edges, bisect */

int cutsubpolygon (struct subpolygon *sp)
{
	int n, j, i0, i1, result, orn;
	short orns[MAXPED];
	char message[MAXLINE];
	struct subedge *se, *se0, *se1;
	struct subedge *edges[MAXPED];
	
	if (sp == NULL) return (1);
	n = 0; i0 = -1; i1 = -1;
	for (j = 0; j < sp -> n_edge; j++) {
		se = sp -> edges[j];
		orn = sp -> orn[j];
		se0 = se -> child[0];
		se1 = se -> child[1];
		if (se0 == NULL || se1 == NULL) {
			/* non-split: just copy */
			if (n >= MAXPED) {
				set_error1 ("cutsubpolygon: MAXPED exceeded");
				return (0);
			}
			edges[n] = sp -> edges[j];
			orns[n] = sp -> orn[j];
			n++;
		}
		else {
			/* split edge */
			if (n >= MAXPED) {
				set_error1 ("cutsubpolygon: MAXPED exceeded");
				return (0);
			}
			if (i0 < 0) i0 = n+1;
			else if (i1 < 0) i1 = n+1;
			else {
				sprintf (message, "cutsubpolygon: > 2 cuts, positions: %3d %3d %3d of %3d",
					i0, i1, n, sp -> n_edge);
				set_error1 (message);
				return (0);
			}
			edges[n] = orn ? se1 : se0;
			orns[n] = (short) orn;
			n++;
			if (n >= MAXPED) {
				set_error1 ("cutsubpolygon: MAXPED exceeded");
				return (0);
			}
			edges[n] = orn ? se0 : se1;
			orns[n] = (short) orn;
			n++;
		}
	}
	if (i0 < 0 || i1 < 0) return (1);
	if (i0 >= 0 && i1 < 0) {
		set_error1 ("cutsubpolygon: only one cut");
		return (0);
	}
	for (j = 0; j < n; j++) {
		sp -> edges[j] = edges[j];
		sp -> orn[j] = orns[j];
	}
	sp -> n_edge = (short) n;
	result = bisectPolygon (sp, i0, i1);
	if (!result) {
		set_error1 ("cutsubpolygon: bisectPolygon fails");
		return (0);
	}
	return (1);
}

int bisectPolygon (struct subpolygon *sp, int i0, int i1)
{
	int m, j, orn0, orn1;
	long v0, v1, ne, ne0, ne1;
	struct subpolygon *sp0, *sp1;
	struct subedge *se;
	char message[MAXLINE];
	
	if (i0 >= i1) {
		sprintf (message, "bisectPolygon: i0 (%2d) >= i1 (%2d)", i0, i1);
		set_error1 (message);
		return (0);
	}
	sp0 = allocate_subpolygon ();
	if (sp0 == NULL) {
		return (0);
	}
	sp1 = allocate_subpolygon ();
	if (sp1 == NULL) {
		return (0);
	}
	/* create new edge */
	se = allocate_subedge ();
	if (se == NULL) {
		return (0);
	}
	sp -> se = se;
	orn0 = sp -> orn[i0];
	orn1 = sp -> orn[i1];
	v0 = sp -> edges[i0] -> vns[orn0];
	v1 = sp -> edges[i1] -> vns[orn1];
	se -> vns[0] = v0;
	se -> vns[1] = v1;
	ne = sp -> n_edge;
	ne0 = i0 + (ne - i1);
	ne1 = i1 - i0;
	/* first polygon */
	sp0 -> n_edge = ne0+1;
	for (j = 0; j < ne0; j++) {
		m = i1 + j;
		if (m >= ne) m -= ne;
		sp0 -> orn[j] = sp -> orn[m];
		sp0 -> edges[j] = sp -> edges[m];
	}
	/* bisecting edge (positive orientation) */
	sp0 -> orn[ne0] = 0;
	sp0 -> edges[ne0] = se;
	/* second polygon */
	sp1 -> n_edge = ne1+1;
	for (j = 0; j < ne1; j++) {
		m = i0 + j;
		sp1 -> orn[j] = sp -> orn[m];
		sp1 -> edges[j] = sp -> edges[m];
	}
	/* bisecting edge (reverse orientation) */
	sp1 -> orn[ne1] = 1;
	sp1 -> edges[ne1] = se;
	sp -> child[0] = sp0;
	sp -> child[1] = sp1;
	sp0 -> parent = sp;
	sp1 -> parent = sp;
	return (1);
}

/* finish subdividing polygons with > 3 vertices */

long finishST (struct Plex *plex)
{
	long p, n, np;
	struct subpolygon *sp;
	
	np = 0;
	for (p = 0; p < plex -> mtriangle; p++) {
		sp = *(plex -> subpolygons + p);
		n = finishPST (plex, sp);
		if (n == 0) {
			set_error1 ("finishST: finishPST return zero value");
			return (0);
		}
		np += n;
	}
	return (np);
}

long finishPST (struct Plex *plex, struct subpolygon *sp)
{
	long n0, n1, n;
	struct subpolygon *sp0, *sp1;

	n = 0;
	sp0 = sp -> child[0];
	sp1 = sp -> child[1];
	if (sp0 != NULL && sp1 != NULL) {
		n0 = finishPST (plex, sp0);
		if (n0 == 0) return (0L);
		n1 = finishPST (plex, sp1);
		if (n1 == 0) return (0L);
		n += n0 + n1;
	}
	else {
		n0 = finishPolygon (plex, sp);
		if (n0 == 0) return (0L);
		n += n0;
	}
	return (n);
}

long finishPolygon (struct Plex *plex, struct subpolygon *sp)
{
	int result;
	long n, n0, n1;
	struct subpolygon *sp0, *sp1;
	
	if (sp == NULL) return (0L);
	if (sp -> n_edge == 3) return (1L);	/* done */
	sp0 = sp -> child[0];
	sp1 = sp -> child[1];
	if (sp0 == NULL && sp1 == NULL) {
		result = bisectPolygon (sp, 0, 2);
		if (!result) {
			set_error1 ("finishPolygon: bisectPolygon fails");
			return (0L);
		}
		sp0 = sp -> child[0];
		sp1 = sp -> child[1];
	}
	n0 = finishPolygon (plex, sp0);
	n1 = finishPolygon (plex, sp1);
	n = n0 + n1;
	return (n);
}

int fetchST (struct Plex *plex)
{
	long ne, np;

	ne = fetchSubEdges (plex);
	if (ne == 0) return (0);
	np = fetchSubPolygons (plex);
	if (np == 0) return (0);
	return (1);
}

long fetchSubEdges (struct Plex *plex)
{
	long e, ne, n, p;
	struct subedge *se;
	struct subpolygon *sp;
	
	plex -> n_edge = 0;
	ne = 0;
	for (e = 0; e < plex -> medge; e++) {
		se = *(plex -> subedges + e);
		n = fetchETree (plex, se);
		if (n <= 0) return (0);
		ne +=n;
	}
	for (p = 0; p < plex -> mtriangle; p++) {
		sp = *(plex -> subpolygons + p);
		n = fetchPETree (plex, sp);
		if (n < 0) return (0);
		ne += n;
	}
	return (ne);
}

long fetchETree (struct Plex *plex, struct subedge *se)
{
	int j;
	long n0, n1;
	struct subedge *se0, *se1;
	struct PlexEdge *pe;
	char message[MAXLINE];
	
	if (se == NULL) return (0);
	if (se -> child[0] == NULL) {
		if (plex -> n_edge >= plex -> maxedge) {
			sprintf (message, "fetchETree: too many edges (%6ld)",
				plex -> maxedge);
			set_error1 (message);
			return (0);
		}
		pe = plex -> plexedges + plex -> n_edge;
		plex -> n_edge++;
		se -> number = plex -> n_edge;
		for (j = 0; j < 2; j++) {
			if (se -> vns[j] == 0) {
				sprintf (message, "edge %6ld with null vertex", se -> number);
				set_error1 (message);
				return (0);
			}
			pe -> vns[j] = se -> vns[j];
		}
		return (1L);
	}
	else {
		se0 = se -> child[0];
		se1 = se -> child[1];
		n0 = fetchETree (plex, se0);
		n1 = fetchETree (plex, se1);
		return (n0 + n1);
	}
}

long fetchPETree (struct Plex *plex, struct subpolygon *sp)
{
	long ne, n0, n1, n;
	struct subedge *se;
	struct subpolygon *sp0, *sp1;

	n = 0;
	sp0 = sp -> child[0];
	sp1 = sp -> child[1];
	se = sp -> se;
	if (se != NULL) {
		ne = fetchETree (plex, se);
		if (ne <= 0) return (0L);
		n += ne;
	}
	if (sp0 != NULL && sp1 != NULL) {
		n0 = fetchPETree (plex, sp0);
		n1 = fetchPETree (plex, sp1);
		n += n0 + n1;
	}
	return (n);
}

long fetchSubPolygons (struct Plex *plex)
{
	long p, np, n;
	struct subpolygon *sp;
	
	plex -> n_triangle = 0;
	np = 0;
	for (p = 0; p < plex -> mtriangle; p++) {
		sp = *(plex -> subpolygons + p);
		n = fetchPTree (plex, sp);
		if (n <= 0) return (0);
		np += n;
	}
	return (np);
}

long fetchPTree (struct Plex *plex, struct subpolygon *sp)
{
	int j, orn;
	long n0, n1;
	struct PlexTriangle *pt;
	struct subpolygon *sp0, *sp1;
	struct subedge *se;
	char message[MAXLINE];
	
	if (sp == NULL) return (0L);
	sp0 = sp -> child[0];
	sp1 = sp -> child[1];
	if (sp0 == NULL || sp1 == NULL) {
		if (sp -> n_edge != 3) return (0L);
		pt = plex -> plextriangles + plex -> n_triangle;
		if (plex -> n_triangle > plex -> maxtriangle) {
			set_error1 ("fetchPTree: too many triangles");
			return (0);
		}
		plex -> n_triangle++;
		for (j = 0; j < 3; j++) {
			if (sp -> edges[j] == NULL) {
				sprintf (message, "triangle %6ld with null edge", plex -> n_triangle);
				set_error1 (message);
				return (0);
			}
			se = sp -> edges[j];
			orn = sp -> orn[j];
			if (se -> number == 0) {
				sprintf (message, "triangle %8ld with 0 edge number", plex -> n_triangle);
				set_error1 (message);
				sprintf (message, "edges numbers: %8ld %8ld %8ld",
					sp -> edges[0] -> number, sp -> edges[1] -> number,
					sp -> edges[2] -> number);
				set_error2 (message);
				return (0);
			}
			pt -> ens[j] = (1-2*orn) * se -> number;
			pt -> vns[j] = se -> vns[orn];
		}
		return (1L);
	}
	else {
		n0 = fetchPTree (plex, sp0);
		if (n0 == 0) return (0L);
		n1 = fetchPTree (plex, sp1);
		if (n1 == 0) return (0L);
		return (n0 + n1);
	}
}

/* free subdivision memory */

int freeST (struct Plex *plex)
{
	long ne, np;
	char message[MAXLINE];

	ne = freeSubEdges (plex);
	if (ne == 0) return (0);
	sprintf (message, "%8ld subedges freed", ne);
	inform (message);
	free_pointers (SUBEDGE, (short *) (plex -> subedges));
	plex -> subedges = NULL;
	np = freeSubPolygons (plex);
	if (np == 0) return (0);
	sprintf (message, "%8ld subpolygons freed", np);
	inform (message);
	free_pointers (SUBPOLYGON, (short *) (plex -> subpolygons));
	plex -> subpolygons = NULL;
	free_cache (SUBEDGE);
	free_cache (SUBPOLYGON);
	return (1);
}

long freeSubEdges (struct Plex *plex)
{
	long e, ne, n, p;
	struct subedge *se;
	struct subpolygon *sp;
	
	ne = 0;
	for (e = 0; e < plex -> medge; e++) {
		se = *(plex -> subedges + e);
		n = freeETree (plex, se);
		if (n <= 0) return (0);
		ne +=n;
	}
	for (p = 0; p < plex -> mtriangle; p++) {
		sp = *(plex -> subpolygons + p);
		n = freePETree (plex, sp);
		if (n < 0) return (0);
		ne += n;
	}
	return (ne);
}

long freeETree (struct Plex *plex, struct subedge *se)
{
	long n0, n1;
	struct subedge *se0, *se1;
	
	if (se == NULL) return (0);
	if (se -> child[0] == NULL) {
		free_subedge (se);
		return (1L);
	}
	else {
		se0 = se -> child[0];
		se1 = se -> child[1];
		n0 = freeETree (plex, se0);
		n1 = freeETree (plex, se1);
		free_subedge (se);
		return (n0 + n1 + 1);
	}
}

long freePETree (struct Plex *plex, struct subpolygon *sp)
{
	long ne, n0, n1, n;
	struct subedge *se;
	struct subpolygon *sp0, *sp1;

	n = 0;
	sp0 = sp -> child[0];
	sp1 = sp -> child[1];
	se = sp -> se;
	if (se != NULL) {
		ne = freeETree (plex, se);
		if (ne <= 0) return (0L);
		n += ne;
	}
	if (sp0 != NULL && sp1 != NULL) {
		n0 = freePETree (plex, sp0);
		n1 = freePETree (plex, sp1);
		n += n0 + n1;
	}
	return (n);
}

long freeSubPolygons (struct Plex *plex)
{
	long p, np, n;
	struct subpolygon *sp;
	
	np = 0;
	for (p = 0; p < plex -> mtriangle; p++) {
		sp = *(plex -> subpolygons + p);
		n = freePTree (plex, sp);
		if (n <= 0) return (0);
		np += n;
	}
	return (np);
}

long freePTree (struct Plex *plex, struct subpolygon *sp)
{
	long n0, n1;
	struct subpolygon *sp0, *sp1;
	
	if (sp == NULL) return (0L);
	sp0 = sp -> child[0];
	sp1 = sp -> child[1];
	if (sp0 == NULL || sp1 == NULL) {
		if (sp -> n_edge != 3) return (0L);
		free_subpolygon (sp);
		return (1L);
	}
	else {
		n0 = freePTree (plex, sp0);
		if (n0 == 0) return (0L);
		n1 = freePTree (plex, sp1);
		if (n1 == 0) return (0L);
		free_subpolygon (sp);
		return (n0 + n1 + 1);
	}
}

int storePolyAfter (struct Plex *plex)
{
	int k;
	long t, e, v, v0, v1, v2, e0, e1, e2;
	long idnumbers[4];
	double measure;
	double center[3], center0[3], center1[3];
	struct PlexVertex *plexvtx, *plexvtx0, *plexvtx1, *plexvtx2;
	struct PlexEdge *plexedg;
	struct PlexTriangle *plextri;
	struct Glass *glass;
	
	/* vertex part of hashing table */
	for (v = 0; v < plex -> n_vertex; v++) {
		plexvtx = plex -> plexvertices + v;
		idnumbers[0] = v + 1;
		for (k = 1; k < 4; k++) idnumbers[k] = 0;
		for (k = 0; k < 3; k++)
			center[k] = plexvtx -> center[k];
		glass = storeGlass (plex, plexvtx -> type, idnumbers, center, 0.0);
		if (glass == NULL) {
			set_error1 ("storePolyAfter: storeGlass returns null");
			return (0);
		}
	}
	/* edge part of hashing table */
	for (e = 0; e < plex -> n_edge; e++) {
		plexedg = plex -> plexedges + e;
		v0 = plexedg -> vns[0];
		v1 = plexedg -> vns[1];
		plexvtx0 = plex -> plexvertices + (v0 - 1);
		plexvtx1 = plex -> plexvertices + (v1 - 1);
		for (k = 0; k < 3; k++)
			center[k] = (plexvtx0 -> center[k] + plexvtx1 -> center[k]) / 2;
		for (k = 0; k < 3; k++) {
			center0[k] = plexvtx0 -> center[k];
			center1[k] = plexvtx1 -> center[k];
		}
		measure = distance (center0, center1);
		for (k = 2; k < 4; k++) idnumbers[k] = 0;
		idnumbers[0] = v0; idnumbers[1] = v1;
		glass = storeGlass (plex, PolyEdge, idnumbers, center, measure);
		if (glass == NULL) {
			set_error1 ("storePolyAfter: storeGlass returns null");
			return (0);
		}
	}
	if (plex -> dimension == 2) return (1);
	/* triangle part of hashing table */
	for (t = 0; t < plex -> n_triangle; t++) {
		plextri = plex -> plextriangles + t;
		v0 = plextri -> vns[0];
		v1 = plextri -> vns[1];
		v2 = plextri -> vns[2];
		e0 = plextri -> ens[0];
		e1 = plextri -> ens[1];
		e2 = plextri -> ens[2];
		if (v0 < 1 || v1 < 1 || v2 < 1) {
			set_error1 ("storePolygonAfter: non-positive vertex numbers");
			return (0);
		}
		if (e0 == 0 || e1 == 0 || e2 == 0) {
			set_error1 ("storePolygonAfter: zero edge numbers");
			return (0);
		}
		e0 = abs(e0);
		e1 = abs(e1);
		e2 = abs(e2);
		plexvtx0 = plex -> plexvertices + (v0 - 1);
		plexvtx1 = plex -> plexvertices + (v1 - 1);
		plexvtx2 = plex -> plexvertices + (v2 - 1);
		idnumbers[0] = v0;
		idnumbers[1] = v1;
		idnumbers[2] = v2;
		idnumbers[3] = 0;
		/* compute center */
		for (k = 0; k < 3; k++)
			center[k] = 0.0;
		for (k = 0; k < 3; k++)
			center[k] += plexvtx0 -> center[k];
		for (k = 0; k < 3; k++)
			center[k] += plexvtx1 -> center[k];
		for (k = 0; k < 3; k++)
			center[k] += plexvtx2 -> center[k];
		for (k = 0; k < 3; k++)
			center[k] /= 3;
		/* later: use Heron's formula for area */
		glass = storeGlass (plex, PolyTriangle, idnumbers, center, 0.0);
		if (glass == NULL) {
			set_error1 ("storePolyAfter: storeGlass returns null");
			return (0);
		}
	}
	return (1);
}

/* low level memory allocation routines: */

struct subedge *allocate_subedge ()
{
	struct subedge *sed;

	/* allocate memory */
	sed = (struct subedge *) allocate_object (SUBEDGE);
	if (error ()) return (NULL);
	if (sed == NULL) {
		set_error1 ("allocate_subedge: ran out of memory");
		return(NULL);
	}
	return (sed);
}

void free_subedge (struct subedge *sed)
{
	free_object (SUBEDGE, (short *) sed);
	if (error()) return;
}


struct subpolygon *allocate_subpolygon ()
{
	struct subpolygon *spg;

	/* allocate memory */
	spg = (struct subpolygon *) allocate_object (SUBPOLYGON);
	if (error ()) return (NULL);
	if (spg == NULL) {
		set_error1 ("allocate_subpolygon: ran out of memory");
		return(NULL);
	}
	return (spg);
}

void free_subpolygon (struct subpolygon *spg)
{
	free_object (SUBPOLYGON, (short *) spg);
	if (error()) return;
}



/* Molecular Surface Package */
/* Copyright 1995 by Michael L. Connolly */

