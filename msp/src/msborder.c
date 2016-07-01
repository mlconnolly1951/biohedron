/* MSP: Polyhedron-Density */
/* Copyright 1995 by Michael L. Connolly */
/* June 5, 1998 */

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"


/* BORDERS    *    *     BORDERS     *      *     BORDERS */

long makeBorderVertex (struct Plex *plex, long x, long y, long z, long pN0, long pN1, long pN2, long pN3, long pN4, long pN5, long pN6, long pN7)
{
	int k;
	long idnumbers[4];
	enum PlexType type;
	double center[3];
	struct PlexVertex *plexvtx;
	struct Glass *glass;

	type = BorderVertex;
	if (plex -> n_vertex >= plex -> maxvertex) return (0);
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
	glass -> pN[0] = pN0;
	glass -> pN[1] = pN1;
	glass -> pN[2] = pN2;
	glass -> pN[3] = pN3;
	glass -> pN[4] = pN4;
	glass -> pN[5] = pN5;
	glass -> pN[6] = pN6;
	glass -> pN[7] = pN7;
	glass -> npn = (plex -> dimension == 2) ? 4 : 8;
	plexvtx -> glass = glass;
	return (plex -> n_vertex);
}

/* lookup or add border edge */

struct Glass *makeBorderEdge (struct Plex *plex, long v0, long v1, long pN0, long pN1, long pN2, long pN3)
{
	int k;
	long idnumbers[4];
	double measure;
	double center[3], center0[3], center1[3];
	struct PlexVertex *plexvtx0, *plexvtx1;
	struct Glass *eglass;

	idnumbers[0] = v0;
	idnumbers[1] = v1;
	idnumbers[2] = 0;
	idnumbers[3] = 0;
	plexvtx0 = plex -> plexvertices + (v0 - 1);
	plexvtx1 = plex -> plexvertices + (v1 - 1);
	for (k = 0; k < 3; k++)
		center[k] = (plexvtx0 -> center[k] + plexvtx1 -> center[k]) / 2;
	for (k = 0; k < 3; k++) {
		center0[k] = plexvtx0 -> center[k];
		center1[k] = plexvtx1 -> center[k];
	}
	measure = distance (center0, center1);
	eglass = storeGlass (plex, BorderEdge, idnumbers, center, measure);
	if (eglass == NULL) return (NULL);
	eglass -> pN[0] = pN0;
	eglass -> pN[1] = pN1;
	eglass -> pN[2] = pN2;
	eglass -> pN[3] = pN3;
	eglass -> npn = (plex -> dimension == 2) ? 2 : 4;
	return (eglass);
}

struct Glass *makeBorderTriangle (struct Plex *plex, long v0, long v1, long v2, long pN0, long pN1, enum PlexType subtype)
{
	int k;
	long idnumbers[4];
	double measure;
	double center[3], center0[3], center1[3], center2[3];
	struct PlexVertex *plexvtx0, *plexvtx1, *plexvtx2;
	struct Glass *tglass;

	idnumbers[0] = v0;
	idnumbers[1] = v1;
	idnumbers[2] = v2;
	idnumbers[3] = 0;
	plexvtx0 = plex -> plexvertices + (v0 - 1);
	plexvtx1 = plex -> plexvertices + (v1 - 1);
	plexvtx2 = plex -> plexvertices + (v2 - 1);
	for (k = 0; k < 3; k++)
		center[k] = (plexvtx0 -> center[k] + plexvtx1 -> center[k] + plexvtx2 -> center[k]) / 3;
	for (k = 0; k < 3; k++) {
		center0[k] = plexvtx0 -> center[k];
		center1[k] = plexvtx1 -> center[k];
		center2[k] = plexvtx2 -> center[k];
	}
	measure = distance (center0, center1) * distance (center1, center2);
	tglass = storeGlass (plex, BorderTriangle, idnumbers, center, measure);
	if (tglass == NULL) return (NULL);
	tglass -> pN[0] = pN0;
	tglass -> pN[1] = pN1;
	tglass -> pN[2] = 0;
	tglass -> pN[3] = 0;
	tglass -> npn = 2;
	tglass -> subtype = subtype;
	return (tglass);
}

int createBorderVertices (struct Plex *plex)
{
	int dim;
	long x, y, z, vN, nx, ny, nz;
	long nborder, idx;
	long pN0, pN1, pN2, pN3, pN4, pN5, pN6, pN7;
	long indices[3];
	double center[3];
	struct PlexCube *plexcubes, *cube;
	struct PlexCube *cube0, *cube1, *cube2, *cube3, *cube4, *cube5, *cube6, *cube7;
	char message[MAXLINE];

	nx = plex -> dimensions[0]; ny = plex -> dimensions[1]; nz = plex -> dimensions[2];
	plexcubes = plex -> plexcubes;
	nborder = 0; dim = plex -> dimension;
	for (z = (dim == 3) ? 1 : 0; z < nz; z++) {
		for (y = 1; y < ny; y++) {
			for (x = 1; x < nx; x++) {
				idx = y * nx + x;
				if (dim == 3) idx += z * nx * ny;
				cube = plexcubes + idx;
				cube0 = cube;
				cube1 = cube-1;
				cube2 = cube-nx;
				cube3 = cube-nx-1;
				pN0 = cube0 -> provinceNumber;
				pN1 = cube1 -> provinceNumber;
				pN2 = cube2 -> provinceNumber;
				pN3 = cube3 -> provinceNumber;
				if (pN0 == 0) continue;
				if (pN1 == 0) continue;
				if (pN2 == 0) continue;
				if (pN3 == 0) continue;
				if (dim == 3) {
					cube4 = cube0-nx*ny;
					cube5 = cube1-nx*ny;
					cube6 = cube2-nx*ny;
					cube7 = cube3-nx*ny;
					pN4 = cube4 -> provinceNumber;
					pN5 = cube5 -> provinceNumber;
					pN6 = cube6 -> provinceNumber;
					pN7 = cube7 -> provinceNumber;
					if (pN4 == 0) continue;
					if (pN5 == 0) continue;
					if (pN6 == 0) continue;
					if (pN7 == 0) continue;
				}
				else {
					pN4 = 0;
					pN5 = 0;
					pN6 = 0;
					pN7 = 0;
				}
				indices[0] = x; indices[1] = y; indices[2] = z;
				indices2center (plex, indices, center, 0);
				if (cube -> type == EmptyCube) continue;
				if (cube -> type == PartialCube) {
					if (!insidePoly (plex, center)) continue;
				}
				vN = makeBorderVertex (plex, x, y, z, pN0, pN1, pN2, pN3, pN4, pN5, pN6, pN7);
				if (vN == 0) return (0);
				nborder++;
			}
		}
	}
	sprintf (message, "%8ld border vertices created", nborder);
	inform (message);
	return (1);
}

int createBorderEdges (struct Plex *plex)
{
	long x, y, z, nx, ny, nz, n,ncreated;
	long pN000, pN001, pN010, pN011, pN100, pN101, pN110;
	long vns[MAXVRT];
	int k, result; double vector[3];
	struct PlexCube *plexcubes;
	struct PlexCube *cube000, *cube001, *cube010, *cube011;
	struct PlexCube *cube100, *cube101, *cube110;
	char message[MAXLINE];
	
	plexcubes = plex -> plexcubes; 
	ncreated = 0;
	nx = plex -> dimensions[0]; ny = plex -> dimensions[1]; nz = plex -> dimensions[2];
	if (plex -> dimension == 2) {
		for (x = 0; x < nx; x++)
			for (y = 0; y < ny-1; y++) {
				for (k = 0; k < 3; k++)
					vector[k] = -(k == 0);
				cube000 = plexcubes + y * nx + x;
				pN000 = cube000 -> provinceNumber;
				cube010 = plexcubes + (y+1) * nx + x;
				pN010 = cube010 -> provinceNumber;
				if (pN000 != 0 && pN010 != 0) {
					n = collectBorderVertices (plex, vns, pN000, pN010, 0, 0);
					if (n < 0) return (0);
					if (n > MAXVRT) {
						set_error1 ("createBorderEdges: MAXVRT too small");
						return (0);
					}
					if (n > 0) {
						result = sortBorderVertices (plex, vector, n, vns);
						if (!result) return (0);
						ncreated += processBorderVertices (plex, n, vns, pN000, pN010, 0, 0);
					}
				}
			}
		for (y = 0; y < ny; y++)
			for (x = 0; x < nx-1; x++) {
				for (k = 0; k < 3; k++)
					vector[k] = -(k == 1);
				cube001 = plexcubes + y * nx + (x+1);
				cube000 = plexcubes + y * nx + x;
				pN001 = cube001 -> provinceNumber;
				pN000 = cube000 -> provinceNumber;
				if (pN001 != 0 && pN000 != 0) {
					n = collectBorderVertices (plex, vns, pN001, pN000, 0, 0);
					if (n < 0) return (0);
					if (n > MAXVRT) {
						set_error1 ("createBorderEdges: MAXVRT too small");
						return (0);
					}
					if (n > 0) {
						result = sortBorderVertices (plex, vector, n, vns);
						if (!result) return (0);
						ncreated += processBorderVertices (plex, n, vns, pN001, pN000, 0, 0);
					}
				}
			}
	}
	else if (plex -> dimension == 3) {
		for (x = 0; x < nx; x++)
			for (y = 0; y < ny-1; y++)
				for (z = 0; z < nz-1; z++) {
					for (k = 0; k < 3; k++)
						vector[k] = -(k == 0);
					cube000 = plexcubes + z * nx * ny + y * nx + x;
					pN000 = cube000 -> provinceNumber;
					cube010 = plexcubes + z * nx * ny + (y+1) * nx + x;
					pN010 = cube010 -> provinceNumber;
					cube100 = plexcubes + (z+1) * nx * ny + y * nx + x;
					pN100 = cube100 -> provinceNumber;
					cube110 = plexcubes + (z+1) * nx * ny + (y+1) * nx + x;
					pN110 = cube110 -> provinceNumber;
					if (pN000 != 0 && pN010 != 0 && pN100 != 0 && pN110 != 0) {
						n = collectBorderVertices (plex, vns, pN000, pN010, pN110, pN100);
						if (n < 0) return (0);
						if (n > MAXVRT) {
							set_error1 ("createBorderEdges: MAXVRT too small");
							return (0);
						}
						if (n > 0) {
							result = sortBorderVertices (plex, vector, n, vns);
							if (!result) return (0);
							ncreated += processBorderVertices (plex, n, vns, pN000, pN010, pN110, pN100);
						}
					}
				}
		for (y = 0; y < ny; y++)
			for (x = 0; x < nx-1; x++)
				for (z = 0; z < nz-1; z++) {
					for (k = 0; k < 3; k++)
						vector[k] = -(k == 1);
					cube001 = plexcubes + z * nx * ny + y * nx + (x+1);
					pN001 = cube001 -> provinceNumber;
					cube000 = plexcubes + z * nx * ny + y * nx + x;
					pN000 = cube000 -> provinceNumber;
					cube101 = plexcubes + (z+1) * nx * ny + y * nx + (x+1);
					pN101 = cube101 -> provinceNumber;
					cube100 = plexcubes + (z+1) * nx * ny + y * nx + x;
					pN100 = cube100 -> provinceNumber;
					if (pN001 != 0 && pN000 != 0 && pN101 != 0 && pN100 != 0) {
						n = collectBorderVertices (plex, vns, pN000, pN100, pN101, pN001);
						if (n < 0) return (0);
						if (n > MAXVRT) {
							set_error1 ("createBorderEdges: MAXVRT too small");
							return (0);
						}
						if (n > 0) {
							result = sortBorderVertices (plex, vector, n, vns);
							if (!result) return (0);
							ncreated += processBorderVertices (plex, n, vns, pN000, pN100, pN101, pN001);
						}
					}
				}
		for (z = 0; z < nz; z++)
			for (x = 0; x < nx-1; x++)
				for (y = 0; y < ny-1; y++) {
					for (k = 0; k < 3; k++)
						vector[k] = -(k == 2);
					cube001 = plexcubes + z * nx * ny + y * nx + (x+1);
					pN001 = cube001 -> provinceNumber;
					cube000 = plexcubes + z * nx * ny + y * nx + x;
					pN000 = cube000 -> provinceNumber;
					cube011 = plexcubes + z * nx * ny + (y+1) * nx + (x+1);
					pN011 = cube011 -> provinceNumber;
					cube010 = plexcubes + z * nx * ny + (y+1) * nx + x;
					pN010 = cube010 -> provinceNumber;
					if (pN001 != 0 && pN000 != 0 && pN011 != 0 && pN010 != 0) {
						n = collectBorderVertices (plex, vns, pN000, pN001, pN011, pN010);
						if (n < 0) return (0);
						if (n > MAXVRT) {
							set_error1 ("createBorderEdges: MAXVRT too small");
							return (0);
						}
						if (n > 0) {
							result = sortBorderVertices (plex, vector, n, vns);
							if (!result) return (0);
							ncreated += processBorderVertices (plex, n, vns, pN000, pN001, pN011, pN010);
						}
					}
				}
	}
	sprintf (message, "%8ld border edges created", ncreated);
	inform (message);
	return (1);
}

long collectBorderVertices (struct Plex *plex, long vns[], long pN0, long pN1, long pN2, long pN3)
{
	long n, maxhash, t, h, j, nok, dim;
	enum PlexType vertextypes[] = {PolyVertex, DenVertex, BorderVertex};
	enum PlexType type;
	struct Plexi *plexi;
	struct Glass *glass;
	
	n = 0; maxhash = plex -> maxhash; dim = plex -> dimension;
	for (t = 0; t < 3; t++) {
		type = vertextypes[t];
		plexi = plex -> plexis + (long) type * maxhash;
		for (h = 0; h < maxhash; h++) {
			for (glass = (plexi + h) -> head; glass != NULL; glass = glass -> next) {
				if ((dim == 2 && glass -> npn < 2) || (dim == 3 && glass -> npn < 4)) continue;
				nok = 0;
				for (j = 0; j < 8; j++) {
					if      (pN0 == glass -> pN[j]) nok++;
					else if (pN1 == glass -> pN[j]) nok++;
					else if (pN2 > 0 && pN2 == glass -> pN[j]) nok++;
					else if (pN3 > 0 && pN3 == glass -> pN[j]) nok++;
				}
				if ((dim == 2 && nok == 2) || (dim == 3 && nok == 4)) {
					if (n >= MAXVRT) {
						set_error1 ("collectBorderVertices: MAXVRT too small");
						return (-1);
					}
					vns[n++] = glass -> vertexnumbers[0];
				}
			}
		}
	}
	return (n);
}

int sortBorderVertices (struct Plex *plex, double vector[3], long n, long vns[])
{
	int k;
	double center[3];
	double lowval, dns[MAXVRT]; int used[MAXVRT];
	long wns[MAXVRT], i, idx, low;
	struct PlexVertex *plexvtx, *plexvertices;
	
	if (n > MAXVRT) {
		set_error1 ("sortBorderVertices: MAXVRT too small");
		return (0);
	}
	plexvertices = plex -> plexvertices;
	for (i = 0; i < n; i++) {
		wns[i] = vns[i];
		plexvtx = plexvertices + (vns[i] - 1);
		for (k = 0; k < 3; k++)
			center[k] = plexvtx -> center[k];
		dns[i] = dot_product (vector, center);
		used[i] = 0;
	}
	idx = 0;
	for (;;) {
		low = -1;
		lowval = 1000000.0;
		for (i = 0; i < n; i++) {
			if (used[i]) continue;
			if (dns[i] < lowval) {
				lowval = dns[i];
				low = i;
			}
		}
		if (low < 0) break;
		used[low] = 1;
		vns[idx++] = wns[low];
	}
	if (idx != n) return (0);
	return (1);
}

int processBorderVertices (struct Plex *plex, long n, long vns[], long pN0, long pN1, long pN2, long pN3)
{
	int i, ncreated;
	long v0, v1;
	struct Glass *eglass;
	
	ncreated = 0;
	for (i = 0; i < n - 1; i++) {
		v0 = vns[i]; v1 = vns[i+1];
		eglass = makeBorderEdge (plex, v0, v1, pN0, pN1, pN2, pN3);
		if (eglass == NULL) return (0);
		ncreated++;
	}
	return (ncreated);
}

int createBorderTriangles (struct Plex *plex)
{
	int k, result;
	long x, y, z, nx, ny, nz, n, ncreated;
	long pN000, pN001, pN010, pN100;
	long indices[3];
	double center[3], vector[3], center000[3], center001[3], center010[3], center100[3];
	struct PlexCube *plexcubes;
	struct PlexCube *cube000, *cube001, *cube010;
	struct PlexCube *cube100;
	int orns[MAXVRT];
	struct Glass *edges[MAXVRT];
	enum PlexType subtypes[MAXVRT];
	char message[MAXLINE];
	
	plexcubes = plex -> plexcubes; 
	ncreated = 0;
	nx = plex -> dimensions[0]; ny = plex -> dimensions[1]; nz = plex -> dimensions[2];
	for (x = 0; x < nx-1; x++)
		for (y = 0; y < ny-1; y++)
			for (z = 0; z < nz-1; z++) {
				cube000 = plexcubes + z * nx * ny + y * nx + x;
				pN000 = cube000 -> provinceNumber;
				indices[0] = x; indices[1] = y; indices[2] = z;
				indices2center (plex, indices, center000, 0);
				cube001 = plexcubes + z * nx * ny + y * nx + (x+1);
				pN001 = cube001 -> provinceNumber;
				indices[0] = x+1; indices[1] = y; indices[2] = z;
				indices2center (plex, indices, center001, 0);
				cube010 = plexcubes + z * nx * ny + (y+1) * nx + x;
				pN010 = cube010 -> provinceNumber;
				indices[0] = x; indices[1] = y+1; indices[2] = z;
				indices2center (plex, indices, center010, 0);
				cube100 = plexcubes + (z+1) * nx * ny + y * nx + x;
				pN100 = cube100 -> provinceNumber;
				indices[0] = x; indices[1] = y; indices[2] = z+1;
				indices2center (plex, indices, center100, 0);
				if (pN000 != 0 && pN001 != 0) {
					center[0] =  center001[0];
					center[1] = (center000[1] + center010[1])/2;
					center[2] = (center000[2] + center100[2])/2;
					for (k = 0; k < 3; k++)
						vector[k] = -(k == 0);
					n = collectBorderEdges (plex, edges, orns, subtypes, pN000, pN001, center, vector);
					if (n < 0) return (0);
					if (n > MAXVRT) {
						set_error1 ("createBorderTriangles: MAXVRT too small");
						return (0);
					}
					if (n > 0) {
						result = groupBorderEdges (plex, n, edges, orns, subtypes);
						if (!result) return (0);
						ncreated += processBorderEdges (plex, n, edges, orns, subtypes, pN000, pN001);
					}
				}
				if (pN000 != 0 && pN010 != 0) {
					center[0] = (center000[0] + center001[0])/2;
					center[1] =  center010[1];
					center[2] = (center000[2] + center100[2])/2;
					for (k = 0; k < 3; k++)
						vector[k] = -(k == 1);
					n = collectBorderEdges (plex, edges, orns, subtypes, pN000, pN010, center, vector);
					if (n < 0) return (0);
					if (n > MAXVRT) {
						set_error1 ("createBorderTriangles: MAXVRT too small");
						return (0);
					}
					if (n > 0) {
						result = groupBorderEdges (plex, n, edges, orns, subtypes);
						if (!result) return (0);
						ncreated += processBorderEdges (plex, n, edges, orns, subtypes, pN000, pN010);
					}
				}
				if (pN000 != 0 && pN100 != 0) {
					center[0] = (center000[0] + center001[0])/2;
					center[1] = (center000[1] + center010[1])/2;
					center[2] =  center100[2];
					for (k = 0; k < 3; k++)
						vector[k] = -(k == 2);
					n = collectBorderEdges (plex, edges, orns, subtypes, pN000, pN100, center, vector);
					if (n < 0) return (0);
					if (n > 0) {
						result = groupBorderEdges (plex, n, edges, orns, subtypes);
						if (!result) return (0);
						ncreated += processBorderEdges (plex, n, edges, orns, subtypes, pN000, pN100);
					}
				}
			}
	sprintf (message, "%8ld border triangles created", ncreated);
	inform (message);
	return (1);
}

long collectBorderEdges (struct Plex *plex, struct Glass *edges[], int orns[], enum PlexType subtypes[], long pN0, long pN1, double center[3], double vector[3])
{
	int pos0, pos1, j, k, t;
	long n, maxhash, h, nok, v0, v1;
	enum PlexType edgetypes[] = {PolyEdge, DenEdge, BorderEdge};
	enum PlexType type;
	double axis[3], radial[3], midpnt[3], sgn;
	struct Plexi *plexi;
	struct Glass *glass;
	struct PlexVertex *plexvtx0, *plexvtx1, *plexvertices;
	
	plexvertices = plex -> plexvertices;
	n = 0; maxhash = plex -> maxhash;
	for (t = 0; t < 3; t++) {
		type = edgetypes[t];
		plexi = plex -> plexis + (long) type * maxhash;
		for (h = 0; h < maxhash; h++) {
			for (glass = (plexi + h) -> head; glass != NULL; glass = glass -> next) {
				if (glass -> npn < 2) continue;
				nok = 0; pos0 = 0; pos1 = 0;
				for (j = 0; j < 4; j++) {
					if      (pN0 == glass -> pN[j]) { nok++; pos0 = j;}
					else if (pN1 == glass -> pN[j]) { nok++; pos1 = j;}
				}
				if (nok != 2) continue;
				if (n >= MAXVRT) {
					set_error1 ("collectBorderEdges: MAXVRT too small");
					return (-1);
				}
				edges[n] = glass;
				subtypes[n] = glass -> type;
				v0 = glass -> vertexnumbers[0];
				v1 = glass -> vertexnumbers[1];
				plexvtx0 = plexvertices + (v0 - 1);
				plexvtx1 = plexvertices + (v1 - 1);
				for (k = 0; k < 3; k++) {
					axis[k] = plexvtx1 -> center[k] - plexvtx0 -> center[k];
					midpnt[k] = (plexvtx0 -> center[k] + plexvtx1 -> center[k])/2;
					radial[k] = midpnt[k] - center[k];
				}
				sgn = triple_product (axis, radial, vector);
				if (glass -> type == PolyEdge) orns[n] = (pos0 > pos1);
				if (glass -> type == BorderEdge) orns[n] = (sgn > 0);
				if (glass -> type == DenEdge) orns[n] = (sgn > 0);
				n++;
			}
		}
	}
	return (n);
}

/* do this later, to make zone picture prettier; should have no effect on density */
int groupBorderEdges (struct Plex *plex, long n, struct Glass *edges[], int orns[], enum PlexType subtypes[])
{
	if (plex == NULL) return (0);
	if (n == 0L) return (0);
	if (edges == NULL) return (0);
	if (orns == NULL) return (0);
	if (subtypes == NULL) return (0);
	return (1);
}

int processBorderEdges (struct Plex *plex, long n, struct Glass *edges[], int orns[], enum PlexType subtypes[], long pN0, long pN1)
{
	int orn, ncreated;
	long e, v0, v1, v2;
	long vns[MAXVRT][2];
	struct Glass *tglass;
	
	if (n > MAXVRT) {
		set_error1 ("processBorderEdges: MAXVRT too small");
		return (0);
	}
	ncreated = 0;
	for (e = 0; e < n; e++) {
		orn = orns[e];
		vns[e][0] = edges[e] -> vertexnumbers[orn];
		vns[e][1] = edges[e] -> vertexnumbers[1-orn];
	}
	v0 = vns[0][0];
	for (e = 1; e < n; e++) {
		v1 = vns[e][0]; v2 = vns[e][1];
		tglass = makeBorderTriangle (plex, v0, v1, v2, pN1, pN0, subtypes[e]);
		if (tglass == NULL) return (0);
		ncreated++;
	}
	return (ncreated);
}

