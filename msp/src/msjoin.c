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

/* PROVINCES   * *  PROVINCES   * *   PROVINCES  */

/* create all provinces, one for each surface cube */

int createProvinces (struct Plex *plex)
{
	long n, idx, n_cube, nx, ny, nz, x, y, z;
	long indices[3];
	struct Province *provinces;
	struct PlexCube *cube;
	char message[MAXLINE];
	
	n_cube = plex -> n_cube;
	n = 0;
	for (idx = 0; idx < n_cube; idx++) {
		cube = plex -> plexcubes + idx;
		if (cube -> type == PartialCube) n++;
	}
	provinces = (struct Province *) allocate_objects (PROVINCE, n);
	if (provinces == NULL)
		return (0);
	plex -> provinces = provinces;
	plex -> n_surface = n;
	nx = plex -> dimensions[0];
	ny = plex -> dimensions[1];
	nz = plex -> dimensions[2];
	for (z = 0; z < nz; z++) 
		for (y = 0; y < ny; y++)
			for (x = 0; x < nx; x++) {
				idx = z * (nx * ny) + y * nx + x;
				cube = plex -> plexcubes + idx;
				if (cube -> type != PartialCube) continue;
				indices[0] = x; indices[1] = y; indices[2] = z;
				cube -> provinceNumber = makeProvince (plex, indices);
				if (cube -> provinceNumber == 0) return (0);
			}
	sprintf (message, "%8ld provinces created", n);
	inform (message);
	return (1);
}


/* recompute each province center as centroid of vertices belonging to province */

int provinceCenter (struct Plex *plex)
{
	int j, k;
	long h, maxhash, iprov, pN;
	struct Glass *glass;
	struct Plexi *plexi;
	struct Province *prov;
	struct PlexVertex *plexvtx;
	char message[MAXLINE];

	maxhash = plex -> maxhash;
	for (iprov = 0; iprov < plex -> n_province; iprov++) {
		prov = plex -> provinces + iprov;
		/* initialize */
		for (k = 0; k < 3; k++)
			prov -> center[k] = 0.0;
		prov -> n_vertex = 0;
		prov -> measure = 0.0;
	}
	for (h = 0; h < maxhash; h++) {
		plexi = plex -> plexis + ((int) PolyEdge) * maxhash + h;
		for (glass = plexi -> head; glass != NULL; glass = glass -> next) {
			for (j = 0; j < glass -> npn; j++) {
				pN = glass -> pN[j];
				if (pN == 0) continue;
				prov = plex -> provinces + (pN-1);
				prov -> n_vertex++;
				prov -> measure += glass -> measure;
				for (k = 0; k < 3; k++)
					prov -> center[k] += glass -> measure * glass -> center[k];
			}
		}
	}
	for (h = 0; h < maxhash; h++) {
		plexi = plex -> plexis + ((int) BorderEdge) * maxhash + h;
		for (glass = plexi -> head; glass != NULL; glass = glass -> next) {
			for (j = 0; j < glass -> npn; j++) {
				pN = glass -> pN[j];
				if (pN == 0) continue;
				prov = plex -> provinces + (pN-1);
				prov -> n_vertex++;
				prov -> measure += glass -> measure;
				for (k = 0; k < 3; k++)
					prov -> center[k] += glass -> measure * glass -> center[k];
			}
		}
	}
	for (h = 0; h < maxhash; h++) {
		plexi = plex -> plexis + ((int) DenEdge) * maxhash + h;
		for (glass = plexi -> head; glass != NULL; glass = glass -> next) {
			for (j = 0; j < glass -> npn; j++) {
				pN = glass -> pN[j];
				if (pN == 0) continue;
				prov = plex -> provinces + (pN-1);
				prov -> n_vertex++;
				prov -> measure += glass -> measure;
				for (k = 0; k < 3; k++)
					prov -> center[k] += glass -> measure * glass -> center[k];
			}
		}
	}
	for (iprov = 0; iprov < plex -> n_province; iprov++) {
		prov = plex -> provinces + iprov;
		if (prov -> measure == 0.0) continue;
		for (k = 0; k < 3; k++)
			prov -> center[k] /= prov -> measure;
		/* change vertex also */
		plexvtx = plex -> plexvertices + (prov -> vertexNumber - 1);
		for (k = 0; k < 3; k++)
			plexvtx -> center[k] = prov -> center[k];
	}
	sprintf (message, "%8ld province centers recomputed", plex -> n_province);
	inform (message);
	return (1);
}

/*   POLYHEDRON  PROVINCE     *   *   *   *   *     POLYHEDRON  PROVINCE    */

/* make one province (part of surface belonging to a cube) */

long makeProvince (struct Plex *plex, long indices[3])
{
	int k;
	long pN, idns[4];
	double center[3];
	struct Province *prov;

	indices2center (plex, indices, center, 1);
	prov = plex -> provinces + plex -> n_province;
	for (k = 0; k < 3; k++) {
		prov -> cubeIndices[k] = indices[k];
		prov -> center[k] = center[k];
	}
	prov -> vertexNumber = makeVertex(plex, ProvVertex, center);
	if (prov -> vertexNumber == 0) return (0L);
	for (k = 0; k < 4; k++) idns[k] = 0;
	idns[0] = prov -> vertexNumber;
	plex -> n_province++;
	pN = plex -> n_province;
	return (pN);
}

/* assign polyhedron elements to provinces */

int polyhedronProvince (struct Plex *plex)
{
	/* square province numbers already set */
	if (plex -> dimension == 3) {
		if (!triangleProvince (plex)) return (0);
	}
	/* now handle edges */
	if (!edgeProvince (plex)) return (0);
	/* now handle vertices */
	if (!vertexProvince (plex)) return (0);
	/* free memory for PlexLink objects */
	free_links (plex, PolyVertex);
	return (1);
}

/* add polyhedron triangles to province of surface cube they lie in */

int triangleProvince (struct Plex *plex)
{
	int k, result;
	long h, maxhash, nx, ny, idx, indices[3], pN, n;
	double center[3];
	struct Plexi *plexi;
	struct Glass *glass;
	struct PlexCube *cube;
	char message[MAXLINE];
	
	nx = plex -> dimensions[0];
	ny = plex -> dimensions[1];
	maxhash = plex -> maxhash;
	n = 0;
	for (h = 0; h < maxhash; h++) {
		plexi = plex -> plexis + ((int) PolyTriangle) * maxhash + h;
		for (glass = plexi -> head; glass != NULL; glass = glass -> next) {
			for (k = 0; k < 3; k++)
				center[k] = glass -> center[k];
			result = center2indices (plex, center, indices);
			if (!result) return(0);
			idx = indices[2] * nx * ny + indices[1] * nx + indices[0];
			cube = plex -> plexcubes + idx;
			pN = cube -> provinceNumber;
			if (pN == 0) {
				set_error1 ("triangleProvince: 0 province number for triangle");
				return (0);
			}
			glass -> pN[0] = pN;
			glass -> npn = 1;
			n++;
		}
	}
	sprintf (message, "%8ld polyhedron triangles added to provinces", n);
	inform (message);
	return (1);
}

/* add polyhedron edges to provinces they lie in or border */
int edgeProvince (struct Plex *plex)
{
	int orn, e, j, k, result;
	long h, maxhash, pN, n, nx, idx, indices[3];
	long idnumbers[4], vns[4], ednumbers[4];
	double center[3];
	struct Plexi *plexi;
	struct Glass *glass, *eglass;
	enum PlexType type;
	struct PlexCube *cube;
	char message[MAXLINE];
	
	nx = plex -> dimensions[0];
	maxhash = plex -> maxhash;
	n = 0;
	if (plex -> dimension == 2) {
		for (h = 0; h < maxhash; h++) {
			plexi = plex -> plexis + ((int) PolyEdge) * maxhash + h;
			for (glass = plexi -> head; glass != NULL; glass = glass -> next) {
				for (k = 0; k < 3; k++)
					center[k] = glass -> center[k];
				result = center2indices (plex, center, indices);
				if (!result) return(0);
				idx = indices[1] * nx + indices[0];
				cube = plex -> plexcubes + idx;
				pN = cube -> provinceNumber;
				glass -> pN[0] = pN;
				glass -> npn = 1;
			}
		}
		return (1);
	}
	for (h = 0; h < maxhash; h++) {
		plexi = plex -> plexis + ((int) PolyTriangle) * maxhash + h;
		for (glass = plexi -> head; glass != NULL; glass = glass -> next) {
			pN = glass -> pN[0];
			type = PolyEdge;
			/* gather polyhedron triangle vertex numbers */
			vns[0] = glass -> vertexnumbers[0];
			vns[1] = glass -> vertexnumbers[1];
			vns[2] = glass -> vertexnumbers[2];
			vns[3] = 0;
			idnumbers[2] = 0; idnumbers[3] = 0;
			for (e = 0; e < 3; e++) {
				for (j = 0; j < 2; j++) {
					k = j + e; if (k >= 3) k -= 3;
					idnumbers[j] = vns[k];
				}
				eglass = lookupGlass (plex, type, idnumbers);
				if (eglass == NULL) {
					sprintf (message, "edgeProvince: lookupGlass fails");
					set_error1 (message);
					sprintf (message, "idnumbers: %ld %ld", idnumbers[0], idnumbers[1]);
					set_error2 (message);
					return (0);
				}
				for (j = 0; j < 4; j++)
					ednumbers[j] = eglass -> idnumbers[j];
				orn = idmatch (type, idnumbers, ednumbers);
				if (orn == 1) eglass -> pN[0] = pN;
				else if (orn == -1) eglass -> pN[1] = pN;
				n++;
			}
		}
	}
	for (h = 0; h < maxhash; h++) {
		plexi = plex -> plexis + ((int) PolyEdge) * maxhash + h;
		for (eglass = plexi -> head; eglass != NULL; eglass = eglass -> next) {
			eglass -> npn = (eglass -> pN[0] == eglass -> pN[1]) ? 1 : 2;
		}
	}
	sprintf (message, "%8ld polyhedron edges added to provinces", n);
	inform (message);
	return (1);
}

/* add vertices to provinces */

int vertexProvince (struct Plex *plex)
{
	if (plex -> dimension == 2) {
		if (!vertexProvince0 (plex, PolyEdge)) return (0);
		return (1);
	}

	if (!vertexProvince1 (plex, PolyTriangle)) return (0);
	if (!vertexProvince2 (plex, PolyVertex)) return (0);
	
	return (1);
}

/* add vertices to provinces */

int vertexProvince0 (struct Plex *plex, enum PlexType type)
{
	long maxhash, h, pN, vns[2], idns[4];
	struct Plexi *plexi;
	struct Glass *eglass, *vglass, *vglass0, *vglass1;
	
	maxhash = plex -> maxhash;
	idns[1] = 0; idns[2] = 0; idns[3] = 0;
	for (h = 0; h < maxhash; h++) {
		plexi = plex -> plexis + ((int) type) * maxhash + h;
		for (eglass = plexi -> head; eglass != NULL; eglass = eglass -> next) {
			pN = eglass -> pN[0];
			vns[0] = eglass -> vertexnumbers[0];
			vns[1] = eglass -> vertexnumbers[1];
			idns[0] = vns[0];
			vglass0 = lookupGlass (plex, PolyVertex, idns);
			if (vglass0 == NULL)
				return (0);
			idns[0] = vns[1];
			vglass1 = lookupGlass (plex, PolyVertex, idns);
			if (vglass1 == NULL)
				return (0);
			/* province number of edge assigned to second array slot of first vertex*/
			vglass0 -> pN[1] = pN; vglass0 -> npn++;
			/* province number of edge assigned to first array slot of second vertex*/
			vglass1 -> pN[0] = pN; vglass1 -> npn++;
		}
	}
	if (type == PolyEdge) {	/* should always be true, since param is now superfl. */
		idns[1] = 0; idns[2] = 0; idns[3] = 0;
		for (h = 0; h < maxhash; h++) {
			plexi = plex -> plexis + ((int) PolyVertex) * maxhash + h;
			for (vglass = plexi -> head; vglass != NULL; vglass = vglass -> next) {
				if (vglass -> npn == 2 && vglass -> pN[0] == vglass -> pN[1]) {
					/* collapse duplication */
					vglass -> pN[1] = 0;
					vglass -> npn = 1;
				}
			}
		}
		
	}
	return (1);
}

/* add vertices to provinces; also add (homology theory) link edges to vertex */

int vertexProvince1 (struct Plex *plex, enum PlexType type)
{
	int result;
	long maxhash, h, pN, n, vns[4];
	struct Plexi *plexi;
	struct Glass *glass;
	char message[MAXLINE];
	
	n = 0;
	maxhash = plex -> maxhash;
	for (h = 0; h < maxhash; h++) {
		plexi = plex -> plexis + ((int) type) * maxhash + h;
		for (glass = plexi -> head; glass != NULL; glass = glass -> next) {
			pN = glass -> pN[0];
			vns[0] = glass -> vertexnumbers[0];
			vns[1] = glass -> vertexnumbers[1];
			vns[2] = glass -> vertexnumbers[2];
			vns[3] = glass -> vertexnumbers[3];
			result = addLink (plex, pN, vns[0], vns[1], vns[2]);
			if (result == 0) return (0);
			result = addLink (plex, pN, vns[1], vns[2], vns[0]);
			if (result == 0) return (0);
			result = addLink (plex, pN, vns[2], vns[0], vns[1]);
			if (result == 0) return (0);
			n++;
		}
	}
	sprintf (message, "%8ld polyhedron vertices added to provinces", n);
	inform (message);
	return (1);
}

int vertexProvince2 (struct Plex *plex, enum PlexType type)
{
	long maxhash, h, vn;
	struct PlexVertex *plexvtx;
	struct Plexi *plexi;
	struct Glass *glass;
	
	maxhash = plex -> maxhash;
	for (h = 0, plexi = plex -> plexis + (int) type * maxhash;
		h < maxhash; h++, plexi++) {
		for (glass = plexi -> head; glass != NULL; glass = glass -> next) {
			vn = glass -> vertexnumbers[0];
			plexvtx = plex -> plexvertices + (vn - 1);
			if (!sortLink (plexvtx, glass)) {
				set_error1 ("vertexProvince2: sortLink failed");
				return (0);
			}
		}
	}
	return (1);
}

int free_links (struct Plex *plex, enum PlexType type)
{
	long maxhash, h, vn;
	struct PlexVertex *plexvtx;
	struct Plexi *plexi;
	struct Glass *glass;
	struct PlexLink *next, *plnk;
	
	maxhash = plex -> maxhash;
	for (h = 0, plexi = plex -> plexis + (int) type * maxhash;
		h < maxhash; h++, plexi++) {
		for (glass = plexi -> head; glass != NULL; glass = glass -> next) {
			vn = glass -> vertexnumbers[0];
			plexvtx = plex -> plexvertices + (vn - 1);
			/* free PlexLink memory for this Plex vertex */
			for (plnk = plexvtx -> head; plnk != NULL; plnk = next) {
				next = plnk -> next;
				free_PlexLink (plnk);
			}
			plexvtx -> head = (struct PlexLink *) NULL;
		}
	}
	free_cache (PLEXLINK);
	return (1);
}

/* add edge to linked list of edges opposite central vertex */

int addLink (struct Plex *plex, long pN, long vbefore, long vat, long vafter)
{
	struct PlexVertex *plexvtx;
	struct PlexLink *plnk;

	plexvtx = plex -> plexvertices + (vat - 1);
	plnk = allocate_PlexLink ();
	if (plnk == NULL) {
		return (0);
	}
	/* where in linked list added, not important */
	plnk -> next = plexvtx -> head;
	plexvtx -> head = plnk;
	plnk -> pN = pN;
	plnk -> vns[0] = vbefore;
	plnk -> vns[1] = vafter;
	return (1);
}

/* sort edges so that the vertices match up,
   so we can identify which provinces the edges belong to,
   and store province numbers in the vertex's hash entry */
   
/* not a sort anymore, just unique partition numbers */
int sortLink (struct PlexVertex *plexvtx, struct Glass *glass)
{
	int found, i, pnidx, npn;
	long linkpns[8];					/* province numbers; should be at most 8 */
	struct PlexLink *plnk;
	
	for (pnidx = 0; pnidx < 8; pnidx++)
		linkpns[pnidx] = 0;
	pnidx = 0;
	
	for (plnk = plexvtx -> head; plnk != NULL; plnk = plnk -> next) {
		found = 0;
		for (i = 0; i < pnidx; i++) {
			/* check for partition number already in our local table */
			if (plnk -> pN == linkpns[i]) {
				found = 1;
				break;
			}
		}
		plnk -> used = 1;
		if (!found) {
			/* new province number */
			if (pnidx > 7) {
				set_error1 ("sortLink: partition number overflow");
				return (0);
			}
			linkpns[pnidx] = plnk -> pN;
			pnidx++;
		}
	}
	npn = pnidx;
	if (npn > 8) {
		set_error1 ("sortLink: partition number overflow");
		return (0);
	}
	glass -> npn = npn;
	for (pnidx = 0; pnidx < npn; pnidx++)
		glass -> pN[pnidx] = linkpns[pnidx];
	return (1);
}

/* INNER DENSITY BOUNDARIES & PROVINCES  * * * * * *  INNER DENSITY BOUNDARIES & PROVINCES  */

/* routines for converting the boundary of the full-density to squares in hash table */

/* find dual Province vertex, edge, triangle or square(s) */

int densityProvince (struct Plex *plex)
{
	int result;
	if (plex -> dimension == 2)
		result = densityProvince2 (plex);
	else if (plex -> dimension == 3)
		result = densityProvince3 (plex);
	else result = 0;
	return (result);
}

int densityProvince2 (struct Plex *plex)
{
	int j, nempty, nsurf, nfull, npn;
	long pN, dVN, pVN,xpN, ypN;
	long x0, y0;
	long xp, yp, x, y, nx, ny, dx, dy;
	long ex, ey, ex0, ex1, ey0, ey1;
	/* clockwise -- because we are subtracting the squares/cubes */
	static long elist[4][4] = { 0, 0, 0, 2, 0, 2, 2, 2, 2, 2, 2, 0, 2, 0, 0, 0};
	struct PlexCube *plexcubes, *cubec, *cubep;
	struct Glass *eglass, *vglass;
	struct Province *provinces;
	struct locube locubes[3][3];


	plexcubes = plex -> plexcubes;
	provinces = plex -> provinces;
	nx = plex -> dimensions[0];
	ny = plex -> dimensions[1];
	/* add each edge of boundary of inner (full) cubes (squares) to hash table */
	for (y = 1; y < ny-1; y++)
		for (x = 1; x < nx-1; x++) {
			cubec = plexcubes + y * nx + x;
			if (cubec -> type == EmptyCube || cubec -> type == PartialCube) continue;
			/* first pass: count surface cubes touching this cube */
			nfull = 0; nsurf = 0; nempty = 0;
			for (dy = -1; dy <= 1; dy++)
				for (dx = -1; dx <= 1; dx++) {
					xp = x + dx; yp = y + dy;
					cubep = plexcubes + yp * nx + xp;
					if (cubep -> type == PartialCube) nsurf++;
					else if (cubep -> type == EmptyCube) nempty++;
					else if (cubep -> type == FullCube) nfull++;
				}
			if (nsurf == 0) continue;	/* too far interior */
			/* second pass: the central cube is solid interior,
			   and at least one of its neighbors is surface;
			   initialize local array */
			for (dy = -1; dy <= 1; dy++)
				for (dx = -1; dx <= 1; dx++) {
					ex = dx + 1; ey = dy + 1;
					xp = x + dx; yp = y + dy;
					cubep = plexcubes + yp * nx + xp;
					if (cubep -> type == PartialCube) {
						pN = cubep -> provinceNumber;
						pVN = (provinces + (pN - 1)) -> vertexNumber;
					}
					else { pN = 0; pVN = 0; }
					locubes[ex][ey].cube = cubep;
					locubes[ex][ey].glass = NULL;
					locubes[ex][ey].dVN = 0;
					locubes[ex][ey].pN = pN;
					locubes[ex][ey].pVN = pVN;
				}
			/* third pass: set up density vertices */
			for (dy = -1; dy <= 1; dy++)
				for (dx = -1; dx <= 1; dx++) {
					ex = dx + 1; ey = dy + 1;
					xp = x + dx; yp = y + dy;
					if (dx == 0 || dy == 0) continue;
					if (locubes[ex][ey].pN > 0 || locubes[1][ey].pN > 0 || locubes[ex][1].pN > 0) {
						x0 = x + (dx+1)/2;
						y0 = y + (dy+1)/2;
						dVN = loraDenVertex (plex, x0, y0, 0L);
						if (dVN == 0) return (0);
						locubes[ex][ey].dVN = dVN;
						vglass = (plex -> plexvertices + (dVN - 1)) -> glass;
						if (vglass == NULL) return (0);
						locubes[ex][ey].glass = vglass;
					}
				}
			/* fourth pass: density edges */
			for (j = 0; j < 4; j++) {
				ex0 = elist[j][0];
				ey0 = elist[j][1];
				ex1 = elist[j][2];
				ey1 = elist[j][3];
				ex = (ex0+ex1)/2;
				ey = (ey0+ey1)/2;
				if (locubes[ex][ey].pN != 0) {
					if (locubes[ex0][ey0].dVN == 0 || locubes[ex1][ey1].dVN == 0) return (0); 
					eglass = loraDenEdge (plex, locubes[ex0][ey0].dVN, locubes[ex1][ey1].dVN);
					if (eglass == NULL) return (0);
					locubes[ex][ey].glass = eglass;
				}
			}
			/* fifth pass: provinces */
			for (dy = -1; dy <= 1; dy++)
				for (dx = -1; dx <= 1; dx++) {
					ex = dx + 1; ey = dy + 1;
					xp = x + dx; yp = y + dy;
					if (dx == 0 && dy == 0) continue;	/* no vertex or edge in center */
					/* corners */
					if (dx != 0 && dy != 0 && (vglass = locubes[ex][ey].glass) != NULL && vglass -> npn == 0) {
						dVN = locubes[ex][ey].dVN;
						xpN = locubes[ex][1].pN;
						pVN = locubes[ex][ey].pVN; pN = locubes[ex][ey].pN;
						ypN = locubes[1][ey].pN;
						npn = 0;	/* may be zeroes inbetween non-zeroes */
						if (dx != dy) {
							if (xpN > 0) vglass -> pN[0] = xpN; npn++;
							if (pN  > 0) vglass -> pN[1] =  pN; npn++;
							if (ypN > 0) vglass -> pN[2] = ypN; npn++;
						}
						else {
							if (xpN > 0) vglass -> pN[2] = xpN; npn++;
							if (pN  > 0) vglass -> pN[1] =  pN; npn++;
							if (ypN > 0) vglass -> pN[0] = ypN; npn++;
						}
						vglass -> npn = npn;
					}
					/* edges */
					if ((dx == 0 || dy == 0) && (eglass = locubes[ex][ey].glass) != NULL && eglass -> npn == 0) {
						pVN = locubes[ex][ey].pVN;
						if (pVN > 0) {
							eglass -> pN[0] = locubes[ex][ey].pN; eglass -> npn = 1;
						}
					}
				}
		}
	return (1);
}

int densityProvince3 (struct Plex *plex)
{
	int j;
	int ndzero, nempty, nsurf, nfull, npn;
	long nsquare;
	long x1, y1, z1, x, y, z, idx, nx, ny, nz;
	long pN, dVN, pVN;
	long xpN, ypN, zpN, xypN, xzpN, yzpN, xyzpN;
	long x0, y0, z0, xp, yp, zp, dx, dy, dz;
	long ex, ey, ez, ex0, ex1, ey0, ey1, ez0, ez1;
	long vns[4], idns[4];
	static long elist[12][6] = {
		0, 0, 0, 0, 2, 0,
		0, 2, 0, 2, 2, 0,
		2, 2, 0, 2, 0, 0,
		2, 0, 0, 0, 0, 0,
		0, 0, 2, 0, 2, 2,
		0, 2, 2, 2, 2, 2,
		2, 2, 2, 2, 0, 2,
		2, 0, 2, 0, 0, 2,
		0, 0, 0, 0, 0, 2,
		2, 0, 2, 2, 0, 0,
		0, 2, 0, 0, 2, 2,
		2, 2, 2, 2, 2, 0
	};
	struct PlexCube *plexcubes, *cubec, *cubex, *cubey, *cubez, *cubep;
	struct Glass *glass, *sglass, *eglass, *vglass;
	struct Province *provinces;
	struct locube locubes[3][3][3];
	char message[MAXLINE];

	provinces = plex -> provinces;
	plexcubes = plex -> plexcubes;
	nsquare = 0;
	nx = plex -> dimensions[0];
	ny = plex -> dimensions[1];
	nz = plex -> dimensions[2];
	vns[0] = 0; vns[1] = 0; vns[2] = 0; vns[3] = 0;
	idns[0] = 0; idns[1] = 0; idns[2] = 0; idns[3] = 0;

	/* add each square of boundary of inner (full) cubes to hash table */
	for (z = 0; z < nz; z++)
		for (y = 0; y < ny; y++)
			for (x = 0; x < nx; x++) {
				idx = z * (nx * ny) + y * nx + x;
				x1 = x + 1; y1 = y + 1; z1 = z + 1;
				cubec = plexcubes + idx;
				if (x < nx - 1) {
					cubex = cubec + 1;
					if (cubex -> type == PartialCube && cubec -> type == FullCube) {
						pN = cubex -> provinceNumber;
						glass = addSquare (plex, x1, y, z, x1, y, z1, x1, y1, z1, x1, y1, z);
						if (glass == NULL) {
							set_error1 ("densityProvince3: addSquare fails for xy");
							return (0);
						}
						glass -> pN[0] = pN;
						glass -> npn = 1;
						nsquare++;
					}
					else if (cubec -> type == PartialCube && cubex -> type == FullCube) {
						pN = cubec -> provinceNumber;
						glass = addSquare (plex, x1, y, z, x1, y1, z, x1, y1, z1, x1, y, z1);
						if (glass == NULL) {
							set_error1 ("densityProvince3: addSquare fails for xy");
							return (0);
						}
						glass -> pN[0] = pN;
						glass -> npn = 1;
						nsquare++;
					}
				}
				if (y < ny - 1) {
					cubey = cubec + nx;
					if (cubey -> type == PartialCube && cubec -> type == FullCube) {
						pN = cubey -> provinceNumber;
						glass = addSquare (plex, x, y1, z, x1, y1, z, x1, y1, z1, x, y1, z1);
						if (glass == NULL) {
							set_error1 ("densityProvince3: addSquare fails for xy");
							return (0);
						}
						glass -> pN[0] = pN;
						glass -> npn = 1;
						nsquare++;
					}
					else if (cubec -> type == PartialCube && cubey -> type == FullCube) {
						pN = cubec -> provinceNumber;
						glass = addSquare (plex, x, y1, z, x, y1, z1, x1, y1, z1, x1, y1, z);
						if (glass == NULL) {
							set_error1 ("densityProvince3: addSquare fails for xy");
							return (0);
						}
						glass -> pN[0] = pN;
						glass -> npn = 1;
						nsquare++;
					}
				}
				if (z < nz - 1) {
					cubez = cubec + (nx * ny);
					if (cubec -> type == PartialCube && cubez -> type == FullCube) {
						pN = cubec -> provinceNumber;
						glass = addSquare (plex, x, y, z1, x1, y, z1, x1, y1, z1, x, y1, z1);
						if (glass == NULL) {
							set_error1 ("densityProvince3: addSquare fails for xy");
							return (0);
						}
						glass -> pN[0] = pN;
						glass -> npn = 1;
						nsquare++;
					}
					else if (cubez -> type == PartialCube && cubec -> type == FullCube) {
						pN = cubez -> provinceNumber;
						glass = addSquare (plex, x, y, z1, x, y1, z1, x1, y1, z1, x1, y, z1);
						if (glass == NULL) {
							set_error1 ("densityProvince3: addSquare fails for xy");
							return (0);
						}
						glass -> pN[0] = pN;
						glass -> npn = 1;
						nsquare++;
					}
				}
			}
	for (z = 1; z < nz-1; z++)
		for (y = 1; y < ny-1; y++)
			for (x = 1; x < nx-1; x++) {
				idx = z * (nx * ny) + y * nx + x;
				cubec = plexcubes + idx;
				if (cubec -> type == EmptyCube || cubec -> type == PartialCube) continue;
				/* first pass: count surface cubes touching this cube */
				nfull = 0; nsurf = 0; nempty = 0;
				for (dz = -1; dz <= 1; dz++)
					for (dy = -1; dy <= 1; dy++)
						for (dx = -1; dx <= 1; dx++) {
							xp = x + dx; yp = y + dy; zp = z + dz;
							cubep = plexcubes + zp * ny * nx + yp * nx + xp;
							if (cubep -> type == PartialCube) nsurf++;
							else if (cubep -> type == EmptyCube) nempty++;
							else if (cubep -> type == FullCube) nfull++;
						}
				if (nsurf == 0) continue;	/* too far interior */
				/* second pass: the central cube is solid interior,
				   and at least one of its neighbors is surface;
				   initialize local array */
				for (dz = -1; dz <= 1; dz++)
					for (dy = -1; dy <= 1; dy++)
						for (dx = -1; dx <= 1; dx++) {
							ex = dx + 1; ey = dy + 1; ez = dz + 1;
							xp = x + dx; yp = y + dy; zp = z + dz;
							cubep = plexcubes + zp * ny * nx + yp * nx + xp;
							if (cubep -> type == PartialCube) {
								pN = cubep -> provinceNumber;
								pVN = (provinces + (pN - 1)) -> vertexNumber;
							}
							else { pN = 0; pVN = 0; }
							locubes[ex][ey][ez].cube = cubep;
							locubes[ex][ey][ez].glass = NULL;
							locubes[ex][ey][ez].dVN = 0;
							locubes[ex][ey][ez].pN = pN;
							locubes[ex][ey][ez].pVN = pVN;
						}
				/* third pass: set up density vertices */
				for (dz = -1; dz <= 1; dz++)
					for (dy = -1; dy <= 1; dy++)
						for (dx = -1; dx <= 1; dx++) {
							ex = dx + 1; ey = dy + 1; ez = dz + 1;
							xp = x + dx; yp = y + dy; zp = z + dz;
							if (dx == 0 || dy == 0 || dz == 0) continue;
							if (locubes[ex][ey][ez].pN > 0 || locubes[1][ey][ez].pN > 0 ||
								locubes[ex][1][ez].pN > 0 || locubes[ex][ey][1].pN > 0) {
								x0 = x + (dx+1)/2;
								y0 = y + (dy+1)/2;
								z0 = z + (dz+1)/2;
								dVN = loraDenVertex (plex, x0, y0, z0);
								if (dVN == 0) return (0);
								locubes[ex][ey][ez].dVN = dVN;
								vglass = (plex -> plexvertices + (dVN - 1)) -> glass;
								if (vglass == NULL) return (0);
								locubes[ex][ey][ez].glass = vglass;
							}
						}
				/* fourth pass: density edges */
				for (j = 0; j < 12; j++) {
					ex0 = elist[j][0];
					ey0 = elist[j][1];
					ez0 = elist[j][2];
					ex1 = elist[j][3];
					ey1 = elist[j][4];
					ez1 = elist[j][5];
					ex = (ex0+ex1)/2;
					ey = (ey0+ey1)/2;
					ez = (ez0+ez1)/2;
					if (locubes[ex][ey][ez].pN != 0) {
						if (locubes[ex0][ey0][ez0].dVN == 0 ||
						locubes[ex1][ey1][ez1].dVN == 0) return (0); 
						eglass = loraDenEdge (plex, locubes[ex0][ey0][ez0].dVN,
							locubes[ex1][ey1][ez1].dVN);
						if (eglass == NULL) return (0);
						locubes[ex][ey][ez].glass = eglass;
					}
				}
				/* fifth pass: provinces */
				for (dz = -1; dz <= 1; dz++)
					for (dy = -1; dy <= 1; dy++)
						for (dx = -1; dx <= 1; dx++) {
							ex = dx + 1; ey = dy + 1; ez = dz + 1;
							xp = x + dx; yp = y + dy;
							ndzero = 0;
							if (dx == 0) ndzero++;
							if (dy == 0) ndzero++;
							if (dz == 0) ndzero++;
							/* no vertex or edge in center */
							if (ndzero == 3) continue;
							/* corners */
							if (ndzero == 0 &&
								(vglass = locubes[ex][ey][ez].glass) != NULL
									&& vglass -> npn == 0) {
								dVN = locubes[ex][ey][ez].dVN;
								xpN = locubes[ex][1][1].pN;
								xypN = locubes[ex][ey][1].pN;
								ypN = locubes[1][ey][1].pN;
								zpN = locubes[1][1][ez].pN;
								xzpN = locubes[ex][1][ez].pN;
								xyzpN = locubes[ex][ey][ez].pN;
								yzpN = locubes[1][ey][ez].pN;
								npn = 0;	/* may be zeroes inbetween non-zeroes */
								if (dx * dy * dz != 1) {
									if (xpN   > 0) vglass -> pN[0] = xpN;   npn++;
									if (xypN  > 0) vglass -> pN[1] = xypN;  npn++;
									if (ypN   > 0) vglass -> pN[2] = ypN;   npn++;
									if (zpN   > 0) vglass -> pN[3] = zpN;   npn++;
									if (xzpN  > 0) vglass -> pN[4] = xzpN;  npn++;
									if (xyzpN > 0) vglass -> pN[5] = xyzpN; npn++;
									if (yzpN  > 0) vglass -> pN[6] = yzpN;  npn++;
								}
								else {
									if (xpN   > 0) vglass -> pN[6] = xpN;   npn++;
									if (xypN  > 0) vglass -> pN[5] = xypN;  npn++;
									if (ypN   > 0) vglass -> pN[4] = ypN;   npn++;
									if (zpN   > 0) vglass -> pN[3] = zpN;   npn++;
									if (xzpN  > 0) vglass -> pN[2] = xzpN;  npn++;
									if (xyzpN > 0) vglass -> pN[1] = xyzpN; npn++;
									if (yzpN  > 0) vglass -> pN[0] = yzpN;  npn++;
								}
								vglass -> npn = npn;
							}
							/* edges */
							if (ndzero == 1 &&
								(eglass = locubes[ex][ey][ez].glass) != NULL
								&& eglass -> npn == 0) {
								npn = 0;	/* may be zeroes inbetween non-zeroes */
								if (dz == 0) {
									xpN = locubes[ex][1][1].pN;
									pN = locubes[ex][ey][1].pN;
									ypN = locubes[1][ey][1].pN;
									if (dx != dy) {
										if (xpN > 0) eglass -> pN[0] = xpN; npn++;
										if (pN  > 0) eglass -> pN[1] =  pN; npn++;
										if (ypN > 0) eglass -> pN[2] = ypN; npn++;
									}
									else {
										if (xpN > 0) eglass -> pN[2] = xpN; npn++;
										if (pN  > 0) eglass -> pN[1] =  pN; npn++;
										if (ypN > 0) eglass -> pN[0] = ypN; npn++;
									}
								}
								if (dy == 0) {
									xpN = locubes[ex][1][1].pN;
									pN = locubes[ex][1][ez].pN;
									zpN = locubes[1][1][ez].pN;
									if (dx == dz) {
										if (xpN > 0) eglass -> pN[0] = xpN; npn++;
										if (pN  > 0) eglass -> pN[1] =  pN; npn++;
										if (zpN > 0) eglass -> pN[2] = zpN; npn++;
									}
									else {
										if (xpN > 0) eglass -> pN[2] = xpN; npn++;
										if (pN  > 0) eglass -> pN[1] =  pN; npn++;
										if (zpN > 0) eglass -> pN[0] = zpN; npn++;
									}
								}
								if (dx == 0) {
									zpN = locubes[1][1][ez].pN;
									pN = locubes[1][ey][ez].pN;
									ypN = locubes[1][ey][1].pN;
									if (dz != dy) {
										if (zpN > 0) eglass -> pN[0] = zpN; npn++;
										if (pN  > 0) eglass -> pN[1] =  pN; npn++;
										if (ypN > 0) eglass -> pN[2] = ypN; npn++;
									}
									else {
										if (zpN > 0) eglass -> pN[2] = zpN; npn++;
										if (pN  > 0) eglass -> pN[1] =  pN; npn++;
										if (ypN > 0) eglass -> pN[0] = ypN; npn++;
									}
								}
								eglass -> npn = npn;
							}
							/* square faces */
							if (ndzero == 2 &&
								(sglass = locubes[ex][ey][ez].glass) != NULL
								&& sglass -> npn == 0) {
								pVN = locubes[ex][ey][ez].pVN;
								if (pVN > 0) {
									sglass -> pN[0] = locubes[ex][ey][ez].pN; sglass -> npn = 1;
								}
							}
						}
			}
	sprintf (message, "%8ld squares added to boundary of inner cubes", nsquare);
	inform (message);
	return (1);
}

/*  JOINS     *      JOINS      *      JOINS  */

int computeJoins (struct Plex *plex)
{
	long nitem;
	char message[MAXLINE];

	nitem = plex -> n_join;
	plex -> joins = (struct PlexJoin *) allocate_objects (PLEXJOIN, nitem);
	if (plex -> joins == NULL)
		return (0);
	if (plex -> dimension == 2) {
		if (!joinType (plex, PolyEdge)) return (0);
		if (!joinType (plex, DenEdge)) return (0);
		if (!joinType (plex, BorderEdge)) return (0);
	}
	else if (plex -> dimension == 3) {
		if (!joinType (plex, PolyTriangle)) return (0);
		if (!joinType (plex, DenSquare)) return (0);
		if (!joinType (plex, BorderTriangle)) return (0);
	}
	if (plex -> ijoin != plex -> n_join) {
		set_error1 ("computeJoins: join count inconsistency");
		sprintf (message, "%8ld != %8ld", plex -> ijoin, plex -> n_join);
		return (0);
	}
	sprintf (message, "%8ld joins computed", plex -> n_join);
	inform (message);
	return (1);
}

int joinType (struct Plex *plex, enum PlexType type)
{
	int j, dim1, dim2, orn, result, iapex;
	long h, maxhash, apexVN;
	long vns[4];
	struct Glass *glass;
	struct Plexi *plexi;
	
	maxhash = plex -> maxhash;
	for (h = 0; h < maxhash; h++) {
		plexi = plex -> plexis + ((int) type) * maxhash + h;
		for (glass = plexi -> head; glass != NULL; glass = glass -> next) {
			for (iapex = 0; iapex < 2; iapex++) {
				apexVN = glass -> apex[iapex];
				if (apexVN == 0) continue;
				dim1 = typeDim (glass -> type);
				dim2 = typeDim (ProvVertex);
				if (dim1 + dim2 != plex -> dimension + 1) {
					set_error1 ("joinType: dimensions don't add up");
					return (0);
				}
				orn = (1 - 2 * iapex) * 1; /* glass -> orn; */
				for (j = 0; j < 4; j++)
					vns[j] = glass -> vertexnumbers[j];
				result = makeJoin (plex, orn, vns, apexVN, type, glass -> subtype);
				if (!result)
					return (0);
			}
		}
	}
	return (1);
}

int makeJoin (struct Plex *plex, int orn, long vns[4], long apex, enum PlexType type, enum PlexType subtype)
{
	int k, i;
	char message[MAXLINE];
	struct PlexJoin *join;
	struct PlexVertex *plexvtx;
	
	join = plex -> joins + plex -> ijoin;
	plex -> ijoin++;
	if (plex -> ijoin > plex -> n_join) {
		set_error1 ("makeJoin: join count inconsistency");
		return (0);
	}
	if (apex <= 0 || apex > plex -> n_vertex) {
		sprintf (message, "makeJoin: invalid apex %8ld", apex);
		set_error1 (message);
		return (0);
	}
	plexvtx = plex -> plexvertices + apex - 1;
	for (k = 0; k < 3; k++) {
		if (plexvtx -> center[k] < -1000.0 || plexvtx -> center[k] > 1000.0) {
			sprintf (message, "makeJoin: invalid vertex coordinate: %12.3f",
				plexvtx -> center[k]);
			return (0);
		}
		join -> apex[k] = plexvtx -> center[k];
	}
	for (i = 0; i < 4; i++)
		join -> vertices[i] = vns[i];
	join -> orn = orn;
	join -> type = type;
	join -> subtype = subtype;
	return (1);
}

/* Intermediate Level Math */

int insidePoly (struct Plex *plex, double center[3])
{
	int inside;
	
	if (plex -> dimension == 2)
		inside = insidePolygon (plex, center);
	else if (plex -> dimension == 3)
		inside = insidePolyhedron (plex, center);
	else inside = 0;
	return (inside);
}

/* compute winding number */

int insidePolygon (struct Plex *plex, double pnt[3])
{
	int k, inside;
	long h, v0, v1, maxhash;
	double angle, winding;
	double vect0[3], vect1[3], zaxis[3];
	struct Plexi *plexi;
	struct Glass *glass;
	struct PlexVertex *vtx0, *vtx1;

	/* get axis */
	for (k = 0; k < 3; k++)
		zaxis[k] = (k == 2);
	winding = 0.0;
	maxhash = plex -> maxhash;
	for (h = 0; h < maxhash; h++) {
		plexi = plex -> plexis + ((int) PolyEdge) * maxhash + h;
		for (glass = plexi -> head; glass != NULL; glass = glass -> next) {
	        /* compute radial vectors */
			v0 = glass -> vertexnumbers[0];
			v1 = glass -> vertexnumbers[1];
			vtx0 = plex -> plexvertices + (v0 - 1);
			vtx1 = plex -> plexvertices + (v1 - 1);
			for (k = 0; k < 3; k++) {
				vect0[k] = vtx0 -> center[k] - pnt[k];
				vect1[k] = vtx1 -> center[k] - pnt[k];
			}
			if (!normalize (vect0)) return (0);
			if (!normalize (vect1)) return (0);
			angle = triple_angle (vect0, vect1, zaxis);
			winding += angle;
		}
	}
	inside = (winding > PI);
	return (inside);
}

int insidePolyhedron (struct Plex *plex, double pnt[3])
{
	int j, k, t;
	int is_inside;
	long vn[3];
	double liquid, solid, delta, dt, circumference, area, approx, dist;
	double centers[3][3], center[3];
	double vector[3], axis[3], vect1[3], vect2[3];
	double *vs[3];

	/* initialize solid angle subtended by polyhedron as seen by point */
	solid = 0.0;
	liquid = 0.0;

	for (t = 0; t < plex -> mtriangle; t++) {
		for (j = 0; j < 3; j++) {
			vn[j] = *(plex -> triangles + 6 * t + 3 + j);
			vs[j] = plex -> vertices + 3 * (vn[j]-1);
			for (k = 0; k < 3; k++)
				centers[j][k] = *(vs[j] + k);
		}
		/* skip very small, possibly pathological, triangles */
		circumference = distance (vs[0], vs[1]) + distance (vs[1], vs[2]) + distance (vs[2], vs[0]);
		if (circumference < 0.01) continue;
		for (k = 0; k < 3; k++) {
			center[k] = (centers[0][k]+centers[1][k]+centers[2][k])/3.0;
			vector[k] = center[k] - pnt[k];
		}
		if (!normalize(vector)) continue;
		for (k = 0; k < 3; k++) {
			vect1[k] = centers[1][k] - centers[0][k];
			vect2[k] = centers[2][k] - centers[0][k];
		}
		cross (vect1, vect2, axis);
		area = norm (axis) / 2.0;
		if (area <= 0.0) continue;
		if (!normalize(axis)) continue;
		dt = dot_product (axis, vector);
		dist = distance (pnt, center);
		approx = (area * dt / (dist * dist));
		liquid += approx;
		/* compute solid angle of triangles as seen by point */
		delta = tetra_solid_angle (pnt, vs[0], vs[1], vs[2]);
		solid += delta;
	}
	is_inside = (liquid > 2 * PI);
	return (is_inside);
}

int computeVolume (struct Plex *plex)
{
	int j, k;
	long e, t, vn;
	double volume, simplex_vol;
	double vect[3][3];
	double fg[3];
	double *f;
	char message[MAXLINE];

	volume = 0.0;
	if (plex -> dimension == 2) {
		for (e = 0; e < plex -> medge; e++) {
			for (j = 0; j < 2; j++) {
				vn = *(plex -> edges + 2 * e + j);
				f = plex -> vertices + 3 * (vn - 1);
				for (k = 0; k < 2; k++)
					vect[j][k] = *(f + k);
			}
			simplex_vol = (vect[0][0] * vect[1][1] - vect[0][1] * vect[1][0]) / 2.0;
			volume += simplex_vol;
		}
	}
	else if (plex -> dimension == 3) {
		for (t = 0; t < plex -> mtriangle; t++) {
			for (j = 0; j < 3; j++) {
				vn = *(plex -> triangles + 6 * t + 3 + j);
				f = plex -> vertices + 3 * (vn - 1);
				for (k = 0; k < 3; k++)
					vect[j][k] = *(f + k);
			}
			fg[0] = vect[0][1] * vect[1][2] - vect[0][2] * vect[1][1];
			fg[1] = vect[0][2] * vect[1][0] - vect[0][0] * vect[1][2];
			fg[2] = vect[0][0] * vect[1][1] - vect[0][1] * vect[1][0];
			simplex_vol = (vect[2][0] * fg[0] + vect[2][1] * fg[1] + vect[2][2] * fg[2]) / 6.0;
			volume += simplex_vol;
		}
	}
	plex -> polyvolume = volume;
	sprintf (message, "%8.3f polyhedron volume", volume);
	inform (message);
	return (1);
}


/*  APICES    *   *   *     APICES     */

int allApices (struct Plex *plex)
{
	int dim;
	char message[MAXLINE];

	dim = plex -> dimension;

	if (dim == 3) {
		if (!setApices (plex, PolyTriangle, 1))
			return (0);
		if (!setApices (plex, DenSquare, 1))
			return (0);
		if (!setApices (plex, BorderTriangle, 2))
			return (0);
		sprintf (message, "%8ld apex joins created", plex -> n_join);
		inform (message);
	}
	else if (dim == 2) {
		if (!setApices (plex, PolyEdge, 1))
			return (0);
		if (!setApices (plex, DenEdge, 1))
			return (0);
		if (!setApices (plex, BorderEdge, 2))
			return (0);
	}
	return (1);
}

int setApices (struct Plex *plex, enum PlexType type, int napex)
{
	int i;
	long h, maxhash, pN;
	struct Plexi *plexi;
	struct Province *prov;
	struct Glass *glass;
	char message[MAXLINE];

	maxhash = plex -> maxhash;
	plexi = plex -> plexis + maxhash * type;
	for (h = 0; h < maxhash; h++) {
		for (glass = (plexi + h) -> head; glass != NULL; glass = glass -> next) {
			/* glass -> orn = 1; */
			for (i = 0; i < napex; i++) {
				pN = glass -> pN[i];
				if (pN == 0) continue;
				prov = plex -> provinces + (pN - 1);
				if (prov -> vertexNumber <= 0 ||
					prov -> vertexNumber > plex -> n_vertex) {
					sprintf (message,
						"setApices: invalid province vertex number %8ld",
						prov -> vertexNumber);
					set_error1 (message);
					return (0);
				}
				glass -> apex[i] = prov -> vertexNumber;
				plex -> n_join++;
			}
		}
	}
	return (1);
}


struct PlexLink *allocate_PlexLink ()
{
	struct PlexLink *plk;

	/* allocate memory */
	plk = (struct PlexLink *) allocate_object (PLEXLINK);
	if (error ()) return (NULL);
	if (plk == NULL) {
		set_error1 ("allocate_PlexLink: ran out of memory");
		return(NULL);
	}
	return (plk);
}

void free_PlexLink (struct PlexLink *plk)
{
	free_object (PLEXLINK, (short *) plk);
	if (error()) return;
}



