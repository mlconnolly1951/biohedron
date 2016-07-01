/* MSP: Piecewise-Linear complEX */
/* Copyright 1995 by Michael L. Connolly */
/* March 7, 2000 */

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

/* conversion from density to polyhedron (polygon) */

struct msdata *doXelp (long dimension, long dimensions[3], double bounds[2][3], double ctrlev, struct vanity *vanities)
{
	int j, k, dim;
	long nx, ny, nz;
	double xwidth, ywidth, zwidth;
	struct Plex *plex;
	struct msdata *msd;

	plex = (struct Plex *) allocate_object (PLEX);
	if (plex == NULL) {
		set_error1 ("doXelp: not enough memory for Plex");
		return (0);
	}
	plex -> dimension = dimension;
	dim = plex -> dimension;
	for (k = 0; k < 3; k++)
		plex -> dimensions[k] = dimensions[k];
	if (dimension == 2) plex -> dimensions[2] = 1;
	nx = (plex -> dimensions[0]);
	ny = (plex -> dimensions[1]);
	nz = (plex -> dimensions[2]);
	plex -> n_cube = (nx * ny * nz);
	if (plex -> n_cube == 0) return (0);
	if (vanities == NULL) return (NULL);
	plex -> ctrlev = ctrlev;
	for (j = 0; j < 2; j++)
		for (k = 0; k < 3; k++)
			plex -> bounds[j][k] = bounds[j][k];
	plex -> maxhash = (integer) 256L;
	plex -> maxvertex = (256L + 4 * plex -> n_cube);
	plex -> maxedge = plex -> maxvertex;
	plex -> maxtriangle =  plex -> maxedge;
	if (sizeof (integer) == 2) {
		if (plex -> n_cube > 32767) {
			set_error1 ("doXlep: n_cube > max short integer");
			return (NULL);
		}
		if (plex -> maxvertex > 32767) {
			set_error1 ("doXlep: maxvertex > max short integer");
			return (NULL);
		}
		if (plex -> maxedge > 32767) {
			set_error1 ("doXlep: maxedge > max short integer");
			return (NULL);
		}
		if (plex -> maxtriangle > 32767) {
			set_error1 ("doXlep: maxtriangle > max short integer");
			return (NULL);
		}
	}
	plex -> glassblock = 1024L;
	xwidth = (plex -> bounds[1][0] - plex -> bounds[0][0]);
	ywidth = (plex -> bounds[1][1] - plex -> bounds[0][1]);
	zwidth = (plex -> bounds[1][2] - plex -> bounds[0][2]);
	if (xwidth <= 0.0 || ywidth <= 0.0) return (NULL);
	if (plex -> dimension == 2)
		plex -> cubeVolume = xwidth * ywidth;
	else plex -> cubeVolume = xwidth * ywidth * zwidth;
	plex -> cubeVolume /= plex -> n_cube;
	if (plex -> cubeVolume <= 0.0) return (NULL);

	if (!initializePlex (plex)) return (NULL);

	if (!setOccupancies (plex, vanities)) return (NULL);

	if (!identifyContourCubes (plex)) return (NULL);

	if (!createCubeVertices (plex)) return (NULL);

	if (!createCubeEdges (plex)) return (NULL);

	if (!createCubeTriangles (plex)) return (NULL);

	if (dim == 3) {
		if (!createCubeTetrahedra (plex)) return (NULL);
	}

	if (!createContourVertices (plex)) return (NULL);

	if (!createContourEdges (plex)) return (NULL);

	if (dim == 3) {
		if (!createContourTriangles (plex)) return (NULL);
		if (!computeBoundary (plex)) return (NULL);
	}

	msd = transferContour (plex);
	if (msd == NULL) return (NULL);

	if (!freePlex (plex)) return (NULL);

	return (msd);
}

int setOccupancies (struct Plex *plex, struct vanity *vanities)
{
	int k;
	long i;
	struct vanity *van;
	struct PlexCube *cube, *plexcubes;
	
	plexcubes = plex -> plexcubes;
	for (i = 0; i < plex -> n_cube; i++) {
		cube = plexcubes + i;
		van = vanities + i;
		cube -> occupancy = van -> scalar;
		for (k = 0; k < 3; k++)
			cube -> vector[k] = van -> vector[k];
	}
	return (1);
}

int identifyContourCubes (struct Plex *plex)
{
	int nzero, atbelow, atabove, atedge;
	int dim;
	long nabove, nbelow, njustabove, njustbelow, nat, nouterabove, nouterbelow;
	long x, y, z, idx, nx, ny, nz;
	long xp, yp, zp, dx, dy, dz;
	char message[MAXLINE];
	struct PlexCube *plexcubes, *cube;
	struct PlexCube *cube0, *cubep, *cubec;
	
	dim = plex -> dimension;
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
				if (cube0 -> occupancy > plex -> ctrlev)
					cube0 -> type = AboveCube;
				else if (cube0 -> occupancy < plex -> ctrlev)
					cube0 -> type = BelowCube;
				else cube0 -> type = AtCube;
			}
	if (dim == 2) {
		for (y = 1; y < ny-1; y++)
			for (x = 1; x < nx-1; x++) {
				idx = y * nx + x;
				cubec = plexcubes + idx;
				if (cubec -> type == BelowCube ||
					cubec -> type == JustBelow) {
					for (dy = -1; dy <= 1; dy++)
						for (dx = -1; dx <= 1; dx++) {
							nzero = 0;
							if (dx == 0) nzero++;
							if (dy == 0) nzero++;
							if (nzero == dim) continue;
							xp = x + dx; yp = y + dy;
							cubep = plexcubes + yp * nx + xp;
							if (cubep -> type == AboveCube ||
								cubep -> type == JustAbove) {
								cubep -> type = JustAbove;
								cubec -> type = JustBelow;
							}
						}
				}
			}
	}
	else {
		for (z = 0; z < nz; z++)
			for (y = 0; y < ny; y++)
				for (x = 0; x < nx; x++) {
					idx = z * (nx * ny) + y * nx + x;
					cubec = plexcubes + idx;
					if (cubec -> type == BelowCube || cubec -> type == JustBelow) {
						for (dz = -1; dz <= 1; dz++)
							for (dy = -1; dy <= 1; dy++)
								for (dx = -1; dx <= 1; dx++) {
									nzero = 0;
									if (dx == 0) nzero++;
									if (dy == 0) nzero++;
									if (dz == 0) nzero++;
									if (nzero == dim) continue;
									xp = x + dx; if (xp < 0 || xp >= nx) continue;
									yp = y + dy; if (yp < 0 || yp >= ny) continue;
									zp = z + dz; if (zp < 0 || zp >= nz) continue;
									cubep = plexcubes + zp * ny * nx + yp * nx + xp;
									if (cubep -> type == AboveCube ||
										cubep -> type == JustAbove) {
										cubep -> type = JustAbove;
										cubec -> type = JustBelow;
									}
								}
					}
				}
	}
	nabove = nbelow = njustabove = njustbelow = nat = 0;
	nouterabove = nouterbelow = 0;
	for (z = 0; z < nz; z++) 
		for (y = 0; y < ny; y++)
			for (x = 0; x < nx; x++) {
				idx = z * (nx * ny) + y * nx + x;
				cube = plex -> plexcubes + idx;
				if (cube -> type == JustBelow) njustbelow++;
				else if (cube -> type == JustAbove) njustabove++;
				else if (cube -> type == BelowCube) nbelow++;
				else if (cube -> type == AboveCube) nabove++;
				else if (cube -> type == AtCube) nat++;
				atedge = (x==0 || x==nx-1 || y==0 || y==ny-1 || z==0 || z==nz-1);
				atabove = (cube -> type == JustAbove ||
					cube -> type == AboveCube || cube -> type == AtCube);
				atbelow = (cube -> type == JustBelow ||
					cube -> type == BelowCube || cube -> type == AtCube);
				if (atedge && atbelow) nouterbelow++;
				if (atedge && atabove) nouterabove++;
			}
	if (plex -> dimension == 3) {
		sprintf (message, "%8ld cubes below contour level", nbelow);
		inform (message);
		sprintf (message, "%8ld cubes above contour level", nabove);
		inform (message);
		sprintf (message, "%8ld cubes just below", njustbelow);
		inform (message);
		sprintf (message, "%8ld cubes just above", njustabove);
		inform (message);
		sprintf (message, "%8ld cubes at contour level", nat);
		inform (message);
		sprintf (message, "%8ld outer boundary cubes at or above contour level", nouterabove);
		inform (message);
		sprintf (message, "%8ld outer boundary cubes at or below contour level", nouterbelow);
		inform (message);
	}
	return (1);
}
	
int createCubeVertices (struct Plex *plex)
{
	int k;
	long n;
	long vN;
	long x, y, z, idx, nx, ny, nz;
	long pN;
	long indices[3];
	char message[MAXLINE];
	struct PlexVertex *plexvertices, *plexvtx;
	struct Province *provinces, *prov;
	struct PlexCube *cube;
	
	plexvertices = plex -> plexvertices;
	nx = plex -> dimensions[0];
	ny = plex -> dimensions[1];
	nz = plex -> dimensions[2];
	n = 0;
	for (z = 0; z < nz; z++) 
		for (y = 0; y < ny; y++)
			for (x = 0; x < nx; x++) {
				idx = z * (nx * ny) + y * nx + x;
				cube = plex -> plexcubes + idx;
				if (cube -> type == JustBelow ||
					cube -> type == JustAbove) {
					n++;
				}
			}
	if (n <= 0) return (0);
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
				if (cube -> type == JustBelow ||
					cube -> type == JustAbove) {
					indices[0] = x; indices[1] = y; indices[2] = z;
					pN = makeProvince (plex, indices);
					if (pN == 0) return (0);
					cube -> provinceNumber = pN;
					prov = provinces + (pN-1);
					vN = prov -> vertexNumber;
					plexvtx = plexvertices + (vN-1);
					plexvtx -> value = cube -> occupancy;
					for (k = 0; k < 3; k++)
						plexvtx -> vector[k] = cube -> vector[k];
				}
			}
	if (plex -> dimension == 3) {
		sprintf (message, "%8ld scaffolding vertices", plex -> nglass[ProvVertex]);
		inform (message);
	}
	return (1);
}

/* connect a pair of vertices only if
   they lie on opposite sides of the contour level */

int createCubeEdges (struct Plex *plex)
{
	long x, y, z, nx, ny, nz, ncreated;
	long pN000, pN001, pN010, pN100;
	long vN000, vN001, vN010, vN100;
	char message[MAXLINE];
	struct PlexCube *plexcubes;
	struct PlexCube *cube000, *cube001, *cube010;
	struct PlexCube *cube100;
	struct Province *prov000, *prov001, *prov010;
	struct Province *prov100;
	struct Province *provinces;
	struct Glass *eglass;
	
	plexcubes = plex -> plexcubes;
	provinces = plex -> provinces;
	ncreated = 0;
	nx = plex -> dimensions[0];
	ny = plex -> dimensions[1];
	nz = plex -> dimensions[2];
	if (plex -> dimension == 2) {
		for (x = 0; x < nx; x++)
			for (y = 0; y < ny-1; y++) {
				cube000 = plexcubes + y * nx + x;
				cube010 = plexcubes + (y+1) * nx + x;
				pN000 = cube000 -> provinceNumber;
				pN010 = cube010 -> provinceNumber;
				if (pN000 != 0 && pN010 != 0) {
					prov000 = (provinces + pN000 - 1);
					prov010 = (provinces + pN010 - 1);
					vN000 = prov000 -> vertexNumber;
					vN010 = prov010 -> vertexNumber;
					eglass = makeCubeEdge (plex, vN000, vN010);
					if (eglass == NULL) return (0);
					ncreated++;
				}
			}
		for (y = 0; y < ny; y++)
			for (x = 0; x < nx-1; x++) {
				cube001 = plexcubes + y * nx + (x+1);
				cube000 = plexcubes + y * nx + x;
				pN001 = cube001 -> provinceNumber;
				pN000 = cube000 -> provinceNumber;
				if (pN001 != 0 && pN000 != 0) {
					prov000 = (provinces + pN000 - 1);
					prov001 = (provinces + pN001 - 1);
					vN000 = prov000 -> vertexNumber;
					vN001 = prov001 -> vertexNumber;
					eglass = makeCubeEdge (plex, vN000, vN001);
					if (eglass == NULL) return (0);
					ncreated++;
				}
			}
	}
	else if (plex -> dimension == 3) {
		for (z = 0; z < nz; z++)
			for (x = 0; x < nx; x++)
				for (y = 0; y < ny-1; y++) {
					cube000 = plexcubes + z * ny * nx + y * nx + x;
					cube010 = plexcubes + z * ny * nx + (y+1) * nx + x;
					pN000 = cube000 -> provinceNumber;
					pN010 = cube010 -> provinceNumber;
					if (pN000 != 0 && pN010 != 0) {
						prov000 = (provinces + pN000 - 1);
						prov010 = (provinces + pN010 - 1);
						vN000 = prov000 -> vertexNumber;
						vN010 = prov010 -> vertexNumber;
						eglass = makeCubeEdge (plex, vN000, vN010);
						if (eglass == NULL) return (0);
						ncreated++;
					}
				}
		for (z = 0; z < nz; z++)
			for (y = 0; y < ny; y++)
				for (x = 0; x < nx-1; x++) {
					cube000 = plexcubes + z * ny * nx + y * nx + x;
					cube001 = plexcubes + z * ny * nx + y * nx + (x+1);
					pN000 = cube000 -> provinceNumber;
					pN001 = cube001 -> provinceNumber;
					if (pN001 != 0 && pN000 != 0) {
						prov000 = (provinces + pN000 - 1);
						prov001 = (provinces + pN001 - 1);
						vN000 = prov000 -> vertexNumber;
						vN001 = prov001 -> vertexNumber;
						eglass = makeCubeEdge (plex, vN000, vN001);
						if (eglass == NULL) return (0);
						ncreated++;
					}
				}
		for (x = 0; x < nx; x++)
			for (y = 0; y < ny; y++)
				for (z = 0; z < nz-1; z++) {
					cube000 = plexcubes + z * ny * nx + y * nx + x;
					cube100 = plexcubes + (z+1) * ny * nx + y * nx + x;
					pN000 = cube000 -> provinceNumber;
					pN100 = cube100 -> provinceNumber;
					if (pN100 != 0 && pN000 != 0) {
						prov000 = (provinces + pN000 - 1);
						prov100 = (provinces + pN100 - 1);
						vN000 = prov000 -> vertexNumber;
						vN100 = prov100 -> vertexNumber;
						eglass = makeCubeEdge (plex, vN000, vN100);
						if (eglass == NULL) return (0);
						ncreated++;
					}
				}
	}
	if (plex -> dimension == 3) {
		sprintf (message, "%8ld scaffolding edges", ncreated);
		inform (message);
	}
	return (1);
}

int createCubeTriangles (struct Plex *plex)
{
	int nabove, nbelow, odd;
	long x, y, z, nx, ny, nz, ncreated, ecreated;
	long pN000, pN001, pN010, pN011, pN100, pN101, pN110;
	long vN000, vN001, vN010, vN011, vN100, vN101, vN110;
	char message[MAXLINE];
	struct PlexCube *plexcubes;
	struct PlexCube *cube000, *cube001, *cube010, *cube011;
	struct PlexCube *cube100, *cube101, *cube110;
	struct Province *prov000, *prov001, *prov010, *prov011;
	struct Province *prov100, *prov101, *prov110;
	struct Province *provinces;
	struct Glass *eg000, *eg010;
	struct Glass *eg100;
	struct Glass *eglass, *tglass;
	
	plexcubes = plex -> plexcubes;
	provinces = plex -> provinces;
	ncreated = 0; ecreated = 0;;
	nx = plex -> dimensions[0];
	ny = plex -> dimensions[1];
	nz = plex -> dimensions[2];
	z = 0;
	if (plex -> dimension == 2) {
		for (x = 0; x < nx-1; x++)
			for (y = 0; y < ny-1; y++) {
				odd = (x + y + z) % 2;
				cube000 = cube001 = cube010 = cube011 = NULL;
				prov000 = prov001 = prov010 = prov011 = NULL;
				pN000 = pN001 = pN010 = pN011 = 0;
				vN000 = vN001 = vN010 = vN011 = 0;
				cube000 = plexcubes + y * nx + x;
				cube010 = plexcubes + (y+1) * nx + x;
				cube011 = plexcubes + (y+1) * nx + (x+1);
				cube001 = plexcubes + y * nx + (x+1);
				pN000 = cube000 -> provinceNumber;
				pN010 = cube010 -> provinceNumber;
				pN011 = cube011 -> provinceNumber;
				pN001 = cube001 -> provinceNumber;
				nabove = 0; nbelow = 0;
				if (pN000 != 0) {
					prov000 = (provinces + pN000 - 1);
					vN000 = prov000 -> vertexNumber;
					if (cube000 -> type == JustBelow) nbelow++;
					else if (cube000 -> type == JustAbove) nabove++;
				}
				else prov000 = NULL;
				if (pN010 != 0) {
					prov010 = (provinces + pN010 - 1);
					vN010 = prov010 -> vertexNumber;
					if (cube010 -> type == JustBelow) nbelow++;
					else if (cube010 -> type == JustAbove) nabove++;
				}
				else prov010 = NULL;
				if (pN011 != 0) {
					prov011 = (provinces + pN011 - 1);
					vN011 = prov011 -> vertexNumber;
					if (cube011 -> type == JustBelow) nbelow++;
					else if (cube011 -> type == JustAbove) nabove++;
				}
				else prov011 = NULL;
				if (pN001 != 0) {
					prov001 = (provinces + pN001 - 1);
					vN001 = prov001 -> vertexNumber;
					if (cube001 -> type == JustBelow) nbelow++;
					else if (cube001 -> type == JustAbove) nabove++;
				}
				else prov001 = NULL;
				if (nbelow == 0 || nabove == 0) continue;
				if (nabove + nbelow < 3) continue;
				if (nabove + nbelow == 4) {
					/* create 1 new edge and 2 triangles */
					if (!odd) {
						eg000 = makeCubeEdge (plex, vN000, vN011);
						if (eg000 == NULL) return (0);
						ecreated++;
						tglass = makeCubeTriangle (plex, vN000, vN011, vN010);
						if (tglass == NULL) return (0);
						ncreated++;
						tglass = makeCubeTriangle (plex, vN011, vN000, vN001);
						if (tglass == NULL) return (0);
						ncreated++;
					}
					else {
						eg010 = makeCubeEdge (plex, vN001, vN010);
						if (eg010 == NULL) return (0);
						ecreated++;
						tglass = makeCubeTriangle (plex, vN001, vN010, vN000);
						if (tglass == NULL) return (0);
						ncreated++;
						tglass = makeCubeTriangle (plex, vN010, vN001, vN011);
						if (tglass == NULL) return (0);
						ncreated++;
					}
				}
				else if (vN000 != 0 && vN010 != 0 && vN011 != 0) {
					/* create 1 triangle */
					eglass = makeCubeEdge (plex, vN000, vN011);
					if (eglass == NULL) return (0);
					ecreated++;
					tglass = makeCubeTriangle (plex, vN000, vN010, vN011);
					if (tglass == NULL) return (0);
					ncreated++;
				}
				else if (vN010 != 0 && vN011 != 0 && vN001 != 0) {
					/* create 1 triangle */
					eglass = makeCubeEdge (plex, vN010, vN001);
					if (eglass == NULL) return (0);
					ecreated++;
					tglass = makeCubeTriangle (plex, vN010, vN011, vN001);
					if (tglass == NULL) return (0);
					ncreated++;
				}
				else if (vN011 != 0 && vN001 != 0 && vN000 != 0) {
					/* create 1 triangle */
					eglass = makeCubeEdge (plex, vN000, vN011);
					if (eglass == NULL) return (0);
					ecreated++;
					tglass = makeCubeTriangle (plex, vN011, vN001, vN000);
					if (tglass == NULL) return (0);
					ncreated++;
				}
				else if (vN001 != 0 && vN000 != 0 && vN010 != 0) {
					/* create 1 triangle */
					eglass = makeCubeEdge (plex, vN001, vN010);
					if (eglass == NULL) return (0);
					ecreated++;
					tglass = makeCubeTriangle (plex, vN001, vN000, vN010);
					if (tglass == NULL) return (0);
					ncreated++;
				}
			}
	}
	else if (plex -> dimension == 3) {
		/* XY squares */
		for (z = 0; z < nz; z++)
			for (x = 0; x < nx-1; x++)
				for (y = 0; y < ny-1; y++) {
					odd = (x + y + z) % 2;
					cube000 = cube001 = cube010 = cube011 = NULL;
					prov000 = prov001 = prov010 = prov011 = NULL;
					pN000 = pN001 = pN010 = pN011 = 0;
					vN000 = vN001 = vN010 = vN011 = 0;
					cube000 = plexcubes + z * ny * nx + y * nx + x;
					cube010 = plexcubes + z * ny * nx + (y+1) * nx + x;
					cube011 = plexcubes + z * ny * nx + (y+1) * nx + (x+1);
					cube001 = plexcubes + z * ny * nx + y * nx + (x+1);
					pN000 = cube000 -> provinceNumber;
					pN010 = cube010 -> provinceNumber;
					pN011 = cube011 -> provinceNumber;
					pN001 = cube001 -> provinceNumber;
					nabove = 0; nbelow = 0;
					if (pN000 != 0) {
						prov000 = (provinces + pN000 - 1);
						vN000 = prov000 -> vertexNumber;
						if (cube000 -> type == JustBelow) nbelow++;
						else if (cube000 -> type == JustAbove) nabove++;
					}
					else prov000 = NULL;
					if (pN010 != 0) {
						prov010 = (provinces + pN010 - 1);
						vN010 = prov010 -> vertexNumber;
						if (cube010 -> type == JustBelow) nbelow++;
						else if (cube010 -> type == JustAbove) nabove++;
					}
					else prov010 = NULL;
					if (pN011 != 0) {
						prov011 = (provinces + pN011 - 1);
						vN011 = prov011 -> vertexNumber;
						if (cube011 -> type == JustBelow) nbelow++;
						else if (cube011 -> type == JustAbove) nabove++;
					}
					else prov011 = NULL;
					if (pN001 != 0) {
						prov001 = (provinces + pN001 - 1);
						vN001 = prov001 -> vertexNumber;
						if (cube001 -> type == JustBelow) nbelow++;
						else if (cube001 -> type == JustAbove) nabove++;
					}
					else prov001 = NULL;
					if (nabove + nbelow < 4) continue;
					/* create 1 new edge and 2 triangles */
					if (!odd) {
						eg000 = makeCubeEdge (plex, vN000, vN011);
						if (eg000 == NULL) {
							sprintf (message, "%8ld %8ld", vN000, vN011);
							set_error1 (message);
							return (0);
						}
						ecreated++;
						tglass = makeCubeTriangle (plex, vN000, vN011, vN010);
						if (tglass == NULL) {
							sprintf (message, "%8ld %8ld %8ld", vN000, vN011, vN010);
							set_error1 (message);
							return (0);
						}
						ncreated++;
						tglass = makeCubeTriangle (plex, vN011, vN000, vN001);
						if (tglass == NULL) {
							sprintf (message, "%8ld %8ld %8ld", vN011, vN000, vN001);
							set_error1 (message);
							return (0);
						}
						ncreated++;
					}
					else {
						eg010 = makeCubeEdge (plex, vN001, vN010);
						if (eg010 == NULL) {
							sprintf (message, "%8ld %8ld", vN001, vN010);
							set_error1 (message);
							return (0);
						}
						ecreated++;
						tglass = makeCubeTriangle (plex, vN001, vN010, vN000);
						if (tglass == NULL) {
							sprintf (message, "%8ld %8ld %8ld", vN001, vN010, vN000);
							set_error1 (message);
							return (0);
						}
						ncreated++;
						tglass = makeCubeTriangle (plex, vN010, vN001, vN011);
						if (tglass == NULL) {
							sprintf (message, "%8ld %8ld %8ld", vN010, vN001, vN011);
							set_error1 (message);
							return (0);
						}
						ncreated++;
					}
				}
		/* XZ squares */
		for (y = 0; y < ny; y++)
			for (z = 0; z < nz-1; z++)
				for (x = 0; x < nx-1; x++) {
					odd = (x + y + z) % 2;
					cube000 = cube001 = cube100 = cube101 = NULL;
					prov000 = prov001 = prov100 = prov101 = NULL;
					pN000 = pN001 = pN100 = pN101 = 0;
					vN000 = vN001 = vN100 = vN101 = 0;
					cube000 = plexcubes + z * ny * nx + y * nx + x;
					cube100 = plexcubes + (z+1) * ny * nx + y * nx + x;
					cube101 = plexcubes + (z+1) * ny * nx + y * nx + (x+1);
					cube001 = plexcubes + z * ny * nx + y * nx + (x+1);
					pN000 = cube000 -> provinceNumber;
					pN100 = cube100 -> provinceNumber;
					pN101 = cube101 -> provinceNumber;
					pN001 = cube001 -> provinceNumber;
					nabove = 0; nbelow = 0;
					if (pN000 != 0) {
						prov000 = (provinces + pN000 - 1);
						vN000 = prov000 -> vertexNumber;
						if (cube000 -> type == JustBelow) nbelow++;
						else if (cube000 -> type == JustAbove) nabove++;
					}
					else prov000 = NULL;
					if (pN100 != 0) {
						prov100 = (provinces + pN100 - 1);
						vN100 = prov100 -> vertexNumber;
						if (cube100 -> type == JustBelow) nbelow++;
						else if (cube100 -> type == JustAbove) nabove++;
					}
					else prov100 = NULL;
					if (pN101 != 0) {
						prov101 = (provinces + pN101 - 1);
						vN101 = prov101 -> vertexNumber;
						if (cube101 -> type == JustBelow) nbelow++;
						else if (cube101 -> type == JustAbove) nabove++;
					}
					else prov101 = NULL;
					if (pN001 != 0) {
						prov001 = (provinces + pN001 - 1);
						vN001 = prov001 -> vertexNumber;
						if (cube001 -> type == JustBelow) nbelow++;
						else if (cube001 -> type == JustAbove) nabove++;
					}
					else prov001 = NULL;
					if (nabove + nbelow < 4) continue;
					/* create 1 new edge and 2 triangles */
					if (!odd) {
						eg000 = makeCubeEdge (plex, vN000, vN101);
						if (eg000 == NULL) {
							sprintf (message, "%8ld %8ld", vN000, vN101);
							set_error1 (message);
							return (0);
						}
						ecreated++;
						tglass = makeCubeTriangle (plex, vN000, vN101, vN100);
						if (tglass == NULL) {
							sprintf (message, "%8ld %8ld %8ld", vN000, vN101, vN100);
							set_error1 (message);
							return (0);
						}
						ncreated++;
						tglass = makeCubeTriangle (plex, vN101, vN000, vN001);
						if (tglass == NULL) {
							sprintf (message, "%8ld %8ld %8ld", vN101, vN000, vN001);
							set_error1 (message);
							return (0);
						}
						ncreated++;
					}
					else {
						eg100 = makeCubeEdge (plex, vN001, vN100);
						if (eg100 == NULL) {
							sprintf (message, "%8ld %8ld", vN001, vN100);
							set_error1 (message);
							return (0);
						}
						ecreated++;
						tglass = makeCubeTriangle (plex, vN001, vN100, vN000);
						if (tglass == NULL) {
							sprintf (message, "%8ld %8ld %8ld", vN001, vN100, vN000);
							set_error1 (message);
							return (0);
						}
						ncreated++;
						tglass = makeCubeTriangle (plex, vN100, vN001, vN101);
						if (tglass == NULL) {
							sprintf (message, "%8ld %8ld %8ld", vN100, vN001, vN101);
							set_error1 (message);
							return (0);
						}
						ncreated++;
					}
				}
		/* YZ squares */
		for (x = 0; x < nx; x++)
			for (y = 0; y < ny-1; y++)
				for (z = 0; z < nz-1; z++) {
					odd = (x + y + z) % 2;
					cube000 = cube100 = cube010 = cube110 = NULL;
					prov000 = prov100 = prov010 = prov110 = NULL;
					pN000 = pN100 = pN010 = pN110 = 0;
					vN000 = vN100 = vN010 = vN110 = 0;
					cube000 = plexcubes + z * ny * nx + y * nx + x;
					cube010 = plexcubes + z * ny * nx + (y+1) * nx + x;
					cube110 = plexcubes + (z+1) * ny * nx + (y+1) * nx + x;
					cube100 = plexcubes + (z+1) * ny * nx + y * nx + x;
					pN000 = cube000 -> provinceNumber;
					pN010 = cube010 -> provinceNumber;
					pN110 = cube110 -> provinceNumber;
					pN100 = cube100 -> provinceNumber;
					nabove = 0; nbelow = 0;
					if (pN000 != 0) {
						prov000 = (provinces + pN000 - 1);
						vN000 = prov000 -> vertexNumber;
						if (cube000 -> type == JustBelow) nbelow++;
						else if (cube000 -> type == JustAbove) nabove++;
					}
					else prov000 = NULL;
					if (pN010 != 0) {
						prov010 = (provinces + pN010 - 1);
						vN010 = prov010 -> vertexNumber;
						if (cube010 -> type == JustBelow) nbelow++;
						else if (cube010 -> type == JustAbove) nabove++;
					}
					else prov010 = NULL;
					if (pN110 != 0) {
						prov110 = (provinces + pN110 - 1);
						vN110 = prov110 -> vertexNumber;
						if (cube110 -> type == JustBelow) nbelow++;
						else if (cube110 -> type == JustAbove) nabove++;
					}
					else prov110 = NULL;
					if (pN100 != 0) {
						prov100 = (provinces + pN100 - 1);
						vN100 = prov100 -> vertexNumber;
						if (cube100 -> type == JustBelow) nbelow++;
						else if (cube100 -> type == JustAbove) nabove++;
					}
					else prov100 = NULL;
					if (nabove + nbelow < 4) continue;
					/* create 1 new edge and 2 triangles */
					if (!odd) {
						eg000 = makeCubeEdge (plex, vN000, vN110);
						if (eg000 == NULL) {
							sprintf (message, "%8ld %8ld", vN000, vN110);
							set_error1 (message);
							return (0);
						}
						ecreated++;
						tglass = makeCubeTriangle (plex, vN000, vN110, vN010);
						if (tglass == NULL) {
							sprintf (message, "%8ld %8ld %8ld", vN000, vN110, vN010);
							set_error1 (message);
							return (0);
						}
						ncreated++;
						tglass = makeCubeTriangle (plex, vN110, vN000, vN100);
						if (tglass == NULL) {
							sprintf (message, "%8ld %8ld %8ld", vN110, vN000, vN100);
							set_error1 (message);
							return (0);
						}
						ncreated++;
					}
					else {
						eg010 = makeCubeEdge (plex, vN100, vN010);
						if (eg010 == NULL) {
							sprintf (message, "%8ld %8ld", vN100, vN010);
							set_error1 (message);
							return (0);
						}
						ecreated++;
						tglass = makeCubeTriangle (plex, vN100, vN010, vN000);
						if (tglass == NULL) {
							sprintf (message, "%8ld %8ld %8ld", vN100, vN010, vN000);
							set_error1 (message);
							return (0);
						}
						ncreated++;
						tglass = makeCubeTriangle (plex, vN010, vN100, vN110);
						if (tglass == NULL) {
							sprintf (message, "%8ld %8ld %8ld", vN010, vN100, vN110);
							set_error1 (message);
							return (0);
						}
						ncreated++;
					}
				}
	}
	if (plex -> dimension == 3) {
		sprintf (message, "%8ld scaffolding edges more", ecreated);
		inform (message);
		sprintf (message, "%8ld scaffolding triangles", ncreated);
		inform (message);
	}
	return (1);
}

int createCubeTetrahedra (struct Plex *plex)
{
	int nabove, nbelow, odd;
	long x, y, z, nx, ny, nz, ncreated;
	long pN000, pN001, pN010, pN011, pN100, pN101, pN110, pN111;
	long vN000, vN001, vN010, vN011, vN100, vN101, vN110, vN111;
	char message[MAXLINE];
	struct PlexCube *plexcubes;
	struct PlexCube *cube000, *cube001, *cube010, *cube011;
	struct PlexCube *cube100, *cube101, *cube110, *cube111;
	struct Province *prov000, *prov001, *prov010, *prov011;
	struct Province *prov100, *prov101, *prov110, *prov111;
	struct Province *provinces;
	struct Glass *tglass;
	
	plexcubes = plex -> plexcubes;
	provinces = plex -> provinces;
	ncreated = 0;
	nx = plex -> dimensions[0];
	ny = plex -> dimensions[1];
	nz = plex -> dimensions[2];
	for (x = 0; x < nx-1; x++)
		for (y = 0; y < ny-1; y++)
			for (z = 0; z < nz-1; z++) {
				odd = (x + y + z) % 2;
				prov000 = prov001 = prov010 = prov011 = 0;
				prov100 = prov101 = prov110 = prov111 = 0;
				vN000 = vN001 = vN010 = vN011 = 0;
				vN100 = vN101 = vN110 = vN111 = 0;
				cube000 = plexcubes + z * ny * nx + y * nx + x;
				cube001 = plexcubes + z * ny * nx + y * nx + (x+1);
				cube010 = plexcubes + z * ny * nx + (y+1) * nx + x;
				cube011 = plexcubes + z * ny * nx + (y+1) * nx + (x+1);
				cube100 = plexcubes + (z+1) * ny * nx + y * nx + x;
				cube101 = plexcubes + (z+1) * ny * nx + y * nx + (x+1);
				cube110 = plexcubes + (z+1) * ny * nx + (y+1) * nx + x;
				cube111 = plexcubes + (z+1) * ny * nx + (y+1) * nx + (x+1);
				pN000 = cube000 -> provinceNumber;
				pN001 = cube001 -> provinceNumber;
				pN010 = cube010 -> provinceNumber;
				pN011 = cube011 -> provinceNumber;
				pN100 = cube100 -> provinceNumber;
				pN101 = cube101 -> provinceNumber;
				pN110 = cube110 -> provinceNumber;
				pN111 = cube111 -> provinceNumber;
				nabove = 0; nbelow = 0;
				if (pN000 != 0) {
					prov000 = (provinces + pN000 - 1);
					vN000 = prov000 -> vertexNumber;
					if (cube000 -> type == JustBelow) nbelow++;
					else if (cube000 -> type == JustAbove) nabove++;
				}
				if (pN001 != 0) {
					prov001 = (provinces + pN001 - 1);
					vN001 = prov001 -> vertexNumber;
					if (cube001 -> type == JustBelow) nbelow++;
					else if (cube001 -> type == JustAbove) nabove++;
				}
				if (pN010 != 0) {
					prov010 = (provinces + pN010 - 1);
					vN010 = prov010 -> vertexNumber;
					if (cube010 -> type == JustBelow) nbelow++;
					else if (cube010 -> type == JustAbove) nabove++;
				}
				if (pN011 != 0) {
					prov011 = (provinces + pN011 - 1);
					vN011 = prov011 -> vertexNumber;
					if (cube011 -> type == JustBelow) nbelow++;
					else if (cube011 -> type == JustAbove) nabove++;
				}
				if (pN100 != 0) {
					prov100 = (provinces + pN100 - 1);
					vN100 = prov100 -> vertexNumber;
					if (cube100 -> type == JustBelow) nbelow++;
					else if (cube100 -> type == JustAbove) nabove++;
				}
				if (pN101 != 0) {
					prov101 = (provinces + pN101 - 1);
					vN101 = prov101 -> vertexNumber;
					if (cube101 -> type == JustBelow) nbelow++;
					else if (cube101 -> type == JustAbove) nabove++;
				}
				if (pN110 != 0) {
					prov110 = (provinces + pN110 - 1);
					vN110 = prov110 -> vertexNumber;
					if (cube110 -> type == JustBelow) nbelow++;
					else if (cube110 -> type == JustAbove) nabove++;
				}
				if (pN111 != 0) {
					prov111 = (provinces + pN111 - 1);
					vN111 = prov111 -> vertexNumber;
					if (cube111 -> type == JustBelow) nbelow++;
					else if (cube111 -> type == JustAbove) nabove++;
				}
				if (nbelow == 0 || nabove == 0) continue;
				if (nabove + nbelow < 8) continue;
				if (odd) {
					/* create 4 new triangles */
					tglass = makeCubeTriangle (plex, vN001, vN010, vN100);
					if (tglass == NULL) {
						sprintf (message, "vN001, vN010, vN100: %8ld %8ld %8ld",
							vN001, vN010, vN100);
						set_error1 (message);
						sprintf (message, "x, y, z: %8ld %8ld %8ld", x, y, z);
						set_error2 (message);
						return (0);
					}
					ncreated++;
					tglass = makeCubeTriangle (plex, vN001, vN010, vN111);
					if (tglass == NULL) {
						sprintf (message, "vN001, vN010, vN111: %8ld %8ld %8ld",
							vN001, vN010, vN111);
						set_error1 (message);
						sprintf (message, "x, y, z: %8ld %8ld %8ld", x, y, z);
						set_error2 (message);
						return (0);
					}
					ncreated++;
					tglass = makeCubeTriangle (plex, vN010, vN100, vN111);
					if (tglass == NULL) {
						sprintf (message, "vN010, vN100, vN111: %8ld %8ld %8ld", vN010, vN100, vN111);
						set_error1 (message);
						sprintf (message, "x, y, z: %8ld %8ld %8ld",
							x, y, z);
						set_error2 (message);
						return (0);
					}
					ncreated++;
					tglass = makeCubeTriangle (plex, vN001, vN100, vN111);
					if (tglass == NULL) {
						sprintf (message, "vN001, vN100, vN111: %8ld %8ld %8ld", vN001, vN100, vN111);
						set_error1 (message);
						sprintf (message, "x, y, z: %8ld %8ld %8ld",
							x, y, z);
						set_error2 (message);
						return (0);
					}
					ncreated++;
					/* create 5 new tetrahedra */
					tglass = makeCubeTetrahedron (plex, vN000, vN001, vN010, vN100);
					if (tglass == NULL) {
						sprintf (message, "%8ld %8ld %8ld %8ld",
							vN000, vN001, vN010, vN100);
						set_error1 (message);
						return (0);
					}
					ncreated++;
					tglass = makeCubeTetrahedron (plex, vN011, vN001, vN010, vN111);
					if (tglass == NULL) {
						sprintf (message, "%8ld %8ld %8ld %8ld",
							vN011, vN001, vN010, vN111);
						set_error1 (message);
						return (0);
					}
					ncreated++;
					tglass = makeCubeTetrahedron (plex, vN110, vN010, vN100, vN111);
					if (tglass == NULL) {
						sprintf (message, "%8ld %8ld %8ld %8ld",
							vN110, vN010, vN100, vN111);
						set_error1 (message);
						return (0);
					}
					ncreated++;
					tglass = makeCubeTetrahedron (plex, vN101, vN001, vN100, vN111);
					if (tglass == NULL) {
						sprintf (message, "%8ld %8ld %8ld %8ld",
							vN101, vN001, vN100, vN111);
						set_error1 (message);
						return (0);
					}
					ncreated++;
					tglass = makeCubeTetrahedron (plex, vN010, vN001, vN100, vN111);
					if (tglass == NULL) {
						sprintf (message, "%8ld %8ld %8ld %8ld",
							vN010, vN001, vN100, vN111);
						set_error1 (message);
						return (0);
					}
					ncreated++;
				}
				else {
					/* create 4 new triangles */
					tglass = makeCubeTriangle (plex, vN110, vN101, vN011);
					if (tglass == NULL) {
						sprintf (message, "vN110, vN101, vN011: %8ld %8ld %8ld", vN110, vN101, vN011);
						set_error1 (message);
						sprintf (message, "x, y, z: %8ld %8ld %8ld", x, y, z);
						set_error2 (message);
						return (0);
					}
					ncreated++;
					tglass = makeCubeTriangle (plex, vN110, vN101, vN000);
					if (tglass == NULL) {
						sprintf (message, "vN110, vN101, vN000: %8ld %8ld %8ld", vN110, vN101, vN000);
						set_error1 (message);
						sprintf (message, "x, y, z: %8ld %8ld %8ld", x, y, z);
						set_error2 (message);
						return (0);
					}
					ncreated++;
					tglass = makeCubeTriangle (plex, vN101, vN011, vN000);
					if (tglass == NULL) {
						sprintf (message, "vN101, vN011, vN000: %8ld %8ld %8ld", vN101, vN011, vN000);
						set_error1 (message);
						sprintf (message, "x, y, z: %8ld %8ld %8ld", x, y, z);
						set_error2 (message);
 						return (0);
					}
					ncreated++;
					tglass = makeCubeTriangle (plex, vN110, vN011, vN000);
					if (tglass == NULL) {
						sprintf (message, "vN110, vN011, vN000: %8ld %8ld %8ld", vN110, vN011, vN000);
						set_error1 (message);
						sprintf (message, "x, y, z: %8ld %8ld %8ld", x, y, z);
						set_error2 (message);
						return (0);
					}
					ncreated++;
					/* create 5 new tetrahedra */
					tglass = makeCubeTetrahedron (plex, vN111, vN110, vN101, vN011);
					if (tglass == NULL) {
						sprintf (message, "%8ld %8ld %8ld %8ld", vN111, vN110, vN101, vN011);
						set_error1 (message);
						return (0);
					}
					ncreated++;
					tglass = makeCubeTetrahedron (plex, vN100, vN110, vN101, vN000);
					if (tglass == NULL) {
						sprintf (message, "%8ld %8ld %8ld %8ld", vN100, vN110, vN101, vN000);
						set_error1 (message);
						return (0);
					}
					ncreated++;
					tglass = makeCubeTetrahedron (plex, vN001, vN101, vN011, vN000);
					if (tglass == NULL) {
						sprintf (message, "%8ld %8ld %8ld %8ld", vN001, vN101, vN011, vN000);
						set_error1 (message);
						return (0);
					}
					ncreated++;
					tglass = makeCubeTetrahedron (plex, vN010, vN110, vN011, vN000);
					if (tglass == NULL) {
						sprintf (message, "%8ld %8ld %8ld %8ld", vN010, vN110, vN011, vN000);
						set_error1 (message);
						return (0);
					}
					ncreated++;
					tglass = makeCubeTetrahedron (plex, vN101, vN110, vN011, vN000);
					if (tglass == NULL) {
						sprintf (message, "%8ld %8ld %8ld %8ld", vN101, vN110, vN011, vN000);
						set_error1 (message);
						return (0);
					}
					ncreated++;
				}
			}
	if (plex -> dimension == 3) {
		sprintf (message, "%8ld scaffolding tetrahedra", ncreated);
		inform (message);
	}
	return (1);
}


struct Glass *makeCubeEdge (struct Plex *plex, long v0, long v1)
{
	int k;
	long idnumbers[4];
	double measure;
	double center[3], center0[3], center1[3];
	struct PlexVertex *plexvtx0, *plexvtx1;
	struct Glass *eglass;

	if (v0 == 0 || v1 == 0) return (0);
	idnumbers[0] = v0;
	idnumbers[1] = v1;
	idnumbers[2] = 0;
	idnumbers[3] = 0;
	eglass = lookupGlass (plex, ProvEdge, idnumbers);
	if (eglass != NULL) return (eglass);
	plexvtx0 = plex -> plexvertices + (v0 - 1);
	plexvtx1 = plex -> plexvertices + (v1 - 1);
	for (k = 0; k < 3; k++)
		center[k] = (plexvtx0 -> center[k] + plexvtx1 -> center[k]) / 2;
	for (k = 0; k < 3; k++) {
		center0[k] = plexvtx0 -> center[k];
		center1[k] = plexvtx1 -> center[k];
	}
	measure = distance (center0, center1);
	eglass = storeGlass (plex, ProvEdge, idnumbers, center, measure);
	if (eglass == NULL) return (NULL);
	return (eglass);
}

struct Glass *makeCubeTriangle (struct Plex *plex, long v0, long v1, long v2)
{
	int k;
	long idnumbers[4], vdns[4], ednumbers[4];
	double measure;
	double center[3], center0[3], center1[3], center2[3];
	char message[MAXLINE];
	struct PlexVertex *plexvtx0, *plexvtx1, *plexvtx2;
	struct Glass *tglass, *eglass, *vglass;

	if (v0 == 0 || v1 == 0 || v2 == 0) return (0);
	idnumbers[0] = v0;
	idnumbers[1] = v1;
	idnumbers[2] = v2;
	idnumbers[3] = 0;
	tglass = lookupGlass (plex, ProvTriangle, idnumbers);
	if (tglass != NULL) return (tglass);
	ednumbers[0] = v0;
	ednumbers[1] = v1;
	ednumbers[2] = 0;
	ednumbers[3] = 0;
	eglass = lookupGlass (plex, ProvEdge, ednumbers);
	if (eglass == NULL) {
		vdns[0] = v0; vdns[1] = 0; vdns[2] = 0;vdns[3] = 0;
		vglass = lookupGlass (plex, ProvVertex, vdns);
		sprintf (message, "makeCubeTriangle: undefined edge vertex %6ld: %8.3f %8.3f %8.3f",
			v0, vglass -> center[0],
				vglass -> center[1],
				vglass -> center[2]);
		set_error1 (message);
		vdns[0] = v1; vdns[1] = 0; vdns[2] = 0;vdns[3] = 0;
		vglass = lookupGlass (plex, ProvVertex, vdns);
		sprintf (message, "makeCubeTriangle: undefined edge vertex %6ld: %8.3f %8.3f %8.3f",
			v1, vglass -> center[0],
				vglass -> center[1],
				vglass -> center[2]);
		set_error2 (message);
		return (NULL);
	}
	ednumbers[0] = v1;
	ednumbers[1] = v2;
	ednumbers[2] = 0;
	ednumbers[3] = 0;
	eglass = lookupGlass (plex, ProvEdge, ednumbers);
	if (eglass == NULL) {
		vdns[0] = v1; vdns[1] = 0; vdns[2] = 0;vdns[3] = 0;
		vglass = lookupGlass (plex, ProvVertex, vdns);
		sprintf (message, "makeCubeTriangle: undefined edge, vertex %6ld: %8.3f %8.3f %8.3f",
			v1, vglass -> center[0],
				vglass -> center[1],
				vglass -> center[2]);
		set_error1 (message);
		vdns[0] = v2; vdns[1] = 0; vdns[2] = 0;vdns[3] = 0;
		vglass = lookupGlass (plex, ProvVertex, vdns);
		sprintf (message, "makeCubeTriangle: undefined edge, vertex %6ld: %8.3f %8.3f %8.3f",
			v2, vglass -> center[0],
				vglass -> center[1],
				vglass -> center[2]);
		set_error2 (message);
		return (NULL);
	}
	ednumbers[0] = v2;
	ednumbers[1] = v0;
	ednumbers[2] = 0;
	ednumbers[3] = 0;
	eglass = lookupGlass (plex, ProvEdge, ednumbers);
	if (eglass == NULL) {
		vdns[0] = v2; vdns[1] = 0; vdns[2] = 0;vdns[3] = 0;
		vglass = lookupGlass (plex, ProvVertex, vdns);
		sprintf (message, "makeCubeTriangle: undefined edge, vertex %6ld: %8.3f %8.3f %8.3f",
			v2, vglass -> center[0],
				vglass -> center[1],
				vglass -> center[2]);
		set_error1 (message);
		vdns[0] = v0; vdns[1] = 0; vdns[2] = 0;vdns[3] = 0;
		vglass = lookupGlass (plex, ProvVertex, vdns);
		sprintf (message, "makeCubeTriangle: undefined edge, vertex %6ld: %8.3f %8.3f %8.3f",
			v0, vglass -> center[0],
				vglass -> center[1],
				vglass -> center[2]);
		set_error2 (message);
		return (NULL);
	}
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
	tglass = storeGlass (plex, ProvTriangle, idnumbers, center, measure);
	if (tglass == NULL) return (NULL);
	return (tglass);
}

struct Glass *makeCubeTetrahedron (struct Plex *plex, long v0, long v1, long v2, long v3)
{
	int k;
	long idnumbers[4];
	double measure;
	double center[3], center0[3], center1[3], center2[3], center3[3];
	struct PlexVertex *plexvtx0, *plexvtx1, *plexvtx2, *plexvtx3;
	struct Glass *tglass;

	if (v0 == 0 || v1 == 0 || v2 == 0 || v3 == 0) return (0);
	idnumbers[0] = v0;
	idnumbers[1] = v1;
	idnumbers[2] = v2;
	idnumbers[3] = v3;
	plexvtx0 = plex -> plexvertices + (v0 - 1);
	plexvtx1 = plex -> plexvertices + (v1 - 1);
	plexvtx2 = plex -> plexvertices + (v2 - 1);
	plexvtx3 = plex -> plexvertices + (v3 - 1);
	for (k = 0; k < 3; k++)
		center[k] = (plexvtx0 -> center[k] +
					plexvtx1 -> center[k] +
					plexvtx2 -> center[k] +
					plexvtx3 -> center[k]) / 4;
	for (k = 0; k < 3; k++) {
		center0[k] = plexvtx0 -> center[k];
		center1[k] = plexvtx1 -> center[k];
		center2[k] = plexvtx2 -> center[k];
		center3[k] = plexvtx3 -> center[k];
	}
	measure = distance (center0, center1) * distance (center1, center2) *
				distance (center2, center3);
	tglass = storeGlass (plex, ProvTetrahedron, idnumbers, center, measure);
	if (tglass == NULL) return (NULL);
	return (tglass);
}

int createContourVertices (struct Plex *plex)
{
	int j, k;
	long h, maxhash, ncreated, vN, vns[2], idns[4];
	double centers[2][3], values[2], fraction, dvalue, center[3];
	double vectors[2][3], vector[3], ctrlev;
	char message[MAXLINE];
	struct Glass *vglass, *eglass;
	struct Plexi *plexi;
	struct PlexVertex *plexvtx;

	ncreated = 0;
	maxhash = plex -> maxhash;
	ctrlev = plex -> ctrlev;
	for (h = 0; h < maxhash; h++) {
		plexi = plex -> plexis + ((int) ProvEdge) * maxhash + h;
		for (eglass = plexi -> head; eglass != NULL; eglass = eglass -> next) {
			for (j = 0; j < 2; j++) {
				vns[j] = eglass -> vertexnumbers[j];
				if (vns[j] == 0) return (0);
				plexvtx = plex -> plexvertices + (vns[j] - 1);
				for (k = 0;k < 3; k++) {
					centers[j][k] = plexvtx -> center[k];
					vectors[j][k] = plexvtx -> vector[k];
				}
				values[j] = plexvtx -> value;
			}
			dvalue = values[1] - values[0];
			if (dvalue == 0.0) continue;
			fraction = (ctrlev - values[0]) / dvalue;
			if (fraction <= 0.0 || fraction >= 1.0) continue;
			if (fraction <= 0.000001) {
				sprintf (message, "createContourVertices: near-zero fraction = %15.12f", fraction);
				inform (message);
			}
			if (fraction >= 0.999999) {
				sprintf (message, "createContourVertices: near-one  fraction = %15.12f", fraction);
				inform (message);
			}
			for (k = 0; k < 3; k++) {
				center[k] = (1.0 - fraction) * centers[0][k] + fraction * centers[1][k];
				vector[k] = (1.0 - fraction) * vectors[0][k] + fraction * vectors[1][k];
			}
			vN = makeContourVertex(plex, center);
			if (vN == 0) {
				sprintf (message, "createContourVertices: fraction = %10.6f", fraction);
				inform (message);
				if (error()) return (0L);
			}
			idns[0] = vN; idns[1] = 0; idns[2] = 0; idns[3] = 0;
			vglass = lookupGlass (plex, PolyVertex, idns);
			if (vglass == NULL) return (0);
			plexvtx = plex -> plexvertices + (vN - 1);
			for (k = 0; k < 3; k++)
				plexvtx -> vector[k] = vector[k];
			/* edge specifies its contour vertex (or 0) */
			eglass -> vertexnumbers[2] = vN;
			ncreated++;
		}
	}
	if (plex -> dimension == 3) {
		sprintf (message, "%8ld contour vertices created", ncreated);
	}
	return (1);
}

int createContourEdges (struct Plex *plex)
{
	int j, k, sgn[3], nc;
	long h, maxhash, ncreated, v0, v1, vns[3], cns[3];
	long idns[4], edns[4];
	char message[MAXLINE];
	struct Glass *glass, *eglass[3], *cglass;
	struct Plexi *plexi;

	ncreated = 0;
	maxhash = plex -> maxhash;
	for (j = 0; j < 4; j++)
		idns[j] = 0;
	for (h = 0; h < maxhash; h++) {
		plexi = plex -> plexis + ((int) ProvTriangle) * maxhash + h;
		for (glass = plexi -> head; glass != NULL; glass = glass -> next) {
			for (j = 0; j < 3; j++) {
				vns[j] = glass -> vertexnumbers[j];
			}
			if (vns[0] == 0 || vns[1] == 0 || vns[2] == 0) {
				sprintf (message, "createContourEdges: zero vertex number %8ld %8ld %8ld", vns[0], vns[1], vns[2]);
				set_error1 (message);
				return (0);
			}
			nc = 0;
			for (j = 0; j < 3; j++) {
				k = (j < 2) ? (j + 1) : 0;
				idns[0] = vns[j];
				idns[1] = vns[k];
				eglass[j] = lookupGlass (plex, ProvEdge, idns);
				if (eglass[j] == NULL) {
					set_error1 ("createContourEdges: cannot find edge");
					return (0);
				}
				edns[0] = eglass[j] -> idnumbers[0];
				edns[1] = eglass[j] -> idnumbers[1];
				edns[2] = eglass[j] -> idnumbers[2];
				edns[3] = eglass[j] -> idnumbers[3];
				sgn[j] = idmatch (ProvEdge, edns, idns);
				if (sgn[j] == 0) {
					set_error1 ("createContourEdges: idmatch inconsistency");
					return (0);
				}
				cns[j] = eglass[j] -> vertexnumbers[2];
				if (cns[j] > 0) nc++;
			}
			if (nc == 0) continue;
			if (nc % 2 == 1) continue;
			/* we have two contour vertices to connect */
			v0 = 0; v1 = 0;
			if (cns[0] != 0 && cns[1] != 0) 	 { v0 = cns[0]; v1 = cns[1]; }
			else if (cns[1] != 0 && cns[2] != 0) { v0 = cns[1]; v1 = cns[2]; }
			else if (cns[2] != 0 && cns[0] != 0) { v0 = cns[2]; v1 = cns[0]; }
			cglass = makeContourEdge (plex, v0, v1);
			if (cglass == NULL) {
				set_error1 ("createContourEdges: no memory for edge");
				return (0);
			}
			glass -> contour = cglass;
			ncreated++;
		}
	}
	if (plex -> dimension == 3) {
		printf (message, "%8ld contour edges created", ncreated);
		inform (message);
	}
	return (1);
}

int createContourTriangles (struct Plex *plex)
{
	int i, j, k, nc, w, nabove, nbelow, okay;
	int try, found, ngot, tries, jsmall;
	int sgn[4], cused[4];
	long nquad, h, maxhash, ncreated, v0, v1, v2, v3;
	long cvns[4], vns[4], idns[4], tdns[4], ucvns[4], jdns[4];
	long triv[4][3] = { 0, 1, 2, 0, 1, 3, 0, 2, 3, 1, 2, 3};  
	double ctrlev, smallest;
	double tvalues[4];
	double tcenters[4][3], ccenters[4][3];
	double cvectors[4][3];
	double ccentroid[3], abovecentroid[3], belowcentroid[3];
	double gradient[3], outgradient[3];
	double cangles[4];
	char message[MAXLINE];
	struct Glass *glass, *cglass, *bglass;
	struct Glass *eglass[4], *tglass[4];
	struct Plexi *plexi;
	struct PlexVertex *plexvtx;

	ncreated = 0; nquad = 0;
	maxhash = plex -> maxhash;
	ctrlev = plex -> ctrlev;
	for (j = 0; j < 4; j++)
		idns[j] = 0;
	for (h = 0; h < maxhash; h++) {
		plexi = plex -> plexis + ((int) ProvTetrahedron) * maxhash + h;
		for (glass = plexi -> head; glass != NULL; glass = glass -> next) {
			for (j = 0; j < 4; j++) {
				vns[j] = glass -> vertexnumbers[j];
				if (vns[j] == 0) {
					set_error1 ("createContourTriangles: vertex number is 0");
					return (0);
				}
			}
			for (j = 0; j < 4; j++) {
				plexvtx = plex -> plexvertices + (vns[j] - 1);
				tvalues[j] = plexvtx -> value;
				for (k = 0; k < 3; k++)
					tcenters[j][k] = plexvtx -> center[k];
			}
			/* compute density gradient */
			for (k = 0; k < 3; k++) {
				abovecentroid[k] = 0.0;
				belowcentroid[k] = 0.0;
			}
			nabove = nbelow = 0;
			for (j = 0; j < 4; j++) {
				if (tvalues[j] > ctrlev) {
					nabove++;
					for (k = 0; k < 3; k++)
						abovecentroid[k] += tcenters[j][k];
				}
				else if (tvalues[j] < ctrlev) {
					nbelow++;
					for (k = 0; k < 3; k++)
						belowcentroid[k] += tcenters[j][k];
				}
			}
			if (nabove == 0 || nbelow == 0) continue;
			for (k = 0; k < 3; k++) {
				abovecentroid[k] /= nabove;
				belowcentroid[k] /= nbelow;
			}
			for (k = 0; k < 3; k++)
				gradient[k] = abovecentroid[k] - belowcentroid[k];
			normalize (gradient);
			for (k = 0; k < 3; k++)
				outgradient[k] = -gradient[k];
			nc = 0;
			for (j = 0; j < 4; j++) {
				for (k = 0; k < 3; k++) {
					w = triv[j][k];
					idns[k] = vns[w];
				}
				tglass[j] = lookupGlass (plex, ProvTriangle, idns);
				if (tglass[j] == NULL) {
					set_error1 ("createContourTriangles: triangle glass not found");
					return (0);
				}
				tdns[0] = tglass[j] -> idnumbers[0];
				tdns[1] = tglass[j] -> idnumbers[1];
				tdns[2] = tglass[j] -> idnumbers[2];
				tdns[3] = tglass[j] -> idnumbers[3];
				sgn[j] = idmatch (ProvTriangle, tdns, idns);
				if (sgn[j] == 0) {
					set_error1 ("createContourTriangles: sign is zero");
					return (0);
				}
				if (tglass[j] -> contour != NULL) {
					if (nc >= 4) {
						set_error1 ("createContourTriangles: too many contour edges for tetrahedron");
						return (0);
					}
					eglass[nc] = tglass[j] -> contour;
					nc++;
				}
			}
			if (nc < 3) continue;
			/* we have three or four contour edges, learn vertices */
			for (j = 0; j < 4; j++) {
				ucvns[j] = 0;
			}
			ngot = 0;
			while (ngot < nc) {
				found = 0;
				for (j = 0; j < nc; j++) {
					try = eglass[j] -> vertexnumbers[0];
					okay = 1;
					for (i = 0; i < ngot; i++)
						if (try == ucvns[i]) okay = 0;
					if (okay) {
						found = 1;
						break;
					}
					try = eglass[j] -> vertexnumbers[1];
					okay = 1;
					for (i = 0; i < ngot; i++)
						if (try == ucvns[i]) okay = 0;
					if (okay) {
						found = 1;
						break;
					}
				}
				if (!found) break;
				ucvns[ngot] = try;
				ngot++;
			}
			if (ngot < nc) {
				set_error1 ("createContourTriangles: repeated vertex");
				return (0);
			}
			for (j = 0; j < nc; j++) {
				plexvtx = plex -> plexvertices + (ucvns[j] - 1);
				for (k = 0;k < 3; k++) {
					ccenters[j][k] = plexvtx -> center[k];
				}
			}
			/* compute centroid of contour vertices */
			for (k = 0; k < 3; k++)
				ccentroid[k] = 0.0;
			for (j = 0; j < nc; j++) {
				for (k = 0; k < 3; k++)
					ccentroid[k] += ccenters[j][k];
			}
			for (k = 0; k < 3; k++)
				ccentroid[k] /= nc;
			/* sort contour vertices */
			for (j = 0; j < nc; j++) {
				for (k = 0; k < 3; k++)
					cvectors[j][k] = ccenters[j][k] - ccentroid[k];
				normalize (cvectors[j]);
			}
			for (j = 0; j < nc; j++) {
				if (j == 0)
					cangles[j] = 0.0;
				else
					cangles[j] = positive_angle (cvectors[0], cvectors[j],  outgradient);
			}
			for (j = 0; j < 4; j++)
				cused[j] = 0;
			ngot = 0;
			for (tries = 0; tries < 4; tries++) {
				smallest = 2 * PI; jsmall = -1;
				for (j = 0; j < nc; j++) {
					if (cused[j]) continue;
					if (cangles[j] < smallest) {
						smallest = cangles[j];
						jsmall = j;
					}
				}
				if (jsmall < 0) {
					set_error1 ("createContourTriangles: sort fails");
					return (0);
				}
				cused[jsmall] = 1;
				cvns[ngot] = ucvns[jsmall];
				ngot++;
				if (ngot >= nc) break;
			}
			if (ngot != nc) {
				set_error1 ("createContourTriangles: sort fails");
				return (0);
			}
			/* recompute for contour vertices now sorted */
			for (j = 0; j < nc; j++) {
				plexvtx = plex -> plexvertices + (cvns[j] - 1);
				for (k = 0;k < 3; k++) {
					ccenters[j][k] = plexvtx -> center[k];
				}
			}
			v0 = cvns[0]; v1 = cvns[1]; v2 = cvns[2]; v3 = cvns[3];
			if (nc == 3) {
				cglass = makeContourTriangle (plex, v0, v1, v2);
				if (cglass == NULL) {
					set_error1 ("createContourTriangles: contour glass is null");
					return (0);
				}
				ncreated++;
			}
			else if (nc == 4) {
				jdns[0] = v0; jdns[1] = v2; jdns[2] = 0; jdns[3] = 0;
				bglass = lookupGlass (plex, PolyEdge, jdns);
				if (bglass != NULL) {
					set_error1 ("createContourTriangles: bisecting edge already in table");
					return (0);
				}
				bglass = makeContourEdge (plex, v0, v2);
				if (bglass == NULL) {
					set_error1 ("createContourTriangles: no memory for edge");
					return (0);
				}
				nquad++;
				cglass = makeContourTriangle (plex, v0, v1, v2);
				if (cglass == NULL) {
					set_error1 ("createContourTriangles: contour glass is null");
					return (0);
				}
				ncreated++;
				cglass = makeContourTriangle (plex, v0, v2, v3);
				if (cglass == NULL) {
					set_error1 ("createContourTriangles: contour glass is null");
					return (0);
				}
				ncreated++;
			}
		}
	}
	if (plex -> dimension == 3) {
		sprintf (message, "%8ld more contour edges created", nquad);
		inform (message);
		sprintf (message, "%8ld contour triangles created", ncreated);
		inform (message);
	}
	return (1);
}
	
long makeContourVertex(struct Plex *plex, double center[3])
{
	int k;
	long vN, v;
	double d, pcenter[3];
	struct PlexVertex *plexvtx;

	/* check for duplicate vertex */
	for (v = 0; v < plex -> n_vertex; v++) {
		plexvtx = plex -> plexvertices + v;
		if (plexvtx -> type != PolyVertex) continue;
		for (k = 0; k < 3; k++)
			pcenter[k] = plexvtx -> center[k];
		if (center[0] == pcenter[0] && center[1] == pcenter[1] && center[2] == pcenter[2]) {
			set_error1 ("makeContourVertex: duplicate coordinates");
			return (0L);
		}
		d = distance (center, pcenter);
		if (d < 0.000001) {
			set_error1 ("makeContourVertex: almost duplicate coordinates");
			return (0L);
		}
	}
	vN = makeVertex(plex, PolyVertex, center);
	return (vN);
}

struct Glass *makeContourEdge (struct Plex *plex, long v0, long v1)
{
	int k;
	long idnumbers[4];
	double measure;
	double center[3], center0[3], center1[3];
	char message[MAXLINE];
	struct PlexVertex *plexvtx0, *plexvtx1;
	struct PlexEdge *plexedg;
	struct Glass *eglass;

	if (v0 == 0 || v1 == 0) return (0);
	idnumbers[0] = v0;
	idnumbers[1] = v1;
	idnumbers[2] = 0;
	idnumbers[3] = 0;
	eglass = lookupGlass (plex, PolyEdge, idnumbers);
	if (eglass != NULL) return (eglass);
	if (plex -> n_edge >= plex -> maxedge) {
		sprintf (message, "makeContourEdge: too many edges (%6ld)",
			plex -> maxedge);
		set_error1 (message);
		return (NULL);
	}
	plexedg = plex -> plexedges + plex -> n_edge;
	plexedg -> type = PolyEdge;
	plexedg -> vns[0] = v0;
	plexedg -> vns[1] = v1;
	plex -> n_edge++;
	plexvtx0 = plex -> plexvertices + (v0 - 1);
	plexvtx1 = plex -> plexvertices + (v1 - 1);
	for (k = 0; k < 3; k++)
		center[k] = (plexvtx0 -> center[k] + plexvtx1 -> center[k]) / 2;
	for (k = 0; k < 3; k++) {
		center0[k] = plexvtx0 -> center[k];
		center1[k] = plexvtx1 -> center[k];
	}
	measure = distance (center0, center1);
	eglass = storeGlass (plex, PolyEdge, idnumbers, center, measure);
	if (eglass == NULL) return (NULL);
	return (eglass);
}

struct Glass *makeContourTriangle (struct Plex *plex, long v0, long v1, long v2)
{
	int k;
	long idnumbers[4];
	double measure;
	double center[3], center0[3], center1[3], center2[3];
	char message[MAXLINE];
	struct PlexVertex *plexvtx0, *plexvtx1, *plexvtx2;
	struct PlexTriangle *plextri;
	struct Glass *tglass;

	if (v0 == 0 || v1 == 0 || v2 == 0) return (0);
	if (plex -> n_triangle >= plex -> maxtriangle) {
		sprintf (message, "makeContourTriangle: too many triangles (%6ld)",
			plex -> maxtriangle);
		return (NULL);
	}
	plextri = plex -> plextriangles + plex -> n_triangle;
	plextri -> type = PolyTriangle;
	plextri -> vns[0] = v0;
	plextri -> vns[1] = v1;
	plextri -> vns[2] = v2;
	plex -> n_triangle++;
	idnumbers[0] = v0;
	idnumbers[1] = v1;
	idnumbers[2] = v2;
	idnumbers[3] = 0;
	tglass = lookupGlass (plex, PolyTriangle, idnumbers);
	if (tglass != NULL) return (tglass);
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
	tglass = storeGlass (plex, PolyTriangle, idnumbers, center, measure);
	if (tglass == NULL) return (NULL);
	return (tglass);
}

int computeBoundary (struct Plex *plex)
{
	int result;
	long h, maxhash, v0, v1, v2, nbdy;
	char message[MAXLINE];
	struct Plexi *plexi;
	struct Glass *glass;
	
	maxhash = plex -> maxhash;
	/* Initialize */
	for (h = 0; h < maxhash; h++) {
		plexi = plex -> plexis + ((int) PolyEdge) * maxhash + h;
		for (glass = plexi -> head; glass != NULL; glass = glass -> next) {
			glass -> coefficient = 0;
		}
	}
	/* go through triangles */
	for (h = 0; h < maxhash; h++) {
		plexi = plex -> plexis + ((int) PolyTriangle) * maxhash + h;
		for (glass = plexi -> head; glass != NULL; glass = glass -> next) {
			v0 = glass -> vertexnumbers[0];
			v1 = glass -> vertexnumbers[1];
			v2 = glass -> vertexnumbers[2];
			result = add_edge (plex, v0, v1); if (!result) return (0);
			result = add_edge (plex, v1, v2); if (!result) return (0);
			result = add_edge (plex, v2, v0); if (!result) return (0);
		}
	}
	nbdy = 0;
	/* print boundary */
	for (h = 0; h < maxhash; h++) {
		plexi = plex -> plexis + ((int) PolyEdge) * maxhash + h;
		for (glass = plexi -> head; glass != NULL; glass = glass -> next) {
			if (glass -> coefficient != 0) {
				v0 = glass -> vertexnumbers[0];
				v1 = glass -> vertexnumbers[1];
				/* non-zero boundary is a sign of a problem only for single-molecule surface */
				if (plex -> ctrlev != 0.0) {
					sprintf (message, "%8ld %8ld vertex numbers, %2d coefficient, boundary edge",
						v0, v1, glass -> coefficient);
					inform (message);
				}
				nbdy++;
			}
		}
	}
	sprintf (message, "%8ld boundary edges for contour level surface", nbdy);
	inform (message);
	return (1);
}

struct msdata *transferContour (struct Plex *plex)
{
	int  k, sgn0, sgn1, sgn2;
	long idx, n;
	long nplexline, count, euler;
	long e0, e1, e2, v0, v1, v2, ue0, ue1, ue2;
	long vcount, ecount, tcount, v, e, t;
	long idns[4], edns[4];
	float *plexlines, *vertices;
	long *edges, *triangles;
	char message[MAXLINE];
	struct PlexVertex *plexvtx;
	struct PlexEdge *plexedg;
	struct PlexTriangle *plextri;
	struct Glass *vglass0, *vglass1, *vglass2;
	struct Glass *eglass0, *eglass1, *eglass2;
	struct msdata *msd;

	/* transfer data to msdata */
	msd = newmsdata ();
	if (msd == NULL) {
		set_error1 ("transferContour: no memory for msd");
		return (NULL);
	}
	count = 0;
	nplexline = plex -> nglass[ProvEdge];
	if (nplexline == 0) {
		set_error1 ("transferContour: no ProvEdge");
		return (NULL);
	}
	plexlines = allocate_floats (nplexline * 6);
	if (plexlines == NULL){
		set_error1 ("transferContour: memory failure for ProvEdge plexlines");
		return (NULL);
	}
	n = transfertype (plex, ProvEdge, plexlines, count);
	if (n != plex -> nglass[ProvEdge]) {
		set_error1 ("transferContour: inconsistent ProvEdge count(1)");
		return (NULL);
	}
	count += n;
	if (count != nplexline) {
		set_error1 ("transferContour: inconsistent ProvEdge count(2)");
		return (NULL);
	}
	if (!addfloats (msd, 6L, nplexline, plexlines)) {
		set_error1 ("transferContour: addfloats fails for ProvEdge plexlines");
		return (NULL);
	}


	vcount = plex -> nglass[PolyVertex];
	ecount = plex -> nglass[PolyEdge];
	tcount = plex -> nglass[PolyTriangle];
	if (vcount == 0 || ecount == 0) return (msd);
	vertices = allocate_floats (vcount * 6);
	if (vertices == NULL) {
		set_error1 ("transferContour: memory failure for PolyVertex");
		return (NULL);
	}
	edges = allocate_longs (ecount * 2, 0, EDGES);
	if (edges == NULL) {
		set_error1 ("transferContour: memory failure for PolyEdge");
		return (NULL);
	}
	if (tcount > 0) {
	triangles = allocate_longs (tcount * 6, 0, TRIANGLES);
		if (triangles == NULL) {
			set_error1 ("transferContour: memory failure for PolyTriangle");
			return (NULL);
		}
	}
	else triangles = NULL;
	idx = 0;
	for (v = 0; v < plex -> n_vertex; v++) {
		plexvtx = plex -> plexvertices + v;
		if (plexvtx -> type != PolyVertex) continue;
		for (k = 0; k < 3; k++) {
			*(vertices + 6 * idx + k) = plexvtx -> center[k];
			*(vertices + 6 * idx + 3 + k) = plexvtx -> vector[k];
		}
		if (plex -> dimension == 2) {
			*(vertices + 6 * idx + 2) = 0.0;
			*(vertices + 6 * idx + 5) = 0.0;
		}
		idx++;
	}
	if (idx != vcount) {
		set_error1 ("transferContour: vertex count inconsistency");
		sprintf (message, "%8ld != %8ld", idx, vcount);
		set_error2 (message);
		return (NULL);
	}
	idx = 0;
	for (e = 0; e < plex -> n_edge; e++) {
		plexedg = plex -> plexedges + e;
		if (plexedg -> type != PolyEdge) continue;
		v0 = plexedg -> vns[0];
		v1 = plexedg -> vns[1];
		idns[0] = v0; idns[1] = 0; idns[2] = 0; idns[3] = 0;
		vglass0 = lookupGlass (plex, PolyVertex, idns);
		idns[0] = v1; idns[1] = 0; idns[2] = 0; idns[3] = 0;
		vglass1 = lookupGlass (plex, PolyVertex, idns);
		*(edges + 2 * idx) = vglass0 -> number;
		*(edges + 2 * idx + 1) = vglass1 -> number;
		idx++;
	}
	if (idx != ecount) {
		set_error1 ("transferContour: edge count inconsistency");
		sprintf (message, "%8ld != %8ld", idx, ecount);
		set_error2 (message);
		return (NULL);
	}
	if (plex -> dimension == 3) {
		idx = 0;
		for (t = 0; t < plex -> n_triangle; t++) {
			plextri = plex -> plextriangles + t;
			if (plextri -> type != PolyTriangle) continue;
			v0 = plextri -> vns[0];
			v1 = plextri -> vns[1];
			v2 = plextri -> vns[2];
			idns[0] = v0; idns[1] = 0; idns[2] = 0; idns[3] = 0;
			vglass0 = lookupGlass (plex, PolyVertex, idns);
			idns[0] = v1; idns[1] = 0; idns[2] = 0; idns[3] = 0;
			vglass1 = lookupGlass (plex, PolyVertex, idns);
			idns[0] = v2; idns[1] = 0; idns[2] = 0; idns[3] = 0;
			vglass2 = lookupGlass (plex, PolyVertex, idns);
			idns[0] = v0; idns[1] = v1; idns[2] = 0; idns[3] = 0;
			eglass0 = lookupGlass (plex, PolyEdge, idns);
			edns[2] = 0; edns[3] = 0;
			edns[0] = eglass0 -> idnumbers[0];
			edns[1] = eglass0 -> idnumbers[1];
			sgn0 = idmatch (PolyEdge, edns, idns);
			ue0 = eglass0 -> number;
			idns[0] = v1; idns[1] = v2; idns[2] = 0; idns[3] = 0;
			eglass1 = lookupGlass (plex, PolyEdge, idns);
			edns[0] = eglass1 -> idnumbers[0];
			edns[1] = eglass1 -> idnumbers[1];
			sgn1 = idmatch (PolyEdge, edns, idns);
			ue1 = eglass1 -> number;
			idns[0] = v2; idns[1] = v0; idns[2] = 0; idns[3] = 0;
			eglass2 = lookupGlass (plex, PolyEdge, idns);
			edns[0] = eglass2 -> idnumbers[0];
			edns[1] = eglass2 -> idnumbers[1];
			sgn2 = idmatch (PolyEdge, edns, idns);
			ue2 = eglass2 -> number;
			e0 = sgn0 * ue0; e1 = sgn1 * ue1; e2 = sgn2 * ue2;
			*(triangles + 6 * idx + 0) = e0;
			*(triangles + 6 * idx + 1) = e1;
			*(triangles + 6 * idx + 2) = e2;
			*(triangles + 6 * idx + 3) = vglass0 -> number;
			*(triangles + 6 * idx + 4) = vglass1 -> number;
			*(triangles + 6 * idx + 5) = vglass2 -> number;
			idx++;
		}
		if (idx != tcount) {
			set_error1 ("transferContour: triangle count inconsistency");
			sprintf (message, "%8ld != %8ld", idx, tcount);
			set_error2 (message);
			return (NULL);
		}
	}
	if (!addfloats (msd, 6L, vcount, vertices)) {
		set_error1 ("transferContour: addfloats fails for PolyVertex");
		return (NULL);
	}
	if (!addlongs (msd, 2L, ecount, edges, EDGES)) {
		set_error1 ("transferContour: addlongs fails for PolyEdge");
		return (NULL);
	}
	if (plex -> dimension == 3) {
		if (!addlongs (msd, 6L, tcount, triangles, TRIANGLES)) {
			set_error1 ("transferContour: addlongs fails for PolyTriangle");
			return (NULL);
		}
	}
	if (plex -> dimension == 3) {
		printf (message, "%8ld polyhedron vertices", vcount);
		inform (message);
		printf (message, "%8ld polyhedron edges", ecount);
		inform (message);
		printf (message, "%8ld polyhedron triangle", tcount);
		inform (message);
		euler = vcount - ecount + tcount;
		inform (message);
		printf (message, "%8ld euler characteristic", euler);
		inform (message);
	}
	return (msd);
}


/* MSP: Piecewise-Linear complEX */
/* Copyright 1995 by Michael L. Connolly */
