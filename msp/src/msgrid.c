#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

/* November 22, 2001 */

/* initialize sphere grid */
struct spheregrid *init_sgrid (long gnspheres, double *gcenters, double *gradii)
{
	int k;
	long s;
	double p1, p2;
	struct spheregrid *sg;

	sg = (struct spheregrid *) allocate_object (SPHEREGRID);
	if (sg == NULL) return (NULL);
	if (gnspheres <= 0) return (NULL);
	if (gradii == (double *) NULL) return (NULL);
	if (gcenters == (double *) NULL) return (NULL);
	sg -> nspheres = gnspheres;
	sg -> centers = allocate_doubles (sg -> nspheres * 3, 0, CENTERS);
	if (sg -> centers == (double *) NULL) return (NULL);
	sg -> radii = allocate_doubles (sg -> nspheres, 0, RADII);
	if (sg -> radii == (double *) NULL) return (NULL);
	for (k = 0; k < 3; k++) {
		sg -> bounds[0][k] = 1000000.0;
		sg -> bounds[1][k] = -1000000.0;
	}
	for (s = 0; s < sg -> nspheres; s++) {
		*(sg -> radii + s) = *(gradii + s);
		for (k = 0; k < 3; k++)
			*(sg -> centers + 3 * s + k) = *(gcenters + 3 * s + k);
	}
	for (s = 0; s < sg -> nspheres; s++) {
		for (k = 0; k < 3; k++) {
			p1 = *(sg -> centers + 3 * s + k) - *(sg -> radii + s);
			if (p1 < sg -> bounds[0][k]) sg -> bounds[0][k] = p1;
			p2 = *(sg -> centers + 3 * s + k) + *(sg -> radii + s);
			if (p2 > sg -> bounds[1][k]) sg -> bounds[1][k] = p2;
		}
	}
	/* enlarge a bit as a safety precaution */
	for (k = 0; k < 3; k++) {
		sg -> bounds[0][k] -= 0.3;
		sg -> bounds[1][k] += 0.3;
		sg -> origin[k] = sg -> bounds[0][k];
	}

	return (sg);
}

/* free sphere grid */
void free_sgrid (struct spheregrid *sg)
{
	free_doubles (sg -> centers, 0, CENTERS);
	free_doubles (sg -> radii, 0, RADII);
	sg -> centers = NULL;
	sg -> radii = NULL;
}


/* cube routines */

/* initialize cube grid */
int init_cgrid (struct spheregrid *sg, double cwidth)
{
	int k, result;
	long s, l;
	char message[MAXLINE];
	
	sg -> cubewidth = cwidth;
	sg -> usnumbers = 0;
	for (k = 0; k < 3; k++) {
		sg -> n_cube[k] = 1 + (sg -> bounds[1][k] - sg -> bounds[0][k])/sg -> cubewidth;
	}
	sg -> ncubexyz = sg -> n_cube[0] * sg -> n_cube[1] * sg -> n_cube[2];
	sg -> scubes = (struct scube *) allocate_objects (SCUBE, sg -> ncubexyz);
	if (sg -> scubes == (struct scube *) NULL)
		return (0);
		
	/* sphere numbers */
	/* first block */
	sg -> first_block = (struct sNumberBlock *) allocate_object (SNUMBERBLOCK);
	if (sg -> first_block == (struct sNumberBlock *) NULL) {
		sprintf (message, "memory allocation for sNumberBlock fails");
		set_error1(message);
		return (0);
	}
	sg -> first_block -> snumbers = (struct snumber *)
		allocate_objects (SNUMBER, (long) SNB_SIZE);
	if (sg -> first_block -> snumbers == (struct snumber *) NULL) {
		sprintf (message, "memory allocation for snumbers fails");
		set_error1(message);
		return (0);
	}
	sg -> nsnumbers = SNB_SIZE;
	/* set up free list */
	for (l = 0; l < sg -> nsnumbers - 1; l++)
		(sg -> first_block -> snumbers + l) -> next =
			sg -> first_block -> snumbers + l + 1;
	sg -> freesn = sg -> first_block -> snumbers;

	sg -> union_numbers = allocate_shorts (sg -> nspheres);
	if (sg -> union_numbers == (short *) NULL) {
		sprintf (message, "memory allocation for union_numbers fails");
		set_error1(message);
		return (0);
	}
	sprintf (message,
		"%8ld by %6ld by %6ld cube grid dimensions; width = %5.2f",
		sg -> n_cube[0],sg -> n_cube[1],sg -> n_cube[2],sg -> cubewidth);
	inform(message);
	
	/* add the spheres to the grid */
	for (s = 1; s <= sg -> nspheres; s++) {
		result = add_sphere_to_grid (sg, s);
		if (!result) {
			sprintf (message,"init_sgrid: cannot add sphere %5ld to grid", s);
			set_error1(message);
			return (0);
		}
	}
	sprintf (message,"%8ld sphere numbers added to grid",
		sg -> usnumbers);
	inform(message);
	return (1);
}

/* add sphere to all cubes it intersects */
int add_sphere_to_grid (struct spheregrid *sg, long snum)
{
	int k, x, y, z, result, res;
	int indices[3], bounds[2][3];
	double scenter[3], sradius, ccenter[3], cwidth;
	char message[MAXLINE];
	
	if (!sphere_to_indices (sg, snum, bounds)) {
		sprintf (message, "sphere_to_indices fails for sphere %5ld", snum);
		set_error1(message);
		return (0);
	}
	cwidth = sg -> cubewidth;
	for (k = 0; k < 3; k++)
		scenter[k] = *(sg -> centers + 3 * (snum - 1) + k);
	sradius = *(sg -> radii + (snum - 1));
	for (z = bounds[0][2]; z <= bounds[1][2]; z++) {
		indices[2] = z;
		for (y = bounds[0][1]; y <= bounds[1][1]; y++) {
			indices[1] = y;
			for (x = bounds[0][0]; x <= bounds[1][0]; x++) {
				indices[0] = x;
				if (!indices_to_center (sg, indices, ccenter)) {
					sprintf (message, "indices_to_center fails");
					set_error1(message);
					return (0);
				}
				result = cube_in_sphere (ccenter, cwidth, scenter, sradius);
				if (result == 1 || result == 2) {
					res = add_sphere_to_cube (sg, snum, indices);
					if (res == 0) {
						sprintf (message,
							"add_sphere_to_grid: add_sphere_to_cube failed");
						set_error1(message);
						return(0);
					}
				}
			}
		}
	}
	return (1);
}

/* add sphere to linked list for this cube */
int add_sphere_to_cube (struct spheregrid *sg, long snum, int indices[3])
{
	struct snumber *sn;
	struct scube *cptr;
	
	cptr = get_scube_ptr (sg, indices);
	if (cptr == NULL) return (0);
	sn = new_snumber (sg);
	if (sn == NULL) return (0);
	sn -> number = snum;
	if (cptr -> first == NULL)
		cptr -> first = sn;
	else cptr -> last -> next = sn;
	cptr -> last = sn;
	cptr -> occupancy = 1;
	return (1);
}

/* determine whether cube lies in sphere, intersects it, or separate */
int cube_in_sphere (double ccenter[3], double cwidth, double scenter[3], double sradius)
{
	double cradius, d;
	int k;
	
	cradius = sqrt (3.0) * cwidth;
	d = distance (ccenter, scenter);
	if (d + cradius < sradius) return (2);
	else if (sradius < d - cradius) return (0);
	for (k = 0; k < 3; k++) {
		if (ccenter[k] + cwidth < scenter[k] - sradius) return (0);
		if (ccenter[k] - cwidth > scenter[k] + sradius) return (0);
	}
	return (1);
}

/* convert integer cube grid indices to angstroms */
int indices_to_center (struct spheregrid *sg, int indices[3], double ccenter[3])
{
	/* indices run from 0 */
	int k;
	
	for (k = 0; k < 3; k++) {
		ccenter[k] = sg -> origin[k] + (indices[k] + 0.5) * sg -> cubewidth;
		if (indices[k] < 0 || indices[k] >= sg -> n_cube[k]) return(0);
	}
	return (1);
}

/* convert angstroms to cube grid indices */
int point_to_indices (struct spheregrid *sg, double pnt[3], int indices[3])
{
	int k;
	double rel[3];
	
	for (k = 0; k < 3; k++) {
		rel[k] = pnt[k] - sg -> origin[k];
		if (rel[k] < 0.0) return (0);
	}
	for (k = 0; k < 3; k++) {
		indices[k] = rel[k] / sg -> cubewidth;
		if (indices[k] < 0) return (0);
		if (indices[k] >= sg -> n_cube[k]) return (0);
	}
	return (1);
}

/* find box that bounds sphere */
int sphere_to_indices (struct spheregrid *sg, long sn, int indices[2][3])
{
	int j, k;
	double *sptr;
	double sradius;
	double extremes[2][3];
	char message[MAXLINE];
    struct cept *ex;
	
	/* sphere numbers run from 1 */
	sptr = sg -> centers + 3 * (sn - 1);
	sradius = *(sg -> radii + (sn - 1));
	for (k = 0; k < 3; k++) {
		extremes[0][k] = *(sptr + k) - sradius;
		extremes[1][k] = *(sptr + k) + sradius;
	}
	for (j = 0; j < 2; j++) {
		if (!point_to_indices (sg, extremes[j], indices[j])) {
			ex = new_cept (GRID_ERROR,  BOUNDS,  FATAL_SEVERITY);
			add_function (ex, "sphere_to_indices");
			add_source (ex, "msgrid.c");
			add_message (ex, "point_to_indices fails");
			sprintf (message, "point = %8.3f %8.3f %8.3f",
				extremes[j][0], extremes[j][1], extremes[j][2]);
            add_message (ex, message);
			return (0);
		}
	}
	return (1);
}

/* return pointer to cube given indices */
struct scube *get_scube_ptr (struct spheregrid *sg, int indices[3])
{
	long nx, ny, nz, ix, iy, iz, idx;
	struct scube *cptr;
	
	nx = sg -> n_cube[0];
	ny = sg -> n_cube[1];
	nz = sg -> n_cube[2];
	ix = indices[0];
	iy = indices[1];
	iz = indices[2];
	idx = iz * nx * ny + iy * nx + ix;
	if (idx < 0 || idx >= nx * ny * nz) return (NULL);
	cptr = sg -> scubes + idx;
	return (cptr);
}

/* take a sphere number struct from the free list */
struct snumber *new_snumber (struct spheregrid *sg)
{
	long l;
	struct snumber *sn;
	struct sNumberBlock *another;
	
	sn = sg -> freesn;
	if (sn == NULL) {
		/* another block */
		another = (struct sNumberBlock *)
			allocate_object (SNUMBERBLOCK);
		if (another == (struct sNumberBlock *) NULL) {
			set_error1 ("memory allocation for sNumberBlock fails");
			return (NULL);
		}
		another -> snumbers = (struct snumber *)
			allocate_objects (SNUMBER, (long) SNB_SIZE);
		if (another -> snumbers == (struct snumber *) NULL) {
			set_error1 ("memory allocation for snumbers fails");
			return (NULL);
		}
		sg -> nsnumbers += SNB_SIZE;
		another -> next = sg -> first_block;
		sg -> first_block = another;
		/* set up free list */
		for (l = 0; l < SNB_SIZE - 1; l++)
			(another -> snumbers + l) -> next =
				another -> snumbers + l + 1;
		sg -> freesn = another -> snumbers;
		sn = sg -> freesn;
	}
	
	sg -> freesn = sn -> next;
	sn -> next = (struct snumber *) NULL;
	sg -> usnumbers++;
	return (sn);
}

/* determine the sphere (atom) numbers lying in the box bounding this sphere */
int sphere_inquiry (struct spheregrid *sg, long sn, int maxsn, long numbers[])
{
	int result;
	int indices[2][3];
    struct cept *ex;
	
	result = sphere_to_indices (sg, sn, indices);
	if (error()) return (0);
	result = box_inquiry (sg, indices, maxsn, numbers);
	if (error()) {
		add_function (tail_cept, "sphere_inquiry");
		add_source (tail_cept, "msgrid.c");
		return(0);
	}
	return (result);
}

/* determine the sphere numbers lying in the given box */
int box_inquiry (struct spheregrid *sg, int indices[2][3], int maxsn, long numbers[])
{
	int i, j, k, count, c;
	long n, lowest, highest;
	int inds[3];
	short *shptr, *shbegin, *shend;
	char message[MAXLINE];
    struct cept *ex;

	lowest = sg -> nspheres + 1;
	highest = 0;
	shbegin = sg -> union_numbers;
	shend = sg -> union_numbers + sg -> nspheres;
	for (shptr = shbegin; shptr < shend; shptr++)
		*shptr = 0;
	for (i = indices[0][0]; i <= indices[1][0]; i++)
		for (j = indices[0][1]; j <= indices[1][1]; j++)
			for (k = indices[0][2]; k <= indices[1][2]; k++) {
				inds[0] = i; inds[1] = j; inds[2] = k;
				count = cube_inquiry (sg, inds, maxsn, numbers);
				if (error()) return (0);
				for (c = 0; c < count; c++) {
					n = numbers[c];
					if (n <= 0 || n > sg -> nspheres) {
						ex = new_cept (ARRAY_ERROR,  BOUNDS,  FATAL_SEVERITY);
						add_object (ex, SPHERE, "sphere numbers");
						add_function (ex, "box_inquiry");
						add_source (ex, "msgrid.c");
						add_long (ex, "sphere number", n);
						return(0);
					}
					*(sg -> union_numbers + n - 1) = 1;
					if (n < lowest) lowest = n;
					if (n > highest) highest = n;
				}
			}
	if (lowest > highest) return (0);
	c = 0;
	/* union of sphere-number lists */
	shbegin = sg -> union_numbers + lowest - 1;
	shend = sg -> union_numbers + highest - 1;
	for (n = lowest, shptr = shbegin; shptr <= shend; shptr++, n++) {
		if (*shptr) {
			numbers[c++] = n;
		}
	}
	return (c);
}

/* return the sphere numbers belonging to this cube in an array */
int cube_inquiry (struct spheregrid *sg, int indices[3], int maxsn, long numbers[])
{
	int n;
	struct scube *cptr;
	struct snumber *sn;
    struct cept *ex;

	cptr = get_scube_ptr (sg, indices);
	if (cptr == (struct scube *) NULL) return (0);
	if (cptr -> occupancy == 0) return (0);
	/* get sphere numbers for this cube from its linked list */
	n = 0;
	for (sn = cptr -> first; sn != (struct snumber *) NULL; sn = sn -> next) {
		if (n >= maxsn) {
			ex = new_cept (ARRAY_ERROR,  MSOVERFLOW,  FATAL_SEVERITY);
			add_object (ex, SPHERE, "numbers");
			add_function (ex, "cube_inquiry");
			add_source (ex, "msgrid.c");
			return (0);
		}
		numbers[n] = sn -> number;
		n++;
	}
	return (n);
}

/* free cube grid memory */
void free_cgrid (struct spheregrid *sg)
{
	long c;
	struct scube *cptr;
	struct sNumberBlock *bptr, *next;

	for (c = 0; c < sg -> ncubexyz; c++) {
		cptr = sg -> scubes + c;
		free_cube (sg, cptr);
	}
	next = NULL;
	for (bptr = sg -> first_block; bptr != NULL; bptr = next) {
		next = bptr -> next;
		free_objects (SNUMBER, (short *) (bptr -> snumbers));
		free_object (SNUMBERBLOCK, (short *) bptr);
	}
	free_shorts (sg -> union_numbers);
	free_objects (SCUBE, (short *) (sg -> scubes));
}

/* return sphere number struct to free list */
void free_snumber (struct spheregrid *sg, struct snumber *sn)
{
	sn -> next = sg -> freesn;
	sg -> freesn = sn;
	sg -> usnumbers--;
	if (sg -> usnumbers < 0) {
		set_error1 ("free_snumber: sphere number underflow");
		return;
	}
}

/* free sphere number linked list for this cube */
void free_cube (struct spheregrid *sg, struct scube *cptr)
{
	struct snumber *sn, *next;
	
	next = NULL;
	for (sn = cptr -> first; sn != NULL; sn = next) {
		next = sn -> next;
		free_snumber (sg, sn);
	}
	cptr -> first = (struct snumber *) NULL;
	cptr -> last = (struct snumber *) NULL;
}

/* bit routines */

/* initialize bit grid */
int init_bgrid (struct spheregrid *sg, double gwidth)
{
	int k, b, bit1, bit2;
	long l, word1, word2;
	char message[MAXLINE];
	
	if (gwidth <= 0.0) return (0);
	sg -> bitwidth = gwidth;
	for (k = 0; k < 3; k++) {
		sg -> nbit[k] = 1 + (sg -> bounds[1][k] - sg -> bounds[0][k])/sg -> bitwidth;
	}
	sprintf (message,
		"%8ld by %6ld by %6ld bit  grid dimensions; width = %5.2f",
		sg -> nbit[0], sg -> nbit[1], sg -> nbit[2], sg -> bitwidth);
	inform(message);
	sg -> nbitxyz = sg -> nbit[0] * sg -> nbit[1] * sg -> nbit[2];
	sg -> bx = sg -> nbit[0];
	sg -> by = sg -> nbit[1];
	sg -> bz = sg -> nbit[2];
	sg -> bxy = sg -> bx * sg -> by;
	sg -> nbitlong = sg -> nbitxyz/32;
	if (sg -> nbitlong % 32 != 0) sg -> nbitlong++;
	/* allocate memory for bitmaps */
	/* two bit grids, because three possibilities (free, full, partial) */
	sg -> bit1_grid = (unsigned long *) allocate_longs (sg -> nbitlong, 0, BIT1_GRID);
	if (sg -> bit1_grid == (unsigned long *) NULL) return (0);
	sg -> bit2_grid = (unsigned long *) allocate_longs (sg -> nbitlong, 0, BIT2_GRID);
	if (sg -> bit2_grid == (unsigned long *) NULL) return (0);
	/* add the spheres to the grid */
	add_spheres (sg);
	sg -> ninterior = 0;
	sg -> nsurface = 0;
	sg -> nexterior = 0;
	for (l = 0; l < sg -> nbitlong; l++) {
		word1 = *(sg -> bit1_grid + l);
		word2 = *(sg -> bit2_grid + l);
		for (b = 0; b < 32; b++) {
			bit1 = word1 & 1;
			word1 >>= 1;
			bit2 = word2 & 1;
			word2 >>= 1;
			if (bit2) sg -> ninterior++;
			else if (bit1) sg -> nsurface++;
			else sg -> nexterior++;
			
		}
	}
	sprintf (message,
		"%8ld interior bits, %8ld surface bits, %8ld exterior bits",
		sg -> ninterior, sg -> nsurface, sg -> nexterior);
	inform(message);
	return (1);
}

/* free bit grid */
void free_bgrid (struct spheregrid *sg)
{
	free_longs ((long *) (sg -> bit1_grid), 0, BIT1_GRID);
	free_longs ((long *) (sg -> bit2_grid), 0, BIT2_GRID);
	sg -> bit1_grid = NULL;
	sg -> bit2_grid = NULL;
}

/* add spheres to bit grid */
void add_spheres (struct spheregrid *sg)
{
	long s;
	
	for (s = 0; s < sg -> nspheres; s++)
		add_sphere(sg, sg -> centers+3*s, *(sg -> radii+s));
}

/* add one sphere to bit grid */
void add_sphere(struct spheregrid *sg, double center[3], double radius)
{
	int k;
	long ix, iy, iz, idx, idx1, idx2;
	double x, y, z, x2, y2, z2, c2, d2, rc, rc2, d;
	double planar_radius, linear_radius, cradius, yr, zr;
	double ztop, zbottom, ytop, ybottom;
	double extremes[2][3], planar[2][3], linear[2][3];
	int ifloor[3], iceiling[3];
	int pfloor[3], pceiling[3];
	int lfloor[3], lceiling[3];
	unsigned long mask;
	unsigned long *bit1_ptr, *bit2_ptr;
	
	c2 = 0.75 * sg -> bitwidth * sg -> bitwidth;
	cradius = sqrt (c2);
	rc = radius + cradius;
	rc2 = rc * rc;
	for (k = 0; k < 3; k++) {
		extremes[0][k] = center[k] - radius;
		extremes[1][k] = center[k] + radius;
	}
	bit_floor (sg, extremes[0], ifloor, 3);
	bit_ceiling (sg, extremes[1], iceiling, 3);
	
	for (iz=ifloor[2];iz<iceiling[2];iz++) {
		z = sg -> bounds[0][2] + (iz+0.5) * sg -> bitwidth;
		zr = z - center[2];
		ztop = zr + 0.5 * sg -> bitwidth;
		zbottom = zr - 0.5 * sg -> bitwidth;
		if (ztop * zbottom <= 0.0) zr = 0.0;
		else {
			ztop = fabs(ztop);
			zbottom = fabs(zbottom);
			if (ztop < zbottom) zr = ztop; else zr = zbottom;
		}
		planar_radius = radius*radius - zr*zr;
		if (planar_radius < 0.0) planar_radius = 0.0;
		planar_radius = sqrt (planar_radius);
		for (k = 0; k < 2; k++) {
			planar[0][k] = center[k] - planar_radius;
			planar[1][k] = center[k] + planar_radius;
		}
		bit_floor (sg, planar[0], pfloor, 2);
		bit_ceiling (sg, planar[1], pceiling, 2);
		for (iy=pfloor[1];iy<pceiling[1];iy++) {
			y = sg -> bounds[0][1] + (iy+0.5) * sg -> bitwidth;
			yr = y - center[1];
			ytop = yr + 0.5 * sg -> bitwidth;
			ybottom = yr - 0.5 * sg -> bitwidth;
			if (ytop * ybottom <= 0.0) yr = 0.0;
			else {
				ytop = fabs(ytop);
				ybottom = fabs(ybottom);
				if (ytop < ybottom) yr = ytop; else yr = ybottom;
			}
			linear_radius = radius*radius - zr*zr - yr*yr;
			if (linear_radius < 0.0) linear_radius = 0.0;
			linear_radius = sqrt (linear_radius);
			linear[0][0] = center[0] - linear_radius;
			linear[1][0] = center[0] + linear_radius;
			bit_floor (sg, linear[0], lfloor, 1);
			bit_ceiling (sg, linear[1], lceiling, 1);
			for (ix=lfloor[0];ix<lceiling[0];ix++) {
				x = sg -> bounds[0][0] + (ix+0.5) * sg -> bitwidth;
				idx = iz * sg -> bxy + iy * sg -> bx + ix;
				idx1 = idx / 32;
				idx2 = idx % 32;
				bit1_ptr = sg -> bit1_grid + idx1;
				bit2_ptr = sg -> bit2_grid + idx1;
				mask = 1 << idx2;
				if (*bit2_ptr & mask) continue;
				x2 = x - center[0]; x2 = x2 * x2;
				y2 = y - center[1]; y2 = y2 * y2;
				z2 = z - center[2]; z2 = z2 * z2;
				d2 = x2 + y2 + z2;
				if (d2 >= rc2) continue;
				*bit1_ptr |= mask;
				d = sqrt (d2);
				if (cradius + d < radius) *bit2_ptr |= mask;
			}
		}
	}
}

/* integer indices below given point */
void bit_floor (struct spheregrid *sg, double pnt[3], int indices[3], int klimit)
{
	int k;
	double rel;
	
	for (k = 0; k < klimit; k++) {
		rel = pnt[k] - sg -> bounds[0][k];
		indices[k] = rel / sg -> bitwidth;
		if (indices[k] < 0) indices[k] = 0;
		if (indices[k] >= sg -> nbit[k]) indices[k] = sg -> nbit[k]-1;
	}
}

/* integer indices above given point */
void bit_ceiling (struct spheregrid *sg, double pnt[3], int indices[3], int klimit)
{
	int k;
	double rel;
	
	for (k = 0; k < klimit; k++) {
		rel = pnt[k] - sg -> bounds[0][k];
		indices[k] = 1 + rel / sg -> bitwidth;
		if (indices[k] < 0) indices[k] = 0;
		if (indices[k] > sg -> nbit[k]) indices[k] = sg -> nbit[k];
	}
}

/* does given point (in angstroms) lie in molecule interior, surface or exterior? */
int point_inquiry (struct spheregrid *sg, double pnt[3])
{
	int indx, indy, indz;
	long idx, idx1, idx2;
	unsigned long bit_mask, *bit1_ptr, *bit2_ptr;

	/* bit_floor (pnt, indices, 3); */
	indx = (pnt[0] - sg -> origin[0]) / sg -> bitwidth;
	if (indx < 0 || indx >= sg -> nbit[0]) return (0);	/* outside grid: empty */
	indy = (pnt[1] - sg -> origin[1]) / sg -> bitwidth;
	if (indy < 0 || indy >= sg -> nbit[1]) return (0);	/* outside grid: empty */
	indz = (pnt[2] - sg -> origin[2]) / sg -> bitwidth;
	if (indz < 0 || indz >= sg -> nbit[2]) return (0);	/* outside grid: empty */
	idx = indz * sg -> bxy + indy * sg -> bx + indx;
	idx1 = idx / 32;
	idx2 = idx % 32;
	bit_mask = 1 << idx2;
	bit1_ptr = sg -> bit1_grid + idx1;
	bit2_ptr = sg -> bit2_grid + idx1;
	if (*bit2_ptr & bit_mask) return (2);		/* inside some sphere */
	else if (*bit1_ptr & bit_mask) return (1);  /* surface (indeterminate) */
	else return (0);							/* outside all spheres */
}

