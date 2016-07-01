/*
 * MSDraw
 * Copyright 1986 by Michael L. Connolly
 * All Rights Reserved
 * January 7, 2002
 */

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

/* Depth (Z) Buffer */


void init_zbuffer (struct msscene *ms)
{
	double pw2;
    struct cept *ex;

	ms -> horizontal = (ms -> viewport[1][0] - ms -> viewport[0][0]);
	ms -> vertical = (ms -> viewport[1][1] - ms -> viewport[0][1]);

	/* set up z part of viewport */
	ms -> viewport[0][2] = 0;
	ms -> viewport[1][2] = 16384;
	ms -> zview = ms -> viewport[1][2] - ms -> viewport[0][2];
	ms -> size = (unsigned long) (ms -> horizontal) * (unsigned long) (ms -> vertical);
	if (ms -> size == (unsigned) 0) {
        ex = new_cept (MEMORY_ERROR, INVALID_VALUE, FATAL_SEVERITY);
        add_function (ex, "init_zbuffer");
        add_source (ex, "msrender.c");
        add_message (ex, "picture too small");
		return;
	}
	ms -> pixel_width = (ms -> window[1][0] - ms -> window[0][0]) / ms -> horizontal;
	pw2 = (ms -> window[1][1] - ms -> window[0][1]) / ms -> vertical;
	if (fabs (ms -> pixel_width - pw2) > EPSILON) {
        ex = new_cept (PARAMETER_ERROR, INCONSISTENCY, FATAL_SEVERITY);
        add_function (ex, "init_zbuffer");
        add_source (ex, "msrender.c");
        add_message (ex, "different scaling in x and y");
        add_double (ex, "x pixel width", ms -> pixel_width);
        add_double (ex, "y pixel width", pw2);
		return;
	}
	ms -> db = allocate_buffer (ms, 0);
	if (error()) return;
}

struct depth_buffer *allocate_buffer (struct msscene *ms, int partial)
{
	unsigned nelem, elemsize;
	struct depth_buffer *db;
    struct cept *ex;

	db = (struct depth_buffer *) allocate_object (DEPTH_BUFFER);
	if (db == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_function (ex, "allocate_buffer");
		add_source (ex, "msrender.c");
		return (NULL);
	}

	/* NOTE: On the Mac, one must consider
	   the number of bits in an integer */
	/* this will work with 16-bit or 32-bit integers */

	nelem = ms -> vertical;
	elemsize = ms -> horizontal * sizeof (unsigned char);
	db -> hues = allocate_bytes (nelem * elemsize);

	if (db -> hues == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
		add_function (ex, "allocate_buffer");
		add_source (ex, "msrender.c");
        add_message (ex, "hues");
		add_long (ex, "nelem * elemsize", nelem * elemsize);
		return (NULL);
	}

	elemsize = ms -> horizontal * sizeof (unsigned char);
	db -> shades = allocate_bytes (nelem * elemsize);

	if (db -> shades == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
		add_function (ex, "allocate_buffer");
		add_source (ex, "msrender.c");
        add_message (ex, "shades");
		add_long (ex, "nelem * elemsize", nelem * elemsize);
		return (NULL);
	}

	elemsize = ms -> horizontal * sizeof (unsigned char);
	db -> inners = allocate_bytes (nelem * elemsize);

	if (db -> inners == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
		add_function (ex, "allocate_buffer");
		add_source (ex, "msrender.c");
        add_message (ex, "inners");
		add_long (ex, "nelem * elemsize", nelem * elemsize);
		return (NULL);
	}

	if (partial) return (db);

	elemsize = ms -> horizontal * sizeof (short);
	db -> heights = allocate_shorts (nelem * elemsize);

	if (db -> heights == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
		add_function (ex, "allocate_buffer");
		add_source (ex, "msrender.c");
        add_message (ex, "heights");
		add_long (ex, "nelem * elemsize", nelem * elemsize);
		return (NULL);
	}

	elemsize = ms -> horizontal * sizeof (unsigned char);
	db -> alphas = allocate_bytes (nelem * elemsize);

	if (db -> alphas == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
		add_function (ex, "allocate_buffer");
		add_source (ex, "msrender.c");
        add_message (ex, "alphas");
		add_long (ex, "nelem * elemsize", nelem * elemsize);
		return (NULL);
	}

	return (db);
}

/* put pixel into depth buffer */
void putpix (struct msscene *ms, struct surface *current_surface, int inner, int shade, int hue, double opacity, int x, int y, int z)
{
	unsigned long xl, yl, zl;
	long idx;

	if (/* screendoor && */ !opaque (x, y, opacity)) return;
	/* later: do it right */

	/* clip pixel to viewport */
	if (x <  ms -> viewport[0][0]) return;
	if (x >= ms -> viewport[1][0]) return;
	if (y <  ms -> viewport[0][1]) return;
	if (y >= ms -> viewport[1][1]) return;
	if (z <  ms -> viewport[0][2]) return;
	if (z >= ms -> viewport[1][2]) return;
	xl = x;
	yl = y;
	zl = z;

	/* compute index into depth buffer arrays */
	idx = (xl - ms -> viewport[0][0]) + (ms -> viewport[1][0] - ms -> viewport[0][0]) *
		(yl - ms -> viewport[0][1]);

	if (idx < (long) 0) return;
	if (idx >= ms -> size) return;

	if (current_surface -> db != NULL) {
		if (zl > *(current_surface -> db -> heights + idx)) {
			*(current_surface -> db -> heights+idx) = (short) zl;
			*(current_surface -> db -> shades+idx) = (unsigned char) shade;
			*(current_surface -> db -> hues+idx) = (unsigned char) hue;
			*(current_surface -> db -> inners+idx) = (unsigned char) inner;
			if (opacity > 0.0)
				*(current_surface -> db -> alphas+idx) = 255;
		}
	}

	/* throw away if below current occupant */
	if (zl <= *(ms -> db -> heights + idx)) return;

	/* store info */
	*(ms -> db -> heights+idx) = (short) zl;
	*(ms -> db -> shades+idx) = (unsigned char) shade;
	*(ms -> db -> hues+idx) = (unsigned char) hue;
	*(ms -> db -> inners+idx) = (unsigned char) inner;
	if (opacity > 0.0)
		*(ms -> db -> alphas+idx) = 255;	/* later: do it right */
		
	return;
}

/* feature added January 2, 2002 */
void merge_buffers (struct msscene *ms)
{
	int inner, ninner;
	int this_hue, this_shade, this_inner;
	int that_height, that_inner;
	int inner_hue, inner_shade;
	int solid_hue, solid_shade;
	int nsolid, noverlap, same_solid;
	long total_overlap, overlap_removed;
	long i, j, d, e, pidx, qidx, zl;
	double zc, xf, yf, zf;
	double pnt[3];
	char message[MAXLINE];
	struct molecule *mol;
	struct surface *srf;
	struct depth_buffer *db, *mdb, *ndb;

	total_overlap = 0;
	overlap_removed = 0;
	/* allocate temporary buffers */
	mdb = allocate_buffer (ms, 1);
	if (error()) return;
	ndb = allocate_buffer (ms, 1);
	if (error()) return;

	/* copy solid and overlap pixels to the temporary buffer */

	for (i = 0; i < ms -> vertical; i++) {
		for (j = 0; j < ms -> horizontal; j++) {
			/* compute index into depth buffer arrays */
			pidx = i * ms -> horizontal + j;
			if (*(ms -> db -> shades + pidx) <= 0) continue;
			ninner = 0;
			for (mol = ms -> head_molecule; mol != NULL; mol = mol -> next) {
				for (srf = mol -> head_surface; srf != NULL; srf = srf -> next) {
					db = srf -> db;
					if (db == NULL) continue;
					inner = *(db -> inners + pidx);
					if (*(db -> shades + pidx) <= 0) continue;
					if (inner && srf -> solid_shade > 0) {
						ninner++;
						if (ninner > 1) break;
						inner_shade = srf -> solid_shade;
						inner_hue = *(db -> hues + pidx);
					}
				}
				if (ninner > 1) break;
			}
			if (ninner <= 0) continue;
			else if (ninner == 1) {
				*(mdb -> hues + pidx) = inner_hue;
				*(mdb -> shades + pidx) = inner_shade;
				*(mdb -> inners + pidx) = ninner;
			}
			else if (ninner > 1) {
				*(mdb -> hues + pidx) = ms -> overlap_hue;
				*(mdb -> shades + pidx) = 255;
				*(mdb -> inners + pidx) = ninner;
			}
		}
	}

	/* cleanup of glitches */

	for (i = 0; i < ms -> vertical; i++) {
		for (j = 0; j < ms -> horizontal; j++) {
			/* compute index into depth buffer arrays */
			pidx = i * ms -> horizontal + j;
			if (*(mdb -> shades + pidx) <= 0) continue;
			noverlap = 0;
			nsolid = 0;
			solid_hue = 0;
			solid_shade = 0;
			same_solid = 1;
			if (*(mdb -> inners + pidx) > 1) {
				for (d = i - 2; d <= i + 2; d++) {
					if (d < 0 || d > ms -> vertical - 1) continue;
					for (e = j - 2; e <= j + 2; e++) {
						if (e < 0 || e > ms -> horizontal - 1) continue;
						if (d == i && e == j) continue;
						qidx = d * ms -> horizontal + e;
						if (*(mdb -> shades + qidx) <= 0) continue;
						this_inner = *(mdb -> inners + qidx);
						if (this_inner == 1) {
							nsolid++;
							if (solid_hue == 0) {
								solid_hue = *(mdb -> hues + qidx);
								solid_shade = *(mdb -> shades + qidx);
							}
							else if (solid_hue != *(mdb -> hues + qidx)) {
								same_solid = 0;
							}
						}
						else if (this_inner > 1) noverlap++;
					}
				}
			}
			if (nsolid + noverlap == 24 && same_solid && nsolid >= 18) {
				*(ndb -> hues + pidx) = solid_hue;
				*(ndb -> shades + pidx) = solid_shade;
				*(ndb -> inners + pidx) = 1;
				overlap_removed++;
			}
			else {
				*(ndb -> hues + pidx) = *(mdb -> hues + pidx);
				*(ndb -> shades + pidx) = *(mdb -> shades + pidx);
				*(ndb -> inners + pidx) = *(mdb -> inners + pidx);
			}
			if (*(ndb -> inners + pidx) > 1) total_overlap++;
		}
	}

	/* copy solid and overlap pixels to the real buffer */

	for (i = 0; i < ms -> vertical; i++) {
		xf = itof(ms, i, 0);
		for (j = 0; j < ms -> horizontal; j++) {
			yf = itof(ms, j, 1);
			/* compute index into depth buffer arrays */
			pidx = i * ms -> horizontal + j;
			this_hue = *(ndb -> hues + pidx);
			this_shade = *(ndb -> shades + pidx);
			if (this_shade <= 0) continue;
			this_inner = *(ndb -> inners + pidx);
			if (this_inner <= 0) continue;
			that_inner = *(ms -> db -> inners + pidx);
			that_height = *(ms -> db -> heights + pidx);
			pnt[0] = xf; pnt[1] = yf; pnt[2] = 0.0; /* z does not matter */
			zc = zclipped (ms -> clip_center, ms -> clip_axis, pnt);
			zl = ftoi (ms, zc, 2);
			if (that_height > zl && that_inner == 0) continue;
			*(ms -> db -> hues + pidx) = this_hue;
			*(ms -> db -> shades + pidx) = this_shade;
			*(ms -> db -> inners + pidx) = this_inner;
		}
	}
	free_buffer (mdb);
	if (error()) return;
	free_buffer (ndb);
	if (error()) return;
	sprintf (message, "%8ld overlap pixels", total_overlap);
	if (total_overlap > 0) inform (message);
	sprintf (message, "%8ld overlap pixels removed", overlap_removed);
	if (overlap_removed > 0) inform (message);
}

void free_buffers (struct msscene *ms)
{
	struct molecule *mol;
	struct surface *srf;
	struct depth_buffer *db;

	free_buffer (ms -> db);
	for (mol = ms -> head_molecule; mol != NULL; mol = mol -> next) {
		for (srf = mol -> head_surface; srf != NULL; srf = srf -> next) {
			db = srf -> db;
			if (db != NULL) {
				free_buffer (db);
				if (error()) return;
			}
		}
	}
	free_cache (DEPTH_BUFFER);
}

void free_buffer (struct depth_buffer *db)
{
	if (db -> heights != NULL) free_shorts (db -> heights);
	if (db -> shades != NULL) free_bytes (db -> shades);
	if (db -> hues != NULL) free_bytes (db -> hues);
	if (db -> alphas != NULL) free_bytes (db -> alphas);
	if (db -> inners != NULL) free_bytes (db -> inners);
	free_object (DEPTH_BUFFER, (short *) db);
}

/* double to int */

int ftoi (struct msscene *ms, double r, int k)
{
	int i;
	double f;

	f = (r - ms -> window[0][k]) / (ms -> window[1][k] - ms -> window[0][k]);
	i = ms -> viewport[0][k] + f * (ms -> viewport[1][k] - ms -> viewport[0][k]);
	return (i);
}

/* int to double */

double itof (struct msscene *ms, int i, int k)
{
	double f, d;

	f = (double) (i - ms -> viewport[0][k] + 0.5) / (double) (ms -> viewport[1][k] - ms -> viewport[0][k]);
	d = ms -> window[0][k] + f * (ms -> window[1][k] - ms -> window[0][k]);
	return (d);
}

/*
	MSDraw
	Copyright 1986 by Michael L. Connolly
	All Rights Reserved
*/
