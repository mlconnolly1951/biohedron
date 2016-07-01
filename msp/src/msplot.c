/*
 * Molecular Surface Package
 * Copyright 1993 by Michael L. Connolly
 * All Rights Reserved

 * February 16, 2006
 */

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

/* CLIPPING */

void do_clipping (struct msscene *ms)
{
	clip_polyhedra (ms);
}

void clip_polyhedra (struct msscene *ms)
{
	int k, n_poly_clipped;
	double clip_center[3], clip_axis[3];
	struct molecule *mol;
	struct surface *phn;
	
	for (k = 0; k < 3; k++) {
		clip_center[k] = ms -> clip_center[k];
		clip_axis[k] = ms -> clip_axis[k];
	}
	n_poly_clipped = 0;
	for (mol = ms -> head_molecule; mol != NULL; mol = mol -> next) {
		for (phn = mol -> head_surface; phn != NULL; phn = phn -> next) {
			if (!phn -> clipping) continue;
			n_poly_clipped++;
			clip_bunch (phn, clip_center, clip_axis);
			if (error()) return;
		}
	}
}

void clip_bunch (struct surface *obj, double clip_center[3], double clip_axis[3])
{
	struct polygon *tri;
	struct phnedg *edg;

	for (edg = obj -> head_phnedg; edg != NULL; edg = edg -> next) {
		clip_edge (obj, edg, clip_center, clip_axis);
		if (error()) return;
	}
	for (tri = obj -> head_polygon; tri != NULL; tri = tri -> next) {
		clip_polygon (obj, tri, clip_center, clip_axis);
		if (error()) return;
	}
	redo_handles (obj);
}


void redo_handles (struct surface *obj)
{
	int j, orn;
	long v, e, p;
	long n_phnvtx, n_phnedg, n_polygon;
	long vnew, enew, pnew;
	char message[MAXLINE];
	struct phnvtx *pv, **vhandles;
	struct phnedg *pe, **ehandles;
	struct polygon *pp, **phandles;

	/* new counts after clipping */
	n_phnvtx = 0; n_phnedg = 0; n_polygon = 0;
	for (pv = obj -> head_phnvtx; pv != NULL; pv = pv -> next)
		n_phnvtx++;
	for (pe = obj -> head_phnedg; pe != NULL; pe = pe -> next)
		n_phnedg++;
	for (pp = obj -> head_polygon; pp != NULL; pp = pp -> next)
		n_polygon++;
	vnew = n_phnvtx - obj -> n_phnvtx;
	enew = n_phnedg - obj -> n_phnedg;
	pnew = n_polygon - obj -> n_polygon;
	if (vnew <= 0 || enew <= 0 || pnew <= 0) {
		/* nothing clipped - be careful not to free same handle array twice */
		return;
	}
	sprintf (message, "%8ld new vertices after clipping", vnew);
	inform (message);
	sprintf (message, "%8ld new edges after clipping", enew);
	inform (message);
	sprintf (message, "%8ld new polygons after clipping", pnew);
	inform (message);
	/* save old phnvtx, phnedg and polygon handle arrays */
	obj -> o_phnvtx = obj -> n_phnvtx;
	obj -> original_phnvtxs = obj -> phnvtx_handles;
	obj -> o_phnedg = obj -> n_phnedg;
	obj -> original_phnedgs = obj -> phnedg_handles;
	obj -> o_polygon = obj -> n_polygon;
	obj -> original_polygons = obj -> polygon_handles;
	/* allocate memory for phnvtx, phnedg, and polygon pointers */
	vhandles = (struct phnvtx **)
		allocate_pointers (PHNVTX, n_phnvtx);
	if (vhandles == NULL) {
		set_error1 ("(redo_handles): memory full");
		return;
	}
	ehandles = (struct phnedg **)
		allocate_pointers (PHNEDG, n_phnedg);
	if (ehandles == NULL) {
		set_error1 ("(redo_handles): memory full");
		return;
	}
	phandles = (struct polygon **)
		allocate_pointers (POLYGON, n_polygon);
	if (phandles == NULL) {
		set_error1 ("(redo_handles): memory full");
		return;
	}
	/* from old handle array */
	for (v = 0; v < obj -> o_phnvtx; v++) {
		pv = *(obj -> original_phnvtxs + v);
		*(vhandles + v) = pv;
		pv -> number = v + 1;
	}
	/* new vertices were inserted at head of linked list */
	for (v = obj -> o_phnvtx, pv = obj -> head_phnvtx;
		v < n_phnvtx; v++, pv = pv -> next) {
		*(vhandles + v) = pv;
		pv -> number = v + 1;
	}
	/* from old handle array */
	for (e = 0; e < obj -> o_phnedg; e++) {
		pe = *(obj -> original_phnedgs + e);
		*(ehandles + e) = pe;
		pe -> number = e + 1;
	}
	/* new edges were inserted at head of linked list */
	for (e = obj -> o_phnedg, pe = obj -> head_phnedg;
		e < n_phnedg; e++, pe = pe -> next) {
		*(ehandles + e) = pe;
		pe -> number = e + 1;
	}
	/* from old handle array */
	for (p = 0; p < obj -> o_polygon; p++) {
		pp = *(obj -> original_polygons + p);
		*(phandles + p) = pp;
	}
	/* new polygons were inserted at head of linked list */
	for (p = obj -> o_polygon, pp = obj -> head_polygon;
		p < n_polygon; p++, pp = pp -> next) {
		*(phandles + p) = pp;
	}
	/* set up edge vertex indices */
	for (e = 0; e < n_phnedg; e++) {
		pe = *(ehandles + e);
		for (j = 0; j < 2; j++) {
			pe -> vtxnum[j] = pe -> pvt[j] -> number;
		}
	}
	/* set up polygon vertex indices */
	for (p = 0; p < n_polygon; p++) {
		pp = *(phandles + p);
		for (j = 0; j < pp -> n_side; j++) {
			pe = pp -> edg[j]; orn = pp -> orn[j];
			pv = pe -> pvt[orn];
			pp -> vertex_index[j] = pv -> number - 1;
		}
	}
	/* store local variables in objects */
	obj -> n_phnvtx = n_phnvtx;
	obj -> n_phnedg = n_phnedg;
	obj -> n_polygon = n_polygon;
	obj -> phnvtx_handles = vhandles;
	obj -> phnedg_handles = ehandles;
	obj -> polygon_handles = phandles;
}


void clip_edge (struct surface *pb, struct phnedg *edg, double clip_center[3], double clip_axis[3])
{
	int j, k, clipped0, clipped1;
	double zdiff, zclip, factor;
	double vector0[3], vector1[3];
	struct phnvtx *vtx0, *vtx1, *vtx2;
	struct phnedg *edg0, *edg1;

	vtx0 = edg -> pvt[0];
	vtx1 = edg -> pvt[1];
	clipped0 = pclipped (clip_center, clip_axis, vtx0 -> center);
	clipped1 = pclipped (clip_center, clip_axis, vtx1 -> center);
	if (!clipped0 && !clipped1) return;
	if (clipped0 && clipped1) {
		edg -> clipped = 2;	/* completely clipped */
		return;
	}
	/* edge crossed clipping plane */
	
	/* if normal vector, mark it completely clipped anyway (looks better) */
	if (pb -> type == NML_SURFACE) {
		edg -> clipped = 2;
		return;
	}
	
	/* compute vertex coordinates */
	for (k = 0; k < 3; k++)
		vector0[k] = vtx1 -> center[k] - vtx0 -> center[k];
	zdiff = dot_product (vector0, clip_axis);
	if (fabs (zdiff) == 0.0) { /* edge lies in clipping plane */
		return;
	}
	for (k = 0; k < 3; k++)
		vector1[k] = clip_center[k] - vtx0 -> center[k];
	zclip = dot_product (vector1, clip_axis);
	factor = zclip / zdiff;
	
	edg -> clipped = 1;		/* partially clipped */
	vtx2 = (struct phnvtx *) allocate_phnvtx ();
	if (vtx2 == (struct phnvtx *) NULL) return;
	/* link into head of list to prevent recursion */
	vtx2 -> next = pb -> head_phnvtx;
	pb -> head_phnvtx = vtx2;
	for (k = 0; k < 3; k++) {
		vtx2 -> center[k] = vtx0 -> center[k] + 
			factor * (vtx1 -> center[k] - vtx0 -> center[k]);
		vtx2 -> outward[k] = vtx0 -> outward[k] + 
			factor * (vtx1 -> outward[k] - vtx0 -> outward[k]);
		vtx2 -> values[k] = vtx0 -> values[k] + 
			factor * (vtx1 -> values[k] - vtx0 -> values[k]);
	}
	normalize (vtx2 -> outward);

	edg0 = (struct phnedg *) allocate_phnedg ();
	if (edg0 == (struct phnedg *) NULL) return;
	edg1 = (struct phnedg *) allocate_phnedg ();
	if (edg1 == (struct phnedg *) NULL) return;
	/* link into head of list to prevent recursion */
	edg1 -> next = pb -> head_phnedg;
	edg0 -> next = edg1;
	pb -> head_phnedg = edg0;
	
	edg0 -> hue = edg -> hue; edg1 -> hue = edg -> hue;
	edg0 -> type = edg -> type; edg1 -> type = edg -> type;
	edg0 -> comp = edg -> comp; edg1 -> comp = edg -> comp;
	edg0 -> atm = edg -> atm; edg1 -> atm = edg -> atm;
	for (j = 0; j < 2; j++) {
		edg0 -> on[j] = edg -> on[j];
		edg1 -> on[j] = edg -> on[j];
	}
	edg0 -> pvt[0] = vtx0; edg0 -> pvt[1] = vtx2;
	edg1 -> pvt[0] = vtx2; edg1 -> pvt[1] = vtx1;
	if (vtx0 -> center[2] > vtx1 -> center[2]) {
		edg0 -> clipped = 1; edg1 -> clipped = 0;
	}
	else {
		edg0 -> clipped = 0; edg1 -> clipped = 1;
	}
	edg -> split[0] = edg0;
	edg -> split[1] = edg1;
}

void clip_polygon (struct surface *pb, struct polygon *tri, double clip_center[3], double clip_axis[3])
{
	int j, n;
	char message[MAXLINE];
	
	n = 0;
	for (j = 0; j < 3; j++)
		if (tri -> edg[j] -> clipped == 1) n++;
	if (n % 2 != 0) {
		sprintf (message, "(clip_polygon) odd number of clipped edges: %d",n);
		set_error1 (message);
		return;
	}
	if (n <= 0) {	/* triangle all under or all over clipping plane */
		tri -> clipped = pclipped (clip_center, clip_axis, tri -> center);
		return;
	}
	tri -> clipped = 1;
	if (tri -> edg[0] -> clipped == 1 && tri -> edg[1] -> clipped == 1)
		clip_tri (pb, tri, 0, 1, 2, clip_center, clip_axis);
	else if (tri -> edg[1] -> clipped == 1 && tri -> edg[2] -> clipped == 1)
		clip_tri (pb, tri, 1, 2, 0, clip_center, clip_axis);
	else if (tri -> edg[2] -> clipped == 1 && tri -> edg[0] -> clipped == 1)
		clip_tri (pb, tri, 2, 0, 1, clip_center, clip_axis);
}

void clip_tri (struct surface *pb, struct polygon *tri, int e0, int e1, int e2, double clip_center[3], double clip_axis[3])
{
	int orn0, orn1, orn2;
	struct phnedg *edg;
	struct polygon *pgn0, *pgn1;

	edg = (struct phnedg *) allocate_phnedg ();
	if (edg == (struct phnedg *) NULL) return;
	/* link into head of list to prevent recursion */
	edg -> next = pb -> head_phnedg;
	pb -> head_phnedg = edg;

	edg -> hue = tri -> hue; edg -> comp = tri -> comp; edg -> atm = tri -> atm;
	edg -> type = BOUNDARY; edg -> clipped = 0;
	pgn0 = allocate_polygon();
	if (pgn0 == (struct polygon *) NULL) return;
	pgn0 -> n_side = 3;
	pgn0 -> hue = tri -> hue;
	pgn0 -> comp = tri -> comp;
	pgn0 -> atm = tri -> atm;
	pgn1 = allocate_polygon();
	if (pgn1 == (struct polygon *) NULL) return;
	pgn1 -> n_side = 4;
	pgn1 -> hue = tri -> hue;
	pgn1 -> comp = tri -> comp;
	pgn1 -> atm = tri -> atm;
	pgn1 -> next = pb -> head_polygon;
	pgn0 -> next = pgn1;
	pb -> head_polygon = pgn0;
	edg -> on[0] = pgn0; edg -> on[1] = pgn1;
	orn0 = tri -> orn[e0];
	orn1 = tri -> orn[e1];
	orn2 = tri -> orn[e2];
	if (tri -> edg[e0] -> pvt[1-orn0] != tri -> edg[e1] -> pvt[orn1]) {
		set_error1("(clip_tri): inconsistent edge orientations");
		return;
	}
	
	pgn0 -> edg[0] = tri -> edg[e0] -> split[1-orn0];
	pgn0 -> orn[0] = (short) orn0;
	pgn0 -> edg[0] -> on[orn0] = pgn0;
	pgn0 -> edg[1] = tri -> edg[e1] -> split[orn1];
	pgn0 -> orn[1] = (short) orn1;
	pgn0 -> edg[1] -> on[orn1] = pgn0;
	edg -> pvt[0] = tri -> edg[e1] -> split[orn1] -> pvt[1-orn1];
	edg -> pvt[1] = tri -> edg[e0] -> split[1-orn0] -> pvt[orn0];
	pgn0 -> edg[2] = edg;
	pgn0 -> orn[2] = (short) 0;
	edg -> on[0] = pgn0;
	pgn1 -> edg[0] = tri -> edg[e0] -> split[orn0];
	pgn1 -> orn[0] = (short) orn0;
	pgn1 -> edg[0] -> on[orn0] = pgn1;
	pgn1 -> edg[1] = edg;
	pgn1 -> orn[1] = (short) 1;
	edg -> on[1] = pgn1;
	pgn1 -> edg[2] = tri -> edg[e1] -> split[1-orn1];
	pgn1 -> orn[2] = (short) orn1;
	pgn1 -> edg[2] -> on[orn1] = pgn1;
	pgn1 -> edg[3] = tri -> edg[e2];
	pgn1 -> orn[3] = (short) orn2;
	pgn1 -> edg[3] -> on[orn2] = pgn1;
	
	
	
	/* set up other fields of polygons */
	do_pgn_axis (pgn0);
	if (error()) return;
	do_pgn_bound (pgn0);
	if (error()) return;
	pgn0 -> clipped = pclipped(clip_center, clip_axis, pgn0 -> center);
	do_pgn_axis (pgn1);
	if (error()) return;
	do_pgn_bound (pgn1);
	if (error()) return;
	pgn1 -> clipped = pclipped(clip_center, clip_axis, pgn1 -> center);
	if (pgn0 -> clipped == pgn1 -> clipped) {
		set_error1 ("(clip_tri): inconsistent triangle clipping");
		return;
	}
}


/* Hidden-line Elimination Section */

/* Squares */

void setup_squares (struct msscene *ms)
{
	long i, j, n;
	double xmin, xmax, ymin, ymax, xsquare, ysquare;
	long averagen;
	char message[MAXLINE];
	struct square *sq;

	/* current polygons has its own space during this routine */
	ms -> current_npolygon = 0;
	ms -> current_polygons = (struct polygon **)
		allocate_pointers (POLYGON, MAX_TIS);
	if (ms ->current_polygons == (struct polygon **) NULL) {
		set_error1 ("msdraw: cannot allocate memory for polygon pointers");
		return;
	}
	
	ms -> n_square = ms -> n_row * ms -> n_column;
	if (ms -> n_square < 1 || ms -> n_square > 1024) {
		sprintf (message, "msdraw: nsquare = %ld", ms -> n_square);
		set_error1(message); return;
	}
	ms -> squares = (struct square *) allocate_objects (SQUARE, ms -> n_square);
	if (ms -> squares == (struct square *) NULL) {
		set_error1 ("msdraw: no memory for squares");
		return;
	}
	xsquare = ms -> plot_width / ms -> n_column;
	ysquare = ms -> plot_width / ms -> n_row;
	averagen = 0;
	for (i = 0; i < ms -> n_row; i++) {
		ymin = ms -> window[0][1] + i * ysquare;
		ymax = ms -> window[0][1] + (i + 1) * ysquare;
		for (j = 0; j < ms -> n_column; j++) {
			xmin = ms -> window[0][0] + j * xsquare;
			xmax = ms -> window[0][0] + (j + 1) * xsquare;
			sq = ms -> squares + (i * ms -> n_column + j);
			sq -> xmin = xmin; sq -> xmax = xmax;
			sq -> ymin = ymin; sq -> ymax = ymax;
			n = setup_square (ms, sq);
			averagen += n;
			if (error()) return;
			if (n == 0) continue;
		}
	}
	averagen /= (ms -> n_row * ms -> n_column);
	sprintf (message, "%8ld triangles in square on average", averagen);
	inform (message);
	
	/* initialize to no current square or polygons */
	ms -> current_npolygon = 0;
	if (ms -> current_polygons != (struct polygon **) NULL)
		free_pointers (POLYGON, ms -> current_polygons);
	ms -> current_polygons = (struct polygon **) NULL;
	ms -> current_square = (struct square *) NULL;
}

void free_squares (struct msscene *ms)
{
	long i, j, n;
	char message[MAXLINE];
	struct square *sq;
	struct linseg *ls,*nextls;
	struct molecule *mol;
	struct surface *phn;

	for (mol = ms -> head_molecule; mol != NULL; mol = mol -> next) {
		for (phn = mol -> head_surface; phn != NULL; phn = phn -> next) {
			nextls = NULL;
			for (ls = phn -> head_linseg; ls != NULL; ls = nextls) {
				nextls = ls -> next;
				free_linseg(ls);
			}
			phn -> head_linseg = NULL;
			phn -> tail_linseg = NULL;
		}
	}

	for (i = 0; i < ms -> n_row; i++) {
		for (j = 0; j < ms -> n_column; j++) {
			sq = ms -> squares + (i * ms -> n_column + j);
			if (sq -> polygons != (struct polygon **) NULL)
				free_pointers (POLYGON, sq -> polygons);
			sq -> polygons = (struct polygon **) NULL;
		}
	}
	free_objects (SQUARE, (short *) (ms -> squares));
	ms -> squares = NULL;
	free_cache (LINSEG);
	free_cache (SQUARE);
}

int setup_square (struct msscene *ms, struct square *sq)
{
	int comp, atm, inner, shape;
	double opacity;
	double xmin, xmax, ymin, ymax;
	struct polygon *tri;
	struct molecule *mol;
	struct surface *po;
	
	xmin = sq -> xmin;
	xmax = sq -> xmax;
	ymin = sq -> ymin;
	ymax = sq -> ymax;
	
	ms -> current_npolygon = 0;

	
	for (mol = ms -> head_molecule; mol != NULL; mol = mol -> next) {
		po = polyhedron_bunch (mol);
		if (po == NULL) {
			informd ("polyhedral bunch not found");
			continue;
		}
		if (po -> type != PHN_SURFACE) {
			informd ("bunch not polyhedral");
			continue;
		}
		for (tri = po -> head_polygon; tri != NULL; tri = tri -> next) {
			if (tri -> clipped) continue;
			comp = tri -> comp;
			atm = tri -> atm;
			inner = 0;
			shape = CONVEX;
			opacity = detopac (atm, inner, comp, shape, po -> scheme, (double) 1.0);
			if (error()) return (0);
			if (opacity < 0.5) continue;
			if (ms -> current_npolygon >= MAX_TIS) {
				set_error1 ("msdraw: MAX_TIS is too small");
				return (0);
			}
			if (tri -> xmax < xmin) continue;
			if (tri -> xmin > xmax) continue;
			if (tri -> ymax < ymin) continue;
			if (tri -> ymin > ymax) continue;
			*(ms -> current_polygons + ms -> current_npolygon++) = tri;
		}
	}
	
	store_square (ms, sq);
	if (error()) return (0);
	
	return (ms -> current_npolygon);
	
}

void store_square (struct msscene *ms, struct square *sq)
{
	long i;
	struct polygon **polygons;
	
	sq -> npolygon = ms -> current_npolygon;
	if (sq -> npolygon <= 0) {
		sq -> polygons = (struct polygon **) NULL;
		return;
	}
	/* allocate memory for this squares polygon pointer list */
	polygons = (struct polygon **) allocate_pointers (POLYGON, (unsigned long) (ms -> current_npolygon));
	if (polygons == (struct polygon **) NULL) {
		set_error1 ("msdraw: cannot allocate memory for squares polygon pointers");
		return;
	}
	sq -> polygons = polygons;
	/* transfer data from global array to square */
	for (i = 0; i < ms -> current_npolygon; i++)
		*(polygons + i) = *(ms -> current_polygons + i);
}


int fetch_square (struct msscene *ms, double pnt[3])
{
	long i, j;
	double x, y, fj, fi;
	struct square *sq;
	
	x = pnt[0] - ms -> window[0][0];
	y = pnt[1] - ms -> window[0][1];
	fj = x / (ms -> plot_width / ms -> n_column);
	fi = y / (ms -> plot_width / ms -> n_row);
	j = floor (fj);
	i = floor (fi);
	if (j < 0 || j > ms -> n_column - 1 || i < 0 || i > ms -> n_row - 1) {
		return (0);
	}
	sq = ms -> squares + (i * ms -> n_column + j);
	if (sq == ms -> current_square) return (1);
	ms -> current_square = sq;
	ms -> current_npolygon = sq -> npolygon;
	ms -> current_polygons = sq -> polygons;
	return (1);
}

int under_square (struct msscene *ms, double pnt[3], struct polygon *tri0, struct polygon *tri1)
{
	int result;
	long i;
	struct polygon *tri;
	
	/* check all polyhedron polygons intersecting this square */
	for (i = 0; i < ms -> current_npolygon; i++) {
		tri = *(ms -> current_polygons + i);
		if (tri == (struct polygon *) NULL) {
			set_error1 ("msdraw: null polygon pointer");
			return(0);
		}
		/* skip polygons bordering edge */
		if (tri == tri0) continue;
		if (tri == tri1) continue;
		/* skip backward-facing polygons (unless clipping) */
		if (! ms -> clipping && ! ms -> translucency && tri -> axis[2] <= 0.0) continue;
		result = put (pnt, tri);
		if (error()) return(0);
		if (result) return (1);
	}
	return (0);
}


void create_linsegs (struct msscene *ms)
{
	struct surface *head_surface, *b;
	for (ms -> current_molecule = ms -> head_molecule;
		ms -> current_molecule != NULL;
		ms -> current_molecule = ms -> current_molecule -> next) {
		head_surface = ms -> current_molecule -> head_surface;
		for (b = head_surface; b != NULL; b = b -> next) {
			if (b -> type == PQMS_SURFACE) continue;
			if (ms -> current_molecule -> blank && b -> type == PHN_SURFACE) continue;
			bunch_linsegs (b);
			if (error()) return;
		}
	}
}

/* create linseg structs for each edge of molecule or contour */

void bunch_linsegs (struct surface *obj)
{
	int result;
	long nls;
	char message[MAXLINE], typename[25];
	struct polygon *tri0, *tri1;
	struct phnedg *edg;

	obj -> head_linseg = (struct linseg *) NULL;
	obj -> tail_linseg = (struct linseg *) NULL;
	
	nls = 0;
	/* convert edges to line segments */
	for (edg = obj -> head_phnedg; edg != NULL; edg = edg -> next) {
		if (edg -> clipped) continue;
		tri0 = edg -> on[0];
		tri1 = edg -> on[1];
		if (tri0 != NULL && tri1 != NULL) {
			if (tri0 -> clipped && tri1 -> clipped) continue;
		}
		result = edge_to_linseg (obj, edg);
		if (error()) return;
		if (result) nls++;
	}
	get_bunch_name (obj -> type, typename);
	sprintf (message, "%8ld %s edges converted to line segments",
		nls, typename);
	inform(message);
}

int edge_to_linseg (struct surface *bnc, struct phnedg *edg)
{
	int j, k, hue, comp, atm, inner, shape;
	double ends[2][2], height[2];
	double values[2], value, opacity;
	struct phnvtx *vtx;
	
	for (j = 0; j < 2; j++) {
		vtx = edg -> pvt[j];
		for (k = 0; k < 2; k++)
			ends[j][k] = vtx -> center[k];
		height[j] = vtx -> center[2];
		values[j] = get_function (bnc, vtx);
		if (error()) return (0);
	}
	value = (values[0] + values[1]) / 2;
	comp = edg -> comp;
	atm = edg -> atm;
	hue = edg -> hue;
	inner = 0;
	shape = CONVEX;
	opacity = detopac (atm, inner, comp, shape, bnc -> scheme, (double) 1.0);
	if (opacity < 0.5) return (0);
	if (error()) return (0);
	make_linseg (bnc, ends, height, comp, atm, hue, value, edg);
	if (error()) return (0);
	return (1);
}

void make_linseg (struct surface *bnc, double pvts[2][2], double height[2], int comp, int atm, int hue, double value, struct phnedg *edg)
{
	int j, k;
	struct linseg *ls;
	
	ls = allocate_linseg ();
	if (ls == (struct linseg *) NULL) {
		set_error1 ("msdraw: no memory for linseg");
		return;
	}
	/* link into list */
	if (bnc -> head_linseg == (struct linseg *) NULL)
		bnc -> head_linseg = ls;
	else bnc -> tail_linseg -> next = ls;
	bnc -> tail_linseg = ls;
	ls -> next = (struct linseg *) NULL;
	
	/* set up fields */
	ls -> comp = (unsigned short) comp;
	ls -> atm = (unsigned short) atm;
	ls -> hue = (short) hue;
	ls -> value = value;
	ls -> edg = edg;		/* back pointer */
	ls -> hidden = 0;
	for (j = 0; j < 2; j++) {
		for (k = 0; k < 2; k++)
			ls -> ends[j][k] = pvts[j][k];
		ls -> height[j] = height[j];
	}
}

long count_linseg (struct linseg *first)
{
	long n;
	struct linseg *ls;
	
	n = 0;
	for (ls = first; ls != (struct linseg *) NULL; ls = ls -> next)
		n++;
	return (n);
}


void mark_molecules (struct msscene *ms)
{
	struct molecule *mol;
	struct surface *b;

	for (mol = ms -> head_molecule; mol != NULL; mol = mol -> next) {
		for (b = mol -> head_surface; b != NULL; b = b -> next) {
			if (b -> type == PQMS_SURFACE) continue;
			mark_bunch (b);
			if (error()) return;
		}
	}
}

void mark_bunch (struct surface *b)
{
	long n_boundary, n_silhouette;
	long n_front, n_back;
	char message[MAXLINE];
	struct phnedg *edg;
	struct polygon *tri0, *tri1;
	
	n_boundary = 0;
	n_silhouette = 0;
	n_back = 0;
	n_front = 0;
	for (edg = b -> head_phnedg; edg != NULL; edg = edg -> next) {
		tri0 = edg -> on[0];
		tri1 = edg -> on[1];
		switch (b -> type) {
		case PHN_SURFACE:
			if (tri0 == (struct polygon *) NULL &&
				tri1 == (struct polygon *) NULL) {
				informd("warning: polyhedron edge lying on no polygon");
				edg -> type = BOUNDARY;
				n_boundary++;
				continue;
			}
			if (tri0 == (struct polygon *) NULL) {
				edg -> type = BOUNDARY;
				n_boundary++;
				continue;
			}
			if (tri1 == (struct polygon *) NULL) {
				edg -> type = BOUNDARY;
				n_boundary++;
				continue;
			}
			if (tri0 -> clipped != tri1 -> clipped) {
				edg -> type = BOUNDARY;
				n_boundary++;
				continue;
			}
			/*  backward-facing polygon */
			if (tri0 -> axis[2] <= 0.0 &&
				tri1 -> axis[2] <= 0.0) {
				edg -> type = BACKFACING;
				continue;
			}
			/* forward-facing polygons */
			if (tri0 -> axis[2] >= 0.0 &&
				tri1 -> axis[2] >= 0.0) {
				edg -> type = FRONTFACING;
				n_front++;
				continue;
			}
			edg -> type = SILHOUETTE;
			n_silhouette++;
			break;
		case CTR_SURFACE:
			if (tri0 == (struct polygon *) NULL) {
				set_error1("contour edge lying on no polygon");
				return;
			}
			if (tri0 -> axis[2] < 0.0) {
				edg -> type = BACKFACING;
				n_back++;
			}
			else {
				edg -> type = FRONTFACING;
				n_front++;
			}
			break;
		case NML_SURFACE:
		case BAS_SURFACE:
		default:
			edg -> type = FRONTFACING;
			n_front++;
			break;
		}
	}
	if (b -> type == PHN_SURFACE) {
		if (n_boundary > 0) {
			sprintf (message,"%8ld boundary edges for polyhedron", n_boundary);
			inform(message);
		}
		sprintf (message,"%8ld silhouette edges for polyhedron", n_silhouette);
		inform(message);
	}
	if (b -> type == CTR_SURFACE) {
		sprintf (message,"%8ld front-facing edges for contours", n_front);
		if (n_front > 0) inform(message);
		sprintf (message,"%8ld back-facing edges for contours", n_back);
		if (n_back > 0) inform(message);
	}
	if (b -> type == NML_SURFACE) {
		sprintf (message,"%8ld front-facing edges for normals", n_front);
		if (n_front) inform(message);
	}
	if (b -> type == BAS_SURFACE) {
		sprintf (message,"%8ld front-facing edges for sticks", n_front);
		if (n_front) inform(message);
	}
}

void hle_bunch (struct msscene *ms, struct surface *b)
{
	if (ms -> clipping) {
		cut_bunch (ms, b);
		if (error()) return;
		back_face_bunch (ms, b);
		if (error()) return;
	}
	else {
		back_face_bunch (ms, b);
		if (error()) return;
		cut_bunch (ms, b);
		if (error()) return;
	}
	put_bunch (ms, b);
	if (error()) return;
}

void back_face_bunch (struct msscene *ms, struct surface *b)
{
	int k;
	long nbf, nin;
	double dot0, dot1;
	double avcen[3], vec0[3], vec1[3];
	char message[MAXLINE];
	struct linseg *ls;
	struct phnedg *edg;
	struct polygon *tri0, *tri1;

	nbf = 0; nin = 0;
	for (ls = b -> head_linseg; ls != NULL; ls = ls -> next) {
		edg = ls -> edg;
		if (edg == NULL) continue;
		switch (edg -> type) {
		case BOUNDARY:
			tri0 = edg -> on[0];
			tri1 = edg -> on[1];
			if (tri0 != NULL) {
				if (tri0 -> axis[2] < 0.0) {
					ls -> inner = 1;
					nin++;
				}
			}
			else if (tri1 != NULL) {
				if (tri1 -> axis[2] < 0.0) {
					ls -> inner = 1;
					nin++;
				}
			}
			break;
		case SILHOUETTE:
			/* check for concave dihedral angle */
			tri0 = edg -> on[0];
			tri1 = edg -> on[1];
			if (tri0 != NULL && tri1 != NULL) {
				for (k = 0; k < 3; k++)
					avcen[k] = (tri0 -> center[k] + tri1 -> center[k])/2;
				for (k = 0; k < 3; k++) {
					vec0[k] = avcen[k] - tri0 -> center[k];
					vec1[k] = avcen[k] - tri1 -> center[k];
				}
				dot0 = dot_product (vec0, tri0 -> axis);
				dot1 = dot_product (vec1, tri1 -> axis);
				if (dot0 > 0.0 && dot1 > 0.0) {
					if (ms -> clipping || ms -> translucency) {
						ls -> inner = 1;
						nin++;
					}
					else {
						ls -> hidden = 1;
						nbf++;
					}
				}
			}
			break;
		case BACKFACING:
			if ((ms -> clipping || ms -> translucency) && b -> type == PHN_SURFACE) {
				ls -> inner = 1;
				nin++;
			}
			else {
				ls -> hidden = 1;
				nbf++;
			}
			break;
		case FRONTFACING:
		default:
			break;
		}
	}
	if (nbf > 0) {
		sprintf (message,
		"%8ld line segments removed because back-facing", nbf);
		inform(message);
	}
	if (nin > 0) {
		sprintf (message,
		"%8ld line segments inner because back-facing", nin);
		inform(message);
	}
}

void cut_bunch (struct msscene *ms, struct surface *obj)
{
	int j, k, result, okay, finished;
	long npue;
	double avedg, avls, xmax, xmin, ymax, ymin, d;
	double ends[2][2], height[2];
	char message[MAXLINE];
	struct linseg *ls, *last;
	struct phnedg *edg, *lsedg;
	struct phnvtx *vtx, *lsvtx, *edgvtx;
	struct molecule *mol;
	struct surface *poly_bunch;

	npue = 0;
	for (mol = ms -> head_molecule; mol != NULL; mol = mol -> next) {
		poly_bunch = polyhedron_bunch (mol);
		if (poly_bunch == NULL) continue;
		/* check all molecule edges in silhouette */
		for (edg = poly_bunch -> head_phnedg; edg != NULL; edg = edg -> next) {
			if (edg -> clipped) continue;
			if (edg -> type == FRONTFACING) continue;
			if (edg -> type == BACKFACING) continue;
			/* must be silhouette or boundary */
			xmax = -MS_INFINITY;
			xmin = MS_INFINITY;
			ymax = -MS_INFINITY;
			ymin = MS_INFINITY;
			for (j = 0; j < 2; j++) {
				vtx = edg -> pvt[j];
				for (k = 0; k < 2; k++)
					ends[j][k] = vtx -> center[k];
				height[j] = vtx -> center[2];
				if (xmax < ends[j][0]) xmax = ends[j][0];
				if (xmin > ends[j][0]) xmin = ends[j][0];
				if (ymax < ends[j][1]) ymax = ends[j][1];
				if (ymin > ends[j][1]) ymin = ends[j][1];
			}
			avedg = (height[0] + height[1]) / 2;
			/* check for line segment under molecule edge */
			last = obj -> tail_linseg;	/* save initial end of list */
			finished = 0;
			for (ls = obj -> head_linseg; ls != NULL; ls = ls -> next) {
				if (finished) break;
				/* don't cut new line segments added by this edge */
				if (ls == last) finished = 1;
				if (ls -> hidden) continue;
				/* check for line segment higher than edge */
				avls = (ls -> height[0] + ls -> height[1]) / 2;
				if (avls >= avedg) continue;
				if (ls -> ends[0][0] < xmin && ls -> ends[1][0] < xmin) continue;
				if (ls -> ends[0][0] > xmax && ls -> ends[1][0] > xmax) continue;
				if (ls -> ends[0][1] < ymin && ls -> ends[1][1] < ymin) continue;
				if (ls -> ends[0][1] > ymax && ls -> ends[1][1] > ymax) continue;
				lsedg = ls -> edg;
				if (lsedg != (struct phnedg *) NULL) {
					if (!(ms -> clipping || ms -> translucency) && lsedg -> type == BACKFACING) continue;
					if (edg == lsedg) continue;
					/* skip edges touching edge */
					okay = 1;
					for (j = 0; j < 2; j++) {
						lsvtx = lsedg -> pvt[j];
						/* contour edge check */
						if (lsvtx -> on == edg) okay = 0;
						if (!okay) break;
						for (k = 0; k < 2; k++) {
							/* molecule edge check */
							edgvtx = edg -> pvt[k];
							if (lsvtx == edgvtx) okay = 0;
							if (!okay) break;
							/* normal vector edge check */
							d = distance (lsvtx -> center, edgvtx -> center);
							if (d <= 0.0) okay = 0;
							if (!okay) break;
						}
					}
					if (!okay) continue;
				}
				/* check for line segment under edge */
				result = lines_intersect (ls -> ends[0], ls -> ends[1], ends[0], ends[1]);
				if (error()) return;
				if (!result) continue;
				/* line segment and edge intersect */
				cut_ls (obj, ls, edg);
				if (error()) return;
				npue++;
			}
		}
	}
	if (npue > 0) {
		sprintf (message,"%8ld line segments cut because under edge", npue);
		inform(message);
	}
}


void cut_ls (struct surface *obj, struct linseg *ls, struct phnedg *edg)
{
	int j, k, comp, atm, hue;
	double cens[2][3], ends[2][2], height[2];
	double ratio, value;
	double new_ends[2][2], new_height[2];
	struct phnvtx *vtx;
	
	for (j = 0; j < 2; j++) {
		vtx = edg -> pvt[j];
		for (k = 0; k < 3; k++)
			cens[j][k] = vtx -> center[k];
		for (k = 0; k < 2; k++)
			ends[j][k] = cens[j][k];
		height[j] = cens[j][2];
	}
	ratio = line_intersection (ls -> ends[0], ls -> ends[1], ends[0], ends[1]);
	/* create new line segment */
	for (k = 0; k < 2; k++) {
		new_ends[0][k] = (1.0 - ratio) * ls -> ends[0][k] +
			ratio * ls -> ends[1][k];
		new_ends[1][k] = ls -> ends[1][k];
	}
	new_height[0] = (1.0 - ratio) * ls -> height[0] +
		ratio * ls -> height[1];
	new_height[1] = ls -> height[1];
	comp = ls -> comp; atm = ls -> atm; hue = ls -> hue; value = ls -> value;
	make_linseg (obj, new_ends, new_height, comp, atm, hue, value, ls -> edg);
	if (error()) return;
	/* truncate old line segment */
	for (k = 0; k < 2; k++) {
		ls -> ends[1][k] = new_ends[0][k];
	}
	ls -> height[1] = new_height[0];
}


void put_bunch (struct msscene *ms, struct surface *obj)
{
	int k, result;
	long nput, nout;
	double pnt[3];
	char message[MAXLINE];
	struct linseg *ls;
	struct phnedg *edg;
	struct polygon *tri0, *tri1;

	nput = 0; nout = 0;
	for (ls = obj -> head_linseg; ls != NULL; ls = ls -> next) {
		if (ls -> hidden) continue;
		edg = ls -> edg;
		tri0 = edg -> on[0];
		tri1 = edg -> on[1];
		for (k = 0; k < 2; k++)
			pnt[k] = (ls -> ends[0][k] + ls -> ends[1][k]) / 2;
		pnt[2] = (ls -> height[0] + ls -> height[1]) / 2;
		if (pnt[0] < ms -> window[0][0] || pnt[0] > ms -> window[1][0]) {
			ls -> hidden = 1;
			nout++;
			continue;
		}
		if (pnt[1] < ms -> window[0][1] || pnt[1] > ms -> window[1][1]) {
			ls -> hidden = 1;
			nout++;
			continue;
		}
		if (!fetch_square (ms, pnt)) {
			if (error()) return;
			ls -> hidden = 1;
			nout++;
			continue;
		}
		if (error()) return;
		result = under_square (ms, pnt, tri0, tri1);
		if (error()) return;
		if (result) {
			ls -> hidden = 1;
			nput++;
		}
	}
	if (nput >= 0) {
		sprintf (message,
			"%8ld line segments removed because linseg under polygon", nput);
		inform(message);
	}
	if (nout > 0) {
		sprintf (message,
			"%8ld line segments removed because outside bounds", nout);
		inform(message);
	}
}

int put (double pnt[3], struct polygon *tri)
{
	int result, up;
	double x, y, d2, r2;
	char message[MAXLINE];

	if (tri -> clipped) return (0);
	if (pnt[2] >= tri -> center[2]) return (0);
	x = (pnt[0] - tri -> center[0]);
	y = (pnt[1] - tri -> center[1]);
	d2 = x * x + y * y;
	r2 = tri -> radius * tri -> radius;
	if (d2 >= r2) return (0);
	up = (tri -> axis[2] > 0.0);
	if (tri -> n_side == 3)
		result = pit (pnt, up, tri -> vc[0], tri -> vc[1], tri -> vc[2]);
	else if (tri -> n_side == 4)
		result = piq (pnt, up, tri -> vc[0], tri -> vc[1], tri -> vc[2], tri -> vc[3]);
	else {
		sprintf (message,
			"msdraw: polygon has wrong # of sides: %d", tri -> n_side);
		set_error1 (message);
		return (0);
	}
	if (error()) return(0);
	return (result);
}


struct linseg *allocate_linseg ()
{
	struct linseg *ls;

	ls = (struct linseg *) allocate_object (LINSEG);
	if (ls == NULL) {
		set_error1 ("(allocate_linseg): mem alloc fails");
		return(NULL);
	}
	return (ls);
}


void free_linseg (struct linseg *ls)
{
	free_object (LINSEG, (short *) ls);
}

