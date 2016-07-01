/* Molecular Surface Package Copyright 1995 Michael L. Connolly */
/* January 5, 2002 */

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

struct phnvtx *num2phnvtx (struct surface *phn, long num)
{
	struct phnvtx *p;
	if (phn -> phnvtx_handles == NULL) return (NULL);
	if (num <= 0) return (NULL);
	if (num > phn -> n_phnvtx) return (NULL);
	p = *(phn -> phnvtx_handles + num - 1);
	return (p);
}

struct phnedg *num2phnedg (struct surface *phn, long num)
{
	struct phnedg *p;
	if (phn -> phnedg_handles == NULL) return (NULL);
	if (num <= 0) return (NULL);
	if (num > phn -> n_phnedg) return (NULL);
	p = *(phn -> phnedg_handles + num - 1);
	return (p);
}

struct phntri *num2phntri (struct surface *phn, long num)
{
	struct phntri *p;
	if (phn -> phntri_handles == NULL) return (NULL);
	if (num <= 0) return (NULL);
	if (num > phn -> n_phntri) return (NULL);
	p = *(phn -> phntri_handles + num - 1);
	return (p);
}

struct polygon *num2polygon (struct surface *phn, long num)
{
	struct polygon *p;
	if (phn -> polygon_handles == NULL) return (NULL);
	if (num <= 0) return (NULL);
	if (num > phn -> n_polygon) return (NULL);
	p = *(phn -> polygon_handles + num - 1);
	return (p);
}


struct surface *polyhedron_bunch (struct molecule *mol)
{
	struct surface *poly_bunch, *h;
	/* find molecule's polyhedron bunch (assume only one) */
	poly_bunch = NULL;
	for (h = mol -> head_surface; h != NULL; h = h -> next)  {
		if (h -> type == PHN_SURFACE) {
			poly_bunch = h;
			break;
		}
	}
	return (poly_bunch);
}


int polygonize_polyhedron (struct surface *phn)
{
	long n_polygon;

	n_polygon = collect_phn_pgn (phn);
	if (n_polygon <= 0) {
		return(0);
	}
	do_polyhedron_bounds (phn);
	if (error ()) return (0);
	return (1);
}

long collect_phn_pgn (struct surface *phn)
{
	int j;
	int orn0, orn1, orn2;
	long idx;
	long n_polygon;
	long en0, en1, en2;
	unsigned long uen0, uen1, uen2;
	char message[MAXLINE];
	struct phnedg *edg0, *edg1, *edg2;
	struct phntri *ptri;
	struct polygon **ph, *pgn;
	struct polygon **polygon_handles;
    struct cept *ex;

	phn -> n_polygon = phn -> n_phntri;
	n_polygon = phn -> n_phntri;
	sprintf (message,"%8ld vertices; %8ld edges; %8ld polygons",
		phn -> n_phnvtx, phn -> n_phnedg, phn -> n_polygon);
	informd(message);
	/* allocate memory for polygon pointers */
	polygon_handles = (struct polygon **)
		allocate_pointers (POLYGON, phn -> n_polygon);
	if (polygon_handles == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_source (ex, "msphn.c");
        add_function (ex, "collect_phn_pgn");
        add_object (ex, POLYGON, "polygon_handles");
		return (0L);
	}

	/* transfer polygons, one by one */
	for (ph = polygon_handles, idx = 0; idx < n_polygon; ph++, idx++) {
		ptri = num2phntri (phn, idx + 1);
		if (ptri == NULL) {
			ex = new_cept (MEMORY_ERROR, INVALID_VALUE,  FATAL_SEVERITY);
			add_source (ex, "msphn.c");
			add_function (ex, "collect_phn_pgn");
			add_object (ex, PHNTRI, "tri");
            add_long (ex, "idx", idx + 1);
			return (0L);
		}
		pgn = allocate_polygon ();
		if (pgn == (struct polygon *) NULL) {
			add_function (tail_cept, "collect_phn_pgn");
			return (0L);
		}
		*(polygon_handles + idx) = pgn;
		*ph = pgn;
		pgn -> n_side = 3;
		pgn -> material_index = ptri -> hue;
		for (j = 0; j < 3; j++)
			pgn -> vertex_index[j] = ptri -> vtxnum[j]-1;
		en0 = ptri -> edgnum[0];
		en1 = ptri -> edgnum[1];
		en2 = ptri -> edgnum[2];
		pgn -> comp = ptri -> comp;
		pgn -> atm = ptri -> atm;
		pgn -> hue = ptri -> hue;
		uen0 = abs (en0);
		uen1 = abs (en1);
		uen2 = abs (en2);
		orn0 = (en0 < 0);
		orn1 = (en1 < 0);
		orn2 = (en2 < 0);
		if (uen0 < 1 || uen0 > phn -> n_phnedg) {
			ex = new_cept (MEMORY_ERROR, INVALID_VALUE,  FATAL_SEVERITY);
			add_source (ex, "msphn.c");
			add_function (ex, "collect_phn_pgn");
            add_long (ex, "uen0", uen0);
            add_long (ex, "n_phnedg", phn -> n_phnedg);
			return (0L);
		}
		if (uen1 < 1 || uen1 > phn -> n_phnedg) {
			ex = new_cept (MEMORY_ERROR, INVALID_VALUE,  FATAL_SEVERITY);
			add_source (ex, "msphn.c");
			add_function (ex, "collect_phn_pgn");
            add_long (ex, "uen1", uen1);
            add_long (ex, "n_phnedg", phn -> n_phnedg);
			return (0L);
		}
		if (uen2 < 1 || uen2 > phn -> n_phnedg) {
			ex = new_cept (MEMORY_ERROR, INVALID_VALUE,  FATAL_SEVERITY);
			add_source (ex, "msphn.c");
			add_function (ex, "collect_phn_pgn");
            add_long (ex, "uen2", uen2);
            add_long (ex, "n_phnedg", phn -> n_phnedg);
			return (0L);
		}
		edg0 = num2phnedg (phn, uen0);
		if (edg0 == NULL) {
			ex = new_cept (MEMORY_ERROR, INVALID_VALUE,  FATAL_SEVERITY);
			add_source (ex, "msphn.c");
			add_function (ex, "collect_phn_pgn");
			add_object (ex, PHNEDG, "edg");
            add_long (ex, "uen0", uen0);
			return (0L);
		}
		edg1 = num2phnedg (phn, uen1);
		if (edg1 == NULL) {
			ex = new_cept (MEMORY_ERROR, INVALID_VALUE,  FATAL_SEVERITY);
			add_source (ex, "msphn.c");
			add_function (ex, "collect_phn_pgn");
			add_object (ex, PHNEDG, "edg");
            add_long (ex, "uen1", uen1);
			return (0L);
		}
		edg2 = num2phnedg (phn, uen2);
		if (edg2 == NULL) {
			ex = new_cept (MEMORY_ERROR, INVALID_VALUE,  FATAL_SEVERITY);
			add_source (ex, "msphn.c");
			add_function (ex, "collect_phn_pgn");
			add_object (ex, PHNEDG, "edg");
            add_long (ex, "uen2", uen2);
			return (0L);
		}
		pgn -> edg[0] = edg0;
		pgn -> edg[1] = edg1;
		pgn -> edg[2] = edg2;
		pgn -> orn[0] = (short) orn0;
		pgn -> orn[1] = (short) orn1;
		pgn -> orn[2] = (short) orn2;
		/* pointers from edge to 2 polygons */
		edg0 -> on[orn0] = pgn;
		edg1 -> on[orn1] = pgn;
		edg2 -> on[orn2] = pgn;
	}

	/* set up linked list */
	phn -> head_polygon = *polygon_handles;
	for (idx = 0; idx < phn -> n_polygon - 1; idx++) {
		pgn = *(polygon_handles + idx);
		pgn -> next = *(polygon_handles + idx + 1);
	}
	phn -> tail_polygon = *(polygon_handles + phn -> n_polygon - 1);
	phn -> polygon_handles = polygon_handles;
	return (n_polygon);
}

void get_bunch_name (int type, char typename[25])
{
	if (type == PHN_SURFACE) strcpy (typename, "polyhedron");
	else if (type == CTR_SURFACE) strcpy (typename, "contours");
	else if (type == NML_SURFACE) strcpy (typename, "normals");
	else if (type == BAS_SURFACE) strcpy (typename, "sticks");
	else if (type == DEN_SURFACE) strcpy (typename, "density");
	else strcpy (typename, "unknown");
}

double get_function (struct surface *phn, struct phnvtx *vtx)
{
	double f;
	char function[MAXLINE];
	
	if (phn == (struct surface *) NULL) return (0.0);
	strcpy (function, phn -> function);
	
	switch (function[0]) {
	case 'u':
		f = vtx -> values[0];
		break;
	case 'v':
		f = vtx -> values[1];
		break;
	case 'w':
		f = vtx -> values[2];
		break;
	case 'x':
		f = vtx -> center[0];
		break;
	case 'y':
		f = vtx -> center[1];
		break;
	case 'z':
		f = vtx -> center[2];
		break;
	default:
		f = 0.0;
		break;
	}
	return (f);
}

int make_phn_simple (struct surface *phn)
{
	int j, k;
	long v, e, t;
	double *vc, *vn;
	long *ev, *tev;
	struct phnvtx *vtx;
	struct phnedg *edg;
	struct phntri *tri;
    struct cept *ex;

	if (phn == NULL) return (0);
	if (phn -> phnvtx_handles == NULL) return (0);
	if (phn -> phnedg_handles == NULL) return (0);
	if (phn -> phntri_handles == NULL) return (0);
	if (phn -> n_phnvtx == 0) return (0);
	if (phn -> n_phnedg == 0) return (0);
	if (phn -> n_phntri == 0) return (0);
	if (phn -> vertex_centers != NULL) return (0);
	if (phn -> vertex_normals != NULL) return (0);
	if (phn -> edgvtx != NULL) return (0);
	if (phn -> triedgvtx != NULL) return (0);
	phn -> vertex_centers = allocate_doubles (phn -> n_phnvtx * 3, 0, VERTEX_CENTERS);
	if (phn -> vertex_centers == NULL) return (0);
	phn -> vertex_normals = allocate_doubles (phn -> n_phnvtx * 3, 0, VERTEX_NORMALS);
	if (phn -> vertex_normals == NULL) return (0);
	phn -> edgvtx = allocate_longs (phn -> n_phnedg * 2, 0, EDGVTX);
	if (phn -> edgvtx == NULL) return (0);
	phn -> triedgvtx = allocate_longs (phn -> n_phntri * 6, 0, TRIEDGVTX);
	if (phn -> triedgvtx == NULL) return (0);
	/* transfer information */
	for (v = 0; v < phn -> n_phnvtx; v++) {
		vtx = num2phnvtx (phn, v + 1);
		if (vtx == NULL) {
			ex = new_cept (MEMORY_ERROR,  INVALID_VALUE,  FATAL_SEVERITY);
			add_source (ex, "msphn.c");
			add_function (ex, "make_phn_simple");
			add_object (ex, PHNVTX, "vtx");
            add_long (ex, "v", v + 1);
            add_long (ex, "n_phnvtx", phn -> n_phnvtx);
			return (0);
		}

		vc = phn -> vertex_centers + 3 * v;
		vn = phn -> vertex_normals + 3 * v;
		for (k = 0; k < 3; k++) {
			*(vc + k) = vtx -> center[k];
			*(vn + k) = vtx -> outward[k];
		}
	}
	for (e = 0; e < phn -> n_phnedg; e++) {
		edg = num2phnedg (phn, e + 1);
		if (edg == NULL) {
			ex = new_cept (MEMORY_ERROR,  INVALID_VALUE,  FATAL_SEVERITY);
			add_source (ex, "msphn.c");
			add_function (ex, "make_phn_simple");
			add_object (ex, PHNEDG, "edg");
            add_long (ex, "e", e + 1);
            add_long (ex, "n_phnedg", phn -> n_phnedg);
			return (0);
		}

		ev = phn -> edgvtx + 2 * e;
		for (j = 0; j < 2; j++)
			*(ev + j) = edg -> vtxnum[j];
	}
	for (t = 0; t < phn -> n_phntri; t++) {
		tri = num2phntri (phn, t + 1);
		if (tri == NULL) {
			ex = new_cept (MEMORY_ERROR,  INVALID_VALUE,  FATAL_SEVERITY);
			add_source (ex, "msphn.c");
			add_function (ex, "make_phn_simple");
			add_object (ex, PHNTRI, "tri");
            add_long (ex, "t", t + 1);
            add_long (ex, "n_phntri", phn -> n_phntri);
			return (0);
		}
		tev = phn -> triedgvtx + 6 * t;
			for (j = 0; j < 3; j++) {
				*(tev + j) = tri -> edgnum[j];
				*(tev + 3 + j) = tri -> vtxnum[j];
			}
	}
	return (1);
}

int link_polyhedron (struct surface *phn)
{
	long i;
	struct phnvtx *vtx;
	struct phnedg *edg;

	/* set up linked list of vertices */
	phn -> head_phnvtx = num2phnvtx (phn, 1);
	for (i = 1; i <= phn -> n_phnvtx - 1; i++) {
		vtx = num2phnvtx (phn, i);
		vtx -> next = num2phnvtx (phn, i + 1);
	}
	phn -> tail_phnvtx = num2phnvtx (phn, phn -> n_phnvtx);

	/* set up linked list of edges */
	phn -> head_phnedg = num2phnedg (phn, 1);
	for (i = 1; i <= phn -> n_phnedg - 1; i++) {
		edg = num2phnedg (phn, i);
		edg -> next = num2phnedg (phn, i + 1);
	}
	phn -> tail_phnedg = num2phnedg (phn, phn -> n_phnedg);
	return (1);
}

/* compute bound for each polygon */

void do_polyhedron_bounds (struct surface *phn)
{
	struct polygon *pgn;

	for (pgn = phn -> head_polygon; pgn != NULL; pgn = pgn -> next) {
		do_pgn_axis (pgn);
		if (error()) return;
		do_pgn_bound (pgn);
		if (error()) return;
	}
}



int free_phn (struct surface *phn)
{
	long v, e, t, p;
	struct phnvtx *pv;
	struct phnedg *pe;
	struct phntri *pt;
	struct polygon *pp;
	struct critlink *cl, *clnext;

	if (phn -> phnvtx_handles != NULL) {
		for (v = 1; v <= phn -> n_phnvtx; v++) {
			pv = num2phnvtx (phn, v);
			for (cl = pv -> head; cl != NULL; cl = clnext) {
				clnext = cl -> next;
				free_object (CRITLINK, (short *) cl);
			}
			free_phnvtx (pv);
		}
		free_pointers (PHNVTX, (phn -> phnvtx_handles));
		phn -> phnvtx_handles = NULL;
	}
	if (phn -> phnedg_handles != NULL) {
		for (e = 1; e <= phn -> n_phnedg; e++) {
			pe = num2phnedg (phn, e);
			free_phnedg (pe);
		}
		free_pointers (PHNEDG, (phn -> phnedg_handles));
		phn -> phnedg_handles = NULL;
	}
	if (phn -> phntri_handles != NULL) {
		for (t = 1; t <= phn -> n_phntri; t++) {
			pt = num2phntri (phn, t);
			free_phntri (pt);
		}
		free_pointers (PHNTRI, (phn -> phntri_handles));
		phn -> phntri_handles = NULL;
	}
	if (phn -> polygon_handles != NULL) {
		for (p = 1; p <= phn -> n_polygon; p++) {
			pp = num2polygon (phn, p);
			free_polygon (pp);
		}
		free_pointers (POLYGON, (phn -> polygon_handles));
		phn -> polygon_handles = NULL;
	}
	if (phn -> evaluation_spheres != NULL) {
		free_objects (EVALPNT, (short *) (phn -> evaluation_spheres));
		phn -> evaluation_spheres = NULL;
	}
	free_pointers (PHNTRI, (phn -> heads));
	phn -> heads = NULL;
	free_pointers (PHNTRI, (phn -> tails));
	phn -> tails = NULL;
	return (1);
}


struct surface *init_phn (long nv, long ne, long nt)
{
	int j, k;
	long i;
	char message[MAXLINE];
	struct surface *phn;
	struct phnvtx *pv;
	struct phnedg *pe;
	struct phntri *pt;
    struct cept *ex;

	phn = (struct surface *) allocate_object (SURFACE);
	if (phn == NULL) return (NULL);
	phn -> format = 1;
	/* store counts of vertices and triangles */
	phn -> n_phnvtx = nv;
	phn -> n_phnedg = ne;
	phn -> n_phntri = nt;
	sprintf (message, "%8ld vertices; %8ld edges; %8ld triangles",
		phn -> n_phnvtx, phn -> n_phnedg, phn -> n_phntri);
	inform (message);

	/* initialize bounds */
	for (k = 0; k < 3; k++) {
		phn -> bounds[0][k] =  1000000.0;
		phn -> bounds[1][k] = -1000000.0;
	}
	for (j = 0; j < 3; j++) {
		phn -> minvals[j] =  1000000.0;
		phn -> maxvals[j] = -1000000.0;
	}

	/* allocate memory for polyhedron */

	phn -> phnvtx_handles = (struct phnvtx **)
		allocate_pointers (PHNVTX, phn -> n_phnvtx);
	if (phn -> phnvtx_handles == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_source (ex, "msphn.c");
        add_function (ex, "init_phn");
        add_object (ex, PHNVTX, "phnvtx_handles");
		return(NULL);
	}
	for (i = 0; i < phn -> n_phnvtx; i++) {
		pv = allocate_phnvtx ();
		if (pv == NULL) return (NULL);
		*(phn -> phnvtx_handles + i) = pv;
	}
	
	/* allocate memory for edges */
	phn -> phnedg_handles = (struct phnedg **)
		allocate_pointers (PHNEDG, phn -> n_phnedg);
	if (phn -> phnedg_handles == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_source (ex, "msphn.c");
        add_function (ex, "init_phn");
        add_object (ex, PHNEDG, "phnedg_handles");
		return(NULL);
	}
	for (i = 0; i < phn -> n_phnedg; i++) {
		pe = allocate_phnedg ();
		if (pe == NULL) {
			add_function (tail_cept, "init_phn");
			return (NULL);
		}
		*(phn -> phnedg_handles + i) = pe;
	}

	phn -> phntri_handles = (struct phntri **)
		allocate_pointers (PHNTRI, phn -> n_phntri);
	if (phn -> phntri_handles == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_source (ex, "msphn.c");
        add_function (ex, "init_phn");
        add_object (ex, PHNTRI, "phntri_handles");
		return(NULL);
	}
	for (i = 0; i < phn -> n_phntri; i++) {
		pt = allocate_phntri ();
		if (pt == NULL) return (NULL);
		*(phn -> phntri_handles + i) = pt;
	}

	return (phn);
}

void do_pgn_axis (struct polygon *pgn)
{
	int j, k, orn0, orn1;
	double centers[MAXGON][3];
	double vect1[3], vect2[3], axis[3];
	struct phnedg *edg0, *edg1;
	struct phnvtx *v0, *v1, *v2;
    struct cept *ex;

	edg0 = pgn -> edg[0];
	edg1 = pgn -> edg[1];
	orn0 = pgn -> orn[0];
	orn1 = pgn -> orn[1];
	v0 = edg0 -> pvt[orn0];
	v1 = edg1 -> pvt[orn1];
	v2 = edg1 -> pvt[1-orn1];
	if (v0 == v1 || v0 == v2 || v1 == v2) {
		ex = new_cept (GEOMETRY_ERROR, DEGENERACY,  FATAL_SEVERITY);
		add_source (ex, "msphn.c");
		add_function (ex, "do_pgn_axis");
		add_object (ex, POLYGON, "pgn");
        add_message (ex, "two vertices of polygon are equal");
		return;
	}
	for (j = 0; j < pgn -> n_side; j++)
		for (k = 0; k < 3; k++)
			centers[j][k] = pgn -> edg[j] -> pvt[pgn -> orn[j]] -> center[k];
	for (k = 0; k < 3; k++) {
		vect1[k] = v1 -> center[k] - v0 -> center[k];
		vect2[k] = v2 -> center[k] - v0 -> center[k];
	}
	cross (vect1, vect2, axis);
	if (!normalize (axis)) {
		ex = new_cept (GEOMETRY_ERROR, DEGENERACY,  FATAL_SEVERITY);
		add_source (ex, "msphn.c");
		add_function (ex, "do_pgn_axis");
        add_double (ex, "centers[0][0]", centers[0][0]);
        add_double (ex, "centers[0][1]", centers[0][1]);
        add_double (ex, "centers[0][2]", centers[0][2]);
        add_double (ex, "centers[1][0]", centers[1][0]);
        add_double (ex, "centers[1][1]", centers[1][1]);
        add_double (ex, "centers[1][2]", centers[1][2]);
        add_message (ex, "cannot normalize axis");
		return;
	}
	for (k = 0; k < 3; k++)
		pgn -> axis[k] = axis[k];
}

void do_pgn_bound (struct polygon *pgn)
{
	int j, k, orn;
	double dvc;
	double pt1[3], pt2[3];
	struct phnedg *e;
	struct phnvtx *vtx;

	/* compute polygon center */
	for (k = 0; k < 3; k++)
		pgn -> center[k] = 0.0;
	for (j = 0; j < pgn -> n_side; j++) {
		e = pgn -> edg[j];
		orn = pgn -> orn[j];
		vtx = e -> pvt[orn];
		for (k = 0; k < 3; k++) {
			pgn -> vc[j][k] = vtx -> center[k];
			pgn -> center[k] += vtx -> center[k];
		}
	}

	for (k = 0; k < 3; k++)
		pgn -> center[k] /= pgn -> n_side;

	/* compute radius of circumscribing circle */
	pgn -> radius = 0.0;
	for (j = 0; j < pgn -> n_side; j++) {
		for (k = 0; k < 3; k++) {
			pt1[k] = pgn -> vc[j][k];
			pt2[k] = pgn -> center[k];
		}
		dvc = distance (pt1, pt2);
		if (dvc > pgn -> radius) pgn -> radius = dvc;
	}
}


/* tangents */

struct surface *polyhedron_tangents (struct surface *phn, double tan_length, int srf_type)
{
	int k;
	long v, nnull;
	double tan_factor, btn_factor;
	char message[MAXLINE];
	struct surface *obj;
	struct phnvtx *vtx, *vtx1, *vtx2, *tan_vtx;
	struct phnvtx **tan_vertices;
	struct phnedg *tan_edg;
	struct phnedg **tan_edges;
    struct cept *ex;

	/* check for null tangent vectors */
	nnull = 0;
	for (vtx = phn -> head_phnvtx; vtx != NULL; vtx = vtx -> next) {
		if (norm (vtx -> base) < 0.5) nnull++;
	}
	if (nnull > 0) return (NULL);

	/* allocate tangents bunch */
	obj = new_surface ();
	if (obj == NULL) {
		add_function (tail_cept, "polyhedron_tangents");
		return(NULL);
	}
	obj -> type = srf_type;
	
	if (phn -> n_phnvtx <= 0) {
		/* empty bunch */
		obj -> n_polygon = 0L;
		obj -> n_phnvtx = 0L;
		obj -> n_phnedg = 0L;
		return (obj);
	}

	/* allocate memory for tangent vertices */
	tan_vertices = (struct phnvtx **)
		allocate_pointers (PHNVTX, 2 * phn -> n_phnvtx);
	if (tan_vertices == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_source (ex, "msphn.c");
        add_function (ex, "polyhedron_tangents");
        add_object (ex, PHNVTX, "tan_vertices");
		return(NULL);
	}
	
	/* allocate memory for tangent edges */
	tan_edges = (struct phnedg **)
		allocate_pointers (PHNEDG, phn -> n_phnvtx);
	if (tan_edges == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_source (ex, "msphn.c");
        add_function (ex, "polyhedron_tangents");
        add_object (ex, PHNEDG, "tan_edges");
		return(NULL);
	}
	
	/* create tangent vertices and edges */
	for (vtx = phn -> head_phnvtx, v=0; vtx != NULL; vtx = vtx -> next, v++) {
		tan_factor = vtx -> values[1] * (1.0 - vtx -> values[0]);
		btn_factor = vtx -> values[1] * vtx -> values[0];
		vtx1 = allocate_phnvtx ();
		if (vtx1 ==  NULL) {
			add_function (tail_cept, "polyhedron_tangents");
			return ((struct surface *) NULL);
		}
		*(tan_vertices + 2 * v) = vtx1;
		if (srf_type == TAN_SURFACE)
			for (k = 0; k < 3; k++)
				vtx1 -> center[k] = vtx -> center[k] - tan_factor * tan_length * vtx -> base[k];
		else if (srf_type == BTN_SURFACE)
			for (k = 0; k < 3; k++)
				vtx1 -> center[k] = vtx -> center[k] - btn_factor * tan_length * vtx -> zenith[k];
		else return (NULL);
		vtx1 -> number = 2 * v + 1;
		vtx2 = allocate_phnvtx ();
		if (vtx2 ==  NULL) {
			add_function (tail_cept, "polyhedron_tangents");
			return ((struct surface *) NULL);
		}
		*(tan_vertices + 2 * v + 1) = vtx2;
		if (srf_type == TAN_SURFACE)
			for (k = 0; k < 3; k++)
				vtx2 -> center[k] = vtx -> center[k] + tan_factor * tan_length * vtx -> base[k];
		else if (srf_type == BTN_SURFACE)
			for (k = 0; k < 3; k++)
				vtx2 -> center[k] = vtx -> center[k] + btn_factor * tan_length * vtx -> zenith[k];
		else return (NULL);
		vtx2 -> number = 2 * v + 2;
		tan_edg = allocate_phnedg ();
		if (tan_edg ==  NULL) {
			add_function (tail_cept, "polyhedron_tangents");
			return ((struct surface *) NULL);
		}
		*(tan_edges + v) = tan_edg;
		tan_edg -> pvt[0] = vtx1;
		tan_edg -> pvt[1] = vtx2;
		tan_edg -> vtxnum[0] = vtx1 -> number;
		tan_edg -> vtxnum[1] = vtx2 -> number;
	}
	/* set up linked list */
	obj -> head_phnvtx = *tan_vertices;
	for (v = 0; v < 2 * phn -> n_phnvtx - 1; v++) {
		tan_vtx = *(tan_vertices + v);
		tan_vtx -> next = *(tan_vertices + v + 1);
	}
	obj -> head_phnedg = *tan_edges;
	for (v = 0; v < phn -> n_phnvtx - 1; v++) {
		tan_edg = *(tan_edges + v);
		tan_edg -> next = *(tan_edges + v + 1);
	}
		
	sprintf (message,"%8.3f tangent length;  %8ld tangent vectors",
		tan_length, phn -> n_phnvtx);
	inform(message);
	/* transfer to structure variables */
	obj -> n_polygon = 0L;
	obj -> n_phnvtx = 2 * phn -> n_phnvtx;
	obj -> n_phnedg = phn -> n_phnvtx;
	obj -> head_phnvtx = *tan_vertices;
	obj -> head_phnedg = *tan_edges;
	obj -> head_polygon = (struct polygon *) NULL;
	
	obj -> phnvtx_handles = tan_vertices;
	obj -> phnedg_handles = tan_edges;

	return (obj);
}

/* normals */

struct surface *polyhedron_normals (struct msscene *ms, double length)
{
	int k;
	long v;
	char message[MAXLINE];
	struct surface *obj, *po;
	struct phnvtx *vtx, *vtx1, *vtx2, *nml_vtx;
	struct phnvtx **nml_vertices;
	struct phnedg *nml_edg;
	struct phnedg **nml_edges;
    struct cept *ex;
	
	sprintf (message,"%8.3f length for polyhedron normals", length);
	inform (message);
	if (ms -> current_molecule == (struct molecule *) NULL)
		return ((struct surface *) NULL);
	po = polyhedron_bunch (ms -> current_molecule);
	if (po == (struct surface *) NULL || po -> type != PHN_SURFACE)
		return ((struct surface *) NULL);

	/* allocate normals bunch */
	obj = new_surface ();
	if (obj == NULL) {
		add_function (tail_cept, "polyhedron_normals");
		return(NULL);
	}
	obj -> type = NML_SURFACE;
	
	
	if (po -> n_phnvtx <= 0) {
		/* empty bunch */
		obj -> n_polygon = 0L;
		obj -> n_phnvtx = 0L;
		obj -> n_phnedg = 0L;
		return (obj);
	}
	
	
	
	/* allocate memory for normal vertices */
	nml_vertices = (struct phnvtx **)
		allocate_pointers (PHNVTX, 2 * po -> n_phnvtx);
	if (nml_vertices == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_source (ex, "msphn.c");
        add_function (ex, "polyhedron_normals");
        add_object (ex, PHNVTX, "nml_vertices");
		return(NULL);
	}
	
	/* allocate memory for normal edges */
	nml_edges = (struct phnedg **)
		allocate_pointers (PHNEDG, po -> n_phnvtx);
	if (nml_edges == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_source (ex, "msphn.c");
        add_function (ex, "polyhedron_normals");
        add_object (ex, PHNEDG, "nml_edges");
		return(NULL);
	}
	
	/* create normal vertices and edges */
	for (vtx = po -> head_phnvtx, v=0; vtx != NULL; vtx = vtx -> next, v++) {
		vtx1 = allocate_phnvtx ();
		if (vtx1 ==  NULL) {
			add_function (tail_cept, "polyhedron_normals");
			return ((struct surface *) NULL);
		}
		*(nml_vertices + 2 * v) = vtx1;
		for (k = 0; k < 3; k++)
			vtx1 -> center[k] = vtx -> center[k];
		vtx1 -> number = 2 * v + 1;
		vtx2 = allocate_phnvtx ();
		if (vtx2 ==  NULL) {
			add_function (tail_cept, "polyhedron_normals");
			return ((struct surface *) NULL);
		}
		*(nml_vertices + 2 * v + 1) = vtx2;
		for (k = 0; k < 3; k++)
			vtx2 -> center[k] = vtx -> center[k] + length * vtx -> outward[k];
		vtx2 -> number = 2 * v + 2;
		nml_edg = allocate_phnedg ();
		if (nml_edg ==  NULL) {
			add_function (tail_cept, "polyhedron_normals");
			return ((struct surface *) NULL);
		}
		*(nml_edges + v) = nml_edg;
		nml_edg -> pvt[0] = vtx1;
		nml_edg -> pvt[1] = vtx2;
		nml_edg -> vtxnum[0] = vtx1 -> number;
		nml_edg -> vtxnum[1] = vtx2 -> number;
	}
	/* set up linked list */
	obj -> head_phnvtx = *nml_vertices;
	for (v = 0; v < 2 * po -> n_phnvtx - 1; v++) {
		nml_vtx = *(nml_vertices + v);
		nml_vtx -> next = *(nml_vertices + v + 1);
	}
	obj -> head_phnedg = *nml_edges;
	for (v = 0; v < po -> n_phnvtx - 1; v++) {
		nml_edg = *(nml_edges + v);
		nml_edg -> next = *(nml_edges + v + 1);
	}
		
	sprintf (message,"%8.3f normal length;  %8ld normal vectors",
		length, po -> n_phnvtx);
	inform(message);
	/* transfer to structure variables */
	obj -> n_polygon = 0L;
	obj -> n_phnvtx = 2 * po -> n_phnvtx;
	obj -> n_phnedg = po -> n_phnvtx;
	obj -> head_phnvtx = *nml_vertices;
	obj -> head_phnedg = *nml_edges;
	obj -> head_polygon = (struct polygon *) NULL;
	
	obj -> phnvtx_handles = nml_vertices;
	obj -> phnedg_handles = nml_edges;

	return (obj);
}

struct phnvtx *allocate_phnvtx ()
{
	struct phnvtx *vtx;
    struct cept *ex;

	vtx = (struct phnvtx *) allocate_object (PHNVTX);
	if (vtx == NULL) {
		ex = new_cept ( MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_source (ex, "msphn.c");
        add_function (ex, "allocate_phnvtx");
        add_object (ex, PHNVTX, "vtx");
		return(NULL);
	}
	return (vtx);
}

void free_phnvtx (struct phnvtx *vtx)
{
	free_object (PHNVTX, (short *) vtx);
}

struct phnedg *allocate_phnedg ()
{
	struct phnedg *edg;
    struct cept *ex;

	edg = (struct phnedg *) allocate_object (PHNEDG);
	if (edg == NULL) {
		ex = new_cept ( MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_source (ex, "msphn.c");
        add_function (ex, "allocate_phnedg");
        add_object (ex, PHNEDG, "edg");
		return(NULL);
	}
	return (edg);
}

void free_phnedg (struct phnedg *edg)
{
	free_object (PHNEDG, (short *) edg);
}

struct phntri *allocate_phntri ()
{
	struct phntri *tri;
    struct cept *ex;

	tri = (struct phntri *) allocate_object (PHNTRI);
	if (tri == NULL) {
		ex = new_cept ( MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_source (ex, "msphn.c");
        add_function (ex, "allocate_phntri");
        add_object (ex, PHNTRI, "tri");
		return(NULL);
	}
	return (tri);
}

void free_phntri (struct phntri *tri)
{
	free_object (PHNTRI, (short *) tri);
}



long phnvtx_size ()
{
	return (sizeof (struct phnvtx));
}

long phnedg_size ()
{
	return (sizeof (struct phnedg));
}

long phntri_size ()
{
	return (sizeof (struct phntri));
}

struct polygon *allocate_polygon ()
{
	struct polygon *pgn;
    struct cept *ex;

	pgn = (struct polygon *) allocate_object (POLYGON);
	if (pgn == NULL) {
		ex = new_cept ( MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_source (ex, "msphn.c");
        add_function (ex, "allocate_polygon");
        add_object (ex, POLYGON, "pgn");
		return(NULL);
	}
	return (pgn);
}

void free_polygon (struct polygon *pgn)
{
	free_object (POLYGON, (short *) pgn);
}


/* compute axes normal to each triangle for later use */

int do_axes (struct surface *msphn)
{
	int i, j, k, orn0, orn1;
	double area;
	double axis[3];
	double trico[3][3];
	struct phntri *tri;
	struct phnedg *edg0, *edg1;
	struct phnvtx *v0, *v1, *v2, *vs[3];
    struct cept *ex;

	for (i = 0; i < msphn -> n_phntri; i++) {
		tri = num2phntri (msphn, i + 1);
		if (tri == NULL) {
			ex = new_cept (MEMORY_ERROR, INVALID_VALUE,  FATAL_SEVERITY);
			add_source (ex, "msphn.c");
			add_function (ex, "do_axes");
			add_object (ex, PHNTRI, "tri");
            add_long (ex, "i", i + 1);
			return (0);
		}

		edg0 = tri -> edg[0];
		edg1 = tri -> edg[1];
		orn0 = tri -> orn[0];
		orn1 = tri -> orn[1];
		v0 = edg0 -> pvt[orn0];
		v1 = edg1 -> pvt[orn1];
		v2 = edg1 -> pvt[1-orn1];
		if (v0 == v1 || v0 == v2 || v1 == v2) {
			ex = new_cept (GEOMETRY_ERROR, DEGENERACY,  FATAL_SEVERITY);
			add_source (ex, "msphn.c");
			add_function (ex, "do_pgn_axis");
			add_object (ex, POLYGON, "pgn");
			add_message (ex, "two vertices of triangle are equal");
			return (0);
		}
		vs[0] = v0;
		vs[1] = v1;
		vs[2] = v2;
		for (j = 0; j < 3; j++)
			for (k = 0; k < 3; k++)
				trico[j][k] = vs[j] -> center[k];
		area = do_axis (trico, axis);
		if (area <= 0.0) {
			/* tiny triangle, use average of vertex normals */
			for (k = 0; k < 3; k++)
				axis[k] = 0.0;
			for (j = 0; j < 3; j++)
				for (k = 0; k < 3; k++)
					axis[k] += vs[j] -> outward[k];
			if (!normalize (axis)) {
				ex = new_cept (GEOMETRY_ERROR, DEGENERACY,  FATAL_SEVERITY);
				add_source (ex, "msphn.c");
				add_function (ex, "do_axes");
				add_message (ex, "triangle too tiny to compute normal");
				return (0);
			}
		}
		for (k = 0; k < 3; k++)
			tri -> axis[k] = axis[k];
		tri -> area = area;
	}
	return (1);
}


/* compute bound for each triangle (for quick intersection check) */

void do_bounds (struct surface *msphn)
{
	int i, j, k, orn;
	double dvc;
	struct phntri *tri;
	struct phnedg *e;
	struct phnvtx *vtx;
    struct cept *ex;

	for (i = 0; i < msphn -> n_phntri; i++) {
		tri = num2phntri (msphn, i + 1);
		if (tri == NULL) {
			ex = new_cept (MEMORY_ERROR, INVALID_VALUE,  FATAL_SEVERITY);
			add_source (ex, "msphn.c");
			add_function (ex, "do_bounds");
			add_object (ex, PHNTRI, "tri");
            add_long (ex, "i", i + 1);
			return;
		}

		/* compute triangle center */
		for (k = 0; k < 3; k++)
			tri -> center[k] = 0.0;
		for (j = 0; j < 3; j++) {
			e = tri -> edg[j];
			orn = tri -> orn[j];
			vtx = e -> pvt[orn];
			for (k = 0; k < 3; k++)
				tri -> center[k] += vtx -> center[k];
		}

		for (k = 0; k < 3; k++)
			tri -> center[k] /= 3;

		/* compute radius of circumscribing circle */
		tri -> radius = 0.0;
		for (j = 0; j < 3; j++) {
			e = tri -> edg[j];
			orn = tri -> orn[j];
			vtx = e -> pvt[orn];
			dvc = distance (vtx -> center, tri -> center);
			if (dvc > tri -> radius) tri -> radius = dvc;
		}
	}
}



int do_poly_bounds (struct surface *phn)
{
	int j, k;
	long i;
	char message[MAXLINE];
	struct phnvtx *pv;
    struct cept *ex;

	for (i = 0; i < phn -> n_phnvtx; i++) {
		pv = num2phnvtx (phn, i + 1);
		if (pv == NULL) {
			ex = new_cept (MEMORY_ERROR, INVALID_VALUE,  FATAL_SEVERITY);
			add_source (ex, "msphn.c");
			add_function (ex, "do_poly_bounds");
			add_object (ex, PHNTRI, "pv");
			add_long (ex, "vertex number", i + 1);
			return (0);
		}

		for (k = 0; k < 3; k++) {
			if (pv -> center[k] < phn -> bounds[0][k])
				phn -> bounds[0][k] = pv -> center[k];
			if (pv -> center[k] > phn -> bounds[1][k])
				phn -> bounds[1][k] = pv -> center[k];
		}
		for (j = 0; j < 3; j++) {
			if (pv -> values[j] < phn -> minvals[j])
				phn -> minvals[j] = pv -> values[j];
			if (pv -> values[j] > phn -> maxvals[j])
				phn -> maxvals[j] = pv -> values[j];
		}
	}
	if (phn -> minvals[2] < phn -> maxvals[2]) {
		sprintf (message, "%8.3f minimum w value", phn -> minvals[2]);
		inform(message);
		sprintf (message, "%8.3f maximum w value", phn -> maxvals[2]);
		inform(message);
	}
	for (k = 0; k < 3; k++)
		phn -> center[k] = (phn -> bounds[0][k] + phn -> bounds[1][k]) / 2.0;
	return (1);
}


/* calculate areas and volumes */

int measure_polyhedron (struct surface *srf)
{
	int j, k, orn, shape, n_share;
	int component_number;
	int vtx_number[3];
	long vertex_number, edge_number;
    long abs_edge_number;
	long t;
	double delta_area, delta_volume;
	double outward[3];
	double vtx_co[3][3];
	double vectors[3][3], side[3][3];
	char message[MAXLINE];
	struct phnvtx *pv_ptr[3];
	struct phnedg *pe;
	struct phntri *pt;
	struct component *cmp_ptr;
    struct cept *ex;

	/* intialization */
	srf -> total_area = 0.0;
	srf -> total_volume = 0.0;
	for (k = 0; k < 3; k++)
		srf -> shape_area[k] = 0.0;
	for (t = 0; t < srf -> n_phntri; t++) {
		pt = num2phntri (srf, t + 1);
		if (pt == NULL) {
			ex = new_cept (MEMORY_ERROR, INVALID_VALUE,  FATAL_SEVERITY);
			add_source (ex, "msphn.c");
			add_function (ex, "measure_polyhedron");
			add_object (ex, PHNTRI, "tri");
            add_long (ex, "t", t + 1);
			return (0);
		}
		for (j = 0; j < 3; j++) {
			edge_number = pt -> edgnum[j];
            abs_edge_number = abs (edge_number);
			pe = num2phnedg (srf, abs_edge_number);
			if (pe == NULL) {
				ex = new_cept (MEMORY_ERROR, INVALID_VALUE,  FATAL_SEVERITY);
				add_source (ex, "msphn.c");
				add_function (ex, "measure_polyhedron");
				add_object (ex, PHNEDG, "pe");
				add_long (ex, "abs_edge_number", abs_edge_number);
				return (0);
			}
			orn = edge_number < 0;
			vtx_number[j] = pe -> vtxnum[orn];
			vertex_number = vtx_number[j];
			pv_ptr[j] = num2phnvtx (srf, vertex_number);
			if (pv_ptr[j] == NULL) {
				ex = new_cept (MEMORY_ERROR, INVALID_VALUE,  FATAL_SEVERITY);
				add_source (ex, "msphn.c");
				add_function (ex, "measure_polyhedron");
				add_object (ex, PHNTRI, "pv_ptr");
				add_long (ex, "vertex_number", vertex_number);
				return (0);
			}
			for (k = 0; k < 3; k++)
				vtx_co[j][k] = pv_ptr[j] -> center[k];
		}
		/* component */
		shape = pt -> shape;
		component_number = pt -> comp;
		cmp_ptr = *(srf -> component_handles + component_number - 1);
		/* area and volume computations */
		for (j = 0; j < 3; j++) {
			for (k = 0; k < 3; k++) {
				vectors[j][k] = (vtx_co[j][k]);
				side[j][k] = (vtx_co[j][k] - vtx_co[0][k]);
			}
		}
		cross (side[1], side[2], outward);
		delta_area = norm (outward) / 2.0;
		/* degenerate triangle */
		if (delta_area <= 0.0) {
			sprintf (message, "measure: delta_area = %8.3f", delta_area);
			informd(message);
		}
		srf -> total_area += delta_area;
		/* add to shape areas */
		shape = pt -> shape;
		if (shape == CONVEX)
			n_share = 1;
		else if (shape == SADDLE)
			n_share = 2;
		else if (shape == CONCAVE)
			n_share = 3;
		else if (shape == CYLINDRICAL)
			n_share = 2;
		else {
            n_share = 1;
		}
		srf -> shape_area[n_share-1] += delta_area;
		delta_volume = triple_product (vectors[0], vectors[1], vectors[2]) / 6.0;
		srf -> total_volume += delta_volume;
		/* add to component sums */
		cmp_ptr -> parea += delta_area;
		cmp_ptr -> pvolume += delta_volume;
	}
	return(1);
}
