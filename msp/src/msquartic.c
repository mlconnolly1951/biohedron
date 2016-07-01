/* Molecular Surface Package */
/* Copyright 1995 by Michael L. Connolly */
/* December 20, 2001 */

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"


struct vertex *allocate_vertex ()
{
	struct vertex *vtx;
    struct cept *ex;

	vtx = (struct vertex *) allocate_object (VERTEX);
	if (vtx == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
		add_function (ex, "allocate_vertex");
		add_source (ex, "msquartic.c");
		return(NULL);
	}
	return (vtx);
}

struct vertex *new_vertex (double center[3], struct sphere *atm, struct probe *prb, struct arc *arcptr, struct face *fac)
{
	int k;
	struct vertex *vtx;
    struct cept *ex;

	/* allocate memory */
	vtx = (struct vertex *) allocate_object (VERTEX);
	if (vtx == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
		add_function (ex, "new_vertex");
		add_source (ex, "msquartic.c");
		return(NULL);
	}
	/* set up fields */
	for (k = 0; k < 3; k++)
		vtx -> center[k] = center[k];
	vtx -> atm = atm;
	vtx -> prb = prb;
	/* check for null value for tent vertices */
	if (arcptr != (struct arc *) NULL) {
		vtx -> lfn = arcptr -> lfn;
		vtx -> ofn = arcptr -> ofn;
	}
	else if (fac != (struct face *) NULL) {
		vtx -> lfn = fac -> lfn;
		vtx -> ofn = fac -> ofn;
	}
	return (vtx);
}

int link_vertex (struct surface *srf, struct vertex *vtx)
{
	/* link into list */
	if (srf -> head_vertex == NULL) srf -> head_vertex = vtx;
	else srf -> tail_vertex -> next = vtx;
	srf -> tail_vertex = vtx;
	/* increment number of vertices for molecule */
	srf -> n_vertex++;
	vtx -> number = srf -> n_vertex;
	return (1);
}

struct probe *new_probe (struct surface *this_srf, double center[3], double height, double altitude[3], int side, struct sphere *atm1, struct sphere *atm2, struct sphere *atm3, struct pair *torus12, struct pair *torus23, struct pair *torus13)
{
	int k, reverse_orientation;
    double volume;
	char message[MAXLINE];
	struct probe *prb;
    struct cept *ex;

	/* allocate memory */
	prb = (struct probe *) allocate_object (PROBE);
	if (prb == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
		add_function (ex, "new_probe");
		add_source (ex, "msquartic.c");
		return(NULL);
	}
	link_probe (this_srf, prb);

	/* set up fields */
	for (k = 0; k < 3; k++)
		prb -> center[k] = center[k];

	/* store whether a low probe */
	prb -> low = (height < this_srf -> probe_radius);
	prb -> height = height;
	prb -> unit_altitude_z = side * altitude[2];

	/* set reverse orientation flag */
    volume = tetrahedron_volume (center, atm1 -> center, atm2 -> center, atm3 -> center);
	if (volume == 0.0) {
		ex = new_cept (GEOMETRY_ERROR,  DEGENERACY,  FATAL_SEVERITY);
		add_function (ex, "new_probe");
		add_source (ex, "msquartic.c");
		add_message (ex, "flat tetrahedron");
		sprintf (message, "height = %12.6f, side = %3d, volume = %12.6f, atoms: %5d %5d %5d",
			height, side, volume, atm1 -> number, atm2 -> number, atm3 -> number);
		add_message (ex, message);
		sprintf (message, "altitude = %12.6f %12.6f %12.6f", altitude[0], altitude[1], altitude[2]);
		add_message (ex, message);
		return(NULL);
	}
	reverse_orientation = (volume < 0.0);

	/* set up prb -> a */
	if (!reverse_orientation) {
		prb -> atm[0] = atm1;
		prb -> atm[1] = atm2;
		prb -> atm[2] = atm3;
		prb -> atom[0] = (atomnum) prb -> atm[0] -> number;
		prb -> atom[1] = (atomnum) prb -> atm[1] -> number;
		prb -> atom[2] = (atomnum) prb -> atm[2] -> number;
		prb -> pairs[0] = torus12;
		prb -> pairs[1] = torus23;
		prb -> pairs[2] = torus13;
	}
	else {
		prb -> atm[0] = atm1;
		prb -> atm[1] = atm3;
		prb -> atm[2] = atm2;
		prb -> atom[0] = (atomnum) prb -> atm[0] -> number;
		prb -> atom[1] = (atomnum) prb -> atm[1] -> number;
		prb -> atom[2] = (atomnum) prb -> atm[2] -> number;
		prb -> pairs[0] = torus13;
		prb -> pairs[1] = torus23;
		prb -> pairs[2] = torus12;
	}
	prb -> natom = 3;
	return (prb);
}

void link_probe (struct surface *srf, struct probe *prb)
{
	/* link into list */
	if (srf -> head_probe == NULL) srf -> head_probe = prb;
	else srf -> tail_probe -> next = prb;
	srf -> tail_probe = prb;

	/* increment number of probes for molecule */
	srf -> n_probe++;
	prb -> number = srf -> n_probe;
	prb -> srf = srf;
}

void delink_probe (struct surface *srf, struct probe *prb)
{
	struct probe *previous, *next, *this;
	if (prb -> delinked) return;
	previous = NULL;
	next = prb -> next;
	for (this = srf -> head_probe; this != NULL; this = this -> next) {
		if (this -> next == prb) previous = this;
	}
	if (previous == NULL)
		srf -> head_probe = next;
	else previous -> next = next;
	if (srf -> tail_probe == prb)
		srf -> tail_probe = previous;

	prb -> next = NULL;
	prb -> delinked = 1;

	/* link into delinked linked list */
	if (srf -> head_delinked == NULL) srf -> head_delinked = prb;
	else srf -> tail_delinked -> next = prb;
	srf -> tail_delinked = prb;

	/* decrement number of probes for molecule */
	srf -> n_probe--;
}

struct arc *allocate_arc ()
{
	struct arc *a;
    struct cept *ex;

	a = (struct arc *) allocate_object (ARC);
	if (a == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
		add_function (ex, "allocate_arc");
		add_source (ex, "msquartic.c");
		return(NULL);
	}
	return (a);
}

struct arc *new_arc (struct circle *cir, struct vertex *vtx1, struct vertex *vtx2,
	int shape, int small, double alpha, long lfn, long ofn, long ffn)
{
	int k;
	double vector1[3], vector2[3];
	struct arc *a;

	/* allocate memory */
	a = (struct arc *) allocate_object (ARC);
	if (a == NULL) {
		return (NULL);
	}

	/* set up fields */
	a -> cir = cir;
	a -> shape = shape;
	a -> small = small;
	a -> alpha = alpha;
	a -> vtx[0] = vtx1;
	a -> vtx[1] = vtx2;
	a -> lfn = lfn;
	a -> ofn = ofn;
	a -> ffn = ffn;
	if (cir != NULL) {
		if (cir -> lfn < lfn) cir -> lfn = lfn;
	}

	/* compute phi angle */

	if (cir == NULL || cir -> radius <= 0.0 || shape == STRAIGHT) a -> phi = 0.0;
	else if (vtx1 != NULL && vtx2 != NULL) {
		for (k = 0; k < 3; k++) {
			vector1[k] = (vtx1 -> center[k] - cir -> center[k]) / cir -> radius;
			vector2[k] = (vtx2 -> center[k] - cir -> center[k]) / cir -> radius;
		}
		a -> phi = positive_angle (vector1, vector2, cir -> axis);
	}
	else a -> phi = 2 * PI;
	return (a);
}

void link_arc (struct surface *srf, struct arc *a)
{
	/* increment number of arcs for molecule */
	srf -> n_arc++;
	a -> number = srf -> n_arc;
	/* link into list */
	if (srf -> head_arc == NULL) srf -> head_arc = a;
	else srf -> tail_arc -> next = a;
	srf -> tail_arc = a;
}

struct edge *allocate_edge ()
{
	struct edge *edg;
    struct cept *ex;

	/* allocate memory */
	edg = (struct edge *) allocate_object (EDGE);
	if (edg == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
		add_function (ex, "allocate_edge");
		add_source (ex, "msquartic.c");
		return(NULL);
	}
	return (edg);
}

struct edge *new_edge (struct arc *a, int orn, struct face *fac, struct cusp *csp)
{
	struct edge *edg;

	if (a == NULL) {
		set_error1("new_edge: null arc pointer");
		return (NULL);
	}
	if (orn != 0 && orn != 1) {
		set_error1("new_edge: invalid orientation");
		return (NULL);
	}
	/* allocate memory */
	edg = allocate_edge ();
	if (edg == NULL) {
		return(NULL);
	}
	/* set up fields */
	edg -> arcptr = a;
	edg -> orn = (short) orn;
	if (csp != NULL) {
		a -> csp = csp;
	}
	a -> edg[orn] = edg;
	if (fac != NULL) {
		edg -> fac = fac;
		edg -> lfn = a -> lfn;
	}
	return (edg);
}

struct central *new_central (double center[3], double radius, double axis[3])
{
	int k;
	struct central *ctl;
    struct cept *ex;

	/* allocate memory */
	ctl = (struct central *) allocate_object (CENTRAL);
	if (ctl == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
		add_function (ex, "new_central");
		add_source (ex, "msquartic.c");
		return (NULL);
	}


	/* set up fields */
	ctl -> radius = radius;
	for (k = 0; k < 3; k++) {
		ctl -> center[k] = center[k];
		ctl -> axis[k] = axis[k];
	}
	return (ctl);
}

void free_central (struct central *cen)
{
	free_object (CENTRAL, (short *) cen);
}

void link_central (struct surface *srf, struct central *ctl)
{
	/* link into list */
	if (srf -> head_central == NULL) srf -> head_central = ctl;
	else srf -> tail_central -> next = ctl;
	srf -> tail_central = ctl;
	/* increment number of central circles for molecule */
	srf -> n_central++;
}

void clear_pqms (struct surface *srf)
{
	free_cycle_handles (srf);
	free_vertex_handles (srf);
	free_arc_handles (srf);
	free_edge_handles (srf);
	free_circle_handles (srf);
	free_face_handles (srf);
}

void free_surface (struct surface *srf)
{
	free_cycles (srf);
	free_faces (srf);
	free_ae (srf);
	free_vertices (srf);
	free_circles (srf);
	free_varieties (srf);
	free_variety_handles (srf);
	free_probes (srf);
	free_tori (srf);
	free_chunks (srf);
	free_components (srf);
	if (srf -> atom_centers != NULL)
		free_doubles (srf -> atom_centers, 0, ATOM_CENTERS);
	if (srf -> atom_radii != NULL)
		free_doubles (srf -> atom_radii, 0, ATOM_RADII);
	if (srf -> densities != NULL)
		free_floats (srf -> densities);
}

struct variety *allocate_variety ()
{
	struct variety *vty;
    struct cept *ex;

	vty = (struct variety *) allocate_object (VARIETY);
	if (vty == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
		add_function (ex, "allocate_variety");
		add_source (ex, "msquartic.c");
		return(NULL);
	}
	return (vty);
}

void free_variety (struct variety *vty)
{
	free_object (VARIETY, (short *) vty);
}

struct torus *allocate_torus ()
{
	struct torus *tor;
    struct cept *ex;

	tor = (struct torus *) allocate_object (TORUS);
	if (tor == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
		add_function (ex, "allocate_torus");
		add_source (ex, "msquartic.c");
		return(NULL);
	}
	return (tor);
}

void free_vertex (struct vertex *vtx)
{
	if (vtx -> ce != (struct cusp_extension *) NULL)
		free_object (CUSP_EXTENSION, (short *) (vtx -> ce));
	free_object (VERTEX, (short *) vtx);
}

void free_torus (struct torus *tor)
{
	if (tor -> ce != (struct cusp_extension *) NULL)
		free_object (CUSP_EXTENSION, (short *) (tor -> ce));
	free_object (TORUS, (short *) tor);
}

void free_probe (struct probe *prb)
{
	free_object (PROBE, (short *) prb);
}

void free_arc (struct arc *a)
{
	free_object (ARC, (short *) a);
}

void free_edge (struct edge *edg)
{
	free_object (EDGE, (short *) edg);
}

void free_prbdot (struct surface *srf)
{
	free_doubles (srf -> prbdot, 0, PRBDOT);
	srf -> prbdot = (double *) NULL;
	srf -> n_prbdot = 0;
}

void free_probes (struct surface *srf)
{
	struct probe *prb, *next_prb;

	next_prb = NULL;
	for (prb = srf -> head_probe; prb != NULL; prb = next_prb) {
		next_prb = prb -> next;
		free_probe (prb);
	}
	next_prb = NULL;
	for (prb = srf -> head_delinked; prb != NULL; prb = next_prb) {
		next_prb = prb -> next;
		free_probe (prb);
	}
	srf -> n_variety -= srf -> n_probe;
	srf -> n_probe = 0;
}

void free_tori (struct surface *srf)
{
	struct torus *tor, *next_tor;

	next_tor = NULL;
	for (tor = srf -> head_torus; tor != NULL; tor = next_tor) {
		next_tor =  tor -> next;
		free_torus (tor);
	}
	srf -> n_variety -= srf -> n_tori;
	srf -> n_tori = 0;
}

void free_variety_handles (struct surface *srf)
{
	if (srf -> variety_handles == NULL) return;
	free_pointers (VARIETY, srf -> variety_handles);
	srf -> variety_handles = (struct variety **) NULL;
}

void free_varieties (struct surface *srf)
{
	int i;
	struct variety *vty;

	for (i = 0; i < srf -> n_variety; i++) {
		vty = *(srf -> variety_handles + i);
		if (vty == (struct variety *) NULL) continue;
		free_variety (vty);
	}
	srf -> n_variety = 0;
}

void free_vertex_handles (struct surface *srf)
{
	if (srf -> vertex_handles == NULL) return;
	free_pointers (VERTEX, srf -> vertex_handles);
	srf -> vertex_handles = NULL;
}

void free_vertices (struct surface *srf)
{
	struct vertex *vtx, *next_vtx;

	next_vtx = NULL;
	for (vtx = srf -> head_vertex; vtx != NULL; vtx = next_vtx) {
		next_vtx = vtx -> next;
		free_vertex (vtx);
	}
	srf -> head_vertex = NULL;
	srf -> tail_vertex = NULL;
	srf -> n_vertex = 0;
}

void free_circle_handles (struct surface *srf)
{
	if (srf -> circle_handles == NULL) return;
	free_pointers (CIRCLE, srf -> circle_handles);
	srf -> circle_handles = NULL;
}

void free_circles (struct surface *srf)
{
	long nfreed = 0;
	char message[MAXLINE];
	struct circle *cir, *next_cir;

	if (srf -> head_circle == NULL) return;
	next_cir = NULL;
	for (cir = srf -> head_circle; cir != NULL; cir = next_cir) {
		next_cir = cir -> next;
		free_circle (cir);
		nfreed++;
	}
	srf -> head_circle = NULL;
	sprintf (message, "%8ld circles freed", nfreed);
	informd(message);
	srf -> tail_circle = NULL;
	srf -> n_circle = 0;
}

void free_face_handles (struct surface *srf)
{
	if (srf -> face_handles == NULL) return;
	free_pointers (FACE, srf -> face_handles);
	srf -> face_handles = NULL;
} 

void free_faces (struct surface *srf)
{
	struct face *fac, *next_fac;

	next_fac = NULL;
	for (fac =  srf -> head_face; fac != NULL; fac = next_fac) {
		next_fac = fac -> next;
		free_face (fac);
	}
	srf -> head_face = NULL;
	srf -> tail_face = NULL;
	srf -> n_face = 0;
} 

void free_arc_handles (struct surface *srf)
{
	if (srf -> arc_handles == NULL) return;
	free_pointers (ARC, srf -> arc_handles);
	srf -> arc_handles = (struct arc **) NULL;
} 


void free_ae (struct surface *srf)
{
	int i, j, m;
	struct arc *a, *next_arc;
	struct face *fac;
	struct cycle *cyc;
	struct edge *edg;
	struct edge *e;

	/* free arcs and edges */
	next_arc = NULL;	/* make lint happy */
	for (a = srf -> head_arc; a != NULL; a = next_arc) {
		next_arc = a -> next;
		for (j = 0; j < 2; j++)
			if ((edg = a -> edg[j]) != NULL) free_edge (edg);
		free_arc (a);
	}
	srf -> n_edge = 0;
	srf -> n_arc = 0;
	srf -> head_arc = srf -> tail_arc = (struct arc *) NULL;
}

void free_edge_handles (struct surface *srf)
{
	if (srf -> edge_handles == NULL) return;
	free_pointers (EDGE, srf -> edge_handles);
	srf -> edge_handles = NULL;
}

void count_problem_faces (struct surface *srf)
{
	char message[MAXLINE];
	struct face *f;

	srf -> n_problem_face = 0;
	for (f = srf -> head_face; f != NULL; f = f -> next)
		if (f -> problem) srf -> n_problem_face++;
	sprintf (message,"%8ld problem faces", srf -> n_problem_face);
	if (srf -> n_problem_face) inform(message);
}


/* free arc and vertices (but not circle!) */
void frearc (struct arc *a)
{
	int j;
	struct vertex *v;

	/* do not free circle: used by leaf */

	for (j = 0; j < 2; j++) {
		v = a -> vtx[j];
		if (v != NULL) free_vertex (v);
		a -> vtx[j] = NULL;
	}
	free_arc (a);
}


struct cycle *allocate_cycle ()
{
	struct cycle *this_cyc;
    struct cept *ex;

	/* allocate memory */
	this_cyc = (struct cycle *) allocate_object (CYCLE);
	if (this_cyc == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
		add_function (ex, "allocate_cycle");
		add_source (ex, "msquartic.c");
		return(NULL);
	}
	return (this_cyc);
}

struct cycle *new_cycle (struct face *given_face, struct edge *first_edge)
{
	struct cycle *prev_cyc, *this_cyc, *cyc;

	/* allocate memory */
	this_cyc = allocate_cycle ();
	if (this_cyc == NULL) {
		return(NULL);
	}
	this_cyc -> first_edge = first_edge;
	/* link into face's list */
	if (given_face != NULL) {
		prev_cyc = NULL;
		for (cyc = given_face -> first_cycle; cyc != NULL; cyc = cyc -> next)
			prev_cyc = cyc;
		if (prev_cyc == NULL) given_face -> first_cycle = this_cyc;
		else prev_cyc -> next = this_cyc;
		given_face -> n_cycle++;
	}
	return (this_cyc);
}

void ink_cycle (struct surface *srf)
{
	/* increment number of cycles for molecule */
	srf -> n_cycle++;
}

int free_cycle (struct cycle *cyc)
{
	return free_object (CYCLE, (short *) cyc);
}

void delete_cycle (struct surface *srf, struct face *given_face, struct cycle *given_cyc)
{
	struct cycle *prev_cyc, *cyc;

	/* delink from face's list */
	prev_cyc = NULL;
	for (cyc = given_face -> first_cycle; cyc != NULL; cyc = cyc -> next)
		if (cyc -> next == given_cyc) prev_cyc = cyc;
	if (prev_cyc == NULL) given_face -> first_cycle = given_cyc -> next;
	else prev_cyc -> next = given_cyc -> next;

	/* decrement number of cycles for face */
	given_face -> n_cycle--;
	/* decrement number of cycles for molecule */
	srf -> n_cycle--;
	free_cycle (given_cyc);
}

void free_cycle_handles (struct surface *srf)
{
	int result;
	long c;
	struct face *fac;
	struct cycle *cyc, *next_cycle;
	struct cycle **cychdl, *cycptr;

	if (srf -> cycle_handles == NULL) return;
	free_pointers (CYCLE, srf -> cycle_handles);
	srf -> cycle_handles = NULL;
}

void free_cycles (struct surface *srf)
{
	int result;
	long c;
	struct face *fac;
	struct cycle *cyc, *next_cycle;
	struct cycle **cychdl, *cycptr;

	for (fac =  srf -> head_face; fac != NULL; fac = fac -> next) {
		if (fac -> problem) continue;
		next_cycle = NULL;
		for (cyc = fac -> first_cycle; cyc != (struct cycle *) NULL;
			cyc = next_cycle) {
			next_cycle = cyc -> next;
			free_cycle (cyc);
		}
	}
	srf -> n_cycle = 0;
}

void free_face (struct face *fac)
{
	struct probe *prb;
	
	if (fac -> shape == CONCAVE) {
		prb = fac -> ptr.prb;
		if (prb != NULL && !prb -> low &&
			fac -> ce != (struct cusp_extension *) NULL)
			free_object (CUSP_EXTENSION, (short *) (fac -> ce));
	}
	free_object (FACE, (short *) fac);
	if (error()) {
		return;
	}
}

void deep_free_face (struct face *fac)
{
	struct cycle *cyc, *nxt;
	struct edge *edg, *next;
	struct arc *a;
	int orn;
	struct vertex *vtx;

	if (fac == NULL) return;
	next = NULL;
	nxt = NULL;
	for (cyc = fac -> first_cycle; cyc != NULL; cyc = nxt) {
		for (edg = cyc -> first_edge; edg != NULL; edg = next) {
			a = edg -> arcptr;
			orn = edg -> orn;
			vtx = a -> vtx[orn];
			/* free vertex */
			if (vtx != NULL) free_vertex (vtx);
			/* do not free circle, because more than one arc points at it */
			/* free arc */
			free_arc (a);
			/* get pointer before too late */
			next = edg -> next;
			/* free edge */
			free_edge (edg);
		}
		/* get pointer before too late */
		nxt = cyc -> next;
		/* free cycle */
		free_cycle (cyc);
	}
	/* free face itself */
	free_face (fac);
}


struct face *allocate_face ()
{
	struct face *fac;
    struct cept *ex;

	fac = (struct face *) allocate_object (FACE);
	if (fac == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
		add_function (ex, "allocate_face");
		add_source (ex, "msquartic.c");
		return(NULL);
	}
	return (fac);
}

struct face *new_face (struct variety *vty, int shape)
{
	struct face *fac;

	/* allocate face */
	fac = (struct face *) allocate_object (FACE);
	if (fac == NULL || error ()) {
		return(NULL);
	}
	/* set up fields */
	fac -> vty = vty;
	fac -> shape = (short) shape;
	return (fac);
}

/* duplicate old face to create new face */
struct face *duplicate_face (struct face *oldfac)
{
	int shape;
	struct variety *vty;
	struct face *fac;
	struct surface *srf;

	shape = oldfac -> shape;
	vty = oldfac -> vty;
	srf = oldfac -> srf;
	fac = new_face (vty, shape);
	if (fac == NULL) {
		set_error1 ("duplicate_face: mem alloc fails");
		return (NULL);
	}
	fac -> ptr = oldfac -> ptr;
	fac -> comp = oldfac -> comp;
	fac -> alpha = oldfac -> alpha;
	fac -> lfn = oldfac -> lfn;
	fac -> ofn = oldfac -> ofn;
	/* not sure about this next one */
	fac -> simplified = oldfac -> simplified;
	oldfac -> largeok = 0;			/* remaining face no longer large enough for add-circle */
	if (!link_face (srf, fac)) return (NULL);
	return (fac);
}

int link_face (struct surface *srf, struct face *fac)
{
	/* link into list */
	if (srf -> head_face == NULL) srf -> head_face = fac;
	else srf -> tail_face -> next = fac;
	srf -> tail_face = fac;
	
	fac -> srf = srf;

	/* increment number of faces for molecule */
	srf -> n_face++;
	return (1);
}


struct component *allocate_component ()
{
	struct component *cmp;
    struct cept *ex;
	
	cmp = (struct component *) allocate_object (COMPONENT);
    if (cmp == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
		add_function (ex, "allocate_component");
		add_source (ex, "msquartic.c");
		return (NULL);
    }
	return (cmp);
}

void free_component (struct component *cmp)
{
	free_object (COMPONENT, (short *) cmp);
}

void free_components (struct surface *srf)
{
	int component_number;
	struct component *cmp_ptr;
	
	for (component_number = 1; component_number <= srf -> n_component;
		component_number++) {
		cmp_ptr = get_component_ptr (srf, component_number);
		free_component (cmp_ptr);
	}
	if (srf -> component_handles != (struct component **) NULL)
		free_pointers (COMPONENT, srf -> component_handles);
	srf -> component_handles = (struct component **) NULL;
	srf -> n_component = 0;
}

void print_face (struct face *fac) {
	char message[MAX_STRING];
	int orn;
	struct vertex *vtx0, *vtx1;
	struct arc *a;
	struct edge *edg;
	struct cycle *cyc;

	sprintf(message, "face shape = %d", fac -> shape);
	inform (message);
	if (fac -> first_cycle != NULL) {
		for (cyc = fac -> first_cycle; cyc != NULL; cyc = cyc -> next) {
			for (edg = cyc -> first_edge; edg != NULL; edg = edg -> next) {
				a = edg -> arcptr;
				orn = edg -> orn;
				vtx0 = a -> vtx[orn];
				vtx1 = a -> vtx[1-orn];
				if (vtx0 == NULL || vtx1 == NULL) {
					inform ("circular edge");
					continue;
				}
				sprintf(message, "arc from %8ld to %8ld", vtx0 -> number, vtx1 -> number);
				inform (message);
			}
		}
	}
	else {
		for (edg = fac -> first_edge; edg != NULL; edg = edg -> next) {
			a = edg -> arcptr;
			orn = edg -> orn;
			print_edge (edg);
		}
	}
	return;
}

int print_edge (struct edge *edg) {
	int orn;
	struct vertex *vtx0, *vtx1;
	struct arc *a;
	char message[MAX_STRING];

	a = edg -> arcptr;
	orn = edg -> orn;
	vtx0 = a -> vtx[orn];
	vtx1 = a -> vtx[1-orn];
	if (vtx0 == NULL || vtx1 == NULL) {
		inform ("circular edge");
	}
	else {
		sprintf(message, "arc from %8ld to %8ld (%6ld)",
			vtx0 -> number, vtx1 -> number, a -> ofn);
		inform (message);
	}
	return (1);
}


