/*
	Molecular Surface Package

	based on: PIECEWISE QUARTIC MOLECULAR SURFACE Computer Program

	Copyright 1986, 1989 by Michael L. Connolly
	All rights reserved

	February 16, 2006
*/

/* CUSP TRIMMING */

/* based mainly on algorithms in my 1985 J.A.C. article */

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

/* analytical cusp routines: */

void cusp_trimming (struct surface *this_srf)
{
	char message[MAXLINE];
	initialize_cusp_trimming (this_srf); if (error()) return;
	do_point_cusps (this_srf); if (error()) return;
	do_saddle_cusps (this_srf); if (error()) return;
	create_point_tori (this_srf); if (error()) return;
	create_cusp_circles (this_srf); if (error()) return;
	assign_vertices_to_cusps (this_srf); if (error()) return;
	all_cusp_arcs (this_srf); if (error()) return;
	allocate_extensions (this_srf); if (error()) return;
	concave_face_edges (this_srf); if (error()) return;
	concave_face_cusps (this_srf); if (error()) return;
	axial_cusp (this_srf);
	if (this_srf -> n_eaten_arc != 0) {
		sprintf (message,"%8ld cusp arcs eaten", this_srf -> n_eaten_arc);
		inform(message);
	}
	if (this_srf -> do_cusp_intersections) {
		cusp_intersections (this_srf); if (error()) return;
	}
	make_full_cusp_arcs (this_srf); if (error()) return;
	probe_cycles (this_srf); if (error()) return;
	free_cusps (this_srf); if (error()) return;
}

/* initialize */
void initialize_cusp_trimming (struct surface *this_srf)
{
	long n_low_quartet;
	char message[MAXLINE];
	struct probe *prb;
	struct probe **prb_hdl;
    struct cept *ex;

	n_low_quartet = 0;
	this_srf -> n_low_probe = 0;
	this_srf -> n_cusp_circle = 0;
	this_srf -> n_cusp_arc = 0;
	this_srf -> n_circle_eaten = 0;
	this_srf -> n_point_torus = 0;
	this_srf -> n_eaten_arc = 0;
	this_srf -> n_non_axial = 0;
	this_srf -> n_cusp_resort = 0;
	this_srf -> n_vertex_eaten = 0;

	this_srf -> head_cusp = this_srf -> tail_cusp = NULL;

	/* count low probes */
	if (this_srf -> probe_radius > 0.0) {
		for (prb = this_srf -> head_probe; prb != NULL; prb = prb -> next) {
			if (prb -> low) {
				this_srf -> n_low_probe++;
				if (prb -> natom >= 4)
					n_low_quartet++;
			}
		}
	}

	if (this_srf -> n_low_probe <= 0) {
		/* no edge cusps to worry about */
		this_srf -> low_probe_hdl = NULL;
		return;
	}

	sprintf (message,"%8ld low probes", this_srf -> n_low_probe);
	inform(message);
	if (n_low_quartet > 0) {
		sprintf (message,"%8ld low quartets", n_low_quartet);
		inform(message);
	}

	/* allocate memory for low probe list */
	this_srf -> low_probe_hdl = (struct probe **)
		allocate_pointers (PROBE, this_srf -> n_low_probe);
	if (this_srf -> low_probe_hdl == NULL) {
		ex = new_cept ( MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_object (ex, PROBE, "low probe list");
        add_function (ex, "initialize_cusp_trimming");
        add_source (ex, "mscusp.c");
		return;
	}

	/* store pointers to low probes */
	prb_hdl = this_srf -> low_probe_hdl;			/* initialize pointer */
	for (prb = this_srf -> head_probe; prb != NULL; prb = prb -> next)
		if (prb -> low)
			*prb_hdl++ = prb;
			
	/* initialize numerical cusp stuff */
	makdot (this_srf);
}

/* create point cusp vertices and assign to tori */

void do_point_cusps (struct surface *this_srf)
{
	int k, one_of_two;
	double point_radius_squared, point_radius, torus_radius;
	double point_center[3];
	char message[MAXLINE];
	struct torus *tor;
	struct vertex *vtx;
	struct cusp_extension *vce, *tce;
    struct cept *ex;

	/* loop through tori */
	for (tor = this_srf -> head_torus; tor != NULL; tor = tor -> next) {
		torus_radius = tor -> radius;
		if (torus_radius >= this_srf -> probe_radius)
			continue;		/* check for self-intersection */
		/* cusp points must lie between atom centers */
		if (tor -> cir[0] -> theta <= 0.0) {
			informd ("         torus self-intersects outside saddle");
			continue;
		}
		if (tor -> cir[1] -> theta <= 0.0) {
			informd ("         torus self-intersects outside saddle");
			continue;
		}
		tce = (struct cusp_extension *) allocate_cusp_extension ();
		if (error()) {
			add_function (tail_cept, "do_point_cusps");
			add_object (tail_cept, CUSP_EXTENSION, "tce");
			return;
		}
		tor -> ce = tce;
		/* compute radius of s0 (see 1985 JAC article) */
		point_radius_squared = this_srf -> probe_radius * this_srf -> probe_radius -
			torus_radius * torus_radius;
		if (point_radius_squared <= 0.0) continue;
		point_radius = sqrt (point_radius_squared);

		/* compute point cusp coordinates */
		for (one_of_two = 0; one_of_two < 2; one_of_two++) {
			for (k = 0; k < 3; k++)
				point_center[k] = tor -> center[k] +
					(2 * one_of_two - 1) * point_radius * tor -> axis[k];
			/* create new vertex */
			vtx = new_vertex (point_center, tor -> atm[one_of_two], NULL, NULL, NULL);
			if (vtx == NULL) return;
			link_vertex (this_srf, vtx);
			vce = (struct cusp_extension *) allocate_cusp_extension ();
			if (error()) {
				add_function (tail_cept, "do_point_cusps");
				add_object (tail_cept, CUSP_EXTENSION, "vce");
				return;
			}
			vtx -> ce = vce;
			tce -> vtx[one_of_two] = vtx;
			vce -> tor = tor;
			this_srf -> n_point_cusp_vertices++;
		}
	}
	sprintf (message, "%8ld point cusp vertices", this_srf -> n_point_cusp_vertices);
	inform (message);
}

/* create saddle cusp faces */
void do_saddle_cusps (struct surface *this_srf)
{
	struct face *saddle_fac;
	struct torus *tor;
	

	/* look through face list for saddles with point cusps */
	for (saddle_fac = this_srf -> head_face; saddle_fac != NULL;
		saddle_fac = saddle_fac -> next) {
		if (saddle_fac -> shape != SADDLE) continue;
		tor = saddle_fac -> ptr.tor;
		if (tor -> ce == NULL) continue;	/* no cusp vtx */
		if (saddle_fac -> n_arc == 2)			/* no concave edges */
			saddle_cone (saddle_fac);
		else if (saddle_fac -> n_arc == 4)	/* concave arcs on saddle */
			saddle_triangle (saddle_fac);
		else continue;			/* new face from previous iteration */
	}
}


/* create two cones from self-intersecting hoop */

void saddle_cone (struct face *saddle_fac)
{
	struct arc *arc2;
	struct face *saddle_fac2;
	struct surface *this_srf;
    struct cept *ex;

	this_srf = saddle_fac -> srf;
	if (this_srf == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
		add_object (ex, SURFACE, "this_srf");
		add_function (ex, "saddle_cone");
		add_source (ex, "mscusp.c");
		return;
	}
	arc2 = saddle_fac -> arcsp[1];
    if (arc2 == NULL) {
		ex = new_cept (POINTER_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
		add_object (ex, SURFACE, "this_srf");
		add_function (ex, "saddle_cone");
		add_source (ex, "mscusp.c");
		return;
	}
	/* create second face */
	saddle_fac2 = new_face (NULL, SADDLE);
	if (saddle_fac2 == NULL) return;
	/* copy pointer to torus */
	saddle_fac2 -> ptr.tor = saddle_fac -> ptr.tor;
	link_face (this_srf, saddle_fac2);
	add2face (saddle_fac2, arc2, NULL, NULL, NULL);
	/* modify old face */
	saddle_fac -> arcsp[1] = NULL;
	saddle_fac -> n_arc = 1;
}

/* create two triangles from self-intersecting saddle rectangle */
void saddle_triangle (struct face *saddle_fac)
{
	struct torus *tor;
	struct arc *arc1, *arc2, *arc3, *arc4, *arc5, *arc6;
	struct vertex *point_1, *point_2;
	struct face *saddle_fac2, *concave_fac1, *concave_fac2;
	struct cusp_extension *tce;
	struct surface *this_srf;
    struct cept *ex;

	this_srf = saddle_fac -> srf;
	if (this_srf == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
		add_object (ex, SURFACE, "this_srf");
		add_function (ex, "saddle_triangle");
		add_source (ex, "mscusp.c");
		return;
	}
	/* put info into local variables */
	tor = saddle_fac -> ptr.tor;
	arc1 = saddle_fac -> arcsp[0];
    if (arc1 == NULL) {
		ex = new_cept (POINTER_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
		add_object (ex, SURFACE, "this_srf");
		add_function (ex, "saddle_triangle");
		add_source (ex, "mscusp.c");
		return;
	}
	arc2 = saddle_fac -> arcsp[1];
    if (arc2 == NULL) {
		ex = new_cept (POINTER_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
		add_object (ex, SURFACE, "this_srf");
		add_function (ex, "saddle_triangle");
		add_source (ex, "mscusp.c");
		return;
	}
	arc3 = saddle_fac -> arcsp[2];
    if (arc3 == NULL) {
		ex = new_cept (POINTER_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
		add_object (ex, SURFACE, "this_srf");
		add_function (ex, "saddle_triangle");
		add_source (ex, "mscusp.c");
		return;
	}
	arc4 = saddle_fac -> arcsp[3];
    if (arc4 == NULL) {
		ex = new_cept (POINTER_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
		add_object (ex, SURFACE, "this_srf");
		add_function (ex, "saddle_triangle");
		add_source (ex, "mscusp.c");
		return;
	}
	tce = tor -> ce;
	if (tce == (struct cusp_extension *) NULL) return;
	point_1 = tce -> vtx[0];
	point_2 = tce -> vtx[1];

	/* create two new arcs for 2nd halves of split arcs */
	arc5 = new_arc (arc1 -> cir, point_2, arc1 -> vtx[1], CONCAVE, 0, (double) 0.0, 0L, 0L, 0L);
	this_srf -> n_arc++;
	arc5 -> number = this_srf -> n_arc;
	arc5 -> tor = tor;		/* copy torus pointer */
	arc5 -> subtype = SHORTENED_SUBTYPE;
	arc6 = new_arc (arc3 -> cir, point_1, arc3 -> vtx[1], CONCAVE, 0, (double) 0.0, 0L, 0L, 0L);
	this_srf -> n_arc++;
	arc6 -> number = this_srf -> n_arc;
	arc6 -> tor = tor;		/* copy torus pointer */
	arc6 -> subtype = SHORTENED_SUBTYPE;

	/* modify old arcs (shorten) */

	arc1 -> vtx[1] = point_1;
	arc1 -> phi = arc_ang (arc1);		/* recompute */
	arc1 -> subtype = SHORTENED_SUBTYPE;

	arc3 -> vtx[1] = point_2;
	arc3 -> phi = arc_ang (arc3);		/* recompute */
	arc3 -> subtype = SHORTENED_SUBTYPE;

	/* create new saddle face */
	saddle_fac2 = new_face (NULL, SADDLE);
	saddle_fac2 -> ptr.tor = tor;		/* copy torus pointer */
	link_face (this_srf, saddle_fac2);
	add2face (saddle_fac2, arc5, arc2, arc3, NULL);

	/* modify old saddle face */
	saddle_fac -> n_arc = 3;
	saddle_fac -> arcsp[1] = arc6;
	saddle_fac -> arcsp[2] = arc4;
	saddle_fac -> arcsp[3] = NULL;

	/* add shortened edges to adjacent concave faces */
	concave_fac1 = arc1 -> fac;
	if (concave_fac1 == NULL) {
		ex = new_cept (PARAMETER_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
		add_object (ex, FACE, "concave_fac1");
		add_function (ex, "saddle_triangle");
		add_source (ex, "mscusp.c");
		return;
	}
	concave_fac2 = arc3 -> fac;
	if (concave_fac2 == NULL) {
		ex = new_cept (PARAMETER_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
		add_object (ex, FACE, "concave_fac2");
		add_function (ex, "saddle_triangle");
		add_source (ex, "mscusp.c");
		return;
	}

	append_arc (arc1, arc5, concave_fac1);
	if (error()) return;
	append_arc (arc3, arc6, concave_fac2);
	if (error()) return;
}

/* add new (shortened) concave edge to concave face */
void append_arc (struct arc *old_arc, struct arc *arc_new, struct face *concave_fac)
{
	int arc_index, done;
    struct cept *ex;

	done = 0;

	/* look for old arc,
		add new arc to corresponding position in array */

	for (arc_index = 0; arc_index < MAXPA; arc_index++) {
		if (concave_fac -> arcsp[arc_index] == NULL) continue;
		if (concave_fac -> arcsp[arc_index] == old_arc) {
			if (++concave_fac -> n_arc > MAX_ARC_ARRAY) {
				ex = new_cept (ARRAY_ERROR,  MSOVERFLOW,  FATAL_SEVERITY);
				add_object (ex, ARC, "face arc array");
				add_function (ex, "append_arc");
				add_source (ex, "mscusp.c");
				return;
			}
			concave_fac -> arcsp[arc_index+MAXPA] = arc_new;
			done = 1;
		}
	}

	/* check for success */
	if (!done) {
		ex = new_cept (LOGIC_ERROR,  NOT_FOUND,  FATAL_SEVERITY);
		add_object (ex, ARC, "old arc");
		add_function (ex, "append_arc");
		add_source (ex, "mscusp.c");
		return;
	}
}

/* create cusp circles */
void create_cusp_circles (struct surface *this_srf)
{
	int k, circle_eaten;
	long n_quartet_circle, n_quartet_eaten;
	double circle_radius, probe_separation;
	double probe_axis[3], circle_center[3], circle_axis[3];
	double signed_distance;
	double circle_probe_vector[3];
	char message[MAXLINE];
	struct probe *prb1, *prb2, *prb3;
	int l1, l2, l3;
	struct circle *cir;
	struct cusp *csp;
	struct cusp_link *cusplink1, *cusplink2;
    struct cept *ex;

	if (this_srf -> n_low_probe <= 1) return;	/* no cusp circle possible */
	n_quartet_circle = 0;
	n_quartet_eaten = 0;

	/* consider each pair of low probes */
	for (l1 = 0; l1 < this_srf -> n_low_probe; l1++)
		for (l2 = 0; l2 < this_srf -> n_low_probe; l2++) {
			/* retrieve pointers to probes */
			prb1 = *(this_srf -> low_probe_hdl + l1);
			prb2 = *(this_srf -> low_probe_hdl + l2);
			if (prb1 >= prb2) continue;	/* no duplication */
			/* compute interprobe vector and distance */
			for (k = 0; k < 3; k++)
				probe_axis[k] = prb2 -> center[k] - prb1 -> center[k];
			probe_separation = norm (probe_axis);
			if (probe_separation >= 2 * this_srf -> probe_radius) continue;
			/* no intersection */

			/* compute axis unit vector and midpoint */
			for (k = 0; k < 3; k++) {
				circle_axis[k] = probe_axis[k] / probe_separation;
				circle_center[k] =
					(prb1 -> center[k] + prb2 -> center[k]) / 2;
			}
			/* compute radius of cusp circle (see 1985 JAC article) */
			circle_radius = this_srf -> probe_radius * this_srf -> probe_radius -
				probe_separation * probe_separation / 4;
			if (circle_radius <= 0.0) continue;
			circle_radius = sqrt (circle_radius);

			/* check for circle inside some other probe */
			circle_eaten = 0;
			for (l3 = 0; l3 < this_srf -> n_low_probe; l3++) {
				/* retrieve pointers to probes */
				prb3 = *(this_srf -> low_probe_hdl + l3);
				if (prb3 == prb1) continue;
				if (prb3 == prb2) continue;
				circle_eaten = probe_eats_circle (this_srf -> probe_radius,
					prb3 -> center, circle_center,
					circle_radius, circle_axis);
				if (circle_eaten) break;
			}
			if (circle_eaten) {
				this_srf -> n_circle_eaten++;
				if (prb1 -> natom > 3 || prb2 -> natom > 3)
					n_quartet_circle++;
				continue;
			}

			/* create new circle */
			cir = new_circle (circle_center, circle_radius, circle_axis);
			if (cir == NULL) return;
			link_circle (this_srf, cir);
			cir -> subtype = CUSP_SUBTYPE;
			/* create cusp circle structure */
			csp = allocate_cusp ();
			if (error()) {
				add_function (tail_cept, "create_cusp_circles");
				add_object (tail_cept, CUSP, "csp");
				return;
			}
			this_srf -> n_cusp_circle++;
			if (prb1 -> natom > 3 || prb2 -> natom > 3)
				n_quartet_circle++;
			csp -> order = 1;	/* default value */
			/* set up fields for cusp circle */
			csp -> cir = cir;
			csp -> prb[0] = prb1;
			csp -> prb[1] = prb2;

			/* link into list of cusp circles for molecule */
			if (this_srf -> head_cusp == NULL) this_srf -> head_cusp = csp;
			else this_srf -> tail_cusp -> next = csp;
			this_srf -> tail_cusp = csp;

			/* add cusp circle to probe lists */
			/* first probe: */
			cusplink1 = allocate_cusp_link() ;
			if (error()) {
				add_function (tail_cept, "create_cusp_circles");
				add_object (tail_cept, CUSP_LINK, "cusplink1");
				return;
			}
			cusplink1 -> csp = csp;
			cusplink1 -> next = prb1 -> first_cusp;
			prb1 -> first_cusp = cusplink1;

			/* second probe: */
			cusplink2 = allocate_cusp_link() ;
			if (error()) {
				add_function (tail_cept, "create_cusp_circles");
				add_object (tail_cept, CUSP_LINK, "cusplink2");
				return;
			}
			cusplink2 -> csp = csp;
			cusplink2 -> next = prb2 -> first_cusp;
			prb2 -> first_cusp = cusplink2;

			/* add cusp circle to tori lists */
			append_cusp (this_srf, csp);

			/* compute theta angle (related to geodesic curvature) */
			for (k = 0; k < 3; k++)
				circle_probe_vector[k] =
					circle_center[k] - prb2 -> center[k]; /* or prb1? */
			signed_distance = (-dot_product (circle_axis,
				circle_probe_vector));
			cir -> theta = atan2 (signed_distance, circle_radius);
		}
	if (this_srf -> n_cusp_circle > 0) {
		sprintf (message,"%8ld cusp circles created", this_srf -> n_cusp_circle);
		inform(message);
	}
	if (n_quartet_circle > 0) {
		sprintf (message,"%8ld quartet cusp circles", n_quartet_circle);
		inform(message);
	}
	if (this_srf -> n_circle_eaten > 0) {
		sprintf (message,"%8ld cusp circles eaten", this_srf -> n_circle_eaten);
		inform(message);
	}
	if (n_quartet_eaten > 0) {
		sprintf (message,"%8ld quartet cusp circles eaten", n_quartet_eaten);
		inform(message);
	}
}

int probe_eats_circle (double probe_radius, double probe_center[3], double circle_center[3], double circle_radius, double circle_axis[3])
{
	int k, circle_eaten;
	double dist3, ang, sinang, law;
	double pcvect[3];

	circle_eaten = 0;
	dist3 = distance (circle_center, probe_center);
	if (dist3 + circle_radius <= probe_radius)
		circle_eaten = 1;
	if (circle_eaten) return (1);
	/* more elaborate check depending on circle orientation */
	for (k = 0; k < 3; k++)
		pcvect[k] = circle_center[k] - probe_center[k];
	ang = simple_angle (circle_axis, pcvect);
	sinang = sin (ang);
	sinang = fabs (sinang);		/* probably not needed */
	/* law of cosines (sin of related angle is used) */
	law = dist3 * dist3 + 2 * sinang * dist3 * circle_radius +
		circle_radius * circle_radius;
	if (law < 0.0) law = 0.0;
	law = sqrt (law);
	circle_eaten = (law < probe_radius);
	return (circle_eaten);
}

/* create list of point cusp tori */
void create_point_tori (struct surface *this_srf)
{
	char message[MAXLINE];
	struct torus *tor;
	struct torus ** torus_hdl;
	struct cusp_extension *tce;
    struct cept *ex;

	/* count number of point cusp tori */
	if (this_srf -> probe_radius > 0.0)
		for (tor = this_srf -> head_torus; tor != NULL; tor = tor -> next) {
			tce = tor -> ce;
			if (tce == (struct cusp_extension *) NULL) continue;
			if (tce -> vtx[0] != NULL) this_srf -> n_point_torus++;
		}

	/* if none, clear pointer to list */
	if (this_srf -> n_point_torus <= 0) {
		this_srf -> point_torus_hdl = NULL;
		return;
	}

	sprintf (message,"%8ld point cusp tori", this_srf -> n_point_torus);
	inform(message);

	this_srf -> point_torus_hdl = (struct torus **)
		allocate_pointers (TORUS, this_srf -> n_point_torus);
	if (this_srf -> point_torus_hdl == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
		add_object (ex, TORUS, "point torus handles");
		add_function (ex, "create_point_tori");
		add_source (ex, "mscusp.c");
		return;
	}

	torus_hdl = this_srf -> point_torus_hdl;			/* initialize pointer */
	/* store pointer in list */
	for (tor = this_srf -> head_torus; tor != NULL; tor = tor -> next) {
		tce = tor -> ce;
		if (tce == (struct cusp_extension *) NULL) continue;
		if (tce -> vtx[0] != NULL)
			*torus_hdl++ = tor;
	}
}

/* add cusp circle to tori lists */
void append_cusp (struct surface *this_srf, struct cusp *csp)
{
	int k1, k2, l1, l2;
	struct atom *a1, *a2, *b1, *b2;
	struct probe *prb1, *prb2;
	struct torus *tor;
	int torus_index;
	struct cusp_link *cusplink;
	struct cusp_extension *tce;
    struct cept *ex;

	if (this_srf -> n_point_torus <= 0) return;	/* no point cusp tori */

	/* store pointers to two probes in local variables */
	prb1 = csp -> prb[0];
	prb2 = csp -> prb[1];

	/* look for tori given probes by means of atoms */
	for (k1 = 0; k1 <= MAXPA; k1++) {
		/* atom and next atom */
		l1 = ((k1 < 2) ? k1 + 1 : 0);
		a1 = (struct atom *) (prb1 -> atm[k1]);
		if (a1 == NULL) continue;
		b1 = (struct atom *) (prb1 -> atm[l1]);
		if (b1 == NULL) continue;
		for (k2 = 0; k2 <= MAXPA; k2++) {
			/* atom and next atom */
			l2 = ((k2 < 2) ? k2 + 1 : 0);
			a2 = (struct atom *) (prb2 -> atm[k2]);
			if (a2 == NULL) continue;
			b2 = (struct atom *) (prb2 -> atm[l2]);
			if (b2 == NULL) continue;
			/* look for match */
			if (a1 == b2 && b1 == a2 || a1 == a2 && b1 == b2)
				/* probes share pair of atoms */
				/* look for torus for these 2 atoms */
				for (torus_index = 0; torus_index < this_srf -> n_point_torus;
					torus_index++) {
					tor = *(this_srf -> point_torus_hdl + torus_index);
					tce = tor -> ce;
					if (tce == (struct cusp_extension *) NULL) continue;
					if (a1 == (struct atom *) (tor -> atm[0]) && b1 == (struct atom *) (tor -> atm[1]) ||
						a1 == (struct atom *) (tor -> atm[1]) && b1 == (struct atom *) (tor -> atm[0])) {
						/* allocate cusp circle list entry */
						cusplink = allocate_cusp_link();
						if (error()) {
							add_function (tail_cept, "append_cusp");
							add_object (tail_cept, CUSP_LINK, "cusplink");
							return;
						}
						/* add to linked list */
						cusplink -> csp = csp;
						cusplink -> next = tce -> first_cusp;
						tce -> first_cusp = cusplink;
						cusplink -> orn = (a1 == (struct atom *) (tor -> atm[0]));
					}
				}		/* end of torus search */
		}				/* end of second probe loop */
	}					/* end of first probe loop */
}

/* add cusp vertices to cusp circles */

void assign_vertices_to_cusps (struct surface *this_srf)
{
	int k, orn;
	struct torus *tor;
	int torus_index;
	struct cusp_link *cusplink;
	struct cusp *csp;
	struct vertex *vtx;
	struct cusp_extension *tce;
    struct cept *ex;

	/* search through list of point cusp tori */
	for (torus_index = 0; torus_index < this_srf -> n_point_torus;
		torus_index++) {
		/* retrieve pointer to torus */
		tor = *(this_srf -> point_torus_hdl + torus_index);
		tce = tor -> ce;
		if (tce == (struct cusp_extension *) NULL) continue;
		/* search through list of cusp circles */
		for (cusplink = tce -> first_cusp; cusplink != NULL;
			cusplink = cusplink -> next ) {
			/* retrieve pointer to cusp circle */
			csp = cusplink -> csp;
			/* two vertices per torus */
			for (k = 0; k < 2; k++) {
				vtx = tce -> vtx[k];
				/* check for inconsistency */
				if (vtx == NULL) {
					ex = new_cept (PARAMETER_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
					add_object (ex, VERTEX, "torus cusp edge");
					add_function (ex, "assign_vertices_to_cusps");
					add_source (ex, "mscusp.c");
					return;
				}
				orn = (k == cusplink -> orn);
				add_vertex_to_cusp (csp, vtx, orn);
				if (error()) return;
			}
		}
	}
}

void add_vertex_to_cusp (struct cusp *csp, struct vertex *vtx, int orn)
{
	int n;
    struct cept *ex;

	n = csp -> n_vertex;
	if (n >= MAX_VERTEX_PER_CIRCLE) {
		ex = new_cept (ARRAY_ERROR,  MSOVERFLOW,  FATAL_SEVERITY);
		add_object (ex, VERTEX, "cusp circle array");
		add_function (ex, "add_vertex_to_cusp");
		add_source (ex, "mscusp.c");
		return;
	}
	csp -> vtx[n] = vtx;
	csp -> orn[n] = (short) orn;
	csp -> n_vertex = n + 1;
}

/* create arcs on cusp circles */

void all_cusp_arcs (struct surface *this_srf)
{
	char message[MAXLINE];
	struct cusp *csp;

	for (csp = this_srf -> head_cusp; csp != NULL; csp = csp -> next)
		make_cusp_arcs (this_srf, csp);
	if (this_srf -> n_cusp_arc > 0) {
		sprintf (message,"%8ld cusp arcs created", this_srf -> n_cusp_arc);
		inform(message);
	}
}

/* create arcs for this cusp circle */
void make_cusp_arcs (struct surface *this_srf, struct cusp *csp)
{
	int k, nv, i0, i1, e, iv;
	struct circle *cir;
	struct arc *a;
	struct vertex *v1, *v2, *v;
	struct arc *arcs[MAX_VERTEX_PER_CIRCLE/2];
	struct edge *edge_handles[MAX_VERTEX_PER_CIRCLE/2];
	struct cusp_extension *vce1, *vce2;
	short indices[MAX_VERTEX_PER_CIRCLE];
	double circle_points[MAX_VERTEX_PER_CIRCLE][3];
	double center[3];
	double radius;
	double axis[3];
	char message[MAXLINE];
    struct cept *ex;

	cir = csp -> cir;		/* transfer pointer to local variable */

	nv = csp -> n_vertex;

	if (nv > MAX_VERTEX_PER_CIRCLE) {
		ex = new_cept (ARRAY_ERROR,  MSOVERFLOW,  FATAL_SEVERITY);
		add_object (ex, VERTEX, "cusp circle array");
		add_function (ex, "make_cusp_arcs");
		add_source (ex, "mscusp.c");
		return;
	}

	/* should be even number of vertices */
	if (nv % 2 != 0) {
		ex = new_cept (LOGIC_ERROR,  INVALID_VALUE,  FATAL_SEVERITY);
		add_function (ex, "make_cusp_arcs");
		add_source (ex, "mscusp.c");
        add_message (ex, "total number of vertices not even");
		add_long (ex, "nv", (long) nv);
		return;
	}

	/* special routine for full (no vertices) cusp circles */
	if (nv <= 0) {
		csp -> full = 1;
		/* full_cusp (csp); */
		return;
	}

	for (k = 0; k < 3; k++) {
		center[k] = cir -> center[k];
		axis[k] = cir -> axis[k];
	}
	radius = cir -> radius;
	for (iv = 0; iv < nv; ++iv) {
		v = csp -> vtx[iv];
		for (k = 0; k < 3; k++)
			circle_points[iv][k] = v -> center[k];
	}
	
	
	/* call subroutine to sort circle points */

	sort_points (center, radius, axis, nv, circle_points, csp -> orn, indices);


	/* create cusp arcs from sorted vlist */

	e = 0;					/* initialize array index */
	for (iv = 0; iv < nv; e++, iv += 2) {
		i0 = indices[iv];
		i1 = indices[iv+1];
		v1 = csp -> vtx[i0];
		v2 = csp -> vtx[i1];
		/* create new arc */
		a = new_arc (cir, v1, v2, CONCAVE, 0, (double) 0.0, 0L, 0L, 0L);
		this_srf -> n_arc++;
		a -> number = this_srf -> n_arc;
		a -> subtype = AXIAL_SUBTYPE;
		this_srf -> n_cusp_arc++;
		arcs[e] = a;	/* store new arc */
		/* for cusp intersections, check for both arc vertices on
			same axis */
		vce1 = v1 -> ce;
		if (vce1 == (struct cusp_extension *) NULL) continue;
		vce2 = v2 -> ce;
		if (vce2 == (struct cusp_extension *) NULL) continue;
		if (vce1 -> tor == vce2 -> tor) {
			a -> perm = 1;		/* mark as permitting intersections */
			a -> tor = vce1 -> tor;	/* store pointer to torus */
		}
	}

	e = 0;					/* initialize array index */
	for (iv = 0; iv < nv; e++, iv += 2) {
		a = arcs[e];
		edge_handles[e] = new_edge (a, 0, NULL, csp);	/* store new edge */
		/* add to linked list for cusp circle */
		if (e == 0) csp -> first_edge = edge_handles[e];
		else edge_handles[e-1] -> next = edge_handles[e];
	}
}

void make_full_cusp_arcs (struct surface *this_srf)
{
	struct cusp *csp;

	for (csp = this_srf -> head_cusp; csp != NULL; csp = csp -> next) {
		if (csp -> full) full_cusp (this_srf, csp);
	}
}

/* check full cusp circle without vertices */
int full_cusp_ok (struct cusp *csp)
{
	int k, m, l1, l2;
	int in[2];
	double vect[3];
	struct circle *cir;
	struct probe *prb;
	struct face *fac;

	cir = csp -> cir;		/* retrieve pointer to circle */

	/* check for circle inside each concave face */
	for (m = 0; m < 2; m++) {
		/* retrieve probe and face pointers */
		prb = csp -> prb[m];
		fac = prb -> fac;
		/* vector from probe center to cusp circles center */
		for (k = 0; k < 3; k++)
			vect[k] = cir -> center[k] - prb -> center[k];

		/* cusp circle center inside concave face ? */
		l1 = vector_in_face (vect, fac);
		/* cusp circle intersects or might contain face ? */
		l2 = ciomcfn (cir, fac);
		in[m] = (l1 && !l2);	/* cusp circle inside face */
	}

	/* if prospective cusp circle fails
	containment for either face, return 0 */

	if (!in[0] || !in[1]) return (0);

	return (1);
}


/* full cusp circle without vertices */
void full_cusp (struct surface *this_srf, struct cusp *csp)
{
	int m, ok;
	struct circle *cir;
	struct arc *a;
	struct edge *edg;
	struct probe *prb;
	struct face *fac;

	cir = csp -> cir;		/* retrieve pointer to circle */

	ok = full_cusp_ok (csp);
	if (!ok) return;

	/* new arc */
	a = new_arc (cir, NULL, NULL, CONCAVE, 0, (double) 0.0, 0L, 0L, 0L);
	this_srf -> n_arc++;
	a -> number = this_srf -> n_arc;
	a -> subtype = NONAXIAL_SUBTYPE;

	/* add arc to both concave faces */
	for (m = 0; m < 2; m++) {
		prb = csp -> prb[m];
		fac = prb -> fac;
		edg = new_edge (a, m, NULL, csp);	/* allocate edge */
		/* link into list */
		edg -> next = fac -> first_edge;
		fac -> first_edge = edg;
	}
}

/* vector in concave face ? */
int vector_in_face (double vect[3], struct face *concave_fac)
{
	int m, k, n;
	double t;
	double va[MAXPA][3];
	struct probe *prb;

	/* store probe pointer in local variable */
	prb = concave_fac -> ptr.prb;

	/* create vectors */
	for (m = 0; m < prb -> natom; m++)
		for (k = 0; k < 3; k++)
			va[m][k] = prb -> atm[m] -> center[k] -
				prb -> center[k];

	/* check triple products */
	for (m = 0; m < prb -> natom; m++) {
		n = ((m < prb -> natom - 1) ? m + 1 : 0);
		t = triple_product (va[m], va[n], vect);
		if (t >= 0.0) return (0);		/* outside */
	}

	return (1);					/* inside */
}

int point_in_pyr (double pnt[3], struct face *fac)
{
	int k, inside;
	double vect[3];
	struct probe *prb;

	if (fac -> shape != CONCAVE) return (0);
	prb = fac -> ptr.prb;
	for (k = 0; k < 3; k++)
		vect[k] = pnt[k] - prb -> center[k];
	inside = vector_in_face (vect, fac);
	return (inside);
}

/* circle intersects or might contain negatively curved face  */
int ciomcfn (struct circle *cir, struct face *concave_fac)
{
	int m;
	double xpnt[2][3];
	char message[MAXLINE];
	struct arc *a;
	struct edge *edg;
	struct probe *prb;

	prb = concave_fac -> ptr.prb;
	/* check whether given circle intersects any of the circles of
		the boundary arcs */

	for (m = 0; m < MAX_ARC_ARRAY; m++) {
		a = concave_fac -> arcsp[m];
        if (a == NULL) continue;
		if (circle_circle (cir, a -> cir, xpnt)) return (1);
	}
	if (prb -> low) {
		for (edg = concave_fac -> first_edge; edg != NULL;
			edg = edg -> next) {
			a = edg -> arcptr;
			if (circle_circle (cir, a -> cir, xpnt)) {
				return (1);
			}
		}
	}
	return (0);
}

void allocate_extensions (struct surface *this_srf)
{
	struct face *fac;
	struct probe *prb;
	struct cusp_extension *ce;
	
	/* search for faces of low probes */
	for (fac = this_srf -> head_face; fac != NULL; fac = fac -> next) {
		if (fac -> shape != CONCAVE) continue;
		prb = fac -> ptr.prb;
		if (!prb -> low) continue;
		ce = (struct cusp_extension *) allocate_cusp_extension ();
		if (error()) {
			add_function (tail_cept, "allocate_extensions");
			add_object (tail_cept, CUSP_EXTENSION, "ce");
			return;
		}
		fac -> ce = ce;
	}
}

/* edge list for concave face: part1 */

void concave_face_edges (struct surface *this_srf)
{
	int m;
	struct face *fac;
	struct arc *a;
	struct probe *prb;
	struct edge *edg;

	/* search for faces of low probes */
	for (fac = this_srf -> head_face; fac != NULL; fac = fac -> next) {
		if (fac -> shape != CONCAVE) continue;
		prb = fac -> ptr.prb;
		if (!prb -> low) continue;
		/* add shortened arcs or long arcs */
		/* MAX_ARC_ARRAY originally defined to be 6 */
		for (m = 0; m < MAX_ARC_ARRAY; m++) {
			a = fac -> arcsp[m];
			if (a == NULL) continue;
			/* allocate new edge */
			edg = new_edge (a, 0, NULL, NULL);
			/* link into list for face */
			edg -> next = fac -> first_edge;
			fac -> first_edge = edg;
		}
	}
}

/* edge list for concave face: part 2 */
void concave_face_cusps (struct surface *this_srf)
{
	int orn;
	struct face *fac;
	struct probe *prb;
	struct arc *a;
	struct cusp_link *cusplink;
	struct cusp *csp;
	struct edge *edg, *dup_edg;

	/* search for faces of low probes */
	for (fac = this_srf -> head_face; fac != NULL; fac = fac -> next) {
		if (fac -> shape != CONCAVE) continue;
		prb = fac -> ptr.prb;
		if (!prb -> low) continue;

		/* search through cusp circles of low faces */
		for (cusplink = prb -> first_cusp; cusplink != NULL;
			cusplink = cusplink -> next) {
			csp = cusplink -> csp;
			/* reverse orientation for 2nd probe of cusp circle */
			orn = (csp -> prb[1] == prb);	/* 1 for 2nd probe */

			/* add cusp circle arcs to concave face */
			for (edg = csp -> first_edge; edg != NULL; edg = edg -> next) {
				a = edg -> arcptr;
				/* allocate new edge */
				dup_edg = new_edge (a, orn, NULL, csp);
				/* link into list for face */
				dup_edg -> next = fac -> first_edge;
				fac -> first_edge = dup_edg;
			}
		}
	}
}

void axial_cusp (struct surface *this_srf)		/* axial cusp arc intersections */
{
	struct probe *prb;
	struct face *fac;

	for (fac = this_srf -> head_face; fac != NULL; fac = fac -> next) {
		if (fac -> shape != CONCAVE) continue;
		prb = fac -> ptr.prb;
		if (!prb -> low) continue;
		one_axial (fac);
		if (error ()) return;
	}

	for (fac = this_srf -> head_face; fac != NULL; fac = fac -> next) {
		if (fac -> shape != CONCAVE) continue;
		prb = fac -> ptr.prb;
		if (!prb -> low) continue;
		delete_fully_eaten (this_srf, fac);
		if (error ()) return;
	}

}

/* axial cusp intersections and deletions for face */
void one_axial (struct face *concave_fac)
{
	struct probe *prb;
	struct edge *edg, *next_edg;
	struct arc *arcptr;

	prb = concave_fac -> ptr.prb;

	/* initialization */
	next_edg = NULL;

	/* consider each arc of concave face */
	for (edg = concave_fac -> first_edge; edg != NULL; edg = next_edg) {
		next_edg = edg -> next;		/* pointer to next arc */
		/* function returns true if edge eaten */
		if (edge_eaten (edg, prb)) {
			arcptr = edg -> arcptr;
			/* mark arc half-eaten */
			arcptr -> eaten++;
		}
	}
}

/* is edge eaten? */
int edge_eaten (struct edge *edg, struct probe *prb)
{
	int k, l, m, n;
	double na, nb, angle, engle;
	double base[3], axis[3], avect[3], evect[3];
	struct arc *a;
	struct torus *tor;
	struct probe *aprb, *eprb;
	struct sphere *other, *atm1, *atm2;
	struct cusp_link *ecusplink;
	struct cusp *ecsp, *csp;
	struct cusp_extension *tce;
    struct cept *ex;

	csp = edg -> arcptr -> csp;
	if (csp == NULL) return (0);	/* not a cusp arc */

	a = edg -> arcptr;
	tor = a -> tor;
	if (tor == NULL)	return (0);	/* not on a single torus */
	tce = tor -> ce;
	if (tce == (struct cusp_extension *) NULL) return (0);

	/* check that this is a permissable arc for arc-arc intersections */
	if (!a -> perm) return (0);

	/* find other atom, i.e. atom opposite torus */
	other = NULL;
	for (m = 0; m < prb -> natom; m++)
		if (prb -> atm[m] != tor -> atm[0] &&
			prb -> atm[m] != tor -> atm[1]) {
			other = (prb -> atm[m]);
			break;
		}
	if (other == NULL) {
		ex = new_cept (LOGIC_ERROR,  NOT_FOUND,  FATAL_SEVERITY);
		add_function (ex, "edge_eaten");
		add_source (ex, "mscusp.c");
		add_message (ex, "other atom (opposite torus) not found");
		return(0);
	}

	/* find atom associated with arc */
	l = ((m > 0) ? m - 1 : prb -> natom - 1);
	n = ((m < prb -> natom - 1) ? m + 1 : 0);
	/* these atoms are in order defined by probe, not by torus */
	/* if atm1 and atm2 were in wrong order -- trouble */
	atm1 = (prb -> atm[n]);
	atm2 = (prb -> atm[l]);

	/* base vector from torus to probe */
	for (k = 0; k < 3; k++) {
		base[k] = prb -> center[k] - tor -> center[k];
		/* don't use torus axis because orientation crucial */
		axis[k] = atm2 -> center[k] - atm1 -> center[k];
	}

	/* normalize base and axis vectors */
	nb = norm (base);
	na = norm (axis);
	for (k = 0; k < 3; k++) {
		base[k] /= nb;
		axis[k] /= na;
	}

	/* find the probe that helped generate this arc */
	if (csp -> prb[0] == prb) aprb = csp -> prb[1];
	else if (csp -> prb[1] == prb) aprb = csp -> prb[0];
	else {
		ex = new_cept (LOGIC_ERROR,  NOT_FOUND,  FATAL_SEVERITY);
		add_object (ex, PROBE, "the probe that helped generate this arc");
		add_function (ex, "edge_eaten");
		add_source (ex, "mscusp.c");
		return(0);
	}

	/* compute vector to this probe */
	for (k = 0; k < 3; k++)
		avect[k] = (aprb -> center[k] - tor -> center[k])
			/ tor -> radius;

	/* angle about axis */
	angle = positive_angle (base, avect, axis);

	/* search through cusp circles belonging to torus */
	for (ecusplink = tce -> first_cusp; ecusplink != NULL;
		ecusplink = ecusplink -> next) {
		/* store cusp cir in local variable */
		ecsp = ecusplink -> csp;
		if (ecsp == csp) continue;

		/* find the probe that helped create this cusp circle */
		/* only an arc on this concave face can be a replacement */
		if (ecsp -> prb[0] == prb) eprb = ecsp -> prb[1];
		else if (ecsp -> prb[1] == prb) eprb = ecsp -> prb[0];
		else {
			continue;
		}

		/* compute vector from torus center to this probe */
		for (k = 0; k < 3; k++)
			evect[k] = (eprb -> center[k] - tor -> center[k])
				/ tor -> radius;
		engle = positive_angle (base, evect, axis);

		/* compare angles as eating check */
		if (engle > angle) return (1);		/* eaten */
	}
	return (0);					/* uneaten */
}

void delete_fully_eaten (struct surface *this_srf, struct face *fac)
{
	struct arc *arcptr;
	struct edge *pedg, *edg, *next_edg;

	pedg = NULL;
	next_edg = NULL;
	/* consider each arc of concave face */
	for (edg = fac -> first_edge; edg != NULL; edg = next_edg) {
		next_edg = edg -> next;
		arcptr = edg -> arcptr;
		if (arcptr -> eaten >= 2) {
			/* fix pointers */
			if (pedg == NULL) fac -> first_edge = next_edg;
			else pedg -> next = next_edg;
			free_edge (edg);
			/* mark arc half-eater */
			arcptr -> eater++;
			if (arcptr -> eater >= 2) {
				/* the free_arc was previously commented out */
				free_arc (arcptr);
				this_srf -> n_arc--;
				/* increment number of arcs eaten */
				this_srf -> n_eaten_arc++;
			}
		}
		else pedg = edg; /* store pointer to previous arc */
	}
}

void cusp_intersections (struct surface *this_srf)		/* cusp arc intersections */
{
	long pass, nmore;
	char message[MAXLINE];

	/* the non-axial cusp intersections */
	for (pass = 0; pass < this_srf -> do_cusp_intersections; pass++) {
		nmore = non_axial_cusp (this_srf, pass);
		if (nmore <= 0) break;
		sprintf (message,"%8ld intersections found in pass %2ld",
			nmore, pass + 1);
		if (error ()) return;
	}

	if (this_srf -> n_non_axial != 0) {
		sprintf (message,
			"%8ld intersections for pairs of cusp edges",
			this_srf -> n_non_axial);
		inform(message);
	}
	if (this_srf -> n_cusp_resort != 0) {
		sprintf (message,
			"%8ld re-orderings of vertices along a cusp edge",
			this_srf -> n_cusp_resort);
		inform(message);
	}
	if (this_srf -> n_vertex_eaten > 0) {
		sprintf (message,
			"%8ld cusp vertices eaten", this_srf -> n_vertex_eaten);
		inform(message);
	}
}


/* non-axial cusp arc intersections */
long non_axial_cusp (struct surface *this_srf, long pass)
{
	long n, i, prev_non_axial;
	struct probe *prb;
	struct face *fac;
	struct face **faces, **ptr;
    struct cept *ex;

	n = 0;
	for (fac = this_srf -> head_face; fac != NULL; fac = fac -> next) {
		if (fac -> shape != CONCAVE) continue;
		prb = fac -> ptr.prb;
		if (!prb -> low) continue;
		n++;
	}
	if (n <= 0) return (0L);
	faces = (struct face **)
		allocate_pointers (FACE, n);
	if (faces == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_object (ex, FACE, "faces");
        add_function (ex, "non_axial_cusp");
        add_source (ex, "mscusp.c");
		return (0L);
	}
	for (fac = this_srf -> head_face, ptr = faces; fac != NULL; fac = fac -> next) {
		if (fac -> shape != CONCAVE) continue;
		prb = fac -> ptr.prb;
		if (!prb -> low) continue;
		*ptr++ = fac;
	}
	prev_non_axial = this_srf -> n_non_axial;
	for (i = 0; i < n; i++) {
		fac = *(faces + i);
		if (fac -> notrim) continue;
		if (fac -> problem) continue;
		fac -> changed = FALSE;
		/* first, set up the table of possible arc intersections */
		pre_non_axial (this_srf, fac, pass);
		if (error ()) return (0L);
		one_non_axial (this_srf, fac, pass);
		if (error ()) return (0L);
	}
	free_pointers (FACE, faces);
	return (this_srf -> n_non_axial - prev_non_axial);
}

void pre_non_axial (struct surface *this_srf, struct face *fac, long pass)
{
	int k, result;
	long i_nai, i, j, n;
	struct arc *arc1, *arc2;
	struct edge *edg1, *edg2;
	struct circle *cir1, *cir2;
	struct edge **edges, **ptr;
    struct cept *ex;

	n = 0;
	for (edg1 = fac -> first_edge; edg1 != NULL; edg1 = edg1 -> next) n++;
	if (n <= 0) return;
	edges = (struct edge **) allocate_pointers (EDGE, n);
	if (edges == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_object (ex, EDGE, "edges");
        add_function (ex, "pre_non_axial");
        add_source (ex, "mscusp.c");
		return;
	}
	for (edg1 = fac -> first_edge, ptr = edges; edg1 != NULL; edg1 = edg1 -> next, ptr++) {
		*ptr = edg1;
	}
	
	/* initialize */
	this_srf -> n_nai = 0;
	for (i_nai = 0; i_nai < MAX_NAI; i_nai++) {
		this_srf -> nais[i_nai].nx = (short) 0;
		this_srf -> nais[i_nai].eaten = (short) 0;
		for (k = 0; k < 3; k++)
			this_srf -> nais[i_nai].xpnt[k] = 0.0;
		this_srf -> nais[i_nai].closest_distance = MS_INFINITY;
	}
	/* check for pairs of permissable arcs for intersection */
	for (i = 0; i < n; i++) {
		edg1 = *(edges + i);
		if (fac -> problem) return;	/* maybe set in one_edge_edge */
		/* transfer info to local variables */
		arc1 = edg1 -> arcptr;
		cir1 = arc1 -> cir;
		if (!arc1 -> perm) continue;
		for (j = 0; j < n; j++) {
			edg2 = *(edges + j);
			if (fac -> problem) return;	/* maybe set in one_edge_edge */
			if (fac -> notrim) return;
			if (edg2 == edg1) continue;
			arc2 = edg2 -> arcptr;
			cir2 = arc2 -> cir;
			if (!arc2 -> perm) continue;
			if (cir2 == cir1) continue;
			result = one_edge_edge (this_srf, fac, edg1, edg2, pass);
			if (error ()) return;
			if (!result) continue;
		}
	}
	free_pointers (EDGE, edges);
}

/* return whether any intersection points */
int one_edge_edge (struct surface *this_srf, struct face *fac, struct edge *edg1, struct edge *edg2, int pass)
{
	int k, m, too_close, eatens[2];
	int orn1, swap12, pp, aa, vip1, vip2, vip3;
	double dists[2];
	double xpnt[2][3];
	char message[MAXLINE];
	struct probe *prb, *prb1, *prb2, *close_prb;
	struct arc *arc1, *arc2;
	struct circle *cir1, *cir2, *pcir1, *pcir2;
	struct cusp *csp, *csp1, *csp2;
	struct sphere *atm1, *atm2, *atm3;
	struct face *concave_fac1, *concave_fac2, *concave_fac3;
	struct probe *pprb1, *pprb2, *pprb3;
	struct vertex *vtx;
	struct vertex *old_vtx1, *old_vtx2;
	struct cusp_extension *vce;
    struct cept *ex;

	/* transfer info to local variables */
	prb = fac -> ptr.prb;
	if (!prb -> low) return (0);
	atm1 =  (prb -> atm[0]);
	atm2 =  (prb -> atm[1]);
	atm3 =  (prb -> atm[2]);
	arc1 = edg1 -> arcptr;
	cir1 = arc1 -> cir;
	csp1 = edg1 -> arcptr -> csp;
	arc2 = edg2 -> arcptr;
	cir2 = arc2 -> cir;
	csp2 = edg2 -> arcptr -> csp;
	orn1 = edg1 -> orn;

	/* find neighboring faces and probes */
	if (csp1 -> prb[0] == prb) prb1 = csp1 -> prb[1];
	else if (csp1 -> prb[1] == prb) prb1 = csp1 -> prb[0];
	else {
		ex = new_cept (LOGIC_ERROR,  NOT_FOUND,  FATAL_SEVERITY);
		add_object (ex, PROBE, "the first probe that helped generate this cusp");
		add_function (ex, "one_edge_edge");
		add_source (ex, "mscusp.c");
		return (0);
	}
	if (csp2 -> prb[0] == prb) prb2 = csp2 -> prb[1];
	else if (csp2 -> prb[1] == prb) prb2 = csp2 -> prb[0];
	else {
		ex = new_cept (LOGIC_ERROR,  NOT_FOUND,  FATAL_SEVERITY);
		add_object (ex, PROBE, "the second probe that helped generate this cusp");
		add_function (ex, "one_edge_edge");
		add_source (ex, "mscusp.c");
		return (0);
	}
	concave_fac1 = prb1 -> fac;
	concave_fac2 = prb2 -> fac;
	csp = find_cusp (prb1, prb2);
	if (csp == (struct cusp *) NULL) return (0);
	/* new arc depends on probe order */
	swap12 = (csp -> prb[0] != prb1);
	if (swap12) return (0);
	for (pp = 0; pp < pass - 1; pp++) {
		pcir1 = fac -> ce -> nacirs[pp][0];
		pcir2 = fac -> ce -> nacirs[pp][1];
		if (pcir1 == cir1 && pcir2 == cir2) return (0);
		if (pcir2 == cir1 && pcir1 == cir2) return (0);
	}
	/* clear intersection point coordinate array */
	for (m = 0; m < 2; m++)
		for (k = 0; k < 3; k++)
			xpnt[m][k] = 0.0;

	aa = edg_edg (edg1, edg2, xpnt);
	if (aa <= 0) return (0);		/* no intersection */
	if (aa == 1) {
		eatens[1] = 1;
		informd ("         single edge-edge intersection point");
	}
	/* check whether we already have this triple-probe intersection */
	old_vtx1 = NULL; old_vtx2 = NULL;
	for (vtx = this_srf -> head_vertex; vtx != NULL; vtx = vtx -> next) {
		vce = vtx -> ce;
		if (vce == (struct cusp_extension *) NULL) continue;
		pprb1 = vce -> prbs[0];
		pprb2 = vce -> prbs[1];
		pprb3 = vce -> prbs[2];
		if (pprb1 == (struct probe *) NULL) continue;
		if (pprb2 == (struct probe *) NULL) continue;
		if (pprb3 == (struct probe *) NULL) continue;
		if (pprb1 == prb1 && pprb2 == prb2 && pprb3 == prb) {old_vtx1 = vtx; }
		if (pprb1 == prb1 && pprb3 == prb2 && pprb2 == prb) {old_vtx2 = vtx; }
		if (pprb2 == prb1 && pprb3 == prb2 && pprb1 == prb) {old_vtx1 = vtx; }
		if (pprb2 == prb1 && pprb1 == prb2 && pprb3 == prb) {old_vtx2 = vtx; }
		if (pprb3 == prb1 && pprb1 == prb2 && pprb2 == prb) {old_vtx1 = vtx; }
		if (pprb3 == prb1 && pprb2 == prb2 && pprb1 == prb) {old_vtx2 = vtx; }
	}
	if (old_vtx1 != NULL) {
		sprintf (message, "%8ld old vertex(1) found, aa = %d",
			old_vtx1 -> number, aa);
		informd (message);
	}
	if (old_vtx2 != NULL) {
		sprintf (message, "%8ld old vertex(2) found, aa = %d",
			old_vtx2 -> number, aa);
		informd (message);
	}
	eatens[0] = 0;
	eatens[1] = 0;
	vip1 = 0;
	vip2 = 0;
	vip3 = 0;
	for (m = 0; m < 2; m++) {
		if (!eatens[m]) {
			dists[m] = find_closest_probe (this_srf,
				xpnt[m], prb1, prb2, prb, &close_prb);
			too_close = (dists[m] < this_srf -> probe_radius +
				EPSILON * EPSILON);
			concave_fac3 = close_prb -> fac;
			vip1 = point_in_pyr (xpnt[m], concave_fac1);
			vip2 = point_in_pyr (xpnt[m], concave_fac2);
			vip3 = point_in_pyr (xpnt[m], concave_fac3);
            sprintf (message,
				 "%8.5f (distance) %8.5f  %8.5f  %8.5f (vertex)",
				 dists[m], xpnt[m][0], xpnt[m][1], xpnt[m][2]);
			if (too_close) {
				informd ("         vertex eaten because too close");
				informd(message);
				eatens[m] = 1;
			}
			else if (!vip1) {
				informd ("         vertex eaten because not in pyramid (vip1)");
				informd(message);
				eatens[m] = 1;
			}
			else if (!vip2) {
				informd ("         vertex eaten because not in pyramid (vip2)");
				informd(message);
				eatens[m] = 1;
			}
		}
	}
	if (eatens[0] != eatens[1]) {	/* problems ! */
		sprintf (message,
		"%8ld %8ld %8ld atoms: %8.3f %8.3f %8.3f probe center: cusp not trimmed",
			atm1 -> number, atm2 -> number, atm3 -> number,
			prb -> center[0], prb -> center[1], prb -> center[2]);
		informd(message);
		this_srf -> n_nontrimmed++;
		/* mark as a notrim face */
		fac -> notrim = TRUE;
		return (0);
	}
	for (m = 0; m < 2; m++) {
		/* add  intersection point to table */
		if (this_srf -> n_nai >= MAX_NAI) {
			ex = new_cept (ARRAY_ERROR,  MSOVERFLOW,  FATAL_SEVERITY);
			add_function (ex, "pre_non_axial");
            add_message (ex, "arc intersection point array");
			add_source (ex, "mscusp.c");
			return (0);
		}
		this_srf -> nais[this_srf -> n_nai].closest_distance = dists[m];
		this_srf -> nais[this_srf -> n_nai].eaten = (short) eatens[m];
		this_srf -> nais[this_srf -> n_nai].fac = fac;
		this_srf -> nais[this_srf -> n_nai].prb = prb;
		this_srf -> nais[this_srf -> n_nai].edg1 = edg1;
		this_srf -> nais[this_srf -> n_nai].edg2 = edg2;
		this_srf -> nais[this_srf -> n_nai].arc1 = arc1;
		this_srf -> nais[this_srf -> n_nai].arc2 = arc2;
		this_srf -> nais[this_srf -> n_nai].cir1 = cir1;
		this_srf -> nais[this_srf -> n_nai].cir2 = cir2;
		this_srf -> nais[this_srf -> n_nai].csp1 = csp1;
		this_srf -> nais[this_srf -> n_nai].csp2 = csp2;
		this_srf -> nais[this_srf -> n_nai].prb1 = prb1;
		this_srf -> nais[this_srf -> n_nai].prb2 = prb2;
		this_srf -> nais[this_srf -> n_nai].nx = (short) aa;
		for (k = 0; k < 3; k++)
			this_srf -> nais[this_srf -> n_nai].xpnt[k] = xpnt[m][k];
		this_srf -> n_nai++;
	}
	return (1);
}

/* non-axial cusp intersection for face */

void one_non_axial (struct surface *this_srf, struct face *concave_fac, int pass)
{
	int k, m, i_nai, eaten1, eaten2;
	int in_arc1, in_arc2, rev1, rev2;
	int orn1, orn2, orn;
	double xpnt[2][3];
	struct probe *prb, *prb1, *prb2;
	struct vertex *vtx1, *vtx2;
	struct arc *arc1, *arc2, *na_arc;
	struct edge *etemp;
	struct edge *edg1, *edg2, *new_edg1, *new_edg2;
	struct circle *cir1, *cir2;
	struct cusp *csp, *csp1, *csp2;
	struct sphere *atm1, *atm2, *atm3;
	struct face *concave_fac1, *concave_fac2;
	struct cusp_extension *vce1, *vce2;
    struct cept *ex;

	prb = concave_fac -> ptr.prb;
	atm1 =  (prb -> atm[0]);
	atm2 =  (prb -> atm[1]);
	atm3 =  (prb -> atm[2]);

	for (i_nai = 0; i_nai < this_srf -> n_nai; i_nai += 2) {
		if (this_srf -> nais[i_nai].fac != concave_fac) continue;
		eaten1 = this_srf -> nais[i_nai].eaten;
		eaten2 = this_srf -> nais[i_nai+1].eaten;
		if (eaten1 && eaten2) continue;
		edg1 = this_srf -> nais[i_nai].edg1;
		edg2 = this_srf -> nais[i_nai].edg2;
		orn1 = edg1 -> orn;
		orn2 = edg2 -> orn;
		arc1 = this_srf -> nais[i_nai].arc1;
		arc2 = this_srf -> nais[i_nai].arc2;
		cir1 = this_srf -> nais[i_nai].cir1;
		cir2 = this_srf -> nais[i_nai].cir2;
		csp1 = this_srf -> nais[i_nai].csp1;
		csp2 = this_srf -> nais[i_nai].csp2;
		prb1 = this_srf -> nais[i_nai].prb1;
		prb2 = this_srf -> nais[i_nai].prb2;
		concave_fac1 = prb1 -> fac;
		concave_fac2 = prb2 -> fac;
		/* assume it is the next one */
		for (m = 0; m < 2; m++) {
			for (k = 0; k < 3; k++)
				xpnt[m][k] = this_srf -> nais[i_nai+m].xpnt[k];
		}
		this_srf -> n_non_axial++;
		csp = find_cusp (prb1, prb2);
		csp -> order = 2;	/* cusp-cusp intersection value */
		/* create two new vertices */
		vtx1 = new_vertex (xpnt[0], NULL, NULL, NULL, NULL);
		if (vtx1 == NULL) return;
		link_vertex (this_srf, vtx1);
		vce1 = (struct cusp_extension *) allocate_cusp_extension ();
		if (error()) {
			add_function (tail_cept, "one_non_axial");
			add_object (tail_cept, CUSP_EXTENSION, "vce1");
			return;
		}
		vtx1 -> ce = vce1;
		vtx2 = new_vertex (xpnt[1], NULL, NULL, NULL, NULL);
		if (vtx2 == NULL) return;
		link_vertex (this_srf, vtx2);
		vce2 = (struct cusp_extension *) allocate_cusp_extension ();
		if (error()) {
			add_function (tail_cept, "one_non_axial");
			add_object (tail_cept, CUSP_EXTENSION, "vce2");
			return;
		}
		vtx2 -> ce = vce2;
		orn = orn1; 	/* (empirical) orientation for middle edge */
		/* set up three probe pointers in each vertex structure */
		vce1 = vtx1 -> ce;
		if (vce1 == (struct cusp_extension *) NULL) continue;
		vce1 -> prbs[0] = prb;
		vce1 -> prbs[1] = prb1;
		vce1 -> prbs[2] = prb2;
		vce2 = vtx2 -> ce;
		if (vce2 == (struct cusp_extension *) NULL) continue;
		vce2 -> prbs[0] = prb;
		vce2 -> prbs[1] = prb1;
		vce2 -> prbs[2] = prb2;

		/* call routine to create new arc */
		in_arc1 = point_in_arc (vtx1 -> center, arc1);
		in_arc2 = point_in_arc (vtx1 -> center, arc2);
		if (!in_arc1 || !in_arc2) {
			informd ("         (one_non_axial): point in arc inconsistency");
			continue;
		}
		in_arc1 = point_in_arc (vtx2 -> center, arc1);
		in_arc2 = point_in_arc (vtx2 -> center, arc2);
		if (!in_arc1 || !in_arc2) {
			informd ("         (one_non_axial): point in arc inconsistency");
			continue;
		}
		na_arc = non_axial_arc (this_srf, vtx1, vtx2, csp);
		
		/* allocate new edges */
		new_edg1 = new_edge (na_arc, orn, NULL, csp);
		new_edg2 = new_edge (na_arc, 1 - orn, NULL, csp);
		/* link into lists for concave faces */
		/* store face pointers in local variables */
		new_edg2 -> next = concave_fac1 -> first_edge;
		concave_fac1 -> first_edge = new_edg2;
		new_edg1 -> next = concave_fac2 -> first_edge;
		concave_fac2 -> first_edge = new_edg1;
		/* handle concave faces on both sides */
		rev1 = split_cusp_edge (this_srf, new_edg1, prb, prb1, edg1);
		if (rev1) {
			informd("         reverse edge orientations");
			new_edg1 -> orn = 1 - new_edg1 -> orn;
			new_edg2 -> orn = 1 - new_edg2 -> orn;
			etemp = na_arc -> edg[0];
			na_arc -> edg[0] = na_arc -> edg[1];
			na_arc -> edg[1] = etemp;
		}
		rev2 = split_cusp_edge (this_srf, new_edg2, prb, prb2, edg2);
		if (rev2) {
			informd("         reverse edge orientations");
			new_edg1 -> orn = 1 - new_edg1 -> orn;
			new_edg2 -> orn = 1 - new_edg2 -> orn;
			etemp = na_arc -> edg[0];
			na_arc -> edg[0] = na_arc -> edg[1];
			na_arc -> edg[1] = etemp;
		}		
		/* make vertices point at each other for old volume code */
		vce1 -> partner = vtx2;
		vce2 -> partner = vtx1;
			
		/* it is not clear when this vertex swapping should be done */
		vertex_resort (this_srf, prb, na_arc, csp);

		/* finished with this face */
		concave_fac -> changed = 1;
		concave_fac1 -> changed = 1;
		concave_fac2 -> changed = 1;
		concave_fac -> ce -> nacirs[pass-1][0] = cir1;
		concave_fac -> ce -> nacirs[pass-1][1] = cir2;
		/* we cannot handle another intersection now */
		return;
	}
}


struct cusp *find_cusp (struct probe *prb1, struct probe *prb2)
{
	int found;
	struct cusp *csp;
	struct cusp_link *cusplink;

	/* look for cusp circle for these two probes */
	found = 0;
	for (cusplink = prb1 -> first_cusp; cusplink != NULL;
		cusplink = cusplink -> next) {
		csp = cusplink -> csp;
		if (csp -> prb[0] == prb1 && csp -> prb[1] == prb2
			|| csp -> prb[0] == prb2 && csp -> prb[1] == prb1) {
			found = 1;
			break;
		}
	}
	if (!found)
		return ((struct cusp *) NULL);
	return (csp);
}

int vertex_resort (struct surface *this_srf, struct probe *prb, struct arc *arcptr, struct cusp *csp)
{
	int fiddle;
	double sum1203, sum1302, midist;
	double mid12[3], mid03[3], mid13[3], mid02[3];
	int mid12_eaten, mid13_eaten, mid03_eaten, mid02_eaten;
	struct arc *fiddle_arc;
	struct arc ta1, ta2;
	struct circle *cir;
	struct vertex *vtx0, *vtx1, *vtx2, *vtx3;
	struct cusp_extension *fvce1, *fvce2, *avce1, *avce2;
	char message[MAXLINE];

	/* check for edge list */
	fiddle = 0;			/* initialization */
	vtx1 = arcptr -> vtx[0];
	vtx2 = arcptr -> vtx[1];
	cir = csp -> cir;
	if (csp -> first_edge != NULL) {
		if (csp -> first_edge -> next != NULL) {	/* too much */
			inform ("(fiddle_with_arcs): cusp for 2nd order cusp not free");
			inform (" obscure case not handled yet");
			this_srf -> n_nontrimmed++;
		}
		else {	/* we'll have to fiddle with this arc */
			fiddle = 1;
			fiddle_arc = csp -> first_edge -> arcptr;
		}
	}
	if (!fiddle) return (0);
	
	/* dangerous fiddling for obscure case */

	vtx0 = fiddle_arc -> vtx[0];
	vtx3 = fiddle_arc -> vtx[1];
	sum1203 = arcptr -> phi + fiddle_arc -> phi;
	middle (arcptr, mid12);
	midist = distance (mid12, prb -> center);
	mid12_eaten = (midist < this_srf -> probe_radius);
	if (mid12_eaten)
		informd ("         (fiddle): middle point of arc eaten by probe");
	middle (fiddle_arc, mid03);
	midist = distance (mid03, prb -> center);
	mid03_eaten = (midist < this_srf -> probe_radius);
	if (mid03_eaten)
		informd ("         (fiddle): middle point of arc eaten by probe");


	store_arc (cir, CONCAVE, vtx1, vtx3, &ta1);
	ta1.phi = arc_ang (&ta1);
	store_arc (cir, CONCAVE, vtx0, vtx2, &ta2);
	ta2.phi = arc_ang (&ta2);
	sum1302 = ta1.phi + ta2.phi;
	middle (&ta1, mid13);
	midist = distance (mid13, prb -> center);
	mid13_eaten = (midist < this_srf -> probe_radius);
	middle (&ta2, mid02);
	midist = distance (mid02, prb -> center);
	mid02_eaten = (midist < this_srf -> probe_radius);
	
	/* require stricly less than sum1203 */
	
	if (sum1302 + EPSILON * EPSILON < sum1203  && !mid13_eaten && !mid02_eaten) {
		informd ("         fiddle: swap ends (23)");
		arcptr -> vtx[0] = vtx1;
		arcptr -> vtx[1] = vtx3;
		fiddle_arc -> vtx[0] = vtx0;
		fiddle_arc -> vtx[1] = vtx2;		
	}
	else {
		/* leave alone - do nothing */
		informd ("         fiddle: do no swap");
	}
	this_srf -> n_cusp_resort++;
	sprintf (message, "%8f %8f %8f center of probe for cusp resort",
		prb -> center[0], prb -> center[1], prb -> center[2]);
	informd(message);
	arcptr -> phi = arc_ang (arcptr); 			/* recompute */
	fiddle_arc -> phi = arc_ang (fiddle_arc); 	/* recompute */
	arcptr -> subtype = RESORT_SUBTYPE;
	fiddle_arc -> subtype = RESORT_SUBTYPE;

	fvce1 = fiddle_arc -> vtx[0] -> ce;
	if (fvce1 == (struct cusp_extension *) NULL) return(0);
	fvce2 = fiddle_arc -> vtx[0] -> ce;
	if (fvce2 == (struct cusp_extension *) NULL) return(0);
	/* check whether new arcs belong to particular tori */
	if (fvce1 -> tor == fvce2 -> tor) {
		fiddle_arc -> tor = fvce1 -> tor;
		fiddle_arc -> perm = 1;
	}
	else {
		fiddle_arc -> tor = NULL;
		fiddle_arc -> perm = 1;
	}
	avce1 = arcptr -> vtx[0] -> ce;
	if (avce1 == (struct cusp_extension *) NULL) return(0);
	avce2 = arcptr -> vtx[0] -> ce;
	if (avce2 == (struct cusp_extension *) NULL) return(0);
	if (avce1 -> tor == avce2 -> tor) {
		arcptr -> tor = avce1 -> tor;
		arcptr -> perm = 1;
	}
	else {
		arcptr -> tor = NULL;
		arcptr -> perm = 1;
	}
	arcptr -> phi = arc_ang (arcptr); 	/* recompute */
	fiddle_arc -> phi = arc_ang (fiddle_arc); 	/* recompute */
	return (1);
}


/* non-axial arc intersection -- new arc */

struct arc *non_axial_arc (struct surface *this_srf, struct vertex *vtx1, struct vertex *vtx2, struct cusp *csp)
{
	struct circle *cir;
	struct arc *arcptr;

	/* circle pointer for this cusp circle */
	cir = csp -> cir;

	/* allocate new arc */
	arcptr = new_arc (cir, vtx1, vtx2, CONCAVE, 0, (double) 0.0, 0L, 0L, 0L);
	this_srf -> n_arc++;
	arcptr -> number = this_srf -> n_arc;
	arcptr -> subtype = NONAXIAL_SUBTYPE;
	arcptr -> perm = 1;		/* new arc permitted to intersect */
	return (arcptr);
}

/* split cusp edge into two cusp edges */

int split_cusp_edge (struct surface *this_srf, struct edge *na_edg, struct probe *prb1, struct probe *prb2, struct edge *edg)
{
	int na_orn, orn, arc_reversed;
	double angle1, angle2, angle_before;
	struct face *concave_fac1, *concave_fac2;
	struct arc *na_arc, *arc1, *arc2;
	struct vertex *vtx1, *vtx2, *vtx3;
	struct circle *cir;
	struct cusp *csp;
	struct edge *edg1, *edg2;
	char message[MAXLINE];
    struct cept *ex;

	arc_reversed = 0;
	
	/* store info in local variables */

	na_orn = na_edg -> orn;
	if (na_orn != 0 && na_orn != 1) {
		ex = new_cept (LOGIC_ERROR,  INVALID_VALUE,  FATAL_SEVERITY);
		add_function (ex, "split_cusp_edge");
		add_source (ex, "mscusp.c");
        add_message (ex, "invalid arc orientation");
		add_long (ex, "arc orientation", (long) na_orn);
		return (0);
	}
	na_arc = na_edg -> arcptr;
	if (na_arc == NULL) {
		ex = new_cept (POINTER_ERROR,  NULL_POINTER,  FATAL_SEVERITY);
        add_object (ex, ARC, "na_arc");
        add_function (ex, "split_cusp_edge");
        add_source (ex, "mscusp.c");
		return (0);
	}
	vtx1 = na_arc -> vtx[na_orn];
	if (vtx1 == NULL) {
		ex = new_cept (POINTER_ERROR,  NULL_POINTER,  FATAL_SEVERITY);
        add_object (ex, VERTEX, "vtx1");
        add_function (ex, "split_cusp_edge");
        add_source (ex, "mscusp.c");
		return (0);
	}
	vtx2 = na_arc -> vtx[1-na_orn];
	if (vtx2 == NULL) {
		ex = new_cept (POINTER_ERROR,  NULL_POINTER,  FATAL_SEVERITY);
        add_object (ex, VERTEX, "vtx2");
        add_function (ex, "split_cusp_edge");
        add_source (ex, "mscusp.c");
		return (0);
	}

	orn = edg -> orn;
	arc1 = edg -> arcptr;
	cir = arc1 -> cir;
	concave_fac1 = prb1 -> fac;
	if (concave_fac1 == NULL) {
		ex = new_cept (POINTER_ERROR,  NULL_POINTER,  FATAL_SEVERITY);
        add_object (ex, FACE, "concave_fac1");
        add_function (ex, "split_cusp_edge");
        add_source (ex, "mscusp.c");
		return (0);
	}
	concave_fac2 = prb2 -> fac;
	if (concave_fac2 == NULL) {
		ex = new_cept (POINTER_ERROR,  NULL_POINTER,  FATAL_SEVERITY);
        add_object (ex, FACE, "concave_fac2");
        add_function (ex, "split_cusp_edge");
        add_source (ex, "mscusp.c");
		return (0);
	}
	csp = edg -> arcptr -> csp;
	if (csp == NULL) {
		ex = new_cept (POINTER_ERROR,  NULL_POINTER,  FATAL_SEVERITY);
        add_object (ex, CUSP, "csp");
        add_function (ex, "split_cusp_edge");
        add_source (ex, "mscusp.c");
		return (0);
	}
	vtx3 = arc1 -> vtx[1-orn];
	if (vtx3 == NULL) {
		ex = new_cept (POINTER_ERROR,  NULL_POINTER,  FATAL_SEVERITY);
        add_object (ex, VERTEX, "vtx3");
        add_function (ex, "split_cusp_edge");
        add_source (ex, "mscusp.c");
		return (0);
	}
	angle_before = arc1 -> phi;


	/* allocate new arc */
	if (orn)
	arc2 = new_arc (cir, vtx3, vtx2, CONCAVE, 0, (double) 0.0, 0L, 0L, 0L);
	else
	arc2 = new_arc (cir, vtx2, vtx3, CONCAVE, 0, (double) 0.0, 0L, 0L, 0L);
	if (arc2 == NULL) return(0);
	if (error()) return (0);
	this_srf -> n_arc++;
	arc2 -> number = this_srf -> n_arc;
	arc2 -> perm = 1;
	arc2 -> subtype = arc1 -> subtype;
			

	/* modify old arc */
	arc1 -> vtx[1-orn] = vtx1;
	arc1 -> perm = 1;
	arc1 -> phi = arc_ang (arc1);		/* recompute */

	/* link new edge into lists for both concave faces */
	edg1 = new_edge (arc2, orn, NULL, csp);
	if (edg1 == NULL) return(0);
	if (error()) return (0);
	edg2 = new_edge (arc2, 1-orn, NULL, csp);
	if (edg2 == NULL) return(0);
	if (error()) return (0);
	edg1 -> next = concave_fac1 -> first_edge;
	concave_fac1 -> first_edge = edg1;
	edg2 -> next = concave_fac2 -> first_edge;
	concave_fac2 -> first_edge = edg2;
	angle1 = arc1 -> phi;
	angle2 = arc2 -> phi;
	if (angle1 + angle2 > angle_before) {
		sprintf(message, "         split_cusp_edge: arc reversed, angles: %8.4f %8.4f %8.4f",
			angle1, angle2, angle_before);
		informd(message);
		if (orn) arc2 -> vtx[1] = vtx1;
		else arc2 -> vtx[0] = vtx1;
		arc1 -> vtx[1-orn] = vtx2;
		arc_reversed = 1;
		arc1 -> phi = arc_ang (arc1);		/* recompute */
		arc2 -> phi = arc_ang (arc2);		/* recompute */
		na_arc -> phi = arc_ang (na_arc);	/* recompute */
		angle1 = arc1 -> phi;
		angle2 = arc2 -> phi;
		sprintf(message, "         split_cusp_edge: angles after: %8.4f %8.4f",
			angle1, angle2);
		informd(message);
	}
	return (arc_reversed);
}

/* create and group cycles of arcs on probe spheres */
void probe_cycles (struct surface *this_srf)
{
	int finished;
	long ncp;
	int nc = 0;
	long number1 = 0L;
	long number2 = 0L;
	long number3 = 0L;
	long number4 = 0L;
	struct probe *prb;
	struct face *fac, *old_tail;
	char message[MAXLINE];
	struct sphere *atm1, *atm2, *atm3, *atm4;

	old_tail = this_srf -> tail_face;

	finished = 0; ncp = 0;
	for (fac = this_srf -> head_face; fac != NULL; fac = fac -> next) {
		if (fac == old_tail) finished = 1;
		if (fac -> shape != CONCAVE) continue;
		prb = fac -> ptr.prb;
		if (!prb -> low) continue;

		/* finished if beyond old last face */
		if (finished && fac != old_tail) break;
		number1 = 0L;
		number2 = 0L;
		number3 = 0L;
		number4 = 0L;
		atm1 = prb -> atm[0];
		atm2 = prb -> atm[1];
		atm3 = prb -> atm[2];
		atm4 = prb -> atm[3];

		/* form probe cycles */
		nc = form_cycles (fac, NULL);
		if (nc < 0) {
			sprintf (message,
				"(probe_cycles) cycle does not close for atoms %5ld %5ld %5ld %5ld",
				number1, number2, number3, number4);
			inform(message);
			/* mark as a problem face */
			fac -> problem = TRUE;
			if (atm1 != NULL) {
				number1 = atm1 -> number;
				atm1 -> problem = TRUE;
			}
			if (atm2 != NULL) {
				number2 = atm2 -> number;
				atm2 -> problem = TRUE;
			}
			if (atm3 != NULL) {
				number3 = atm3 -> number;
				atm3 -> problem = TRUE;
			}
			if (atm4 != NULL) {
				number4 = atm4 -> number;
				atm4 -> problem = TRUE;
			}
		}
		else {
			this_srf -> n_cycle += nc;
			if (error()) return;
		}

		if (!fac -> problem) {
			group_cycles (this_srf, fac, NULL);
			if (error()) return;
		}
		else {
			inform("(probe_cycles): problem with concave face");
			ncp++;
		}

		/* don't look at faces newly created by group_probe_cycles */
		if (finished) break;
	}
	if (ncp > 0) {
		sprintf (message, "%8ld concave faces with bad boundaries", ncp);
		inform(message);
	}
	if (this_srf -> n_nontrimmed > 0) {
		sprintf (message, "%8ld non-trimmed cusps", this_srf -> n_nontrimmed);
		inform(message);
	}
}



double find_closest_probe (struct surface *this_srf, double pnt[3], struct probe *prb1, struct probe *prb2, struct probe *prb3, struct probe **prb4)
{
	int l;
	double closest_dist, dist;
	struct probe *prb;
	struct probe *cprb;

	closest_dist = MS_INFINITY;
	cprb = NULL;
	for (l = 0; l < this_srf -> n_low_probe; l++) {
		/* retrieve pointer to probe */
		prb = *(this_srf -> low_probe_hdl + l);
		if (prb == prb1) continue;
		if (prb == prb2) continue;
		if (prb == prb3) continue;
		dist = distance (prb -> center, pnt);
		if (dist < closest_dist) {
			closest_dist = dist;
			cprb = prb;
		}
	}
	if (cprb != NULL) *prb4 = cprb;
	return (closest_dist);
}


struct cusp *allocate_cusp ()
{
	struct cusp *csp;
    struct cept *ex;
	
	csp = (struct cusp *) allocate_object (CUSP);
	if (csp == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_function (ex, "allocate_cusp");
        add_source (ex, "mscusp.c");
		return (NULL);
	}
	return (csp);
}

struct cusp_link *allocate_cusp_link ()
{
	struct cusp_link *clk;
    struct cept *ex;
	
	clk = (struct cusp_link *) allocate_object (CUSP_LINK);
	if (clk == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_function (ex, "allocate_cusp_link");
        add_source (ex, "mscusp.c");
		return (NULL);
	}
	return (clk);
}


struct cusp_extension *allocate_cusp_extension ()
{
	struct cusp_extension *cex;
    struct cept *ex;
	
	cex = (struct cusp_extension *) allocate_object (CUSP_EXTENSION);
	if (cex == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_function (ex, "allocate_cusp_extension");
        add_source (ex, "mscusp.c");
		return (NULL);
	}
	return (cex);
}


void free_cusp (struct cusp *csp)
{
	free_object (CUSP, (short *) csp);
}

void free_cusp_edges (struct cusp *csp)
{
	struct edge *edg, *nxtedg;

	for (edg = csp -> first_edge; edg != NULL; edg = nxtedg) {
		nxtedg = edg -> next;
		free_edge (edg);
	}
}

void free_cusp_link (struct cusp_link *clk)
{
	free_object (CUSP_LINK, (short *) clk);
}

void free_cusp_extension (struct cusp_extension *cex)
{
	free_object (CUSP_EXTENSION, (short *) cex);
}

void free_cusps (struct surface *srf)
{
	long cusps_freed, extensions_freed, links_freed;
	int l, j, torus_index;
	char message[MAXLINE];
	struct face *fac;
	struct probe *prb;
	struct torus *tor;
	struct cusp *csp, *next_csp;
	struct cusp_link *clk, *next_clk;
	struct cusp_extension *fce, *tce, *vce;
	struct vertex *vtx;
	
	cusps_freed = 0;
	extensions_freed = 0;
	links_freed = 0;

	
	/* search for faces of low probes */
	for (fac = srf -> head_face; fac != NULL; fac = fac -> next) {
		if (fac -> shape != CONCAVE) continue;
		prb = fac -> ptr.prb;
		if (!prb -> low) continue;
		fce = fac -> ce;
		if (fce == NULL) continue;
		free_cusp_extension(fce);
		extensions_freed++;
		fac -> ce = NULL;
	}
	/* consider each low probe */
	for (l = 0; l < srf -> n_low_probe; l++) {
			/* retrieve pointer to probe */
			prb = *(srf -> low_probe_hdl + l);
			next_clk = NULL;
			for (clk = prb -> first_cusp; clk != (struct cusp_link *) NULL;
				clk = next_clk) {
				next_clk = clk -> next;
				free_cusp_link (clk);
				links_freed++;
			}
	}
	for (torus_index = 0; torus_index < srf -> n_point_torus;
		torus_index++) {
		tor = *(srf -> point_torus_hdl + torus_index);
		tce = tor -> ce;
		if (tce == (struct cusp_extension *) NULL) continue;
		/* free vertex cusp extension */
		for (j = 0; j < 2; j++) {
			vtx = tce -> vtx[j];
			if (vtx == NULL) continue;
			vce = vtx -> ce;
			if (vce == NULL) continue;
			for (clk = vce -> first_cusp; clk != (struct cusp_link *) NULL;
				clk = next_clk) {
				next_clk = clk -> next;
				free_cusp_link (clk);
				links_freed++;
			}
			free_cusp_extension(vce);
			vtx -> ce = NULL;
		}
		for (clk = tce -> first_cusp; clk != (struct cusp_link *) NULL;
			clk = next_clk) {
			next_clk = clk -> next;
			free_cusp_link (clk);
			links_freed++;
		}
		free_cusp_extension(tce);
		tor -> ce = NULL;
	}
	next_csp = NULL;
	for (csp = srf -> head_cusp; csp != (struct cusp *) NULL; csp = next_csp) {
		next_csp = csp -> next;
		free_cusp_edges (csp);
		free_cusp (csp);
		cusps_freed++;
	}
	if (srf -> point_torus_hdl != (struct torus **) NULL)
		free_pointers (TORUS, srf -> point_torus_hdl);
	srf -> point_torus_hdl = (struct torus **) NULL;
	if (srf -> low_probe_hdl != (struct probe **) NULL)
		free_pointers (PROBE, srf -> low_probe_hdl);
	srf -> low_probe_hdl = (struct probe **) NULL;
	free_prbdot (srf);
	free_cache (CUSP);
	free_cache (CUSP_LINK);
	free_cache (CUSP_EXTENSION);
	sprintf(message, "%8ld cusps freed", cusps_freed);
	informd(message);
	sprintf(message, "%8ld cusp extensions freed", extensions_freed);
	informd(message);
	sprintf(message, "%8ld cusp links freed", links_freed);
	informd(message);
}

