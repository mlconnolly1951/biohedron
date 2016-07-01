#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

/* Copyright 1995 by Michael L. Connolly */
/* Last revised: December 1, 2001 */

void create_centrals (struct surface *this_srf)
{
	int j;
	double center[3];
	double radius;
	double axis[3];
	struct central *central_ptr;
	struct pair *pair_ptr;

	/* create central circles for non-buried temporary tori */
	this_srf -> n_central = 0;
	if (this_srf -> n_pair > 0) {
		for (pair_ptr = this_srf -> pair_array; pair_ptr <
			this_srf -> pair_array + this_srf -> n_pair;
			pair_ptr++) {
			if (pair_ptr -> buried) continue;
			setup_torus_fields (this_srf -> probe_radius, pair_ptr, center, &radius, axis);
			central_ptr = new_central (center, radius, axis);
			link_central(this_srf, central_ptr);
			/* transfer info from local record to central structure */
			central_ptr -> free = pair_ptr -> free;
			for (j = 0; j < 2; j++)
				central_ptr -> sph[j] = pair_ptr -> sph[j];
			central_ptr -> first_arc = pair_ptr -> first_arc;
		}
	}
}

/* create tori */
void create_tori (struct surface *this_srf)
{
	int k;
	long n_pair;
	double root1, root2, distance12, squared_distance12;
	double radius1, radius2, torus_radius;
	double axis[3];
	double *sphere1_center, *sphere2_center;
	char message[MAXLINE];
	struct sphere *sphere1, *sphere2;
	struct neighbor *first_neighbor, *last_neighbor, *neighbor;
	struct neighbor *first_neighbor2, *last_neighbor2, *neighbor2;
	struct pair *pair_ptr;
    struct cept *ex;

	if (this_srf -> n_pair <= 0) return;

	/* initialize to tmp tor array beginning */
	pair_ptr = this_srf -> pair_array;
	n_pair = 0;

	/* loop through sphere list */
	for (sphere1 = (struct sphere *) (this_srf -> head_atom); sphere1 != NULL;
		sphere1 = sphere1 -> next) {
		if ((first_neighbor = sphere1 -> first_neighbor) == NULL) continue;
		last_neighbor = sphere1 -> last_neighbor;
		/* transfer info to local variables */
		sphere1_center = sphere1 -> center;
		radius1 = sphere1 -> radius;
		/* loop through neighbors of this sphere */
		for (neighbor = first_neighbor; neighbor <= last_neighbor;
			neighbor++) {
			sphere2 = neighbor -> sphptr;
			if (sphere1 >= sphere2) continue; /* no duplication */
			/* transfer info to local variables */
			sphere2_center = sphere2 -> center;
			radius2 = sphere2 -> radius;
			/* geometric computations for torus */
			for (k = 0; k < 3; k++)
				axis[k] = *(sphere2_center + k) - *(sphere1_center + k);
			distance12 = norm (axis);
			if (distance12 <= 0.0) {
				ex = new_cept (GEOMETRY_ERROR,  DEGENERACY,  FATAL_SEVERITY);
				add_function (ex, "create_tori");
				add_source (ex, "mstorus.c");
				add_message (ex, "coincident atoms");
				add_atom (ex, sphere1);
				add_atom (ex, sphere2);
				return;
			}
			squared_distance12 = distance12 * distance12;
			if (squared_distance12 <= 0.0) {
				ex = new_cept (GEOMETRY_ERROR,  DEGENERACY,  FATAL_SEVERITY);
				add_function (ex, "create_tori");
				add_source (ex, "mstorus.c");
				add_message (ex, "coincident atoms");
				add_atom (ex, sphere1);
				add_atom (ex, sphere2);
				return;
			}
			root1 = (radius1 + radius2 + 2 * this_srf -> probe_radius) *
				(radius1 + radius2 + 2 * this_srf -> probe_radius) -
				squared_distance12;
			if (root1 < 0.0) continue; /* sphere too far away */
			root1 = sqrt (root1);
			root2 = squared_distance12 - (radius1 - radius2) * (radius1 - radius2);
			if (root2 < 0.0) continue; /* one sphere inside other */
			root2 = sqrt (root2);
			torus_radius = 0.5 * root1 * root2 / distance12;
			if (torus_radius <= 0.0) continue;
			/* store pointer for torus in first spheres list */
			neighbor -> torptr = pair_ptr;
			/* store pointer for torus in second spheres list */
			first_neighbor2 = sphere2 -> first_neighbor;
			last_neighbor2 = sphere2 -> last_neighbor;
			for (neighbor2 = first_neighbor2;
				neighbor2 <= last_neighbor2; neighbor2++)
				if (neighbor2 -> sphptr == sphere1) {
					neighbor2 -> torptr = pair_ptr;
					break;
				}
			pair_ptr -> free = TRUE;
			pair_ptr -> buried = TRUE;
			pair_ptr -> sph[0] = sphere1;
			pair_ptr -> sph[1] = sphere2;
			pair_ptr++;
			n_pair = pair_ptr - this_srf -> pair_array;
			if (n_pair > this_srf -> n_pair) {
				ex = new_cept (LOGIC_ERROR,  INCONSISTENCY,  FATAL_SEVERITY);
				add_function (ex, "create_tori");
				add_source (ex, "mstorus.c");
				add_long (ex, "before", this_srf -> n_pair);
				add_long (ex, "after", n_pair);
				add_message (ex, "inconsistent number of neighbors");
				return;
			}
		}
	}

	this_srf -> n_pair = n_pair;
	sprintf (message, "%8ld neighbor pairs", n_pair);
	inform (message);
	return;
}

/* create permanent tori */
void create_permanent (struct surface *this_srf)
{
	int k;
	long niter;
	char message[MAXLINE];
	struct arc *arcptr;
	struct torus *torus_ptr;
	struct central *central_ptr;
    struct cept *ex;

	/* create permanent tori for central circles */
	this_srf -> n_tori = 0;
	if (this_srf -> n_pair > 0) {
		for (central_ptr = this_srf -> head_central; central_ptr != NULL;
			central_ptr = central_ptr -> next) {
			torus_ptr = allocate_torus ();
			if (torus_ptr == NULL) {
				ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
				add_object (ex, TORUS, "torus");
				add_function (ex, "create_permanent");
				add_source (ex, "mstorus.c");
				return;
			}
			/* store data in structure */
			for (k = 0; k < 3; k++) {
				 torus_ptr -> center[k] = central_ptr -> center[k];
				 torus_ptr -> axis[k] = central_ptr -> axis[k];
			}
			torus_ptr -> radius = central_ptr -> radius;
			torus_ptr -> atm[0] = central_ptr -> sph[0];
			torus_ptr -> atm[1] = central_ptr -> sph[1];
			torus_ptr -> free = central_ptr -> free;
			torus_ptr -> buried = 0;
			torus_ptr -> first_arc = central_ptr -> first_arc;
			/* set up arc -> torus pointers */
			niter = 0;
			for (arcptr = torus_ptr -> first_arc; arcptr != NULL;
				arcptr = arcptr -> next, niter++) {
				arcptr -> tor = torus_ptr;
				if (niter > 2 * this_srf -> n_atom) {
					add_object (ex, ARC, "arcptr");
					add_function (ex, "create_permanent");
					add_source (ex, "mstorus.c");
					add_message (ex, "infinite loop for torus <--> arc");
					return;
				}
			}
			if (this_srf -> head_torus == NULL) this_srf -> head_torus = torus_ptr;
			else this_srf -> tail_torus -> next = torus_ptr;
			this_srf -> tail_torus = torus_ptr;
			this_srf -> n_tori++;
			torus_ptr -> number = this_srf -> n_tori;
		}
	}
	sprintf (message,"%8ld tori", this_srf -> n_tori); inform(message);
}

void free_centrals (struct surface *this_srf)
{
	struct central *central_ptr, *next_ptr;

	if (this_srf -> n_central > 0) {
		next_ptr = NULL;
		for (central_ptr = this_srf -> head_central; central_ptr != NULL;
			central_ptr = next_ptr) {
			next_ptr = central_ptr -> next;
			free_central (central_ptr);
		}
	}
	this_srf -> head_central = (struct central *) NULL;
	this_srf -> tail_central = (struct central *) NULL;
	this_srf -> n_central = 0;
	free_cache(CENTRAL);
}

void free_pairs (struct surface *this_srf)
{
	if (this_srf -> n_pair > 0 && this_srf -> pair_array != NULL) {
		free_objects (PAIR, (short *) (this_srf -> pair_array));
		this_srf -> pair_array = NULL;
		this_srf -> n_pair = 0;
	}
}

/* CONCAVE FACE ROUTINES */


/* calculate torus center for spheres 1 and 2 */

void compute_torus_center (double probe_radius, struct sphere *sphere1,struct sphere *sphere2, double torus_center[3])
{
	int k;
	double radius1, radius2, asymmetry, distance12_squared;
	char message[MAXLINE];
    struct cept *ex;

	distance12_squared =
		distance_squared (sphere1 -> center, sphere2 -> center);
	if (distance12_squared <= 0.0) {
		ex = new_cept (GEOMETRY_ERROR,  DEGENERACY,  FATAL_SEVERITY);
		add_function (ex, "compute_torus_center");
		add_source (ex, "mstorus.c");
		add_message (ex, "coincident atoms");
		add_atom (ex, sphere1);
		add_atom (ex, sphere2);
		return;
	}
	radius1 = sphere1 -> radius;
	radius2 = sphere2 -> radius;
	asymmetry = (radius1 + probe_radius) * (radius1 + probe_radius) -
		(radius2 + probe_radius) * (radius2 + probe_radius);
	for (k = 0; k < 3; k++)
		torus_center[k] = 0.5 * (sphere1 -> center[k] + sphere2 -> center[k])
			+ 0.5 * ( asymmetry / distance12_squared) *
			(sphere2 -> center[k] - sphere1 -> center[k]);
}

/* setup contents of torus fields given pointer to temporary torus */

void setup_torus_fields (double probe_radius, struct pair *pair_ptr, double torus_center[3], double *return_radius, double torus_axis[3])
{
	int k;
	double root1, root2, asymmetry;
	double distance12, distance12_squared;
	double radius1, radius2, torus_radius;
	double *sphere1_center, *sphere2_center;
	char message[MAXLINE];
	struct sphere *sphere1, *sphere2;
    struct cept *ex;

	sphere1 = pair_ptr -> sph[0];
	sphere2 = pair_ptr -> sph[1];
	if (sphere1 >= sphere2) {
		ex = new_cept (LOGIC_ERROR,  INCONSISTENCY,  FATAL_SEVERITY);
		add_function (ex, "setup_torus_fields");
		add_source (ex, "mstorus.c");
		add_atom (ex, sphere1);
		add_atom (ex, sphere2);
		return;
	}
	/* transfer info to local variables */
	sphere1_center = sphere1 -> center;
	radius1 = sphere1 -> radius;
	sphere2_center = sphere2 -> center;
	radius2 = sphere2 -> radius;
	/* geometric computations for torus */
	for (k = 0; k < 3; k++)
		torus_axis[k] = *(sphere2_center + k) - *(sphere1_center + k);
	distance12 = norm (torus_axis);
	if (distance12 <= 0.0) {
		ex = new_cept (GEOMETRY_ERROR,  DEGENERACY,  FATAL_SEVERITY);
		add_function (ex, "setup_torus_fields");
		add_source (ex, "mstorus.c");
		add_message (ex, "coincident atoms");
		add_atom (ex, sphere1);
		add_atom (ex, sphere2);
		return;
	}
	distance12_squared = distance12 * distance12;
	if (distance12_squared <= 0.0) {
		ex = new_cept (GEOMETRY_ERROR,  DEGENERACY,  FATAL_SEVERITY);
		add_function (ex, "setup_torus_fields");
		add_source (ex, "mstorus.c");
		add_message (ex, "coincident atoms");
		add_atom (ex, sphere1);
		add_atom (ex, sphere2);
		return;
	}
	for (k = 0; k < 3; k++) torus_axis[k] /= distance12;
	root1 = (radius1 + radius2 + 2 * probe_radius) *
		(radius1 + radius2 + 2 * probe_radius) - distance12_squared;
	if (root1 < 0.0) {
		ex = new_cept (LOGIC_ERROR,  INCONSISTENCY,  FATAL_SEVERITY);
		add_function (ex, "setup_torus_fields");
		add_source (ex, "mstorus.c");
		add_message (ex, "inconsistent existence of torus");
		add_atom (ex, sphere1);
		add_atom (ex, sphere2);
		add_double (ex, "root1", root1);
		return;
	}
	root1 = sqrt (root1);
	root2 = distance12_squared - (radius1 - radius2) *
		(radius1 - radius2);
	if (root2 < 0.0) {
		ex = new_cept (LOGIC_ERROR,  INCONSISTENCY,  FATAL_SEVERITY);
		add_function (ex, "setup_torus_fields");
		add_source (ex, "mstorus.c");
		add_message (ex, "inconsistent existence of torus");
		add_atom (ex, sphere1);
		add_atom (ex, sphere2);
		add_double (ex, "root2", root2);
		return;
	}
	root2 = sqrt (root2);
	torus_radius = 0.5 * root1 * root2 / distance12;
	if (torus_radius <= 0) {
		ex = new_cept (LOGIC_ERROR,  INCONSISTENCY,  FATAL_SEVERITY);
		add_function (ex, "setup_torus_fields");
		add_source (ex, "mstorus.c");
		add_message (ex, "inconsistent existence of torus");
		add_atom (ex, sphere1);
		add_atom (ex, sphere2);
		add_double (ex, "torus_radius", torus_radius);
		return;
	}
	asymmetry = (radius1 + probe_radius) * (radius1 + probe_radius) -
		(radius2 + probe_radius) * (radius2 + probe_radius);
	for (k = 0; k < 3; k++)
		torus_center[k] = 0.5 * (*(sphere1_center + k) +
			*(sphere2_center + k))
		+ 0.5 * torus_axis[k] * asymmetry / distance12;

	*return_radius = torus_radius;
	return;
}
