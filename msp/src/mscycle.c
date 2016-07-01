#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

/* Copyright 1995 by Michael L. Connolly */
/* Last revised: December 14, 2001  */


/* form preliminary convex face for atom */
struct face *pre_form (struct surface *this_srf, struct sphere *atm)
{
	long n_atom_arcs;
	struct face *convex_fac;
	struct arc *arcptr;
	struct edge *edge_ptr;

	/* new convex face */
	convex_fac = new_face (NULL, CONVEX);
	if (convex_fac == NULL) return (NULL);
	convex_fac -> chi = 2;
	convex_fac -> ptr.atm = (struct sphere *) atm;
	convex_fac -> n_arc = 0;
	convex_fac -> first_cycle = NULL;
	convex_fac -> n_cycle = 0;
	convex_fac -> first_edge = NULL;
	link_face (this_srf, convex_fac);

	/* completely accessible atom */
	if (atm -> first_arc == NULL) return (convex_fac);

	/* count the number of arcs */
	n_atom_arcs = 0;
	for (arcptr = atm -> first_arc; arcptr != NULL;
		arcptr = arcptr -> next)
		n_atom_arcs++;
	convex_fac -> n_arc = (short) n_atom_arcs;

	/* create edges for arcs */
	for (arcptr = atm -> first_arc; arcptr != NULL;
		arcptr = arcptr -> next) {
		edge_ptr = new_edge (arcptr, 0, convex_fac, NULL);
		if (edge_ptr == NULL) return (NULL);
		edge_ptr -> next = convex_fac -> first_edge;
		convex_fac -> first_edge = edge_ptr;
	}
	return (convex_fac);
}

int form_cycles (struct face *given_fac, struct solid_angle *sa)
{
	long fac_nedge;
	long edge_index, n_face_cycles;
	long cyc_nedge = 0;
	long cycle1_index;
	struct edge *head_edge;
	struct vertex *vtx;
	struct edge *edg0, *edg1, *edg2, *previous_edg, *starting_edg;
	struct cycle *cyc, *next_cycle;
	struct cycle **cycle_list = NULL;
	struct edge **edge_list = NULL;
    struct cept *ex;

	fac_nedge = 0;
	if (given_fac != NULL) {
		if (given_fac -> first_edge == NULL) return(0);
		if (given_fac -> problem) return(0);
		given_fac -> chi = 2;			/* initialize Euler characteristic */
		/* count edges */
		for (edg0 = given_fac -> first_edge; edg0 != NULL; edg0 = edg0 -> next)
			fac_nedge++;
		given_fac -> n_arc = (short) fac_nedge;
	}
	else if (sa != NULL) {
		fac_nedge = sa -> n_edge;
		if (fac_nedge == 0) return(0);
	}

	/* allocate memory */
	edge_list = (struct edge **)
		allocate_pointers (EDGE, fac_nedge);
	if (edge_list == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_object (ex, EDGE, "edge_list");
        add_function (ex, "form_cycles");
        add_source (ex, "mscycle.c");
		return(0);
	}
	cycle_list = (struct cycle **) allocate_pointers (CYCLE, fac_nedge);
	if (cycle_list == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_object (ex, CYCLE, "cycle_list");
        add_function (ex, "form_cycles");
        add_source (ex, "mscycle.c");
		return(0);
	}

	/* store pointers to edge list entries */
	if (given_fac != NULL) head_edge = given_fac -> first_edge;
	else if (sa != NULL) head_edge = sa -> head_edge;
	else inform("form_cycles: two nulls");
	for (edg0 = head_edge, edge_index = 0; edge_index < fac_nedge;
		edg0 = edg0 -> next, edge_index++)
		*(edge_list + edge_index) = edg0;

	n_face_cycles = 0;			/* initialize number of cycles */
	for (;;) {
	/* look for starting edge */
		starting_edg = NULL;
		for (edge_index = 0; edge_index < fac_nedge; edge_index++)
			if ((edg2 = *(edge_list+edge_index)) != NULL) {
				starting_edg = edg2;
				*(edge_list + edge_index) = NULL;
				break;
			}
		if (starting_edg == NULL) break;	/* finished forming cycles */
		/* start new cycle */
		cyc = new_cycle (given_fac, starting_edg); if (error()) return(0);
		/* one edge in cycle at this time */
		if (sa != NULL) cyc_nedge = 1;
		/* store pointer to cycle */
		*(cycle_list + n_face_cycles++) = cyc;
		if (sa != NULL) {
			sa -> n_cycle++;
			if (sa -> head_cycle == NULL) sa -> head_cycle = cyc;
			else sa -> tail_cycle -> next = cyc;
			sa -> tail_cycle = cyc;
		}
		if (starting_edg -> arcptr -> vtx[0] == NULL) {
			/* one-edge cycle */
			if (given_fac != NULL) given_fac -> chi--;
			starting_edg -> next = NULL;
			if (sa != NULL)
				cyc -> n_edge = cyc_nedge;
			continue;
		}

		/* multi-edge cycle */
		previous_edg = starting_edg;

		/* loop last until starting edge is reached */
		while (previous_edg -> arcptr -> vtx[1 - previous_edg -> orn]
			!= starting_edg -> arcptr -> vtx[starting_edg -> orn]) {
			/* get vertex */
			vtx = previous_edg -> arcptr ->
				vtx[1 - previous_edg -> orn];
			/* look for next edge */
			edg1 = NULL;
			for (edge_index = 0; edge_index < fac_nedge; edge_index++) {
				/* skip used edges (marked null) */
				if ((edg2 = *(edge_list + edge_index)) == NULL) continue;
				/* edge connectivity found via vertices */
				if (edg2 -> arcptr -> vtx[edg2 -> orn] == vtx) {
					edg1 = edg2;
					/* clear edge from list */
					*(edge_list + edge_index) = NULL;
					break;
				}
			}
			/* if no next edge, we've got problems */
			if (edg1 == NULL) {
				ex = new_cept (GEOMETRY_ERROR,  INCONSISTENCY,  FATAL_SEVERITY);
				add_function (ex, "form_cycles");
				add_source (ex, "mscycle.c");
				add_message(ex, "cycle does not close");
				/* return memory to free list */
				if (edge_list != NULL) {
					free_pointers (EDGE, edge_list);
					edge_list = NULL;
				}
				if (cycle_list != NULL) {
					free_pointers (CYCLE, cycle_list);
					cycle_list = NULL;
				}
				return(-1);
			}
			previous_edg -> next = edg1;		/* set up link */
			previous_edg = edg1;
			/* increment number of edges in cycle */
			if (sa != NULL)
				cyc_nedge++;
		}
		if (given_fac != NULL) given_fac -> chi--;			/* decrement Euler characteristic */
		previous_edg -> next = NULL;		/* mark end of cycle */
		/* store number of edges in cycle */
		if (sa != NULL)
			cyc -> n_edge = cyc_nedge;
	}

	if (given_fac != NULL) {
		/* set up linked list of cycles */
		for (cycle1_index = 0; cycle1_index < n_face_cycles;
			cycle1_index++) {
			/* pointer to the cycle */
			cyc = *(cycle_list + cycle1_index);
			/* pointer to next cycle */
			if (cycle1_index == n_face_cycles - 1) next_cycle = NULL;
			else next_cycle = *(cycle_list + cycle1_index + 1);
			/* link */
			cyc -> next = next_cycle;
		}
	
		given_fac -> first_cycle = *cycle_list;	/* face points to first cycle */
		given_fac -> n_cycle = (short) n_face_cycles;	/* store number of cycles */
	}
	/* return memory to free list */
	if (edge_list != NULL)
		free_pointers (EDGE, edge_list);
	if (cycle_list != NULL)
		free_pointers (CYCLE, cycle_list);
	edge_list = NULL;
	cycle_list = NULL;
	return (n_face_cycles);
}

void group_cycles (struct surface *this_srf, struct face *given_fac, struct solid_angle *sa)
{
	long *contains = NULL;
	long *bound_same = NULL; 
	long *cycle_used = NULL;
	long n_face_cycles, n_squared;
	long index1, index2, index3;
	long index12, index21, index23, index32, index13, index31;
	long contain12, contain21, overlap, n_do, n_dc;
	long n_f, n_face_edges;
	long n_one_edge_cycles, result;
	long number1, number2, number3, number4;
	char message[MAXLINE];
	struct cycle *cyc, *previous_cycle, *head_cycle;
	struct cycle *cyc1, *cyc2;
	struct cycle **cycle_list = NULL;
	struct cycle **cycle_hdl = NULL;
	struct circle *cir1, *cir2;
	struct arc *arc1, *arc2;
	struct sphere *atm1, *atm2, *atm3, *atm4;
	struct face *fac;
	struct probe *prb;
	struct variety *vty;
    struct cept *ex;

	/* check for simple case */
	if (given_fac != NULL) {
		n_face_cycles = given_fac -> n_cycle;
		if (n_face_cycles <= 1) return;
		if (given_fac -> problem) return;
	}
	else if (sa != NULL) {
		n_face_cycles = sa -> n_cycle;
		vty = sa -> vty;
		if (n_face_cycles <= 1) {
			/* create first face */
			fac = new_face (vty, CONVEX);
			add2solid (sa, fac);
			fac -> chi = 2 - n_face_cycles;
			fac -> first_cycle = sa -> head_cycle;
			fac -> n_cycle = (short) n_face_cycles;
			return;
		}
	}
	/* count one-edge cycles */
	n_one_edge_cycles = 0;
	if (given_fac != NULL) {
		for (cyc = given_fac -> first_cycle; cyc != NULL; cyc = cyc -> next)
			if (edges_in_cycle (cyc) == 1) n_one_edge_cycles++;
		}

	/* store atom or probe pointer in local variable */
	atm1 = NULL;
	atm2 = NULL;
	atm3 = NULL;
	atm4 = NULL;
	number1 = 0;
	number2 = 0;
	number3 = 0;
	number4 = 0;
	prb = NULL;
	if (given_fac != NULL) {
		if (given_fac -> shape == CONVEX) {
			atm1 = given_fac -> ptr.atm;
			number1 = atm1 -> number;
		}
		else if (given_fac -> shape == CONCAVE) {
			prb = given_fac -> ptr.prb;
			atm1 = prb -> atm[0];
			number1 = atm1 -> number;
			atm2 = prb -> atm[1];
			number2 = atm2 -> number;
			atm3 = prb -> atm[2];
			number3 = atm3 -> number;
			atm4 = prb -> atm[3];
			if (atm4 != NULL)
				number4 = atm4 -> number;
		}
	}

	if (given_fac != NULL && given_fac -> shape == CONCAVE) {
		/* check for case we cannot haddle */
		if (n_one_edge_cycles > 1) {
			/* make sure nothing weird is going on */
			if (given_fac -> n_cycle != 1 + n_one_edge_cycles) {
				inform ("(group_cycles): concave face too complicated 1+#of1-ecyc!= #cyc");
				given_fac -> problem = TRUE;
				atm1 -> problem = TRUE;
				atm2 -> problem = TRUE;
				atm3 -> problem = TRUE;
				if (atm4 != NULL)
					atm4 -> problem = TRUE;
				return;
			}
			/* very special case: one multiply connected face */
			n_do = 0; n_dc = 0;
			for (cyc1 = given_fac -> first_cycle; cyc1 != NULL; cyc1 = cyc1 -> next) {
				if (edges_in_cycle (cyc1) != 1) continue;
				arc1 = cyc1 -> first_edge -> arcptr; cir1 = arc1 -> cir;
				for (cyc2 = cyc1; cyc2 != NULL; cyc2 = cyc2 -> next) {
					if (cyc1 == cyc2) continue;
					if (edges_in_cycle (cyc2) != 1) continue;
					arc2 = cyc2 -> first_edge -> arcptr; cir2 = arc2 -> cir;
					/* call disk_contain */
					contain12 = disk_contain (prb -> center, this_srf -> probe_radius, cir1 -> center, cir2 -> center,
						cir1 -> radius, cir2 -> radius);
					if (contain12) {
						delete_cycle (this_srf, given_fac, cyc2);
						n_one_edge_cycles--;
						inform ("delete contained disk");
						n_dc++; break;
					}
					/* call disk_contain */
					contain21 = disk_contain (prb -> center, this_srf -> probe_radius,
						cir2 -> center, cir1 -> center, cir2 -> radius, cir1 -> radius);
					if (contain21) {
						delete_cycle (this_srf, given_fac, cyc1);
						n_one_edge_cycles--;
						inform ("delete contained disk");
						n_dc++; break;
					}
					/* call disk_overlap */
					overlap = disk_overlap (prb -> center, this_srf -> probe_radius,
						cir1 -> center, cir2 -> center, cir1 -> radius, cir2 -> radius);
					if (overlap) {
						n_do++;
						sprintf (message,"arcs overlap: %12ld %1d %12ld %1d",
						(long) arc1, arc1 ->perm, (long) arc2, arc2 ->perm);
						informd(message);
					}
				}
			}
			if (n_do > 0 || n_dc > 0) {
				if (n_do > 0)
					inform ("(group_cycles) disks overlap");
				if (n_dc > 0)
					inform ("(group_cycles) disk contained");
				inform ("will attempt to recover from error");
				given_fac -> problem = TRUE;
				atm1 -> problem = TRUE;
				atm2 -> problem = TRUE;
				atm3 -> problem = TRUE;
				if (atm4 != NULL)
					atm4 -> problem = TRUE;
				sprintf (message,"concave face (%12ld), atoms %5ld %5ld %5ld %5ld",
					(long) given_fac, number1, number2, number3, number4);
				inform(message);
				return;
			}
			/* hope we're okay */
			return;
		}
		else if (n_one_edge_cycles == 1) {
			/* one one-edge cycle: */
			/* make sure nothing weird is going on */
			if (given_fac -> n_cycle != 2) {
				inform ("(group_cycles): # of 1-e cyc = 1 and #cyc != 2");
				given_fac -> problem = TRUE;
				atm1 -> problem = TRUE;
				atm2 -> problem = TRUE;
				atm3 -> problem = TRUE;
				atm4 -> problem = TRUE;
				return;
			}
			/* special case: one multiply connected face */
			/* nothing to do */
			return;
		}
	}
	/* more than one cycle */
	n_squared = n_face_cycles * n_face_cycles;
	if (given_fac != NULL && n_squared <= 1) return;
	contains = allocate_longs (n_squared, 0, CONTAINS);
	if (contains == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_variable (ex, CONTAINS, "contains");
        add_function (ex, "group_cycles");
        add_source (ex, "mscycle.c");
		return;
	}
	if (given_fac != NULL)
		head_cycle = given_fac -> first_cycle;
	else if (sa != NULL)
		head_cycle = sa -> head_cycle;
	for (index1 = 0, cyc1 = head_cycle; index1 < n_face_cycles;
		index1++, cyc1 = cyc1 -> next) {
		for (index2 = 0, cyc2 = head_cycle; index2 < n_face_cycles;
			index2++, cyc2 = cyc2 -> next) {
			if (given_fac != NULL) {
				if (given_fac -> shape == CONVEX) {
					result = contain (cyc1, cyc2, atm1 -> center, atm1 -> radius);
					if (result < 0) {
						inform("problem face due to contain returns negative");
						given_fac -> problem = TRUE;
						atm1 -> problem = TRUE;
						if (atm2 != NULL) atm2 -> problem = TRUE;
						if (atm3 != NULL) atm3 -> problem = TRUE;
						if (atm4 != NULL) atm4 -> problem = TRUE;
						if (contains != NULL)
							free_longs (contains, 0, CONTAINS);
						contains = NULL;
						return;
					}
				}
				else {
					result = (cyc1 == cyc2);
					/* for concave faces, assume each cycle is a separate face */
				}
			}
			else if (sa != NULL) {
				result = contain (cyc1, cyc2, vty -> center, vty -> radii[0]);
				if (error()) return;
			}
			*(contains + index1 * n_face_cycles + index2) = result;
		}
	}
	bound_same = allocate_longs (n_squared, 0, BOUND_SAME);
	if (bound_same == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_variable (ex, BOUND_SAME, "bound_same");
        add_function (ex, "group_cycles");
        add_source (ex, "mscycle.c");
		return;
	}

	/* first guess */
	for (index1 = 0; index1 < n_face_cycles;
		index1++)
		for (index2 = 0; index2 < n_face_cycles; index2++) {
			index12 = index1 * n_face_cycles + index2;
			index21 = index2 * n_face_cycles + index1;
			*(bound_same + index12) =
				(*(contains + index12) && *(contains + index21));
		}

	/* check for special case of two cycles bounding one face */
	if (n_face_cycles == 2)
		if (*(bound_same + 1)) {
			if (contains != NULL)
				free_longs (contains, 0, CONTAINS);
			if (bound_same != NULL)
				free_longs (bound_same, 0, BOUND_SAME);
			contains = NULL;
			bound_same = NULL;
			if (sa != NULL) {
				/* create first face */
				fac = new_face (vty, CONVEX);
				add2solid (sa, fac);
				fac -> chi = 2 - n_face_cycles;
				fac -> first_cycle = sa -> head_cycle;
				fac -> n_cycle = (short) n_face_cycles;
			}
			return;
		}

	/* correction */
	for (index1 = 0; index1 < n_face_cycles - 1; index1++)
		for (index2 = index1 + 1; index2 < n_face_cycles; index2++) {
			index12 = index1 * n_face_cycles + index2;
			index21 = index2 * n_face_cycles + index1;
			for (index3 = 0; index3 < n_face_cycles;
				index3++) {
				if (index3 == index1) continue;
				if (index3 == index2) continue;
				index13 = index1 * n_face_cycles + index3;
				index31 = index3 * n_face_cycles + index1;
				index23 = index2 * n_face_cycles + index3;
				index32 = index3 * n_face_cycles + index2;
				if (!*(contains + index13)) continue;
				if (!*(contains + index23)) continue;
				if (*(contains + index31) &&
					*(contains + index32)) continue;
				*(bound_same + index12) = 0;
				*(bound_same + index21) = 0;
			}
		}


	/* form new faces */

	/* list of pointers to cycles */
	cycle_list = (struct cycle **)
		allocate_pointers (CYCLE, n_face_cycles);
	if (cycle_list == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_object (ex, CYCLE, "cycle_list");
        add_function (ex, "group_cycles");
        add_source (ex, "mscycle.c");
		return;
	}
	cycle_hdl = cycle_list;
	if (given_fac != NULL)
		head_cycle = given_fac -> first_cycle;
	else if (sa != NULL)
		head_cycle = sa -> head_cycle;
	for (cyc = head_cycle; cyc != NULL; cyc = cyc -> next)
		*cycle_hdl++ = cyc;
	/* array to mark cycles used for new faces */
	cycle_used = allocate_longs (n_face_cycles, 0, CYCLE_USED);
	if (cycle_used == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_variable (ex, CYCLE_USED, "cycle_used");
        add_function (ex, "group_cycles");
        add_source (ex, "mscycle.c");
		return;
	}

	n_f = 0;
	for (index1 = 0; index1 < n_face_cycles; index1++) {
		cyc1 = *(cycle_list + index1);
		if (*(cycle_used + index1)) continue;
		*(cycle_used + index1) = 1;
		cyc1 -> next = NULL;
		/* start new face */
		if (given_fac != NULL) {
			if (n_f > 0) {
				fac = new_face (NULL, (int) given_fac -> shape);
				if (fac == NULL) return;
				link_face (this_srf, fac);
			}
			else fac = given_fac;
		}
		else if (sa != NULL) {
			fac = new_face (vty, CONVEX);
			add2solid (sa, fac);
		}
		fac -> chi = 1;
		if (given_fac != NULL) {
			n_face_edges = edges_in_cycle (cyc1);
			if (given_fac -> shape == CONVEX)
				fac -> ptr.atm = atm1;
			else if (given_fac -> shape == CONCAVE)
				fac -> ptr.prb = prb;
		}
		else if (sa != NULL) {
			n_face_edges = cyc1 -> n_edge;
		}
		n_f++;
		fac -> first_cycle = cyc1;
		previous_cycle = cyc1;
		/* look for cycles bounding same face */
		for (index2 = 0; index2 < n_face_cycles; index2++) {
			if (*(cycle_used + index2)) continue;
			index12 = index1 * n_face_cycles + index2;
			if (!*(bound_same + index12)) continue;
			*(cycle_used + index2) = 1;
			cyc2 = *(cycle_list + index2);
			previous_cycle -> next = cyc2;
			previous_cycle = cyc2;
			cyc2 -> next = NULL;
			--fac -> chi;
			if (given_fac != NULL)
				n_face_edges += edges_in_cycle (cyc2);
			else if (sa != NULL)
				n_face_edges += cyc2 -> n_edge;
		}
		if (given_fac != NULL) {
			fac -> n_arc = (short) n_face_edges;
			fac -> n_cycle = (short) (2 - fac -> chi);
		}
		else if (sa != NULL) {
			fac -> n_cycle = 2 - fac -> chi;
		}
	}

	if (cycle_used != NULL)
		free_longs (cycle_used, 0, CYCLE_USED);
	if (cycle_list != NULL)
		free_pointers (CYCLE, cycle_list);
	if (bound_same != NULL)
		free_longs (bound_same, 0, BOUND_SAME);
	if (contains != NULL)
		free_longs (contains, 0, CONTAINS);
	cycle_used = NULL;
	cycle_list = NULL;
	bound_same = NULL;
	contains = NULL;
}


/* does cycle a contain cycle b */

int contain (struct cycle *cyc1, struct cycle *cyc2, double cen[3], double rad)
{
	int k;
	long  n_cyc1_edges, vtx_idx, vtx_idx2;
	double south_component, extension, edge_axis_sign;
	double rotation_angle, side_length;
	double vector[3], point_on_cyc2[3], southward_vector[3];
	double *axis_on_cyc2;
	double *projected_vertex_list, *polygon_tangent, *exterior_angle;
	struct edge *edg1, *edg2;
	struct vertex *vtx;
    struct cept *ex;

	n_cyc1_edges = 0;
	if (cyc1 == cyc2) return (1);
	for (edg1 = cyc1 -> first_edge; edg1 != NULL; edg1 = edg1 -> next)
		n_cyc1_edges++;
	if (n_cyc1_edges <= 2) return (1);
	for (edg1 = cyc1 -> first_edge; edg1 != NULL; edg1 = edg1 -> next)
		for (edg2 = cyc2 -> first_edge; edg2 != NULL; edg2 = edg2 -> next)
			if (edg1 -> arcptr -> cir == edg2 -> arcptr -> cir)
				return (0);
	/* find a point on cyc2 */
	vtx = NULL;
	for (edg2 = cyc2 -> first_edge; edg2 != NULL; edg2 = edg2 -> next) {
		if (edg2 -> arcptr -> vtx[edg2 -> orn] != NULL) {
			vtx = edg2 -> arcptr -> vtx[edg2 -> orn];
			break;
		}
		axis_on_cyc2 = edg2 -> arcptr -> cir -> axis;
		edge_axis_sign = ((edg2 -> orn == 0) ? 1.0 : -1.0);
	}
	if (vtx != NULL)
		for (k = 0; k < 3; k++)
			point_on_cyc2[k] = vtx -> center[k];
		else {
		/* move outward to sphere */
		for (k = 0; k < 3; k++)
			point_on_cyc2[k] = cen[k] - rad * edge_axis_sign *
				(*(axis_on_cyc2 + k));
	}
	/* stereographic projection */
	for (k = 0; k < 3; k++)
		southward_vector[k] = (cen[k] - point_on_cyc2[k]) / rad;
	/* allocate memory for projected vertices of a */
	projected_vertex_list = allocate_doubles (n_cyc1_edges * 3, 0, PROJECTED_VERTEX_LIST);
	if (projected_vertex_list == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_variable (ex, PROJECTED_VERTEX_LIST, "projected_vertex_list");
        add_function (ex, "contain");
        add_source (ex, "mscycle.c");
		return(0);
	}
	polygon_tangent  = allocate_doubles (n_cyc1_edges * 3, 0, POLYGON_TANGENT);
	if (polygon_tangent == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_variable (ex, POLYGON_TANGENT, "polygon_tangent");
        add_function (ex, "contain");
        add_source (ex, "mscycle.c");
		return(0);
	}
	exterior_angle = allocate_doubles (n_cyc1_edges, 0, EXTERIOR_ANGLE);
	if (exterior_angle == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_variable (ex, EXTERIOR_ANGLE, "exterior_angle");
        add_function (ex, "contain");
        add_source (ex, "mscycle.c");
		return(0);
	}
	/* project vertices */
	vtx_idx = 0;
	for (edg1 = cyc1 -> first_edge; edg1 != NULL; edg1 = edg1 -> next) {
		vtx = edg1 -> arcptr -> vtx[edg1 -> orn];
		for (k = 0; k < 3; k++)
			vector[k] = vtx -> center[k] - point_on_cyc2[k];
		south_component = dot_product (vector, southward_vector);
		if (south_component <= 0.0) {
			ex = new_cept (GEOMETRY_ERROR,  INVALID_VALUE,  FATAL_SEVERITY);
			add_function (ex, "contain");
			add_source (ex, "mscycle.c");
			add_message(ex, "bad projection of vertex");
			if (projected_vertex_list != NULL)
				free_doubles (projected_vertex_list, 0, PROJECTED_VERTEX_LIST);
			if (polygon_tangent != NULL)
				free_doubles (polygon_tangent, 0, POLYGON_TANGENT);
			if (exterior_angle != NULL)
				free_doubles (exterior_angle, 0, EXTERIOR_ANGLE);
			projected_vertex_list = NULL;
			polygon_tangent = NULL;
			exterior_angle = NULL;
			return(0);
		}
		extension = 2 * rad / south_component;
		for (k = 0; k < 3; k++)
			*(projected_vertex_list + 3 * vtx_idx + k) =
				point_on_cyc2[k] + extension * vector[k];
		vtx_idx++;
	}
	/* compute side vectors */
	for (vtx_idx = 0; vtx_idx < n_cyc1_edges; vtx_idx++) {
		vtx_idx2 = vtx_idx + 1;
		if (vtx_idx2 >= n_cyc1_edges) vtx_idx2 = 0;
		for (k = 0; k < 3; k++)
			*(polygon_tangent + 3 * vtx_idx + k) =
				  *(projected_vertex_list + 3 * vtx_idx2 + k)
				- *(projected_vertex_list + 3 * vtx_idx + k);
		side_length = norm (polygon_tangent + 3 * vtx_idx);
		if (side_length <= 0.0) {
			ex = new_cept (GEOMETRY_ERROR,  INVALID_VALUE,  FATAL_SEVERITY);
			add_function (ex, "contain");
			add_source (ex, "mscycle.c");
			add_message(ex, "bad polygon side");
			if (projected_vertex_list != NULL)
				free_doubles (projected_vertex_list, 0, PROJECTED_VERTEX_LIST);
			if (polygon_tangent != NULL)
				free_doubles (polygon_tangent, 0, POLYGON_TANGENT);
			if (exterior_angle != NULL)
				free_doubles (exterior_angle, 0, EXTERIOR_ANGLE);
			projected_vertex_list = NULL;
			polygon_tangent = NULL;
			exterior_angle = NULL;
			return(0);
		}
		for (k = 0; k < 3; k++)
			*(polygon_tangent + 3 * vtx_idx + k) /= side_length;
	}
	/* compute angles */
	for (vtx_idx = 0; vtx_idx < n_cyc1_edges; vtx_idx++) {
		vtx_idx2 = vtx_idx - 1;
		if (vtx_idx2 < 0) vtx_idx2 = n_cyc1_edges - 1;
		*(exterior_angle + vtx_idx) =
			odd_angle (polygon_tangent + 3 * vtx_idx2,
			polygon_tangent + 3 * vtx_idx, southward_vector, (double) (-1.0));
	}
	/* calculate rotation angle */
	rotation_angle = 0.0;
	for (vtx_idx = 0; vtx_idx < n_cyc1_edges; vtx_idx++)
		rotation_angle += *(exterior_angle + vtx_idx);
	if (projected_vertex_list != NULL)
		free_doubles (projected_vertex_list, 0, PROJECTED_VERTEX_LIST);
	if (polygon_tangent != NULL)
		free_doubles (polygon_tangent, 0, POLYGON_TANGENT);
	if (exterior_angle != NULL)
		free_doubles (exterior_angle, 0, EXTERIOR_ANGLE);
	projected_vertex_list = NULL;
	polygon_tangent = NULL;
	exterior_angle = NULL;
	return ((int) (rotation_angle > 0.0));
}


int disk_overlap (double sphere_center[3], double rad, double disk1_center[3], double disk2_center[3], double rad1, double rad2)
{
	int k;
	double angle1, angle2;
	double axis1[3], axis2[3];
	double angle12;
	double dot12;
	char message[MAXLINE];

	angle1 = asin (rad1 / rad);
	angle2 = asin (rad2 / rad);
	for (k = 0; k < 3; k++)
		axis1[k] =  (disk1_center[k] - sphere_center[k]);
	normalize (axis1);
	for (k = 0; k < 3; k++)
		axis2[k] =  (disk2_center[k] - sphere_center[k]);
	normalize (axis2);
	dot12 = dot_product (axis1, axis2);
	if (dot12 < -1.0) dot12 = -1.0;
	else if (dot12 > 1.0) dot12 = 1.0;
	angle12 = acos (dot12);
	sprintf (message,"angle1 = %8.3f; angle2 = %8.3f; angle12 = %8.3f",
			angle1, angle2, angle12);
	informd(message);
	return (angle1 + angle2 >= angle12);
}

int disk_contain (double sphere_center[3], double rad, double disk1_center[3], double disk2_center[3], double rad1, double rad2)
{
	int k;
	double angle1, angle2;
	double axis1[3], axis2[3];
	double angle12;
	double dot12;
	char message[MAXLINE];

	angle1 = asin (rad1 / rad);
	angle2 = asin (rad2 / rad);
	for (k = 0; k < 3; k++)
		axis1[k] =  (disk1_center[k] - sphere_center[k]);
	normalize (axis1);
	for (k = 0; k < 3; k++)
		axis2[k] =  (disk2_center[k] - sphere_center[k]);
	normalize (axis2);
	dot12 = dot_product (axis1, axis2);
	if (dot12 < -1.0) dot12 = -1.0;
	else if (dot12 > 1.0) dot12 = 1.0;
	angle12 = acos (dot12);
	sprintf (message,"angle1 = %8.3f; angle2 = %8.3f; angle12 = %8.3f",
			angle1, angle2, angle12);
	informd(message);
	return (angle1 >= angle12 + angle2);
}

void add2solid (struct solid_angle *sa, struct face *fac)
{
	if (sa -> head_face == NULL)
		sa -> head_face = fac;
	else sa -> tail_face -> next = fac;
	sa -> tail_face = fac;
	sa -> n_face++;
}


