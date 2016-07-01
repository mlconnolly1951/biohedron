/*

	MSRoll

	Copyright 1986, 1989, 1996 by Michael L. Connolly
	All rights reserved

	Written by Michael L. Connolly.
	March 7, 2000

*/

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"


void sort_faces (struct surface *srf_ptr)
{
	long d;
	char message[MAXLINE];
	struct face *fac;
	struct cycle *cyc;
	struct edge *edg;
	
	d = srf_ptr -> n_face * 0.5 + 1;
	while (d > 0) {
		one_Shell (srf_ptr, d);
		if (error()) return;
		d = d * 0.5;
	}
	/* set up edge --> face pointers */
	for (fac = srf_ptr -> head_face; fac != NULL; fac = fac -> next)
		for (cyc = fac -> first_cycle; cyc != NULL; cyc = cyc -> next)
			for (edg = cyc -> first_edge; edg != NULL; edg = edg -> next)
				edg -> fac = fac;
	sprintf (message, "%8ld faces sorted", srf_ptr -> n_face);
	informd (message);
}

void one_Shell (struct surface *srf, long d)
{
	long i;
	int comparison;
	struct face *fac1;
	struct face *fac2;
	
	for (i = 0; i < srf -> n_face - d; i++) {
		fac1 = *(srf -> face_handles + i);
		fac2 = *(srf -> face_handles + i + d);
		comparison = compare_faces (fac1, fac2);
		if (comparison > 0)
			swap_faces (fac1, fac2);
		if (error()) return;
	}
	srf -> head_face = *(srf -> face_handles);
	srf -> tail_face = *(srf -> face_handles + srf -> n_face - 1);
}

int compare_faces (struct face *fac1, struct face *fac2)
{
	int comparison;
	struct variety *vty1, *vty2;
	
	vty1 = fac1 -> vty;
	vty2 = fac2 -> vty;
	if (vty1 -> center[2] < vty2 -> center[2])
		comparison = -1;
	else if  (vty1 -> center[2] > vty2 -> center[2])
		comparison = 1;
	else comparison = 0;
	return (comparison);
}

void swap_faces (struct face *fac1, struct face *fac2)
{
	unsigned long size;
	struct face temp_face;
	struct face *next1, *next2;
	
	next1 = fac1 -> next;
	next2 = fac2 -> next;
	size = sizeof (struct face);
	copy_bytes ((char *) fac1, (char *) &temp_face, size);
	copy_bytes ((char *) fac2, (char *) fac1, size);
	copy_bytes ((char *) &temp_face, (char *) fac2, size);
	fac1 -> next = next1;
	fac2 -> next = next2;
}


int setup_ffn (struct surface *this_srf)
{
	long j, ffn;
	char message[MAXLINE];
	struct cycle *cyc;
	struct face *fac;
	struct edge *edg;
	struct arc *arcptr;
	
	if (this_srf -> face_handles == NULL) {
		sprintf (message, "null face handles");
		set_error1(message);
		return (0);
	}
	if (this_srf -> arc_handles == NULL || this_srf -> n_arc_original <= 0) {
		return (1);
	}

	/* initialize first face number to n_face + 1 */
	for (j = 0; j < this_srf -> n_arc_original; j++) {
		arcptr = *(this_srf -> arc_handles+j);
		if (arcptr == NULL) {
			sprintf (message, "invalid arc handle (%8ld) for arc %6ld", arcptr, j + 1);
			set_error1(message);
			return (0);
		}
		arcptr -> ffn = this_srf -> n_face + 1;
	}
		
	/* arc first face number based on faces that border it (not kitty-korner) */
	for (ffn = 1; ffn <=  (this_srf -> n_face); ffn++) {
		fac = *(this_srf -> face_handles + ffn - 1);
		if (fac == NULL) {
			sprintf (message, "invalid face handle (%8ld) for face %6ld", fac, ffn);
			set_error1(message);
			return (0);
		}
		/* if (this_srf -> van_der_Waals && fac -> shape != CONVEX) continue; */
		for (cyc = fac -> first_cycle; cyc != (struct cycle *) NULL;
			cyc = cyc -> next) {
			for (edg = cyc -> first_edge; edg != (struct edge *) NULL;
				edg = edg -> next) {
				arcptr = edg -> arcptr;
				/* if (this_srf -> van_der_Waals && arcptr -> shape != CONVEX) continue; */
				if (arcptr -> ffn > ffn) arcptr -> ffn = ffn;
			}
		}
	}
	
	/* check that all have valid values */
	for (j = 0; j < this_srf -> n_arc_original; j++) {
		arcptr = *(this_srf -> arc_handles+j);
		/* if (this_srf -> van_der_Waals && arcptr -> shape != CONVEX) continue; */
		if (arcptr -> ffn >  this_srf -> n_face || arcptr -> ffn < 1) {
			sprintf (message, "invalid first face number (%6ld) for arc %6ld", arcptr -> ffn, j + 1);
			set_error1(message);
			return(0);
		}
	}
	return(1);
}

int original_lfn (struct surface *this_srf)
{
	char message[MAXLINE];
	struct vertex *vtx;
	struct arc *arcptr;
	
	initialize_lfn (this_srf);
	if (error()) return(0);
	iterate_lfn (this_srf);
	if (error()) return(0);
	
	/* check */
	for (arcptr = this_srf -> head_arc; arcptr != (struct arc *) NULL; arcptr = arcptr -> next) {
		/* if (this_srf -> van_der_Waals && arcptr -> shape != CONVEX) continue; */
		if (arcptr -> lfn <= 0) {
			sprintf (message,
				"(original_lfn): check arc lfn = %ld", arcptr -> lfn);
			set_error1(message);
			return(0);
		}
	}
	for (vtx = this_srf -> head_vertex; vtx != (struct vertex *) NULL; vtx = vtx -> next) {
		if (vtx -> lfn <= 0) {
			sprintf (message,
				"(original_lfn): check vertex (%5ld) lfn = %ld",
					vtx -> number, vtx -> lfn);
			set_error1(message);
			return(0);
		}
	}
	return(1);
}

void initialize_lfn (struct surface *this_srf)
{
	long lfn;
	struct face *fac;
	struct variety *vty;
	
	
	/* face last face number (initial value) */
	for (lfn = 1; lfn <= this_srf ->  n_face; lfn++) {
		fac = *(this_srf -> face_handles + lfn - 1);
		fac -> lfn = lfn;
		fac -> vty -> lfn = lfn;
	}

	/* variety last face number - varieties not purged during triangulation */
	for (lfn = 1; lfn <= this_srf ->  n_face; lfn++) {
		fac = *(this_srf -> face_handles + lfn - 1);
		/* if (this_srf ->  van_der_Waals && fac -> shape != CONVEX) continue; */
		vty = fac -> vty;
		if (vty -> lfn < fac -> lfn) vty -> lfn = fac -> lfn;
	}
}

void iterate_lfn (struct surface *this_srf)
{
	int j, scone;
	long lfn;
	char message[MAXLINE];
	struct cycle *cyc;
	struct face *fac, *fac2;
	struct edge *edg;
	struct arc *arcptr;
	struct circle *cir;
	struct vertex *vtx;
	
		
	/* vertex last face number based on face last face number */
	for (lfn = 1; lfn <= this_srf ->  n_face; lfn++) {
		fac = *(this_srf -> face_handles + lfn - 1);
		/* if (this_srf -> van_der_Waals && fac -> shape != CONVEX) continue; */
		for (cyc = fac -> first_cycle; cyc != (struct cycle *) NULL;
			cyc = cyc -> next) {
			for (edg = cyc -> first_edge; edg != (struct edge *) NULL;
				edg = edg -> next) {
				arcptr = edg -> arcptr;
				for (j = 0; j < 2; j++) {
					vtx = arcptr -> vtx[j];
					if (vtx == (struct vertex *) NULL) continue;
					if (vtx -> lfn < lfn) vtx -> lfn = fac -> lfn;
				}
			}
		}
		scone = is_saddle_cone (fac);
		if (scone) {
			vtx = find_cone_vertex (fac);
			if (vtx == (struct vertex *) NULL) continue;
			if (vtx -> lfn < lfn) vtx -> lfn = fac -> lfn;
		}
	}
	
	/* arc last face number based on vertex last face number and face last face number */
	for (lfn = 1; lfn <= this_srf ->  n_face; lfn++) {
		fac = *(this_srf -> face_handles + lfn - 1);
		/* if (this_srf ->  van_der_Waals && fac -> shape != CONVEX) continue; */
		for (cyc = fac -> first_cycle; cyc != (struct cycle *) NULL;
			cyc = cyc -> next) {
			for (edg = cyc -> first_edge; edg != (struct edge *) NULL;
				edg = edg -> next) {
				arcptr = edg -> arcptr;
				for (j = 0; j < 2; j++) {
					vtx = arcptr -> vtx[j];
					if (vtx == (struct vertex *) NULL) continue;
					if (vtx -> lfn > arcptr -> lfn) arcptr -> lfn = vtx -> lfn;
				}
				for (j = 0; j < 2; j++) {
					fac2 = arcptr -> edg[j] -> fac;
					if (fac2 -> lfn > arcptr -> lfn) arcptr -> lfn = fac2 -> lfn;
				}
				if (arcptr -> lfn <= 0) {
					sprintf (message, "(original_lfn): arc lfn = %ld",
						arcptr -> lfn);
					set_error1(message);
					return;
				}
			}
		}
	}
	
	/* edge and circle last face number based on arc last face number */
	for (lfn = 1; lfn <= this_srf ->  n_face; lfn++) {
		fac = *(this_srf -> face_handles + lfn - 1);
		/* if (this_srf ->  van_der_Waals && fac -> shape != CONVEX) continue; */
		for (cyc = fac -> first_cycle; cyc != (struct cycle *) NULL;
			cyc = cyc -> next) {
			for (edg = cyc -> first_edge; edg != (struct edge *) NULL;
				edg = edg -> next) {
				arcptr = edg -> arcptr;
				cir = arcptr -> cir;
				arcptr -> edg[0] -> lfn = arcptr -> lfn;
				arcptr -> edg[1] -> lfn = arcptr -> lfn;
				if (arcptr -> lfn > cir -> lfn) cir -> lfn = arcptr -> lfn;
			}
		}
	}
	
	/* face and cycle last face number based on arc, circle and vertex */
	for (lfn = 1; lfn <= this_srf ->  n_face; lfn++) {
		fac = *(this_srf -> face_handles + lfn - 1);
		/* if (this_srf ->  van_der_Waals && fac -> shape != CONVEX) continue; */
		for (cyc = fac -> first_cycle; cyc != (struct cycle *) NULL;
			cyc = cyc -> next) {
			for (edg = cyc -> first_edge; edg != (struct edge *) NULL;
				edg = edg -> next) {
				arcptr = edg -> arcptr;
				cir = arcptr -> cir;
				if (arcptr -> lfn > fac -> lfn) fac -> lfn = arcptr -> lfn;
				/*
				if (cir -> lfn > fac -> lfn) fac -> lfn = cir -> lfn;
				*/
				for (j = 0; j < 2; j++) {
					vtx = arcptr -> vtx[j];
					if (vtx == (struct vertex *) NULL) continue;
					if (vtx -> lfn > fac -> lfn) fac -> lfn = vtx -> lfn;
				}
			}
		}
	}
}

int original_face_number (struct surface *this_srf)
{
	int j, scone;
	long ofn;
	char message[MAXLINE];
	struct cycle *cyc;
	struct face *fac;
	struct arc *arcptr;
	struct edge *edg;
	struct vertex *vtx;
	
	for (ofn = 1; ofn <= this_srf ->  n_face; ofn++) {
		fac = *(this_srf -> face_handles + ofn - 1);
		fac -> ofn = ofn;
	}
	/* vertex original face number */
	for (ofn = 1; ofn <= this_srf ->  n_face; ofn++) {
		fac = *(this_srf -> face_handles + ofn - 1);
		/* if (this_srf ->  van_der_Waals && fac -> shape != CONVEX) continue; */
		for (cyc = fac -> first_cycle; cyc != (struct cycle *) NULL;
			cyc = cyc -> next) {
			for (edg = cyc -> first_edge; edg != (struct edge *) NULL;
				edg = edg -> next) {
				arcptr = edg -> arcptr;
				for (j = 0; j < 2; j++) {
					vtx = arcptr -> vtx[j];
					if (vtx == (struct vertex *) NULL) continue;
					/* use the first one (most negatively curved) */
					if (vtx -> ofn != 0) continue;
					vtx -> ofn = ofn;
				}
			}
		}
		scone = is_saddle_cone (fac);
		if (scone) {
			vtx = find_cone_vertex (fac);
			if (vtx == (struct vertex *) NULL) continue;
			if (vtx -> ofn < ofn) vtx -> ofn = ofn;
		}
	}
	/* arc original face number */
	for (ofn = 1; ofn <= this_srf ->  n_face; ofn++) {
		fac = *(this_srf -> face_handles + ofn - 1);
		/* if (this_srf -> van_der_Waals && fac -> shape != CONVEX) continue; */
		for (cyc = fac -> first_cycle; cyc != (struct cycle *) NULL;
			cyc = cyc -> next) {
			for (edg = cyc -> first_edge; edg != (struct edge *) NULL;
				edg = edg -> next) {
				arcptr = edg -> arcptr;
				if (arcptr -> ofn != 0) continue;
				arcptr -> ofn = ofn;
			}
		}
	}
	/* check */
	for (arcptr = this_srf ->  head_arc; arcptr != (struct arc *) NULL; arcptr = arcptr -> next) {
		/* if (this_srf -> van_der_Waals && arcptr -> shape != CONVEX) continue; */
		if (arcptr -> ofn <= 0) {
			sprintf (message,
				"(original_face_number): check arc ofn = %ld",
				arcptr -> ofn);
			set_error1(message);
			return(0);
		}
	}
	for (vtx = this_srf ->  head_vertex; vtx != (struct vertex *) NULL; vtx = vtx -> next) {
		if (vtx -> ofn <= 0) {
			sprintf (message,
				"(original_face_number): check vertex (%5ld) ofn = %ld",
					vtx -> number, vtx -> ofn);
			set_error1(message);
			return(0);
		}
	}
	return(1);
}

/* purge every object with a last face number less than the given one */
/* except some original objects */

int purge_surface (struct surface *this_srf, long lfn)
{
	int purge;
	long nv, na, nf, nc, ne;
	char message[MAXLINE];
	struct vertex *prev_vtx, *vtx, *next_vtx;
	struct arc *prev_arc, *a, *next_arc;
	struct edge *edg0, *edg1;
	struct face *fac0, *fac1;
	struct circle *prev_circle, *cir, *next_circle;
	struct face *prev_face, *f, *next_face;
	struct cycle *cyc, *next_cycle;
	

	/* purge vertices */
	nv = 0;
	prev_vtx = (struct vertex *) NULL;
	next_vtx = NULL;
	for (vtx = this_srf -> head_vertex; vtx != NULL; vtx = next_vtx) {
		next_vtx = vtx -> next;
		if (vtx -> lfn <= 0) {
			prev_vtx = vtx;
			continue;
		}
		if (vtx -> converted && vtx -> lfn <= lfn) {
			free_object (VERTEX, (short *) vtx);
			if (prev_vtx == (struct vertex *) NULL)
				this_srf ->  head_vertex = next_vtx;
			else prev_vtx -> next = next_vtx;
			if (vtx == this_srf ->  tail_vertex) this_srf ->  tail_vertex = prev_vtx;
			nv++;
			continue;
		}
		prev_vtx = vtx;
	}

	/* purge arcs and edges */
	na = 0; ne = 0;
	prev_arc = (struct arc *) NULL;
	next_arc = NULL;
	for (a = this_srf -> head_arc; a != NULL; a = next_arc) {
		next_arc = a -> next;
		if (a -> lfn <= 0) {
			prev_arc = a;
			continue;
		}
		if (a -> converted && a -> lfn <= lfn) {
			if (!a -> converted) {
				edg0 = a -> edg[0];
				edg1 = a -> edg[1];
				if (edg0 != NULL) fac0 = edg0 -> fac;
				else fac0 = NULL;
				if (edg1 != NULL) fac1 = edg1 -> fac;
				else fac1 = NULL;
				sprintf (message, "purging before converted, arc number = %8ld; global lfn = %8ld",
					a -> number, lfn);
				set_error1(message);
				if (fac0 != NULL && fac1 != NULL) {
					sprintf (message, "arc ffn = %8ld, ofn = %8ld, lfn = %8ld; ofns: %8ld %8ld",
						a -> ffn, a -> ofn, a -> lfn, fac0 -> ofn, fac1 -> ofn);
					set_error2(message);
				}
				return (0);
			}
			if (a -> original)
				*(this_srf -> arc_handles + a -> number - 1) = (struct arc *) NULL;
			free_object (EDGE, (short *) (a -> edg[0])); ne++;
			free_object (EDGE, (short *) (a -> edg[1])); ne++;
			free_object (ARC, (short *) a); na++;
			if (prev_arc == (struct arc *) NULL)
				this_srf ->  head_arc = next_arc;
			else prev_arc -> next = next_arc;
			if (a == this_srf ->  tail_arc) this_srf ->  tail_arc = prev_arc;
		}
		else prev_arc = a;
	}

	/* purge faces and cycles */
	nf = 0; nc = 0;
	prev_face = (struct face *) NULL;
	next_face = NULL; /* keep lint happy */
	for (f = this_srf -> head_face; f != NULL; f = next_face) {
		next_face = f -> next;
		purge = (! f -> original && f -> lfn <= lfn && f -> converted);
		if (purge) {
			next_cycle = NULL;
			for (cyc = f -> first_cycle; cyc != (struct cycle *) NULL;
				cyc = next_cycle) {
				next_cycle = cyc -> next;
				free_object (CYCLE, (short *) cyc);
				nc++;
			}
			free_object (FACE, (short *) f);
			if (prev_face == (struct face *) NULL)
				this_srf ->  head_face = next_face;
			else prev_face -> next = next_face;
			if (f == this_srf ->  tail_face)
				this_srf ->  tail_face = prev_face;
			nf++;
		}
		else prev_face = f;
	}
	
	/* purge circles */
	nc = 0;
	prev_circle = (struct circle *) NULL;
	next_circle = NULL;
	for (cir = this_srf -> head_circle; cir != NULL; cir = next_circle) {
		next_circle = cir -> next;
		if (cir -> lfn <= 0) {
			prev_circle = cir;
			continue;
		}
		if (cir -> lfn <= lfn) {
			free_object (CIRCLE, (short *) cir);
			if (prev_circle == (struct circle *) NULL)
				this_srf -> head_circle = next_circle;
			else prev_circle -> next = next_circle;
			if (cir == this_srf -> tail_circle)
				this_srf -> tail_circle = prev_circle;
			nc++;
		}
		else prev_circle = cir;
	}
	return(1);
}


/*
   MSRoll
   Copyright 1996 by Michael L. Connolly
   All Rights Reserved

*/
