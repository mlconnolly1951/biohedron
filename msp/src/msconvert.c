/*
	MSRoll
	Copyright 1986, 1989 by Michael L. Connolly
	All rights reserved
	January 8, 2002
*/

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"


/* prepare piecewise polynomial molecular surface for writing */
int cvt_surface (struct surface *srf, double surface_center[3])
{
	int result;
	short *n_edge_per_cycle;
	long i, j, y, e, cycned, fcynum, fednum;
	long *first_cycle_number;
	long *first_edge_number;
	char message[MAXLINE];
	struct edge *edg;
	struct face *fac;
	struct cycle *cyc;

	
	srf -> n_variety = srf -> n_atom + srf -> n_tori + srf -> n_probe;

	if (srf -> n_atom >= MAX_ATOM) {
		sprintf (message,"maximum number of atoms (%d) exceeded: %ld",
			MAX_ATOM, srf -> n_atom);
		set_error1 (message);
		return (0);
	}

	/* converting */
	cvt_header (srf);
	if (error()) return(0);
	srf -> variety_handles = cvt_varieties (srf, surface_center);
	if (error()) return(0);
	srf -> vertex_handles = cvt_vertices (srf);
	if (error()) return(0);
	srf -> circle_handles = cvt_circles (srf);
	if (error()) return(0);
	srf -> arc_handles = cvt_arcs (srf);
	if (error()) return(0);
	if (srf -> n_face > 0) {
		first_cycle_number = allocate_longs (srf -> n_face, 0, FIRST_CYCLE_NUMBER);
		if (first_cycle_number == NULL) {
			print_counts();
			set_error1 ("read_surface: memory allocation failure");
			return (0);
		}
	}
	else {
		first_cycle_number = NULL;
	}
	srf -> face_handles = cvt_faces (srf, first_cycle_number);
	if (error()) return(0);
	if (srf -> n_cycle > 0) {
		n_edge_per_cycle = allocate_shorts (srf -> n_cycle);
		if (n_edge_per_cycle == NULL) {
			print_counts();
			set_error1 ("read_surface: memory allocation failure");
			return (0);
		}
		first_edge_number = allocate_longs (srf -> n_cycle, 0, FIRST_EDGE_NUMBER);
		if (first_edge_number == NULL) {
			print_counts();
			set_error1 ("read_surface: memory allocation failure");
			return (0);
		}
	}
	else {
		n_edge_per_cycle = NULL;
		first_edge_number = NULL;
	}
	srf -> cycle_handles = cvt_cycles (srf, n_edge_per_cycle, first_edge_number);
	if (error()) return(0);
	for (i = 0; i < srf -> n_face; i++) {
		fac = *(srf -> face_handles + i);
		fcynum = *(first_cycle_number + i);
		fac -> first_cycle = (fcynum == 0) ? NULL : *(srf -> cycle_handles + fcynum - 1);
	}
	srf -> edge_handles = cvt_edges (srf, n_edge_per_cycle);
	if (error()) return(0);
	e = 0;
	for (y = 0; y < srf -> n_cycle; y++) {
		cyc = *(srf -> cycle_handles + y);
		fednum = *(first_edge_number + y);
		cyc -> first_edge = *(srf -> edge_handles + fednum - 1);
		cycned = *(n_edge_per_cycle + y);
		for (j = 0; j < cycned; j++, e++) {
			edg = *(srf -> edge_handles + e);
			edg -> next = ((j < cycned - 1) ? *(srf -> edge_handles + e + 1) : NULL);
		}
	}
	result = cvt_components (srf);
	if (result == 0 || error()) return(0);

	/* free temporary memory */
	free_shorts (n_edge_per_cycle);
	free_longs (first_edge_number, 0, FIRST_EDGE_NUMBER);
	free_longs (first_cycle_number, 0, FIRST_CYCLE_NUMBER);
	return (1);
}


void cvt_header (struct surface *srf)
{

	/* initialize surface structure */
	srf -> weight = DEFAULT_WEIGHT;
	srf -> large_omega = DEFAULT_LARGE;
	srf -> surface_type = PQMS_SURFACE;
	srf -> van_der_Waals = (srf -> probe_radius == 0.0);
}

struct variety **cvt_varieties (struct surface *srf, double surface_center[3])
{
	int k;
	long v;
	char message[MAXLINE];
	struct variety **variety_handles;
	struct variety *vty;
	struct sphere *atm;
	struct torus *tor;
	struct probe *prb;

	/* allocate memory for handles */
	variety_handles = (struct variety **)
		allocate_pointers (VARIETY, srf -> n_variety);
	if (variety_handles == (struct variety **) NULL) {
		print_counts();
		set_error1 ("cvt_varieties: memory allocation failure");
		sprintf (message, "%8ld variety handles", srf -> n_variety);
		set_error2 (message);
		return (NULL);
	}
	for (v = 0; v < srf -> n_variety; v++) {
		vty = allocate_variety ();
		if (vty == (struct variety *) NULL) {
			set_error1 ("cvt_varieties: variety allocation failure");
			sprintf (message, "variety %8ld", v+1);
			set_error2 (message);
			return (NULL);
		}
		*(variety_handles + v) = vty;
	}
	/* center of molecule */
	for (k = 0; k < 3; k++)
		surface_center[k] = 0.0;
	v = 0;
	for (atm = srf -> head_atom; atm != NULL; atm = atm -> next, v++) {
		vty = *(variety_handles + v);
		vty -> type = SPHERE;
		for (k = 0; k < 3; k++) {
			vty -> center[k] = atm -> center[k];
			vty -> axis[k] = 0.0;
			surface_center[k] += vty -> center[k];
		}
		for (k = 0; k < MAXPA; k++)
			vty -> atmnum[k] = 0;
		vty -> atmnum[0] = atm -> number;
		vty -> radii[0] = atm -> radius;
		vty -> number = v + 1;
	}
	for (tor = srf -> head_torus; tor != NULL; tor = tor -> next, v++) {
		vty = *(variety_handles + v);
		vty -> type = TORUS;
		for (k = 0; k < 3; k++) {
			vty -> center[k] = tor -> center[k];
			vty -> axis[k] = tor -> axis[k];
		}
		for (k = 0; k < 2; k++)
			vty -> atmnum[k] = tor -> atm[k] -> number;
		vty -> atmnum[2] = 0;
		vty -> radii[0] = tor -> radius;
		vty -> number = v + 1;
	}
	for (prb = srf -> head_probe; prb != NULL; prb = prb -> next, v++) {
		vty = *(variety_handles + v);
		vty -> type = SPHERE;
		vty -> radii[0] = srf -> probe_radius;
		for (k = 0; k < 3; k++) {
			vty -> center[k] = prb -> center[k];
			vty -> axis[k] = 0.0;
		}
		for (k = 0; k < MAXPA; k++) {
			atm = prb -> atm[k];
			if (atm != NULL)
				vty -> atmnum[k] = atm -> number;
		}
		vty -> number = v + 1;
	}
	if (srf -> n_atom <= 0) {
		sprintf (message,
			"cvt_varieties: invalid number of atoms: %ld", srf -> n_atom);
		set_error1(message);
		return (0);
	}
	for (k = 0; k < 3; k++)
		surface_center[k] /= srf -> n_atom;
	sprintf (message,"%8.3f %8.3f %8.3f surface center",
		surface_center[0], surface_center[1], surface_center[2]);
	informd(message);
	sprintf (message,"%8ld varieties converted", srf -> n_variety);
	informd(message);
	return (variety_handles);
}


struct vertex **cvt_vertices (struct surface *srf)
{
	long i;
	char message[MAXLINE];
	struct vertex **vertex_handles;
	struct vertex *vtx;

	if (srf -> n_vertex > 0) {
		vertex_handles = (struct vertex **)
			allocate_pointers (VERTEX, srf -> n_vertex);
		if (vertex_handles == (struct vertex **) NULL) {
			print_counts();
			set_error1 ("read_vertices: memory allocation failure");
			set_error2 ("for vertex handles");
			return (NULL);
		}
		for (i = 0, vtx = srf -> head_vertex; vtx != NULL; vtx = vtx -> next, i++) {
			*(vertex_handles + i) = vtx;
			vtx -> number = i + 1;
			vtx -> cusp = 1;	/* initialize all to cusp */
		}
	}
	else vertex_handles = (struct vertex **) NULL;
	sprintf (message,"%8ld vertices converted", srf -> n_vertex);
	informd(message);
	return (vertex_handles);
}

struct circle **cvt_circles (struct surface *srf)
{
	long i;
	char message[MAXLINE];
	struct circle **circle_handles;
	struct circle *cir;

	if (srf -> n_circle > 0) {
		circle_handles = (struct circle **)
			allocate_pointers (CIRCLE, srf -> n_circle);
		if (circle_handles == (struct circle **) NULL) {
		print_counts();
			set_error1 ("cvt_surface: memory allocation failure");
			set_error2 ("for circle handles");
			return (NULL);
		}
		for (cir = srf -> head_circle, i = 0; cir != NULL; cir = cir -> next, i++) {
			*(circle_handles + i) = cir;
		}
	}
	else circle_handles = (struct circle **) NULL;
	sprintf (message,"%8ld circles converted", srf -> n_circle);
	informd(message);
	return (circle_handles);
}

struct arc **cvt_arcs (struct surface *srf)
{
	int m, j, k;
	long a;
	long cirnum, vtxnum;
	double vector1[3], vector2[3];
	char message[MAXLINE];
	struct arc **arc_handles;
	struct vertex *arc_vtx;
	struct vertex *vtx1, *vtx2;
	struct circle *cir;
	struct arc *ark;
	struct face *fac;
	struct cycle *cyc;
	struct edge *edg;
	struct edge *e;

	if (srf -> n_arc == 0) {
		arc_handles = (struct arc **) NULL;
		return (arc_handles);
	}
	arc_handles = (struct arc **)
		allocate_pointers (ARC, srf -> n_arc);
	if (arc_handles == (struct arc **) NULL) {
		print_counts();
		set_error1 ("cvt_arcs: memory allocation failure");
		set_error2 ("for arc handles");
		return (NULL);
	}
	a = 0;
	for (fac =  srf -> head_face; fac != NULL; fac = fac -> next) {
		switch ((int) fac -> shape) {
		case CONVEX:
			if (fac -> n_cycle <= 0) break;
			for (cyc = fac -> first_cycle; cyc != NULL; cyc = cyc -> next)
				for (e = cyc -> first_edge; e != NULL; e = e -> next) {
					ark = e -> arcptr;
					*(arc_handles + a++) = ark;
					ark -> shape = CONVEX;
					/* index to pointer */
					cirnum = ark -> cir -> number;
					if (cirnum < 1 || cirnum > srf -> n_circle) {
						sprintf (message, "cvt_arcs: invalid arc-circle number %8ld", cirnum);
						set_error1 (message);
						return (NULL);
					}
					for (j = 0; j < 2; j++) {
						vtxnum = ((ark -> vtx[j] == NULL) ? 0 : ark -> vtx[j] -> number);
						if (vtxnum < 0 || vtxnum > srf -> n_vertex) {
							sprintf (message, "cvt_arcs: invalid arc-vertex number %8ld", vtxnum);
							set_error1 (message);
							return (NULL);
						}
					}
				}		/* end of arc loop */
			break;
		case CONCAVE:
			if (fac -> ptr.prb -> low) {
				for (cyc = fac -> first_cycle; cyc != NULL;
					cyc = cyc -> next)
					for (edg = cyc -> first_edge; edg != NULL;
						edg = edg -> next) {
						if (edg -> orn) continue;
						ark = edg -> arcptr;
						*(arc_handles + a++) = ark;
						ark -> shape = CONCAVE;
						/* index to pointer */
						cirnum = ark -> cir -> number;
						if (cirnum < 1 || cirnum > srf -> n_circle) {
							sprintf (message, "cvt_arcs: invalid arc-circle number %8ld", cirnum);
							set_error1 (message);
							return (NULL);
						}
						for (j = 0; j < 2; j++) {
							vtxnum = ((ark -> vtx[j] == NULL) ? 0 : ark -> vtx[j] -> number);
							if (vtxnum < 0 || vtxnum > srf -> n_vertex) {
								sprintf (message, "cvt_arcs: invalid arc-vertex number %8ld", vtxnum);
								set_error1 (message);
								return (NULL);
							}
						}
					}	/* end of arc loop */
			}
			else {
				for (m = 0; m < fac -> n_arc; m++) {
					ark = fac -> arcsp[m];
                    if (ark == NULL) continue;
					*(arc_handles + a++) = ark;
					ark -> shape = CONCAVE;
					/* index to pointer */
					cirnum = ark -> cir -> number;
					if (cirnum < 1 || cirnum > srf -> n_circle) {
						sprintf (message, "cvt_arcs: invalid arc-circle number %8ld", cirnum);
						set_error1 (message);
						return (NULL);
					}
					for (j = 0; j < 2; j++) {
						vtxnum = ((ark -> vtx[m] == NULL) ? 0 : ark -> vtx[j] -> number);
						if (vtxnum < 0 || vtxnum > srf -> n_vertex) {
							sprintf (message, "cvt_arcs: invalid arc-vertex number %8ld", vtxnum);
							set_error1 (message);
							return (NULL);
						}
					}
				}		/* end of arc loop */
			}
			break;
		case SADDLE:
			break;
		default:
			set_error1 ("cvt_arcs: invalid shape");
			return (NULL);
		}
	}
	if (a != srf -> n_arc) {
		sprintf (message, "(cvt_arcs): arc count inconsistency: %8ld %8ld", a, srf -> n_arc);
		set_error1(message);
		return (NULL);
	}
	/* set up head and tail pointers */
	srf -> head_arc = NULL;
	srf -> tail_arc = NULL;
	for (a = 0; a < srf -> n_arc; a++) {
		ark = *(arc_handles + a);
		/* set up head and tail pointers for linked list */
		if (a == 0) srf -> head_arc = ark;
		srf -> tail_arc = ark;
		/* compute phi for omega */
		vtx1 = ark -> vtx[0];
		vtx2 = ark -> vtx[1];
		cir = ark -> cir;
		if (cir -> radius <= 0.0) ark -> phi = 0.0;
		else if (vtx1 != NULL && vtx2 != NULL) {
			for (k = 0; k < 3; k++) {
				vector1[k] = (vtx1 -> center[k] - cir -> center[k]) / cir -> radius;
				vector2[k] = (vtx2 -> center[k] - cir -> center[k]) / cir -> radius;
			}
			ark -> phi = positive_angle (vector1, vector2, cir -> axis);
		}
		else ark -> phi = 2 * PI;
		ark -> next = ((a < srf -> n_arc - 1) ? *(arc_handles + a + 1) : NULL);
		ark -> number = a+1;
		ark -> original = 1;
	}
	/* mark those cusp vertices that show up in no arc */
	/* these vertices are in saddle cones */
	for (ark = srf -> head_arc; ark != NULL; ark = ark -> next) {
		for (k = 0; k < 2; k++) {
			arc_vtx = ark -> vtx[k];
			if (arc_vtx != NULL)
					arc_vtx -> cusp = 0;
		}
	}
	sprintf (message,"%8ld arcs converted", srf -> n_arc);
	informd(message);
	return (arc_handles);
}


struct face **cvt_faces (struct surface *srf, long *first_cycle_number)
{
	long icyc;
	long i, vtynum, fcynum;
	char message[MAXLINE];
	struct face *fac;
	struct face **face_handles;

	face_handles = (struct face **)
		allocate_pointers (FACE, srf -> n_face);
	if (face_handles == (struct face **) NULL) {
		print_counts();
		set_error1 ("cvt_faces: memory allocation failure");
		set_error2 ("for face handles");
		return (NULL);
	}
	icyc = 1;
	for (fac =  srf -> head_face, i = 0; fac != NULL; fac = fac -> next, i++) {
		*(face_handles + i) = fac;
		fac -> srf = srf;
		fac -> original = 1;
		fac -> largeok = 1;			/* ok to use add-circle routine for large solid angle */
		switch ((int) fac -> shape) {
		case CONVEX:
			vtynum = fac -> ptr.atm -> number;
			fcynum = ((fac -> n_cycle > 0) ? icyc : 0);
			icyc += fac -> n_cycle;
			break;
		case SADDLE:
			vtynum = srf -> n_atom + fac -> ptr.tor -> number;
			fcynum = icyc;
			if (fac -> n_arc == 2) icyc += 2;
			else icyc++;
			break;
		case CONCAVE:
			vtynum = srf -> n_atom + srf -> n_tori + fac -> ptr.prb -> number;
			fcynum = icyc;
			if (fac -> ptr.prb -> low)
				icyc += fac -> n_cycle;
			else icyc++;
			break;
		default:
			set_error1("cvt_faces: invalid face shape");
			return (NULL);
		}
		fac -> vty = *(srf -> variety_handles + vtynum - 1);
		*(first_cycle_number + i) = fcynum;
	}
	sprintf (message,"%8ld faces converted", srf -> n_face);
	informd(message);
	return (face_handles);
} 

struct cycle **cvt_cycles (struct surface *srf, short *n_edge_per_cycle, long *first_edge_number)
{
	long  icyc, iedg;
	char message[MAXLINE];
	struct cycle **cycle_handles;
	struct cycle *cyc, *cyc1, *cyc2;
	struct face *fac;
	
	if (srf -> n_cycle == 0) {
		return (NULL);
	}
	cycle_handles = (struct cycle **)
		allocate_pointers (CYCLE, srf -> n_cycle);
	if (cycle_handles == (struct cycle **) NULL) {
		print_counts();
		set_error1 ("cvt_cycles: memory allocation failure");
		set_error2 ("for cycle handles");
		return (NULL);
	}
	icyc = 1;
	iedg = 1;
	for (fac =  srf -> head_face; fac != NULL; fac = fac -> next) {
		if (fac -> n_cycle <= 0) continue;
		for (cyc = fac -> first_cycle; cyc != NULL; cyc = cyc -> next) {
			*(cycle_handles + icyc-1) = cyc;
			*(n_edge_per_cycle + icyc-1) = edges_in_cycle (cyc);
			*(first_edge_number + icyc-1) = iedg;
			iedg += edges_in_cycle (cyc);
			icyc++;
		}
	}
	if (iedg != srf -> n_edge + 1) {
		set_error1 ("(cvt_cycles): edge count inconsistency");
		return (NULL);
	}
	sprintf (message,"%8ld cycles converted", srf -> n_cycle);
	informd(message);
	return (cycle_handles);
}


struct edge **cvt_edges (struct surface *srf, short *n_edge_per_cycle)
{
	int j, orn;
	long e, ne, y;
	char message[MAXLINE];
	struct edge *edg;
	struct edge **edge_handles;
	struct face *fac;
	struct cycle *cyc;
	struct arc *a;

	if (srf -> n_edge == 0) {
		return (NULL);
	}
	edge_handles = (struct edge **)
		allocate_pointers (EDGE, srf -> n_edge);
	if (edge_handles == (struct edge **) NULL) {
		print_counts();
		set_error1 ("cvt_edges: memory allocation failure");
		set_error2 ("for edge handles");
		return (NULL);
	}
	e = 0;
	y = 0;
	for (fac =  srf -> head_face; fac != NULL; fac = fac -> next) {
		switch ((int) fac -> shape) {
		case CONVEX:
			if (fac -> n_cycle <= 0) break;
			for (cyc = fac -> first_cycle; cyc != NULL; cyc = cyc -> next, y++) {
				ne = (*(n_edge_per_cycle+y));
				if (ne == 0) continue;
				for (edg = cyc -> first_edge; edg != NULL; edg = edg -> next, e++) {
					*(edge_handles + e) = edg;
					orn = 0;
					edg -> orn = (short) orn;
					a = edg -> arcptr;
					if (a -> number < 1 || a -> number > srf -> n_arc) {
						set_error1 ("cvt_edges: invalid arc number");
						return (NULL);
					}
					a -> edg[orn] = edg;
				}
			}
			break;
		case SADDLE:
			for (j = fac -> n_arc - 1; j >= 0; j--,  e++) {
				a = fac -> arcsp[j];
                if (a == NULL) continue;
				if (a -> number < 1 || a -> number > srf -> n_arc) {
					set_error1 ("cvt_edges: invalid arc number");
					return (NULL);
				}
				edg = a -> edg[1];
				if (edg == (struct edge *) NULL) {
					set_error1 ("cvt_edges: edge failure");
					return (NULL);
				}
				*(edge_handles + e) = edg;
			}
			if (fac -> n_arc == 2) {
				y += 2;
			}
			else {
				y++;
			}
			break;
		case CONCAVE:
			if (fac -> ptr.prb -> low) {
				for (cyc = fac -> first_cycle; cyc != NULL; cyc = cyc -> next, y++) {
					for (edg = cyc -> first_edge; edg != NULL; edg = edg -> next, e++) {
						*(edge_handles + e) = edg;
						a = edg -> arcptr;
						if (a -> number < 1 || a -> number > srf -> n_arc) {
							set_error1 ("cvt_edges: invalid arc number");
							return (NULL);
						}
						edg -> fac = fac;
					}
				}
			}
			else {
				for (j = 0; j < fac -> n_arc; j++, e++) {
					a = fac -> arcsp[j];
                    if (a == NULL) continue;
					if (a -> number < 1 || a -> number > srf -> n_arc) {
						set_error1 ("cvt_edges: invalid arc number");
						return (NULL);
					}
					edg = a -> edg[0];
					if (edg == (struct edge *) NULL) {
						set_error1 ("cvt_edges: edge failure");
						return (NULL);
					}
					*(edge_handles + e) = edg;
					edg -> fac = fac;
				}
				y++;
			}
			break;
		default:
			break;
		}
	}
	sprintf (message,"%8ld edges converted", srf -> n_edge);
	informd(message);
	if (e != srf -> n_edge) {
		set_error1 ("(cvt_edges): edge count inconsistency");
		return (NULL);
	}
	return (edge_handles);
}


int cvt_components (struct surface *srf)
{
	long comp;
	char message[MAXLINE];
	struct component *cmp_ptr;

	for (comp = 0; comp < srf -> n_component; comp++) {
		cmp_ptr = get_component_ptr (srf, comp+1);
		cmp_ptr -> type = (cmp_ptr -> volume > 0.0) ? OUTER_SUBTYPE : INNER_SUBTYPE;
	}
	sprintf (message,"%8ld components converted", srf -> n_component);
	informd(message);
	return (1);
}

/* POLYHEDRA */

int convert_vet (struct surface *phn)
{
	int j, k, result;
	double n;
	char message[MAXLINE];
	struct vertex *vtx;
	struct vertex *vtx0, *vtx1;
	struct arc *a;
	struct face *f, *fac;
	struct phnvtx *pv;
	struct probe *prb;

	/* check vertices */
	for (vtx = phn -> head_vertex; vtx != NULL; vtx = vtx -> next) {
		if (!vtx -> converted) {
			sprintf (message, "warning: pqms vertex (%8ld) not converted to polyhedron vertex",
				vtx -> number);
			inform (message);
		}
	}
	/* check arcs */
	for (a = phn -> head_arc; a != NULL; a = a -> next) {
		if (a -> eaten > 0) continue;
		if (a -> eater > 0) continue;
		if (!a -> converted) {
			vtx0 = a -> vtx[0]; vtx1 = a -> vtx[1];
			if (vtx0 == NULL && vtx1 == NULL) sprintf (message, "both vertices null");
			else if (vtx0 == NULL) sprintf (message, "null vtx0");
			else if (vtx1 == NULL) sprintf (message, "null vtx1");
			else sprintf (message, "vertices %8ld %8ld", vtx0 -> number, vtx1 -> number);
			inform(message);
		}
	}
	/* check faces */
	for (f = phn -> head_face; f != NULL; f = f -> next) {
		/* skip problem faces */
		if (f -> problem) continue;
		/* skip degenerate faces */
		/* if (phn -> van_der_Waals && f -> shape != CONVEX) continue; */
		if (!f -> converted) {
			if (f -> shape == CONVEX) {
				sprintf (message, "atom = %d", f -> ptr.atm -> number);
				inform (message);
			}
			else if (f -> shape == CONCAVE) {
				prb = f -> ptr.prb;
				if (prb != NULL) {
					sprintf (message, "%8f %8f %8f center of probe",
						prb -> center[0], prb -> center[1], prb -> center[2]);
				}
				else {
					fac = *(phn -> face_handles + f -> ofn - 1);
					if (fac == NULL) return (0);
					prb = fac -> ptr.prb;
					if (prb != NULL)
						sprintf (message, "%8f %8f %8f center of probe",
							prb -> center[0], prb -> center[1], prb -> center[2]);
					else strcpy(message, "probe pointer is null");
				}
				inform(message);
			}
			if (f -> first_cycle == NULL) {
				sprintf (message, "face has 0 cycles");
			}
			else if (f -> first_cycle -> next != NULL) {
				sprintf (message, "face has > 1 cycles");
			}
			else {
				if (f -> first_cycle -> first_edge == NULL)
					sprintf (message, "null first arc");
				else sprintf (message, "arc: %8ld", f -> first_cycle -> first_edge -> arcptr);
				inform (message);
				sprintf (message, "face has one cycle with %d edges", edges_in_cycle (f -> first_cycle));
			}
			inform (message);
		}
	}
	/* average value of properties computed for adjoining faces */
	for (pv = phn -> head_phnvtx; pv != (struct phnvtx *) NULL; pv = pv -> next) {
		for (j = 0; j < 3; j++)
			if (pv -> degree > 0)
				pv -> values[j] /= pv -> degree;
		n = norm (pv -> outward);
		if (fabs (n) > 0.0)
			normalize (pv -> outward);
	}
	result = polyhedron_handles (phn);
	if (result == 0 || error ()) return (0);
	/* initialize bounds */
	for (k = 0; k < 3; k++) {
		phn -> bounds[0][k] =  1000000.0;
		phn -> bounds[1][k] = -1000000.0;
	}
	for (j = 0; j < 3; j++) {
		phn -> minvals[j] =  1000000.0;
		phn -> maxvals[j] = -1000000.0;
	}
	phn -> converted = 1;
	return (1);
}

int polyhedron_handles (struct surface *phn)
{
	int j, orn;
	long v, e, t, atm;
	long vtx_number[3];
	struct phnvtx *pv;
	struct phnedg *pe;
	struct phntri *pt;
	struct phntri **head_phntri;

	/* second conversion */
	phn -> phnvtx_handles = (struct phnvtx **)
		allocate_pointers (PHNVTX, phn -> n_phnvtx);
	if (phn -> phnvtx_handles == NULL) {
		set_error1("(polyhedron_handles): memory full");
		return(0);
	}
	for (pv = phn -> head_phnvtx, v = 0; pv != (struct phnvtx *) NULL;
		pv = pv -> next, v++) {
		*(phn -> phnvtx_handles + v) = pv;
	}
	phn -> phnedg_handles = (struct phnedg **)
		allocate_pointers (PHNEDG, phn -> n_phnedg);
	if (phn -> phnedg_handles == NULL) {
		set_error1 ("(polyhedron_handles): memory full");
		return(0);
	}
	for (pe = phn -> head_phnedg, e = 0; pe != (struct phnedg *) NULL;
		pe = pe -> next, e++) {
		*(phn -> phnedg_handles + e) = pe;
		for (j = 0; j < 2; j++) {
			pe -> vtxnum[j] =  pe -> pvt[j] -> number;
		}
	}
	phn -> phntri_handles = (struct phntri **)
		allocate_pointers (PHNTRI, phn -> n_phntri);
	if (phn -> phntri_handles == NULL) {
		set_error1 ("(polyhedron_handles): memory full");
		return(0);
	}
	
	if (phn -> heads == NULL) {
		set_error1 ("polyhedron_handles: phn -> heads is null");
		return (0);
	}
	for (atm = 0, t = 0; atm < phn -> n_atom; atm++) {
		head_phntri = phn -> heads + atm;
		if (head_phntri == NULL) {
			set_error1 ("polyhedron_handles: head_phntri is null");
			return (0);
		}
		for (pt = *head_phntri; pt != (struct phntri *) NULL;
			pt = pt -> next, t++) {
			*(phn -> phntri_handles + t) = pt;
			for (j = 0; j < 3; j++) {
				orn = pt -> orns[j];
				pe = pt -> edg[j];
				vtx_number[j] = pe -> pvt[orn] -> number;
			}
			for (j = 0; j < 3; j++) {
				pt -> orn[j] = pt -> orns[j];
				pt -> edgnum[j] =
					(1 - 2 * pt -> orns[j]) * pt -> edg[j] -> number;
			}
			for (j = 0; j < 3; j++)
				pt -> vtxnum[j] = vtx_number[j];
		}
	}

	return (1);
}

int convert_face (struct surface *srf, struct face *f)
{
	int j, n_e;
	char message[MAXLINE];
	struct vertex *vtx;
	struct arc *a;
	struct edge *ed;
	struct cycle *cyc;
	
	if (f -> problem) return(0);
	/* skip degenerate faces */
	/* if (srf -> van_der_Waals && f -> shape != CONVEX) return(0); */
	cyc = f -> first_cycle;
	n_e = 0;
	for (ed = cyc -> first_edge; ed != NULL; ed = ed -> next) {
		n_e++;
	}
	if (n_e < 3) {
		sprintf (message, "(convert_face): face with %d sides", n_e);
		set_error1 (message);
		sprintf (message, "ofn = %ld, lfn = %ld", f -> ofn, f -> lfn);
		set_error2 (message);
		return (0);
	}
	if (n_e > 3) {
		sprintf (message, "(convert_face): face with %d sides", n_e);
		set_error1 (message);
		sprintf (message, "ofn = %ld, lfn = %ld", f -> ofn, f -> lfn);
		set_error2 (message);
		return (0);
	}
	n_e = 0;
	for (ed = cyc -> first_edge; ed != NULL; ed = ed -> next) {
		a = ed -> arcptr;
		if (n_e > 2) {
			sprintf (message, "(convert_face): face with %d sides", n_e);
			set_error1 (message);
			sprintf (message, "ofn = %ld, lfn = %ld", f -> ofn, f -> lfn);
			set_error2 (message);
			return (0);
		}
		for (j = 0; j < 2; j++) {
			vtx = a -> vtx[j];
			if (vtx != (struct vertex *) NULL) {
				if (!convert_vertex (srf, vtx)) {
					inform("convert_face: failure to convert vertex");
					return (0);
				}
			}
		}
		if (!convert_edge (srf, a)) {
			inform("convert_face: failure to convert edge");
			return(0);
		}
		n_e++;
	}
	if (!convert_triangle (srf, f)) {
		inform("convert_face: failure to convert triangle");
		return(0);
	}
	return(1);
}

int convert_vertex (struct surface *srf, struct vertex *vtx)
{
	int k;
	struct phnvtx *pv;
	char message[MAXLINE];

	if (vtx -> converted) return(1);	
	sprintf (message, "convert vertex %8ld", vtx -> number);
	informd2(message);
	pv = allocate_phnvtx ();
	if (pv == NULL) {
		print_counts();
		set_error1 ("(convert_vertex): memory allocation failure");
		return (0);
	}
	if (srf -> head_phnvtx == NULL)
		srf -> head_phnvtx = pv;
	else srf -> tail_phnvtx -> next = pv;
	srf -> tail_phnvtx = pv;
	srf -> n_phnvtx++;
	vtx -> pvt = pv;
	pv -> number = srf -> n_phnvtx;
	for (k = 0; k < 3; k++)
		pv -> center[k] = vtx -> center[k];
	vtx -> converted = 1;
	return(1);
}

int convert_edge (struct surface *srf, struct arc *a)
{
	long ofn;
	int comp, atm, shape;
	double di_center[3];
	char message[MAXLINE];
	struct vertex *vtx0, *vtx1;
	struct variety *vty;
	struct face *fac;
	struct phnedg *pe;
	
	if (a -> converted) return(1);	
	pe = allocate_phnedg ();
	if (pe == NULL) {
		print_counts();
		set_error1 ("(convert_edge): memory allocation failure");
		return (0);
	}
	if (srf -> head_phnedg == (struct phnedg *) NULL)
		srf -> head_phnedg = pe;
	else srf -> tail_phnedg -> next = pe;
	srf -> tail_phnedg = pe;
	srf -> n_phnedg++;
	a -> ped = pe;
	pe -> number = srf -> n_phnedg;
	vtx0 = a -> vtx[0];
	vtx1 = a -> vtx[1];
	if (vtx0 == NULL || vtx1 == NULL) {
		sprintf (message, "(convert_edge) null vertex");
		set_error1 (message);
		return (0);
	}
	if (!vtx0 -> converted) {
		sprintf (message, "(convert_edge) vertex not converted %8ld", vtx0 -> number);
		set_error1 (message);
		sprintf (message, "center: %8.3f %8.3f %8.3f; cusp: %d",
			vtx0 -> center[0], vtx0 -> center[1], vtx0 -> center[2], vtx0 -> cusp);
		set_error2 (message);
		return (0);
	}
	if (!vtx1 -> converted) {
		sprintf (message, "(convert_edge) vertex not converted %8ld", vtx1 -> number);
		set_error1 (message);
		sprintf (message, "center: %8.3f %8.3f %8.3f; cusp: %d",
			vtx1 -> center[0], vtx1 -> center[1], vtx1 -> center[2], vtx1 -> cusp);
		set_error2 (message);
		return (0);
	}
	if (vtx0 -> pvt == NULL) {
		sprintf (message, "(convert_edge) null vertex --> phnvtx %8ld", vtx0 -> number);
		set_error1 (message);
		return (0);
	}
	if (vtx1 -> pvt == NULL) {
		sprintf (message, "(convert_edge) null vertex --> phnvtx %8ld", vtx1 -> number);
		set_error1 (message);
		return (0);
	}

	get_di_center (a -> vtx, di_center);
	ofn = a -> ofn;
	if (ofn <= 0) {
		sprintf(message, "convert_edge: ofn = %8ld", ofn);
		set_error1(message);
		return(0);
	}
	fac = *(srf -> face_handles + ofn - 1);
	shape = fac -> shape;
	vty = fac -> vty;
	atm = point_choice (srf, di_center, vty, shape);
	if (error()) return(0);
	comp = fac -> comp;
	pe -> pvt[0] = vtx0 -> pvt;
	pe -> pvt[1] = vtx1 -> pvt;

	pe -> comp = (short) comp;
	pe -> atm = (short) atm;
	atm = point_choice (srf, vtx0 -> center, vty, shape);
	if (error()) return(0);
	vtx0 -> pvt -> comp = (short) comp;
	vtx0 -> pvt -> atm = (short) atm;
	atm = point_choice (srf, vtx1 -> center, vty, shape);
	if (error()) return(0);
	vtx1 -> pvt -> comp = (short) comp;
	vtx1 -> pvt -> atm = (short) atm;
	a -> converted = 1;
	return (1);
}

int convert_triangle (struct surface *srf, struct face *f)
{
	int j, k, n_e, atm, comp, shape;
	int orn;
	int arc_orn[3];
	long ofn;
	long arc_number[3];
	double radius, circumference;
	double tri_center[3], tri_normal[3];
	char message[MAXLINE];
	char shape_string[24];
	struct vertex *tri_vtx[3];
	struct arc *arc_ptrs[3];
	struct arc *a;
	struct edge *ed;
	struct cycle *cyc, *fcyc;
	struct face *fac;
	struct variety *vty;
	struct phnvtx *pv;
	struct phnedg *pe;
	struct phntri *pt;
	struct phntri **head_phntri;
	struct phntri **tail_phntri;
	
	for (k = 0; k < 3; k++)
		tri_normal[k] = 0.0;
	if (f -> converted) return(1);	
	pt = allocate_phntri ();
	if (pt == NULL) {
		print_counts();
		set_error1 ("(convert_triangle): memory allocation failure");
		print_counts();
		return (0);
	}
	comp = f -> comp;
	ofn = f -> ofn;
	if (ofn <= 0) {
		sprintf(message, "convert_triangle: ofn = %8ld", ofn);
		set_error1(message);
		return(0);
	}
	cyc = f -> first_cycle;
	if (cyc == NULL) {
		inform("convert_triangle: no cycles");
		return (0);
	}
	n_e = 0;
	for (ed = cyc -> first_edge; ed != NULL; ed = ed -> next) {
		n_e++;
	}
	if (n_e > 3) {
		sprintf (message, "(convert_triangle): face with %d sides", n_e);
		set_error1 (message);
		sprintf (message, "ofn = %ld, lfn = %ld", f -> ofn, f -> lfn);
		set_error2 (message);
		return (0);
	}
	n_e = 0;
	for (ed = cyc -> first_edge; ed != NULL; ed = ed -> next) {
		a = ed -> arcptr;
		if (a == NULL) {
			set_error1 ("(convert_triangle): null arc");
			return (0);
		}
		if (n_e > 2) {
			sprintf (message, "(convert_triangle): face with %d sides", n_e);
			set_error1 (message);
			sprintf (message, "ofn = %ld, lfn = %ld", f -> ofn, f -> lfn);
			set_error2 (message);
			return (0);
		}
		orn = ed -> orn;
		arc_orn[n_e] = orn;
		arc_number[n_e] = a -> number;
		arc_ptrs[n_e] = a;
		tri_vtx[n_e] = a -> vtx[orn];
		if (tri_vtx[n_e] == NULL) {
			set_error1 ("(convert_triangle): null vertex");
			return (0);
		}
		n_e++;
	}
	fac = *(srf -> face_handles + ofn - 1);
	if (fac == NULL) {
		inform("convert_triangle: face_handles array has null face pointer");
		return (0);
	}
	fcyc = fac -> first_cycle;
	circumference = circum (fcyc);
	shape = fac -> shape;
	vty = fac -> vty;
	if (vty == NULL) {
		inform("convert_triangle: variety pointer is null");
		return (0);
	}
	radius = vty -> radii[0];
	get_tri_center (tri_vtx, tri_center, n_e);
	atm = point_choice (srf, tri_center, vty, shape);
	if (error()) {
		inform("convert_triangle: point_choice fails");
		return(0);
	}
	if (n_e < 3) {
		sprintf (message,
		"atom %d triangle has only %d sides, %hd cycles, %hd arcs",
			atm, n_e, fac -> n_cycle, fac -> n_arc);
		set_error1 (message);
		if (shape == CONVEX)
			strcpy (shape_string, "face shape = convex");
		else if (shape == SADDLE)
			strcpy (shape_string, "face shape = saddle");
		else if (shape == CONCAVE)
			strcpy (shape_string, "face shape = concave");
		else if (shape == FLAT)
			strcpy (shape_string, "face shape = flat");
		else if (shape == CYLINDRICAL)
			strcpy (shape_string, "face shape = cylindrical");
		sprintf (message,
		"%s, omega = %12.6f, circumference = %12.6f",
		shape_string, fac -> area / (radius * radius), circumference);
		set_error2 (message);
		return (0);
	}
	if(!get_normal (tri_vtx, tri_normal)) {
		informd2("warning: degenerate triangle, no normal computed");
	}
	/* store the data in the new structure */
	for (j = 0; j < 3; j++) {
		pt -> axis[j] = tri_normal[j];
		pt -> edg[j] = arc_ptrs[j] -> ped;
		pt -> orns[j] = arc_orn[j];
	}
	pt -> comp = (short) comp;
	pt -> atm = (short) atm;
	pt -> shape = (short) shape;
	/* property of each of three vertices */
	for (j = 0; j < 3; j++) {
		pe = pt -> edg[j];
		orn = pt -> orns[j];
		pv = pe -> pvt[orn];
		pv -> degree++;
		for (k = 0; k < 3; k++)
			pv -> outward[k] += tri_normal[k];
	}
	f -> converted = 1;
	head_phntri = srf -> heads + atm - 1;
	tail_phntri = srf -> tails + atm - 1;
	if (*head_phntri == (struct phntri *) NULL)
		*head_phntri = pt;
	else (*tail_phntri) -> next = pt;
	*tail_phntri = pt;
	srf -> n_phntri++;
	return(1);
}


int point_choice (struct surface *srf, double point[3], struct variety *vty, int shape)
{
	int j, best_of_several, chose_from;
	double minimum_distance;
	double distance_to_atom[3];
	char message[MAXLINE];
	struct variety *atom_variety[MAXPA];

	/* no computations for convex face */
	if (shape == CONVEX)
		return ((int) vty -> atmnum [0]);

	if (shape == CONVEX)
		chose_from = 1;
	else if (shape == SADDLE)
		chose_from = 2;
	else if (shape == CONCAVE)
		chose_from = 3;
	else if (shape == CYLINDRICAL)
		chose_from = 2;
	else {
		sprintf (message, "(point_choice): invalid shape: %d", shape);
		set_error1(message);
		return (0);
	}
	if (vty -> natom > chose_from)
		chose_from = vty -> natom;

	for (j = 0; j < chose_from; j++) {
		atom_variety[j] = *(
			srf -> variety_handles + vty -> atmnum[j] - 1);
		if (atom_variety[j] -> type != SPHERE) {
			set_error1 ("(point_choice): atom variety not sphere");
			return (0);
		}
	}

	/* find the closest atom */

	for (j = 0; j < chose_from; j++)
		distance_to_atom[j] = distance (point, atom_variety[j] -> center) -
			atom_variety[j] -> radii[0];

	minimum_distance = 1000000.0;
	best_of_several = -1;
	for (j = 0; j < chose_from; j++)
		if (distance_to_atom[j] < minimum_distance) {
			minimum_distance = distance_to_atom[j];
			best_of_several = j;
		}

	if (best_of_several < 0 || best_of_several >= chose_from) {
		set_error1 ("(point_choice): best_of_several not found");
		return (0);
	}

	return ( (int) vty -> atmnum[best_of_several]);
}



double get_normal (struct vertex *tri_vtx[3], double tri_normal[3])
{
	int k;
	double area;
	double vect1[3], vect2[3];

	for (k = 0; k < 3; k++) {
		vect1[k] = tri_vtx[1] -> center[k] - tri_vtx[0] -> center[k];
		vect2[k] = tri_vtx[2] -> center[k] - tri_vtx[0] -> center[k];
	}
	cross (vect1, vect2, tri_normal);
	area = norm (tri_normal) / 2;
	area = fabs (area);
	normalize (tri_normal);
	return (area);
}

void get_tri_center (struct vertex *tri_vtx[3], double tri_center[3], int n_e)
{
	int j, k;

	for (k = 0; k < 3; k++)
		tri_center[k] = 0.0;

	for (j = 0; j < n_e; j++)
		for (k = 0; k < 3; k++)
			tri_center[k] += tri_vtx[j] -> center[k]/n_e;
}


void get_di_center (struct vertex *di_vtx[2], double di_center[3])
{
	int j, k;

	for (k = 0; k < 3; k++)
		di_center[k] = 0.0;

	for (j = 0; j < 2; j++)
		for (k = 0; k < 3; k++)
			di_center[k] += di_vtx[j] -> center[k]/2;
}

