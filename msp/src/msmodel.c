/*
 * Molecular Surface Package
 * Copyright 1986 by Michael L. Connolly
 * All Rights Reserved
 * December 13, 2001
 */

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

	/* for rendering and plotting */
struct surface *make_bas (struct molecule *mol)
{
	long n_inverse, atom_set, bond_set, atom_index, bond_index;
	long index1, index2, atom1, atom2;
	long b, i, j, k, nba, is_tuber;
	long atom_number, bond_number;
	long n_phnvtx, n_phnedg;
	long n_cylinder, n_elbow, n_ball, n_variety;
	int atom_numbers[2];
	unsigned long n_atom, n_bond;
	struct array *inverse_index;
	long bnds[MAX_VALENCE], nbrs[MAX_VALENCE];
	double avg_bond_radius, avg_ball_radius;
	double avg_elbow_radii[2];
	double cyl_length, cyl_radius, cyl_dec1, cyl_dec2;
	double cyl_axis[3], cyl_center[3];
	double cyl_end1[3], cyl_end2[3];
	double tcen[3], taxis[3], tradii[2];
	double ccens[2][3], cradii[2];
	double brads[MAX_VALENCE], acens[MAX_VALENCE][3];
	char message[MAXLINE];
	double *atom_centers;
	double *atom_radii;
	struct atom *ball_ptr, *ball_ptr1, *ball_ptr2;
	struct array *atom_array;
	struct bond *stick_ptr;
	struct array *bond_array;
	struct surface *obj;
	struct face *fac;
	struct variety *vty, *vty1, *vty2;
	struct variety *atom_vty, *bond_vty;
	struct variety **variety_handles;
	struct phnvtx *vtx, *vtx0, *vtx1;
	struct phnedg *edg, *edg0, *edg1;
	struct phnvtx **vertices;
	struct phnedg **edges;

	atom_set = mol -> atom_set;
	bond_set = mol -> bond_set;
	inverse_index = NULL;
	
	/* count atoms */
	n_atom = count_set (atom_set);
	if (n_atom == (long) 0) {
		set_error1 ("(make_bas): empty file");
		return(NULL);
	}
	/* count bonds */
	if (bond_set != 0) n_bond = count_set (bond_set);
	else n_bond = 0;

	atom_array = new_array (ATOM, n_atom);
	if (atom_array == (struct array *) NULL) {
		set_error2 ("(make_bas): memory allocation failure");
		return(NULL);
	}
	if (n_bond > 0) {
		bond_array = new_array (BOND, n_bond);
		if (bond_array == (struct array *) NULL) {
			set_error2 ("(make_bas): memory allocation failure");
			return(NULL);
		}
	}
	else bond_array = NULL;

	/* set up inverse index table from atom number (of all molecules)
		to index for this particular molecule */
	n_inverse = 0;
	for (atom_number = init_for (atom_set), atom_index = 0; atom_number != 0;
		atom_number = next_for (atom_set), atom_index++) {
		if (atom_number > n_inverse) n_inverse = atom_number;
	}
	n_inverse++;
	inverse_index = new_array (INTEGER, n_inverse);
	if (inverse_index == (struct array *) NULL) {
		set_error2 ("(make_bas): memory allocation failure for inverse_index");
		return(NULL);
	}
	for (i = 0; i < n_inverse; i++) {
		store_integer (inverse_index, i, -1);
		if (error ()) return (NULL);
	}

	for (atom_number = init_for (atom_set), atom_index = 0; atom_number != 0;
		atom_number = next_for (atom_set), atom_index++) {
		store_integer (inverse_index, atom_number, atom_index);
		if (error ()) return (NULL);
	}

	/* create ball handles and change atomic radii */
	for (atom_number = init_for (atom_set), atom_index = 0; atom_number != 0;
		atom_number = next_for (atom_set), atom_index++) {
		if (atom_index >= n_atom) break;
		ball_ptr = atom_ptr (atom_number);
		store_atom (atom_array, atom_index, ball_ptr);
		if (error ()) return (NULL);
		/* change van der Waals radius to ball radius */
		ball_ptr -> radius = get_atom_ball (atom_number);
	}

	/* create stick handles */
	bond_index = 0;
	if (n_bond > (unsigned long) 0L) {
		for (bond_number = init_for (bond_set), bond_index = 0; bond_number != 0;
			bond_number = next_for (bond_set), bond_index++) {
			if (bond_index >= n_bond) break;
			cyl_radius = mol -> bond_radius;
			if (cyl_radius <= 0.0) cyl_radius = DEFAULT_BOND_RADIUS;
			get_bond_atoms (bond_number, atom_numbers);
			atom1 = atom_numbers[0]; atom2 = atom_numbers[1];
			if (atom1 == atom2) {
				set_error1 ("(make_bas): cannot bond to self");
				return(NULL);
			}
			index1 = fetch_integer (inverse_index, atom1);
			if (error ()) return (NULL);
			if (index1 < 0) {
				set_error1 ("(make_bas): invalid inverse index");
				return(NULL);
			}
			index2 = fetch_integer (inverse_index, atom2);
			if (error () ) return (NULL);
			if (index2 < 0) {
				set_error1 ("(make_bas): invalid inverse index");
				return(NULL);
			}
			stick_ptr = bond_ptr (bond_number);
			store_bond (bond_array, bond_index, stick_ptr);
			if (error ()) return (NULL);
			stick_ptr -> radius = cyl_radius;
			stick_ptr -> atoms[0] = index1 + 1;
			stick_ptr -> atoms[1] = index2 + 1;
		}
	}

	obj = new_surface ();
	if (obj == NULL) {
		set_error1("(make_bas): mem alloc failure");
		return(NULL);
	}
	/* stick figure figuring (2 edges per bond) */
	n_phnvtx = n_atom + n_bond;
	n_phnedg = 2 * n_bond;
	obj -> n_phnvtx = n_phnvtx;
	obj -> n_phnedg = n_phnedg;
	obj -> n_polygon = 0;
	obj -> type = BAS_SURFACE;
	obj -> scheme = define_scheme (UNIFORM_COLORING, 3, -1, 0.0, 0.0);
	n_elbow = 0; n_ball = 0; n_cylinder = 0;
	n_variety = n_atom + n_bond;

	/* allocate memory for vertex pointers */
	vertices = (struct phnvtx **)
		allocate_pointers (PHNVTX, n_phnvtx);
	if (vertices == NULL) {
		set_error1 ("(make_bas): memory full");
		return(NULL);
	}
	/* allocate memory for edge pointers */
	edges = (struct phnedg **)
		allocate_pointers (PHNEDG, n_phnedg);
	if (edges == NULL) {
		set_error1 ("(make_bas): memory full");
		return(NULL);
	}
	atom_centers = allocate_doubles (n_atom * 3, 0, ATOM_CENTERS);
	if (atom_centers == NULL) {
		set_error1("(make_bas) memory allocation failure");
		return(NULL);
	}
	atom_radii = allocate_doubles (n_atom, 0, ATOM_RADII);
	if (atom_radii == NULL) {
		set_error1("(make_bas) memory allocation failure");
		return(NULL);
	}
	
	variety_handles = (struct variety **)
		allocate_pointers (VARIETY, (n_atom + n_bond));
	if (variety_handles == NULL) {
		set_error1("(make_bas) memory allocation failure");
		return(NULL);
	}
	for (i = 0; i < n_variety; i++) {
		vty = allocate_variety ();
		if (vty == (struct variety *) NULL) {
			set_error1 ("(make_bas): variety allocation failure");
			return (0);
		}
		*(variety_handles + i) = vty;
	}

	for (atom_index = 0; atom_index < n_atom; atom_index++) {
		vtx = allocate_phnvtx ();
		if (vtx == (struct phnvtx *) NULL) {
			set_error1 ("(make_bas): vertex allocation failure");
			return(NULL);
		}
		*(vertices + atom_index) = vtx;
		ball_ptr = fetch_atom (atom_array, atom_index);	
		if (error ()) return (NULL);
		vtx -> center[0] = ball_ptr -> center[0];
		vtx -> center[1] = ball_ptr -> center[1];
		vtx -> center[2] = ball_ptr -> center[2];
		vtx -> values[0] = ball_ptr -> radius;
	}
	for (bond_index = 0; bond_index < n_bond; bond_index++) {
		stick_ptr = fetch_bond (bond_array, bond_index);
		if (error ()) return (NULL);
		/* stick_ptr -> radius = cyl_radius; */
		atom1 = stick_ptr -> atoms[0];
		atom2 = stick_ptr -> atoms[1];
		ball_ptr1 = fetch_atom (atom_array, atom1 - 1);
		if (error ()) return (NULL);
		ball_ptr2 = fetch_atom (atom_array, atom2 - 1);
		if (error ()) return (NULL);
		vtx0 = *(vertices + atom1 - 1);
		vtx1 = *(vertices + atom2 - 1);
		for (k = 0; k < 3; k++) {
			cyl_center[k] = (vtx0 -> center[k] + vtx1 -> center[k])/2;
		}
		/* bond midpoint */
		vtx = allocate_phnvtx ();
		if (vtx == (struct phnvtx *) NULL) {
			set_error1 ("(make_bas): vertex allocation failure");
			return(NULL);
		}
		*(vertices + n_atom + bond_index) = vtx;
		vtx -> center[0] = cyl_center[0];
		vtx -> center[1] = cyl_center[1];
		vtx -> center[2] = cyl_center[2];
		vtx -> values[0] = stick_ptr -> radius;
		/* first edge */
		edg0 = allocate_phnedg ();
		if (edg0 == (struct phnedg *) NULL) {
			set_error1 ("(make_bas): edge allocation failure");
			return(NULL);
		}
		*(edges + 2 * bond_index) = edg0;
		edg0 -> pvt[0] = vtx0;
		edg0 -> pvt[1] = vtx;
		edg0 -> vtxnum[0] = atom1;
		edg0 -> vtxnum[1] = n_atom + bond_index + 1;
		edg0 -> comp = 1;
		edg0 -> atm = (unsigned short) atom1;
		edg0 -> hue = 3;
		/* second edge */
		edg1 = allocate_phnedg ();
		if (edg1 == (struct phnedg *) NULL) {
			set_error1 ("(make_bas): edge allocation failure");
			return(NULL);
		}
		*(edges + 2 * bond_index + 1) = edg1;
		edg1 -> pvt[0] = vtx;
		edg1 -> pvt[1] = vtx1;
		edg1 -> vtxnum[0] = n_atom + bond_index + 1;
		edg1 -> vtxnum[1] = atom2;
		edg1 -> comp = 1;
		edg1 -> atm = (unsigned short) atom2;
		edg1 -> hue = 3;
	}
	/* set up linked list */
	obj -> head_phnedg = *edges;
	for (i = 0; i < obj -> n_phnedg - 1; i++) {
		edg = *(edges + i);
		edg -> next = *(edges + i + 1);
	}

	/* set up linked list */
	obj -> head_phnvtx = *vertices;
	for (i = 0; i < obj -> n_phnvtx - 1; i++) {
		vtx = *(vertices + i);
		vtx -> next = *(vertices + i + 1);
	}
	avg_ball_radius = 0.0;
	for (j = 0; j < 2; j++)
		avg_elbow_radii[j] = 0.0;
	for (atom_index = 0; atom_index < n_atom; atom_index++) {
		ball_ptr = fetch_atom (atom_array, atom_index);
		if (error ()) return (NULL);
		/* store atomic coordinates */
		for (k = 0; k < 3; k++)
			*(atom_centers + 3 * (atom_index) + k) = ball_ptr -> center[k];
		*(atom_radii + (atom_index)) = ball_ptr -> radius;
		atom_vty = *(variety_handles + atom_index);
		/* find bonds to this atom */
		nba = 0;
		for (b = 1; b <= n_bond; b++) {
			stick_ptr = fetch_bond (bond_array, b - 1);
			if (error ()) return (NULL);
			atom1 = stick_ptr -> atoms[0];
			atom2 = stick_ptr -> atoms[1];
			if (atom1 == atom_index+1) {
				if (nba >= MAX_VALENCE) break;
				nbrs[nba] = atom2;
				bnds[nba] = b;
				nba++;
			}
			else if (atom2 == atom_index+1) {
				if (nba >= MAX_VALENCE) break;
				nbrs[nba] = atom1;
				bnds[nba] = b;
				nba++;
			}
		}
		is_tuber = 0;
		if (nba >= 2) {
			for (b = 0; b < nba; b++) {
				atom2 = nbrs[b];
				bond_number = bnds[b];
				ball_ptr2 = fetch_atom (atom_array, atom2 - 1);
				if (error ()) return (NULL);
				stick_ptr = fetch_bond (bond_array, bond_number - 1);
				if (error ()) return (NULL);
				brads[b] = stick_ptr -> radius;
				for (k = 0; k < 3; k++)
					acens[b][k] = ball_ptr2 -> center[k];
			}
			/* sort into order of decreasing stick radius */
			sort_nbrs (nba, acens, brads, nbrs, bnds);
			is_tuber = tuber (mol -> elbow, nba, ball_ptr -> center,
				acens, brads, tcen, taxis, tradii, ccens, cradii);
		}
		if (is_tuber) {
			n_elbow++;
			atom_vty -> radii[0] = tradii[0];
			for (k = 0; k < 3; k++)
				atom_vty -> center[k] = tcen[k];
			atom_vty -> type = TORUS;
			atom_vty -> tube = 1;
			for (k = 0; k < 3; k++)
				atom_vty -> axis[k] = taxis[k];
			for (j = 0; j < 2; j++) {
				for (k = 0; k < 3; k++) {
					atom_vty -> ccens[j][k] = ccens[j][k];
				}
				atom_vty -> radii[1] = cradii[j]; /* sic */
			}
			/* fix up sticks */
			for (b = 0; b < 2; b++) {
				bond_number = bnds[b];
				stick_ptr = fetch_bond (bond_array, bond_number - 1);
				if (error ()) return (NULL);
				if (stick_ptr -> atoms[0] == atom_index+1) {
					stick_ptr -> modified1 = 1;
					for (k = 0; k < 3; k++)
						stick_ptr -> end1[k] = ccens[b][k];
				}
				else {
					stick_ptr -> modified2 = 1;
					for (k = 0; k < 3; k++)
						stick_ptr -> end2[k] = ccens[b][k];
				}
			}
			for (j = 0; j < 2; j++)
				avg_elbow_radii[j] += atom_vty -> radii[j];
		}
		else {
			n_ball++;
			atom_vty -> radii[0] = ball_ptr -> radius;
			atom_vty -> center[0] = ball_ptr -> center[0];
			atom_vty -> center[1] = ball_ptr -> center[1];
			atom_vty -> center[2] = ball_ptr -> center[2];
			atom_vty -> type = SPHERE;
			for (k = 0; k < 3; k++) {
				atom_vty -> axis[k] = 0.0;
			}
			for (k = 0; k < MAXPA; k++)
				atom_vty -> atmnum[k] = 0;
		}
		atom_vty -> atmnum[0] = (short) atom_index+1;
		atom_vty -> atmnum[1] = (short) atom_index+1;		/* for toroidal tubes */
		avg_ball_radius += ball_ptr -> radius;

		fac = allocate_face ();
		if (fac == NULL) {
			set_error1("(make_bas): mem alloc failure");
			return(NULL);
		}
		if (obj -> head_face == NULL) obj -> head_face = fac;
		else obj -> tail_face -> next = fac;
		obj -> tail_face = fac;
		fac -> vty = atom_vty;
		fac -> shape = CONVEX;
		fac -> next = NULL;
		fac -> comp = 1;
		fac -> first_cycle = NULL;
	}

	avg_bond_radius = 0.0;
	for (bond_index = 0; bond_index < n_bond; bond_index++) {
		stick_ptr = fetch_bond (bond_array, bond_index);
		if (error ()) return (NULL);
		atom1 = stick_ptr -> atoms[0];
		atom2 = stick_ptr -> atoms[1];
		cyl_radius = stick_ptr -> radius;
		if (atom1 == atom2) {
			set_error1("(make_bas): cannot bond to self");
			return(NULL);
		}
		ball_ptr1 = fetch_atom (atom_array, atom1 - 1);
		if (error ()) return (NULL);
		ball_ptr2 = fetch_atom (atom_array, atom2 - 1);
		if (error ()) return (NULL);
		vty1 = *(variety_handles + atom1 - 1);
		vty2 = *(variety_handles + atom2 - 1);
		bond_vty = *(variety_handles + n_atom + bond_index);
		for (k = 0; k < 3; k++) {
			cyl_center[k] = (ball_ptr1 -> center[k] + ball_ptr2 -> center[k])/2;
			cyl_axis[k] = (ball_ptr2 -> center[k] - ball_ptr1 -> center[k]);
		}
		cyl_length = norm (cyl_axis);
		if (cyl_length <= 0.0) {
			set_error1("(make_bas): zero-length unmodified bond");
			sprintf (message, "atom1 = %ld, atom2 = %ld", atom1, atom2);
			set_error2 (message);
			return(NULL);
		}
		normalize (cyl_axis);
		/* modify so cylinder ends at spheres */
		cyl_dec1 = ball_ptr1 -> radius * ball_ptr1 -> radius -
			cyl_radius * cyl_radius;
		if (cyl_dec1 <= 0.0)
			cyl_dec1 = 0.0;
		else cyl_dec1 = sqrt (cyl_dec1);
		cyl_dec2 = ball_ptr2 -> radius * ball_ptr2 -> radius -
			cyl_radius * cyl_radius;
		if (cyl_dec2 <= 0.0)
			cyl_dec2 = 0.0;
		else cyl_dec2 = sqrt (cyl_dec2);

		/* three possibilities for each end of the bond */
		if (stick_ptr -> modified1)
			for (k = 0; k < 3; k++)
				cyl_end1[k] = stick_ptr -> end1[k];
		else if (vty1 -> type == TORUS)
			for (k = 0; k < 3; k++)
				cyl_end1[k] = ball_ptr1 -> center[k];
		else
			for (k = 0; k < 3; k++)
				cyl_end1[k] = ball_ptr1 -> center[k] + cyl_dec1 * cyl_axis[k];
				
		if (stick_ptr -> modified2)
			for (k = 0; k < 3; k++)
				cyl_end2[k] = stick_ptr -> end2[k];
		else if (vty2 -> type == TORUS)
			for (k = 0; k < 3; k++)
				cyl_end2[k] = ball_ptr2 -> center[k];
		else
			for (k = 0; k < 3; k++)
				cyl_end2[k] = ball_ptr2 -> center[k] - cyl_dec2 * cyl_axis[k];
				
		for (k = 0; k < 3; k++)
			cyl_center[k] = (cyl_end1[k] + cyl_end2[k])/2;
		cyl_length = distance (cyl_end1, cyl_end2);
		if (cyl_length <= 0.0) {
			set_error1("(make_bas): zero-length modified bond");
			sprintf (message, "atom1 = %d, atom2 = %d", atom1, atom2);
			set_error2 (message);
			return(NULL);
		}
		bond_vty -> radii[0] = cyl_radius;
		bond_vty -> length = cyl_length;
		bond_vty -> type = CYLINDER;
		for (k = 0; k < 3; k++) {
			bond_vty -> center[k] = cyl_center[k];
			bond_vty -> axis[k] = cyl_axis[k];
		}
		for (k = 0; k < MAXPA; k++)
			bond_vty -> atmnum[k] = 0;
		bond_vty -> atmnum[0] = (short) atom1;
		bond_vty -> atmnum[1] = (short) atom2;
		avg_bond_radius += bond_vty -> radii[0];

		fac = allocate_face ();
		if (fac == NULL) {
			set_error1("(make_bas): mem alloc failure");
			return(NULL);
		}
		if (obj -> head_face == NULL) obj -> head_face = fac;
		else obj -> tail_face -> next = fac;
		obj -> tail_face = fac;
		fac -> vty = bond_vty;
		fac -> shape = CONVEX;
		fac -> next = NULL;
		fac -> comp = 1;
		fac -> first_cycle = NULL;
		n_cylinder++;
	}

	avg_ball_radius /= (double) (int) n_atom;
	if (n_bond > (long) 0) 
		avg_bond_radius /= (double) (int) n_bond;
	for (j = 0; j < 2; j++) {
			if (n_elbow > (long) 0) 
				avg_elbow_radii[j] /= (double) (int) n_elbow;
			else avg_elbow_radii[j] = 0.0;
	}
	sprintf (message,"%8ld balls with average radius %6.3f", n_ball, avg_ball_radius);
	inform(message);
	sprintf (message,"%8ld elbows with average radii %6.3f %6.3f",
		n_elbow, avg_elbow_radii[0], avg_elbow_radii[1]);
	if (n_elbow > 0) inform(message);
	sprintf (message,"%8ld cylinders with average radius %6.3f", n_cylinder, avg_bond_radius);
	if (n_cylinder > 0) inform(message);

	/* transfer to structure variables */
	obj -> n_atom = n_atom;
	obj -> atom_centers = atom_centers;
	obj -> atom_radii = atom_radii;
	obj -> variety_handles = variety_handles;
	obj -> circle_handles = NULL;
	obj -> vertex_handles = NULL;
	obj -> arc_handles = NULL;
	obj -> face_handles = NULL;
	obj -> cycle_handles = NULL;
	obj -> edge_handles = NULL;
	obj -> n_variety = n_variety;
	obj -> n_circle = 0;
	obj -> n_vertex = 0;
	obj -> n_arc = 0;
	obj -> n_face = 0;
	obj -> n_cycle = 0;
	obj -> n_edge = 0;
	obj -> n_component = 1;
	obj -> clipping = 1;
	obj -> probe_radius = 0.0;
	obj -> phnvtx_handles = vertices;
	obj -> phnedg_handles = edges;
	
	if (atom_array != (struct array *) NULL) free_array (atom_array);
	if (bond_array != (struct array *) NULL) free_array (bond_array);
	if (inverse_index != NULL) free_array (inverse_index);
	return (obj);
}

void sort_nbrs (long nba, double acens[MAX_VALENCE][3], double brads[MAX_VALENCE],
long nbrs[MAX_VALENCE], long bnds[MAX_VALENCE])
{
	int j, k, change, n_iter;
	double temp; int itemp;
	
	
	change = 1;
	n_iter = 0;
	while (change) {
		n_iter++; if (n_iter > 1000) break;
		change = 0;
		for (j = 0; j < nba - 1; j++) {
			n_iter++; if (n_iter > 1000) break;
			if (brads[j] < brads[j+1]) {
				/* swap */
				change = 1;
				temp = brads[j];
				brads[j] = brads[j+1];
				brads[j+1] = temp;
				for (k = 0; k < 3; k++) {
					temp = acens[j][k];
					acens[j][k] = acens[j+1][k];
					acens[j+1][k] = temp;
				}
				itemp = nbrs[j];
				nbrs[j] = nbrs[j+1];
				nbrs[j+1] = itemp;
				itemp = bnds[j];
				bnds[j] = bnds[j+1];
				bnds[j+1] = itemp;
			}
		}
	}
}


/* brads must be provided in decreasing order */


int tuber (double elbow, int nbond, double acen[3], double acens[MAX_VALENCE][3], double brads[MAX_VALENCE],
double tcen[3], double taxis[3], double tradii[2], double ccens[2][3], double cradii[2])
{
	int k;
	double r, rt, lb1, lb2, dot12, beta;
	double thb, chb, check;
	double vect1[3], vect2[3], uvect1[3], uvect2[3];
	double perp[3], inward[3];
	
	/* torus radius must be at least the bond radius */
	if (elbow <= 1.0) return (0);
	/* there must be at least two bonds to the central atom */
	if (nbond < 2) return (0);
	/* and they must have the same radii */
	if (brads[0] != brads[1]) return (0);
	/* if there are more than two bonds to this atom,
	   the first two must have larger bond radii than the remainder */
	if (nbond > 2) {
		if (brads[2] >= brads[1]) return (0);
	}
	for (k = 0; k < 3; k++) {
		vect1[k] = acens[0][k] - acen[k];
		vect2[k] = acens[1][k] - acen[k];
	}
	lb1 = norm (vect1); if (lb1 <= 0.0) return (0);
	lb2 = norm (vect2); if (lb2 <= 0.0) return (0);
	for (k = 0; k < 3; k++) {
		uvect1[k] =  vect1[k] / lb1;
		uvect2[k] =  vect2[k] / lb2;
	}
	/* compute vector perpendicular to bond plane */
	cross (uvect1, uvect2, perp);
	if (norm(perp) <= 0.0) {
		inform ("(tuber): collinear vectors");
		return (0);
	}
	if (!normalize (perp)) return (0);
	cross (perp, uvect1, inward);
	normalize (inward);
	
	/* compute angle between bonds */
	dot12 = dot_product (uvect1, uvect2);
	if (dot12 < -1.0) dot12 = -1.0;
	if (dot12 > 1.0) dot12 = 1.0;
	beta = (double) acos (dot12);
	if (beta <= 0.5) return (0);
	thb = tan (beta/2.0);
	if (thb <= 0.0) return (0);
	chb = 1.0/thb;
	r = brads[0];
	rt = elbow * r;
	if (0.5 * lb1 <= rt * chb) return (0);
	if (0.5 * lb2 <= rt * chb) return (0);
	
	/* set up circles at end of toroidal sector */
	for (k = 0; k < 3; k++) {
		ccens[0][k] = acen[k] + rt * chb * uvect1[k];
		ccens[1][k] = acen[k] + rt * chb * uvect2[k];
	}
	for (k = 0; k < 3; k++)
		tcen[k] = ccens[0][k] + rt * inward[k];
	/* check that distance from torus center to other circle center is torus radius */
	check = distance (tcen, ccens[1]);
	if (fabs (check - rt) > 0.01) {
		set_error1("(tuber): bad torus radius");
		return (0);
	}
	/* set up torus geometry */
	for (k = 0; k < 3; k++)
		taxis[k] = -perp[k];
	tradii[0] = rt;
	tradii[1] = r;
	cradii[0] = r;
	cradii[1] = r;
	/* return true (there is an elbow) */
	return (1);
}


