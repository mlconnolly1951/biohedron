/*

	MSRoll

	Copyright 1986, 1989, 1996 by Michael L. Connolly
	All rights reserved

	Written by Michael L. Connolly.
	January 8, 2002

*/

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"


int setup_alpha (struct molecule *mol)
{
	int atm;
	int atom_set;
	int atom_number;
	double atom_alpha;
	char message[MAXLINE];
    struct cept *ex;

	if (mol -> n_atom < 1) {
		ex = new_cept (LOGIC_ERROR,  INVALID_VALUE,  FATAL_SEVERITY);
		add_function (ex, "setup_alpha");
		add_source (ex, "msalpha.c");
		add_long (ex, "n_atom", mol  -> n_atom);
		return(0);
	}
	atom_set = mol -> atom_set;
	if (atom_set < 1) {
		ex = new_cept (LOGIC_ERROR,  INVALID_VALUE,  FATAL_SEVERITY);
		add_function (ex, "setup_alpha");
		add_source (ex, "msalpha.c");
		add_long (ex, "atom_set", atom_set);
		return(0);
	}
	mol -> atom_alphas = allocate_doubles (mol -> n_atom, 0, ATOM_ALPHAS);
	if (mol -> atom_alphas == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
        add_source (ex, "msalpha.c");
        add_function (ex, "setup_alpha");
        add_variable (ex, ATOM_ALPHAS, "atom_alphas");
        add_long (ex, "n_atom", n_atom);
		return(0);
	}
	for (atom_number = init_for (atom_set), atm = 0; atom_number != 0;
		atom_number = next_for (atom_set), atm++) {
		if (atm >= mol -> n_atom) break;
		atom_alpha = get_atom_angle (atom_number);
		if (error()) return (0);
		sprintf (message, "alpha = %8.3f", atom_alpha);
		informd2 (message);
		if (atom_alpha < MIN_ALPHA) {
            atom_alpha = MIN_ALPHA;
		}
		if (atom_alpha > MAX_ALPHA) {
            atom_alpha = MAX_ALPHA;
		}
		*(mol -> atom_alphas + atm) = atom_alpha;
	}

	return (1);
}

/* modify atom alphas */
void face_alpha (struct molecule *mol, struct surface *srf)
{
	int k;
	int atm, n_share;
	char message[MAXLINE];
	struct face *fac;
	struct variety *vty;
	struct cycle *cyc;
	struct edge *e;
	struct arc *a;
	struct component *cmp_ptr;
    struct cept *ex;

	for (fac = srf -> head_face; fac != NULL; fac = fac -> next) {
		cmp_ptr = *(srf -> component_handles + fac -> comp - 1);
		if (cmp_ptr == NULL) {
			ex = new_cept (LOGIC_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
			add_function (ex, "face_alpha");
			add_source (ex, "msalpha.c");
            add_message (ex, "null component pointer");
			return;
		}
		/* set up face alphas */
		/* this code depends on values of shape being 1, 2, 3 */
		vty = fac -> vty;
		fac -> alpha = 0.0;
		if (fac -> shape == CONVEX)
			n_share = 1;
		else if (fac -> shape == SADDLE)
			n_share = 2;
		else if (fac -> shape == CONCAVE)
			n_share = 3;
		else if (fac -> shape == CYLINDRICAL)
			n_share = 2;
		else {
            n_share = 1;
		}
		if (vty -> natom > n_share) n_share = vty -> natom;
		for (k = 0; k < n_share; k++) {
			atm = vty -> atmnum[k] - 1;
			if (atm < 0 || atm >= mol -> n_atom) {
				ex = new_cept (LOGIC_ERROR,  INVALID_VALUE,  FATAL_SEVERITY);
				add_function (ex, "face_alpha");
				add_source (ex, "msalpha.c");
				add_long (ex, "atom_number", atm + 1);
				add_long (ex, "n_atom", mol  -> n_atom);
				return;
			}
			fac -> alpha += *(mol -> atom_alphas + atm);
		}
		fac -> alpha /= n_share;
		/* override for internal cavities */
		if (mol -> cavity_alpha > 0.0 && cmp_ptr -> volume < 0.0)
			fac -> alpha = mol -> cavity_alpha;
		if (fac -> alpha < MIN_ALPHA) fac -> alpha = MIN_ALPHA;
		if (fac -> alpha > MAX_ALPHA) fac -> alpha = MAX_ALPHA;
		/* set up arc alphas */
		for (cyc = fac -> first_cycle; cyc != NULL; cyc = cyc -> next)
			for (e = cyc -> first_edge; e != NULL; e = e -> next) {
				a = e -> arcptr;
				if (e -> orn) a -> alpha = fac -> alpha;
			}
	}
}


/*
   MSRoll
   Copyright 1996 by Michael L. Connolly
   All Rights Reserved

*/
