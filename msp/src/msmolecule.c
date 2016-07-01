#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

/* MSP Copyright 1995 by Michael L. Connolly */
/* March 6, 2000 */

struct molecule *read_molecule (struct msscene *ms, char *name, char molfilename[], double alpha, FILE *fpm, FILE *fpr, FILE *fpy, FILE *fpf)
{
	int k, nw, result;
	long atom_number;
	long atom_set, n_atom, l, nri;
	long n_duplicate;
	double ball_radius;
	char molname[64], compound_name[64];
	char message[MAXLINE];
	struct array *pdb_lines;
	struct array *radius_lines, *pattern_lines;
	struct array *atom_array;
	struct molecule *mol;
	struct atom *a, *pa, *head_atom, *tail_atom;

	/* read patterns from disk (or else use default) */
	if (fpy != NULL) {
		pattern_lines = read_lines (fpy);
		if (pattern_lines == NULL) return (0);
		fclose (fpy);
		result = read_pattern_file (pattern_lines);
		if (error()) return (NULL);
		free_array (pattern_lines);
		pattern_lines = NULL;
		nri = create_residue_index ();
		if (error ()) return (NULL);
		sprintf (message, "%8ld residue indices for pattern table", nri);
		inform (message);
	}
	/* read radii lines from disk (or else use default) */
	if (fpr != NULL) {
		radius_lines = read_lines (fpr);
		if (radius_lines == NULL) return (NULL);
		fclose (fpr);
		result = read_radii_file (radius_lines);
		if (error()) return (NULL);
		free_array (radius_lines);
	}
	radius_lines = NULL;
	/* read pdb file into string array */
	pdb_lines = read_lines (fpm);
	if (pdb_lines == NULL) return (NULL);
	fclose (fpm);
	if (name == NULL || strlen (name) <= 0) {
		if (molfilename != NULL && strlen (molfilename) > 0)
			setup_default_name (molfilename, molname);
		else strcpy (molname, "unknown");
	}
	else strcpy (molname, name);
	if (strlen(molname) < 1) {
		set_error1 ("read_molecule: null molecule name");
		return (NULL);
	}
	setup_compound_name (pdb_lines, compound_name);
	for (mol = ms -> head_molecule; mol != NULL; mol = mol -> next) {
		if (strcmp (mol -> name, molname) == 0) {
			set_error1 ("duplicate molecule name");
			return (NULL);
		}
	}
	atom_set = named_set (molname, ATOM);

	/* read atoms from string array */
	read_atom_file (pdb_lines, "pdb", molname);
	if (error()) return (NULL);
	free_array (pdb_lines);
	pdb_lines = NULL;

	/* possibly restrict the set of atoms */
	/* also tessellation fineness */
	if (fpf != NULL) {
		read_atom_commands (ms, fpf);
		if (error()) return (NULL);
	}

	/* remove duplicate atoms (same x, y, z) */
	n_duplicate = remove_duplicate_atoms (atom_set);
	if (error ()) return (NULL);
	sprintf (message, "%8ld duplicate atoms removed from molecule", n_duplicate);
	if (n_duplicate > 0) inform (message);

	atom_array = read_atoms (atom_set);
	if (atom_array == NULL) return (NULL);
	n_atom = atom_array -> length;
	if (n_atom <= 0) {
		set_error1 ("read_molecule: no atoms in input file");
		return(NULL);
	}
	
	/* linked list of atoms */
	a = pa = (struct atom *) NULL;
	head_atom = tail_atom = (struct atom *) NULL;
	for (l = 0; l < n_atom; l++) {
		a = fetch_atom (atom_array, l);
		if (head_atom == NULL) head_atom = a;
		else pa -> next = a;
		pa = a;
	}
	tail_atom = pa;
	
	/* new molecule */
	mol = new_molecule (n_atom, 0L, head_atom, tail_atom);
	if (error()) return (NULL);
	if (mol == NULL) {
		set_error1("molecule memory allocation fails");
		return (NULL);
	}
	strcpy(mol -> name, molname);
	mol -> atom_set = atom_set;
	ms -> current_molecule = mol;
	if (ms -> head_molecule == NULL)
		ms -> head_molecule = mol;
	else ms -> tail_molecule -> next = mol;
	ms -> tail_molecule = mol;
	mol -> head_surface = NULL;
	mol -> tail_surface = NULL;
	for (k = 0; k < 3; k++) {
		mol -> center[k] = 0.0;
	}
	/* setup ball radii and subdivision angle parameters */
	ball_radius = mol -> ball_radius;
	setup_ball_radii (ball_radius, atom_set);
	if (alpha > 0.0) {
		for (atom_number = init_for (atom_set); atom_number != 0;
			atom_number = next_for (atom_set)) {
			set_atom_angle (atom_number, alpha);
			if (error()) return (NULL);
		}
	}
	if (alpha > 0.0)
		mol -> alpha = alpha;
	else mol -> alpha = DEFAULT_ALPHA;
	mol -> global_alpha = mol -> alpha;
	mol -> cavity_alpha = DEFAULT_CAVITY_ALPHA;
	mol -> fp_alpha = fpf;
	mol -> n_atom = count_set (atom_set);
	sprintf (message, "%8ld atoms in molecule %-s %-s",
			mol -> n_atom, mol -> name, compound_name);
	inform (message);
	free_array (atom_array);
	return (mol);
}

int setup_compound_name (struct array *pdb_lines, char compound_name[])
{
	int i, m, nw, moloff, npunct;
	long l;
	char id[6];
	char word1[MAXLINE], word2[MAXLINE], word3[MAXLINE];
	char *in_line;
	
	strcpy (compound_name, "");
	/* look for compound record */
	for (l = 0; l < pdb_lines -> length; l++) {
		in_line = fetch_string (pdb_lines, l);
		strcpy (id, "");
		strcpy (word1, "");
		strcpy (word2, "");
		strcpy (word3, "");
		nw = sscanf (in_line, "%s %s %s", word1, word2, word3);
		if (nw < 1) continue;
		moloff = 0;
		if (strcmp (word1, "COMPND") == 0) {
			/* replace pathname by pdb name */
			strncpy (id, in_line + 72, 4);
			id[4] = 0;
			if (isdigit (id[0]) && isalpha(id[1]) &&
				isalpha(id[2]) && isalpha(id[3])) {
				strcpy (compound_name, id);
				compound_name[4] = ' ';
				compound_name[5] = 0;
				moloff = 5;
			}
			if (nw >= 2)
				strncpy (compound_name+moloff, word2, 24);
			if (nw >= 3 && strcmp (id, word3) != 0) {
				moloff = strlen (compound_name);
				compound_name[moloff] = '_';
				moloff++;
				/* look for leading punctuation */
				if (ispunct (word3[0])) m = 1;
				else m = 0;
				npunct = 0;
				/* look for interior punctuation */
				for (i = 1; i < strlen (word3) - 1; i++)
					if (ispunct (word3[i])) npunct++;
				if (npunct == 0)
						strncpy (compound_name+moloff, word3+m, 24);
			}
			moloff = strlen (compound_name);
			/* look for trailing punctuation */
			if (moloff > 4 && ispunct (compound_name[moloff-1]))
				compound_name[moloff-1] = 0;
			moloff = strlen (compound_name);
			for (l = 0; l < moloff; l++) {
				if (isupper (compound_name[l]))
					compound_name[l] = tolower(compound_name[l]);
			}
			return (1);
		}
	}
	return (0);
}


struct array *read_atoms (long atom_set)
{
	long o, n;
	long n_atom;
	struct array *atom_array;
	struct atom *a;
	
	n = count_set (atom_set);
	if (n <= 0) return (NULL);
	atom_array = new_array (ATOM, n);
	if (atom_array == NULL) {
		set_error2 ("read_atoms: memory failure");
		return (NULL);
	}
	/* interpret input atom lines */
	n_atom = 0L;
	for (o = init_for (atom_set); o != 0; o = next_for (atom_set)) {
		/* retrieve PSL structure instance for atomic data */
		a = atom_ptr (o);
		if (a == NULL) {
			set_error1 ("read_atoms: ran out of memory");
			return(NULL);
		}
		/* store handle */
		store_atom (atom_array, n_atom, a);
		n_atom++;
		a -> number = n_atom;
	}
	if (n_atom <= 0) {
		set_error1 ("read_atoms: no atoms in input file");
		return(NULL);
	}
	return (atom_array);
}

void setup_ball_radii (double ball_radius, int atom_set)
{
	int o;

	for (o = init_for (atom_set); o != 0; o = next_for (atom_set)) {
		set_atom_ball (o, ball_radius);
	}
}


/* set atom opacities (opacities) */

void set_atom_opacities (struct molecule *mol, struct surface *obj)
{
	int i, j;
	double *atom_opacities;
	double opacity;
	long atom_set, atom_number;

	obj -> scheme -> opacity_type = ATOM_OPACITY;
	/* allocate memory for atom opacities */
	atom_opacities = allocate_doubles (obj -> n_atom * 4, 0, ATOM_OPACITIES);
	if (atom_opacities == NULL) {
		set_error1 ("(set_atom_opacities): memory allocation failure");
		return;
	}
	obj -> scheme -> atom_opacities = atom_opacities;

	/* one input line for each atom */

	atom_set = mol -> atom_set;
	for (atom_number = init_for (atom_set), i = 0; atom_number != 0;
		atom_number = next_for (atom_set), i++) {
		if (i >= obj -> n_atom) break;
		opacity = get_atom_opacity (atom_number);
		/* contact, reentrant */
		*(atom_opacities+4*i) = opacity;
		*(atom_opacities+4*i+1) = opacity;
		*(atom_opacities+4*i+2) = opacity;
		*(atom_opacities+4*i+3) = opacity;
		/* error checking */
		for (j = 0; j < 4; j++) {
			if (*(atom_opacities+4*i+j) < 0.0) {
				set_error1 ("(set_atom_opacities): invalid atom opacity");
				return;
			}
			if (*(atom_opacities+4*i+j) > 1.0) {
				set_error1 ("(set_atom_opacities): invalid atom opacity");
				return;
			}
		}
	}
}

struct molecule *new_molecule (long n_atom, long n_bond, struct atom *head_atom, struct atom *tail_atom)
{
	struct molecule *mol;

	mol = (struct molecule *) allocate_object (MOLECULE);
	if (mol == NULL) {
		set_error1 ("new_molecule: memory allocation failure");
		return(NULL);
	}
	mol -> n_atom = n_atom;
	mol -> n_bond = n_bond;
	mol -> head_atom = head_atom;
	mol -> tail_atom = tail_atom; 
	mol -> tolerance = DEFAULT_TOLERANCE;
	mol -> elbow = DEFAULT_ELBOW;
	mol -> bond_radius = DEFAULT_BOND_RADIUS;
	mol -> ball_radius = DEFAULT_BALL_RADIUS;
	mol -> outer_width = DEFAULT_OUTER_WIDTH;
	mol -> inner_width = DEFAULT_INNER_WIDTH;
	mol -> cavity_width = DEFAULT_CAVITY_WIDTH;
	mol -> bond_width = DEFAULT_BOND_WIDTH;
	return (mol);
}



