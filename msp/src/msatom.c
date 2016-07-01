/* 
   Molecular Surface Package
   written by Michael L. Connolly
   January 5, 2002
   
*/


#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"


/* ATOMS and BONDS */


/* service routines: */

unsigned long atom_size ()
{
	return (sizeof (struct atom));
}

unsigned long bond_size ()
{
	return (sizeof (struct bond));
}

void get_bond_atoms (int bond_number, int atoms[2])
{
	int j;
	struct bond *bnd_ptr;

	bnd_ptr = bond_ptr (bond_number);
	if (bnd_ptr == (struct bond *) NULL) return;
	for (j = 0; j < 2; j++)
		atoms[j] = bnd_ptr -> atoms[j];
}

void set_bond_atoms (int bond_number, int atoms[2])
{
	int j;
	struct bond *bnd_ptr;

	bnd_ptr = bond_ptr (bond_number);
	if (bnd_ptr == (struct bond *) NULL) return;
	for (j = 0; j < 2; j++)
		bnd_ptr -> atoms[j] = atoms[j];
}

void set_atom (int atom_number, int type, int chemical_element, double center[3], double radius, char group[MAX_ATNAME], char sequence[MAX_ATNAME], char name[MAX_ATNAME])
{
	int k;
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return;
	atm_ptr -> type = (short) type;
	atm_ptr -> chemical_element = (short) chemical_element;
	for (k = 0; k < 3; k++)
		atm_ptr -> center[k] = center[k];
	atm_ptr -> radius = radius;
	strcpy (atm_ptr -> group, group);
	strcpy (atm_ptr -> sequence, sequence);
	strcpy (atm_ptr -> name, name);
}


struct atom *atom_ptr (int number)
{
	return ((struct atom *) generic_ptr (ATOM, number));
}

struct bond *bond_ptr (int number)
{
	return ((struct bond *) generic_ptr (BOND, number));
}

void init_mol ()
{
	define_real ("atom_ratio", DEFAULT_ATOM_RATIO);
	if (error()) return;
	define_real ("bond_ratio", DEFAULT_BOND_RATIO);
	if (error()) return;
	memory_init (ATOM);
	if (error ()) {
			return;
	}
	memory_init (BOND);
	if (error ()) {
			return;
	}
}

/* ATOM ROUTINES */

int put_atom (double center[3], double radius, double ball_radius, double covalent,
char *atom_name, char *residue, char *sequence, char *kind, char *pdb, char *subunit,
int chemical_element, int vdw_type, int anumber, int srn, int mol,
double occupancy, double tfactor, double density)
{
	int atom_number;
	
	/* PSL calls */
	atom_number = allocate_psl (ATOM);
	if (atom_number == 0) return (0);
	set_atom_center (atom_number, center);
	set_atom_pdb (atom_number, pdb);
	set_atom_kind (atom_number, kind);
	set_atom_name (atom_number, atom_name);
	set_atom_group (atom_number, residue);
	set_atom_sequence (atom_number, sequence);
	set_atom_subunit (atom_number, subunit);
	set_atom_anumber (atom_number, anumber);
	set_atom_occupancy (atom_number, occupancy);
	set_atom_tfactor (atom_number, tfactor);
	set_atom_element (atom_number, chemical_element);
	set_atom_type (atom_number, vdw_type);
	set_atom_srn (atom_number, srn);
	set_atom_molecule (atom_number, mol);
	set_atom_color (atom_number, DEFAULT_HUE);
	set_atom_radius (atom_number, radius);
	set_atom_covalent (atom_number, covalent);
	set_atom_ball (atom_number, ball_radius);
	set_atom_angle (atom_number, DEFAULT_ANGLE);
	set_atom_opacity (atom_number, DEFAULT_OPACITY);
	set_atom_density (atom_number, density);
	return (atom_number);
}

void get_atom (int atom_number, double center[3], double *radius,
double *covalent, char atom_name[MAX_ATNAME],
char residue[MAX_ATNAME], char sequence[MAX_ATNAME],
int *chemical_element, int *vdw_type, int *srn, int *mol, int *col,
double *opacity, double *ball, double *angle, double *density)
{
	/* PSL calls */
	get_atom_center (atom_number, center);
	*radius = get_atom_radius (atom_number);
	get_atom_name (atom_number, atom_name);
	get_atom_group (atom_number, residue);
	get_atom_sequence (atom_number, sequence);
	*chemical_element = get_atom_element (atom_number);
	*vdw_type = get_atom_type (atom_number);
	*srn = get_atom_srn (atom_number);
	*mol = get_atom_molecule (atom_number);
	*col = get_atom_color (atom_number);
	*opacity = get_atom_opacity (atom_number);
	*covalent = get_atom_covalent (atom_number);
	*ball = get_atom_ball (atom_number);
	*angle = get_atom_angle (atom_number);
	*density = get_atom_density (atom_number);
}


/* BOND ROUTINES */

/* atom1 < atom2 */
int lookup_bond (int atom1, int atom2)
{
	int b, bond_number;
	struct atom *aptr;
	
	if (atom1 >= atom2) return (0);
	if (atom1 <= 0 || atom2 <= 0) return (0);
	aptr = atom_ptr (atom1);
	if (aptr == NULL) return (0);
	for (b = 0; b < MAX_VALENCE; b++)
		if (aptr -> bonded_to[b] == atom2) {
			bond_number = aptr -> bond_number[b];
			return (bond_number);
		}
	return (0);
}

void add_bond (int atom1, int atom2)
{
	int b, bond_number, mol1, mol2;
	int as[2];
	char message[MAX_STRING];
	struct atom *aptr;

	if (atom1 <= 0 || atom2 <= 0 || atom1 >= atom2) {
		sprintf (message,"invalid atom numbers: %d %d", atom1, atom2);
		set_error1 (message);
		return;
	}
	/* check for same molecule */
	mol1 = get_atom_molecule (atom1);
	if (error()) return;
	mol2 = get_atom_molecule (atom2);
	if (error()) return;
	if (mol1 != mol2) return;
	bond_number = lookup_bond (atom1, atom2);
	if (bond_number > 0) {
		/* already in local array, make sure in PSL set */
		include_member (bonds, BOND, bond_number);
		if (error()) {
			return;
		}
		return;
	}
	aptr = atom_ptr (atom1);
	if (aptr == NULL) return;
	if (error()) return;
	for (b = 0; b < MAX_VALENCE; b++)
		if (aptr -> bonded_to[b] == 0) break;
	if (b >= MAX_VALENCE) return;
	bond_number = allocate_psl(BOND);
	if (bond_number == 0 || error()) {
		sprintf (message, "add_bond: too many bonds");
		set_error1 (message);
		return;
	}
	as[0] = atom1;
	as[1] = atom2;
	set_bond_atoms (bond_number, as);
	aptr -> bond_number[b] = bond_number;
	aptr -> bonded_to[b] = atom2;
	include_member (bonds, BOND, bond_number);
}

void delete_bond (int atom1, int atom2)
{
	int bond_number;
	
	bond_number = lookup_bond (atom1, atom2);
	if (bond_number <= 0) return;
	exclude_member (bonds, BOND, bond_number);
}


void connect_command (char *setname, double tolerance)
{
	int atom_set, nb1, nb2;
	int i1, i2, na, nb, a1, a2, k;
	double d, s, radius1, radius2, max_covalent;
	double center1[3], center2[3];
	struct array *numbers;
	struct array *spheres;
	double   rect[3][2];
	atomnum *jnbr, *inrect;
	char message[MAX_STRING];


	strcpy (function_name, "connect_command");
	atom_set = lookup_set (setname);
	if (atom_set == 0) {
		sprintf (message, "connect_command: invalid set name: %s", setname);
		set_error1(message);
		return;
	}
	nb1 = count_set (bonds);

	na = count_set (atom_set);
	if (na <= 1) return;

	numbers = new_array (INTEGER, na);
	if (numbers == NULL) {
		set_error2("connect_command: numbers array");
		return;
	}
	spheres = new_array (SPHERE, na);
	if (spheres == NULL) {
		set_error2("connect_command: centers array");
		return;
	}
	natom = na;
	atmco = spheres -> centers;
	atmrad = spheres -> radii;
	atmatt = NULL;
	max_covalent = 0.0;
	for (a1 = init_for (atom_set), i1 = 0; a1 != 0;
		a1 = next_for (atom_set), i1++) {
		store_integer (numbers, i1, a1);
		if (error ()) return;
		get_atom_center (a1, center1);
		if (error ()) return;
		radius1 = get_atom_covalent (a1);
		if (error ()) return;
		if (radius1 > max_covalent) max_covalent = radius1;
		store_sphere (spheres, i1, center1, radius1);
		if (error ()) return;
	}

	dsbtree ();     /* create box tree */
    if (errflg) {
       	set_error2 ("connect_command: dsbtree returns error");
        return;
    }

	strcpy (function_name, "connect_command");

	for (i1 = 0; i1 < na - 1; i1++) {
		radius1 = fetch_sphere (spheres, i1, center1);
		if (error ()) return;
        /* set up rectangular range for searching box tree */
        /* make it big enough to include any possible nbr of atom ia */
        for (k = 0; k < 3; k++) {
            rect[k][0] = center1[k] - radius1 - max_covalent - tolerance;
            rect[k][1] = center1[k] + radius1 + max_covalent + tolerance;
        }
        /* get atoms in rectangle */
        inrect = dsgair (rect, rootbox);
        if (errflg) {
        	set_error2 ("dsgair returns error");
            return;
        }
        if (inrect == NULL) {
		continue;
	}
        /* loop through atoms in rectangular range */
        for (jnbr = inrect; *jnbr != EOL; jnbr++) {
            i2 = *jnbr;         /* retrieve atom index from list */
			if (i2 <= i1) continue;
			radius2 = fetch_sphere (spheres, i2, center2);
			if (error ()) return;
			d = distance (center1, center2);
			s = radius1 + radius2;
			if (d > s + tolerance) continue;
			a1 = fetch_integer (numbers, i1);
			if (error ()) return;
			a2 = fetch_integer (numbers, i2);
			if (error ()) return;
			add_bond (a1, a2);
			if (error ()) return;
		}
		/* free list from box tree search */
		free_atomnums (inrect);
	}
	nb2 = count_set (bonds);
	nb = nb2 - nb1;
	sprintf (message, "%8d bonds made for %8d atoms", nb, na);
	inform (message);
	free_array (numbers);
	free_array (spheres);
	dsfbox (rootbox);   /* free box (recursive) */
	free_cache(BOX);
	atmco = NULL; atmrad = NULL;
}

void disconnect_command (char *setname1, char *setname2)
{
	int set1, set2, na1, na2, nb1, nb2, nb;
	int i1, i2, atom1, atom2;
	struct array *numbers1, *numbers2;
	char message[MAX_STRING];


	strcpy (function_name, "disconnect_command");
	set1 = lookup_set (setname1);
	if (set1 == 0) return;
	na1 = count_set (set1);
	set2 = lookup_set (setname2);
	if (set2 == 0) return;
	na2 = count_set (set2);

	if (na1 <= 0 || na2 <= 0) return;
	nb1 = count_set (bonds);

	numbers1 = new_array (INTEGER, na1);
	if (numbers1 == NULL) {
		set_error2("disconnect_command: numbers1 array");
		return;
	}
	numbers2 = new_array (INTEGER, na2);
	if (numbers2 == NULL) {
		set_error1("disconnect_command: numbers2 array");
		return;
	}
	for (atom1 = init_for (set1), i1 = 0; atom1 != 0;
		atom1 = next_for (set1), i1++) {
		store_integer (numbers1, i1, atom1);
		if (error ()) return;
	}
	for (atom2 = init_for (set2), i2 = 0; atom2 != 0;
		atom2 = next_for (set2), i2++) {
		store_integer (numbers2, i2, atom2);
		if (error ()) return;
	}
	
	for (i1 = 0; i1 < na1; i1++) {
		atom1 = fetch_integer (numbers1, i1);
		if (error ()) return;
		for (i2 = 0; i2 < na2; i2++) {
			atom2 = fetch_integer (numbers2, i2);
			if (error ()) return;
			delete_bond (atom1, atom2);
			if (error ()) return;
		}
	}
	nb2 = count_set (bonds);
	nb = nb1 - nb2;
	sprintf (message,"%8d bonds removed", nb);
	inform (message);
	free_array (numbers1);
	free_array (numbers2);
}

/* ATOM SET ROUTINES */

/* defined a named set (atoms or bonds) */

int named_set (char *setname, int type)
{
	int symbol_number;
	int object_number;
	
	symbol_number = new_object_symbol (setname, SET);
	if (error()) return (0);
	object_number = allocate_named (SET, symbol_number);
	if (object_number == 0) return (0);
	set_set_type (object_number, type);
	return (object_number);
}

/* return object number of set */

int lookup_set (char *setname)
{
	int object_number;
	
	object_number = object_number_of (setname);
	if (object_number == 0) return (0);
	return (object_number);
}

/* clear an atom set to the null set */

void clear_command (char *setname)
{
	int set_number;
	char message[MAX_STRING];
	
	set_number = lookup_set (setname);
	if (set_number == 0) {
		sprintf (message, "(clear_command): invalid set");
		set_error1 (message);
		return;
	}
	clear_set (set_number);
}

int protected_set(char *setname)
{
	if (strcmp(setname,"atoms") == 0) return (1);
	else if (strcmp(setname,"these") == 0) return (1);
	else if (strcmp(setname,"unknown") == 0) return (1);
	else if (strcmp(setname,"bonds") == 0) return (1);
	else if (strcmp(setname,"scratch") == 0) return (1);
	else return (0);
}


/* ASSIGNMENT COMMANDS */

void assignment (char *setname, char *field, char *operator, char *value)
{
	int set_number;
	
	set_number = lookup_set (setname);
	if (set_number == 0) {
		set_number = named_set (setname, ATOM);
		if (set_number == 0) {
			return;
		}
		
	} else clear_set (set_number);
	chooser (atoms, set_number, field, operator, value);
}


void times_equals (char *setname, char *field, char *operator, char *value)
{
	int set_number;
	char message[MAX_STRING];
	
	set_number = lookup_set (setname);
	if (set_number == 0) {
		sprintf (message, "invalid setname(%s) in *=",setname);
		set_error1 (message);
		return;
	}
	clear_set (scratch);
	chooser (set_number, scratch, field, operator, value);
	clear_set (set_number);
	set_union (set_number, scratch, set_number);
	clear_set (scratch);
}

void plus_equals (char *setname, char *field, char *operator, char *value)
{
	int set_number;
	char message[MAX_STRING];
	
	set_number = lookup_set (setname);
	if (set_number == 0) {
		sprintf (message, "invalid setname(%s) in +",setname);
		set_error1 (message);
		return;
	}
	chooser (atoms, set_number, field, operator, value);
}

void minus_equals (char *setname, char *field, char *operator, char *value)
{
	int set_number;
	char message[MAX_STRING];
	
	set_number = lookup_set (setname);
	if (set_number == 0) {
		sprintf (message, "invalid setname(%s) in -=",setname);
		set_error1 (message);
		return;
	}
	chooser (set_number, scratch, field, operator, value);
	set_subtraction (set_number, scratch, set_number);
	clear_set (scratch);
}

/* SET OPERATIONS */

void handle_intersection (char *setname1, char *setname2, char *setname3)
{
	int set1, set2, set3;
	char message[MAX_STRING];
	
	set1 = lookup_set (setname1);
	if (set1 == 0) {
		set1 = named_set (setname1, ATOM);
		if (set1 == 0) {
			return;
		}
	}
	set2 = lookup_set (setname2);
	if (set2 == 0) {
		sprintf (message, "invalid setname(%s) in *",setname2);
		set_error1 (message);
		return;
	}
	set3 = lookup_set (setname3);
	if (set3 == 0) {
		sprintf (message, "invalid setname(%s) in *",setname3);
		set_error1 (message);
		return;
	}
	set_intersection (set2, set3, set1);
}

void handle_union (char *setname1, char *setname2, char *setname3)
{
	int set1, set2, set3;
	char message[MAX_STRING];

	set1 = lookup_set (setname1);
	if (set1 == 0) {
		set1 = named_set (setname1, ATOM);
		if (set1 == 0) {
			return;
		}
	}
	set2 = lookup_set (setname2);
	if (set2 == 0) {
		sprintf (message, "invalid setname(%s) in +",setname2);
		set_error1 (message);
		return;
	}
	set3 = lookup_set (setname3);
	if (set3 == 0) {
		sprintf (message, "invalid setname(%s) in +",setname3);
		set_error1 (message);
		return;
	}
	set_union (set2, set3, set1);
}

void handle_subtraction (char *setname1, char *setname2, char *setname3)
{
	int set1, set2, set3;
	char message[MAX_STRING];

	if (strcmp (setname2, setname3) == 0) return;

	set1 = lookup_set (setname1);
	if (set1 == 0) {
		set1 = named_set (setname1, ATOM);
		if (set1 == 0 || error ()) {
			sprintf (message, "problem with set definition");
			set_error1 (message);
			return;
		}
	}
	if (strcmp (setname2, setname3) == 0) {
		clear_set (set1);
		return;
	}
	set2 = lookup_set (setname2);
	if (set2 == 0 || error ()) {
		sprintf (message, "invalid setname(%s) in -",setname2);
		set_error2 (message);
		return;
	}
	set3 = lookup_set (setname3);
	if (set3 == 0 || error ()) {
		sprintf (message, "invalid setname(%s) in -",setname3);
		set_error2 (message);
		return;
	}
	if (error ()) {
		return;
	}
	set_subtraction (set2, set3, set1);
	if (error()) {
		sprintf (message, "problem with set subtraction");
		set_error2 (message);
		return;
	}
}

void simple_set_copy (char *setname1, char *setname2)
{
	int set1, set2, o;
	int count1, count2;
	char message[MAX_STRING];

	if (strcmp (setname1, setname2) == 0) return;
	set1 = lookup_set (setname1);
	if (set1 == 0) {
		set1 = named_set (setname1, ATOM);
		if (set1 == 0 || error ()) {
			sprintf (message, "problem with set definition");
			set_error2 (message);
			return;
		}
	}
	else clear_set (set1);
	
	set2 = lookup_set (setname2);
	if (set2 == 0 || error ()) {
		set_error2 (message);
		sprintf (message, "invalid setname(%s) in -",setname2);
		return;
	}
	count2 = count_set (set2);
	
	/* copy atom numbers from set2 to set1 */
	for (o = init_for (set2); o != 0; o = next_for (set2)) {
		include_member (set1, ATOM, o);
	}
	
	count1 = count_set (set1);
	if (count1 != count2) {
		sprintf (message, "simple_set_copy: inconsistent number of atoms in sets");
		set_error1 (message);
	}
}


/* choose atoms from the from_set and include in the to_set
   based upon satisfying a relation between an atom field
   and a value specified in the command line */
   
void chooser (int from_set, int to_set, char *field, char *operator, char *value)
{
	if (strcmp (field, "atom") == 0)
		atom_field (from_set, to_set, operator, value);
	else if (strcmp (field, "residue") == 0)
		residue_field (from_set, to_set, operator, value);
	else if (strcmp (field, "sequence") == 0)
		sequence_field (from_set, to_set, operator, value);
	else if (strcmp (field, "suffix") == 0)
		suffix_field (from_set, to_set, operator, value);
	else if (strcmp (field, "subunit") == 0)
		subunit_field (from_set, to_set, operator, value);
	else if (strcmp (field, "pdb") == 0)
		pdb_field (from_set, to_set, operator, value);
	else if (strcmp (field, "kind") == 0)
		kind_field (from_set, to_set, operator, value);
	else if (strcmp (field, "center") == 0)
		center_field (from_set, to_set, operator, value);
	else if (strcmp (field, "occupancy") == 0)
		occupancy_field (from_set, to_set, operator, value);
	else if (strcmp (field, "tfactor") == 0)
		tfactor_field (from_set, to_set, operator, value);
	else if (strcmp (field, "element") == 0)
		element_field (from_set, to_set, operator, value);
	else if (strcmp (field, "anumber") == 0)
		anumber_field (from_set, to_set, operator, value);
	else if (strcmp (field, "rnumber") == 0)
		rnumber_field (from_set, to_set, operator, value);
	else if (strcmp (field, "molecule") == 0)
		molecule_field (from_set, to_set, operator, value);
	else if (strcmp (field, "type") == 0)
		type_field (from_set, to_set, operator, value);
}

void atom_field (int from_set, int to_set, char *operator, char *value)
{
	int o, t;
	char atom_name[MAX_ATNAME];
	
	for (o = init_for (from_set); o != 0; o = next_for (from_set)) {
		get_atom_name (o, atom_name);
		t = string_relation (operator, atom_name, value);
		if (t) include_member (to_set, ATOM, o);
	}
}

void residue_field (int from_set, int to_set, char *operator, char *value)
{
	int o, t;
	char residue[MAX_ATNAME];
	
	for (o = init_for (from_set); o != 0; o = next_for (from_set)) {
		get_atom_group (o, residue);
		t = string_relation (operator, residue, value);
		if (t) include_member (to_set, ATOM, o);
	}
}

void sequence_field (int from_set, int to_set, char *operator, char *value)
{
	int o, t;
	char sequence[MAX_ATNAME];
	
	for (o = init_for (from_set); o != 0; o = next_for (from_set)) {
		get_atom_sequence (o, sequence);
		t = string_relation (operator, sequence, value);
		if (t) include_member (to_set, ATOM, o);
	}
}

void subunit_field (int from_set, int to_set, char *operator, char *value)
{
	int o, t;
	char subunit[MAX_ATNAME];
	
	for (o = init_for (from_set); o != 0; o = next_for (from_set)) {
		get_atom_subunit (o, subunit);
		t = string_relation (operator, subunit, value);
		if (t) include_member (to_set, ATOM, o);
	}
}

void suffix_field (int from_set, int to_set, char *operator, char *value)
{
	int o, t, rnumber;
	char suffix[MAX_ATNAME], sequence[MAX_ATNAME];
	
	for (o = init_for (from_set); o != 0; o = next_for (from_set)) {
		get_atom_sequence (o, sequence);
		sscanf (sequence, "%d%s", &rnumber, suffix);
		t = string_relation (operator, suffix, value);
		if (t) include_member (to_set, ATOM, o);
	}
}

void pdb_field (int from_set, int to_set, char *operator, char *value)
{
	int o, t;
	char pdb[MAX_ATNAME];
	
	for (o = init_for (from_set); o != 0; o = next_for (from_set)) {
		get_atom_pdb (o, pdb);
		t = string_relation (operator, pdb, value);
		if (t) include_member (to_set, ATOM, o);
	}
}

void kind_field (int from_set, int to_set, char *operator, char *value)
{
	int o, t;
	char kind[MAX_ATNAME];
	
	for (o = init_for (from_set); o != 0; o = next_for (from_set)) {
		get_atom_kind (o, kind);
		t = string_relation (operator, kind, value);
		if (t) include_member (to_set, ATOM, o);
	}
}

int string_relation (char *operator, char *s1, char *s2)
{
	int t, l1, l2, i;
	char s3[MAX_STRING], s4[MAX_STRING];
	
	if (strcmp (operator, "=") == 0) {
		t = (strcmp (s1, s2) == 0);
	}
	else if (strcmp (operator, "==") == 0) {
		t = (strcmp (s1, s2) == 0);
	}
	else if (strcmp (operator, "!=") == 0) {
		t = (strcmp (s1, s2) != 0);
	}
	else if (strcmp (operator, "matches") == 0) {
		l1 = strlen (s1);
		l2 = strlen (s2);
		/* always ignore case */
		strcpy (s3, s1); strcpy (s4, s2);
		for (i = 0; i < l1; i++)
			s3[i] = toupper (s1[i]);
		for (i = 0; i < l2; i++)
			s4[i] = toupper (s2[i]);
		t = (strncmp (s3, s4, l2) == 0);
	}
	else t = 0;
	return (t);
}

void center_field (int from_set, int to_set, char *operator, char *value)
{
	int region_number, o, t;
	double vrad;
	double center[3], vcen[3], vaxis[3];
	
	region_number = lookup_region (value);
	get_region_center (region_number, vcen);
	vrad = get_region_radius (region_number);
	get_region_direction (region_number, vaxis);
	for (o = init_for (from_set); o != 0; o = next_for (from_set)) {
		get_atom_center (o, center);
		t = space_relation (operator, center, vcen, vrad, vaxis);
		if (t) include_member (to_set, ATOM, o);
	}
}

int space_relation (char *operator, double center[3], double vcen[3], double vrad, double vaxis[3])
{
	int k, t;
	double r, d;
	double vector[3];
	
	for (k = 0; k < 3; k++)
		vector[k] = center[k] - vcen[k];
	r = norm (vector);
	d = dot_product (vector, vaxis);
	if (strcmp (operator, "inside") == 0)
		t = (r <= vrad);
	else if (strcmp (operator, "outside") == 0)
		t = (r > vrad);
	else if (strcmp (operator, "below") == 0)
		t = (d <= 0.0);
	else if (strcmp (operator, "above") == 0)
		t = (d > 0.0);
	else t = 0;
	return (t);
}

void occupancy_field (int from_set, int to_set, char *operator, char *value)
{
	int o, t;
	double fvalue, occupancy;
	
	fvalue = atof (value);
	for (o = init_for (from_set); o != 0; o = next_for (from_set)) {
		occupancy = get_atom_occupancy (o);
		t = real_relation (operator, occupancy, fvalue);
		if (t) include_member (to_set, ATOM, o);
	}
}

void tfactor_field (int from_set, int to_set, char *operator, char *value)
{
	int o, t;
	double fvalue, tfactor;
	
	fvalue = atof (value);
	for (o = init_for (from_set); o != 0; o = next_for (from_set)) {
		tfactor = get_atom_tfactor (o);
		t = real_relation (operator, tfactor, fvalue);
		if (t) include_member (to_set, ATOM, o);
	}
}

int real_relation (char *operator, double avalue, double fvalue)
{
	int t;
	
	if (strcmp (operator, "=") == 0)
		t = (avalue == fvalue);
	else if (strcmp (operator, "==") == 0)
		t = (avalue == fvalue);
	else if (strcmp (operator, "!=") == 0)
		t = (avalue != fvalue);
	else if (strcmp (operator, "<=") == 0)
		t = (avalue <= fvalue);
	else if (strcmp (operator, ">=") == 0)
		t = (avalue >= fvalue);
	else if (strcmp (operator, "<") == 0)
		t = (avalue < fvalue);
	else if (strcmp (operator, ">") == 0)
		t = (avalue > fvalue);
	else t = 0;
	return (t);
}

void element_field (int from_set, int to_set, char *operator, char *value)
{
	int o, t;
	int ivalue, chemical_element;
	
	ivalue = atof (value);
	for (o = init_for (from_set); o != 0; o = next_for (from_set)) {
		chemical_element = get_atom_element (o);
		t = integer_relation (operator, chemical_element, ivalue);
		if (t) include_member (to_set, ATOM, o);
	}
}

void anumber_field (int from_set, int to_set, char *operator, char *value)
{
	int o, t;
	int ivalue, anumber;
	
	ivalue = atof (value);
	for (o = init_for (from_set); o != 0; o = next_for (from_set)) {
		anumber = get_atom_anumber (o);
		t = integer_relation (operator, anumber, ivalue);
		if (t) include_member (to_set, ATOM, o);
	}
}

void rnumber_field (int from_set, int to_set, char *operator, char *value)
{
	int  o, t;
	int ivalue, rnumber;
	char sequence[MAX_ATNAME];
	
	ivalue = atof (value);
	for (o = init_for (from_set); o != 0; o = next_for (from_set)) {
		get_atom_sequence (o, sequence);
		rnumber = atoi (sequence);
		t = integer_relation (operator, rnumber, ivalue);
		if (t) include_member (to_set, ATOM, o);
	}
}

void molecule_field (int from_set, int to_set, char *operator, char *value)
{
	int  o, t;
	int ivalue, mol;
	
	ivalue = atof (value);
	for (o = init_for (from_set); o != 0; o = next_for (from_set)) {
		mol = get_atom_molecule (o);
		t = integer_relation (operator, mol, ivalue);
		if (t) include_member (to_set, ATOM, o);
	}
}

void type_field (int from_set, int to_set, char *operator, char *value)
{
	int  o, t;
	int ivalue, vdw_type;
	
	ivalue = atof (value);
	for (o = init_for (from_set); o != 0; o = next_for (from_set)) {
		vdw_type = get_atom_type (o);
		t = integer_relation (operator, vdw_type, ivalue);
		if (t) include_member (to_set, ATOM, o);
	}
}

int integer_relation (char *operator, int avalue, int ivalue)
{
	int t;
	
	if (strcmp (operator, "=") == 0)
		t = (avalue == ivalue);
	else if (strcmp (operator, "==") == 0)
		t = (avalue == ivalue);
	else if (strcmp (operator, "!=") == 0)
		t = (avalue != ivalue);
	else if (strcmp (operator, "<=") == 0)
		t = (avalue <= ivalue);
	else if (strcmp (operator, ">=") == 0)
		t = (avalue >= ivalue);
	else if (strcmp (operator, "<") == 0)
		t = (avalue < ivalue);
	else if (strcmp (operator, ">") == 0)
		t = (avalue > ivalue);
	else t = 0;
	return (t);
}

/* FIELD ASSIGNMENT */

void field_assignment (char *setname, char *field, char *value)
{
	int set_number, o, n, vdw_type, ubefore, uafter;
	char message[MAX_STRING];
	
	set_number = lookup_set (setname);
	if (set_number == 0) return;
	n = 0;
	for (o = init_for (set_number); o != 0; o = next_for (set_number)) {
		one_field (o, field, value); n++;
	}
	sprintf (message, "%8d assignments of %s to field %s of %s", n, value, field, setname);
	if (n > 0) inform (message);
	/* update unknown vdw type set */
	if (strcmp(field, "type") == 0) {
		ubefore = count_set(unknowns);
		for (o = init_for (atoms); o != 0; o = next_for (atoms)) {
			if (!member_of (o, set_number)) continue;
			vdw_type = get_atom_type (o);
			if (vdw_type == 0) include_member (unknowns, ATOM, o);
			else exclude_member (unknowns, ATOM, o);
		}
		uafter = count_set(unknowns);
		if (uafter > ubefore) sprintf (message, "%d atoms added to unknowns",
			uafter - ubefore);
		else if (uafter < ubefore) printf (message, "%d atoms removed from unknowns",
			ubefore - uafter);
		if (uafter != ubefore) inform (message);
	}
}

void one_field (int atom_number, char *field, char *value)
{
	int type, chemical_element, srn, mol, col;
	double radius, covalent, ball, density, angle, opacity;
	char message[MAX_STRING];
	
	if (strcmp (field, "radius") == 0) {
		radius = atof (value);
		set_atom_radius (atom_number, radius);
	}
	else if (strcmp (field, "covalent") == 0) {
		covalent = atof (value);
		set_atom_covalent (atom_number, covalent);
	}
	else if (strcmp (field, "ball") == 0) {
		ball = atof (value);
		set_atom_ball (atom_number, ball);
	}
	else if (strcmp (field, "density") == 0) {
		density = atof (value);
		set_atom_density (atom_number, density);
	}
	else if (strcmp (field, "opacity") == 0) {
		opacity = atof (value);
		set_atom_opacity (atom_number, opacity);
	}
	else if (strcmp (field, "angle") == 0) {
		angle = atof (value);
		set_atom_angle (atom_number, angle);
	}
	else if (strcmp (field, "element") == 0) {
		chemical_element = atoi (value);
		set_atom_element (atom_number, chemical_element);
	}
	else if (strcmp (field, "type") == 0) {
		type = atoi (value);
		set_atom_type (atom_number, type);
	}
	else if (strcmp (field, "srn") == 0) {
		srn = atoi (value);
		set_atom_srn (atom_number, srn);
	}
	else if (strcmp (field, "molecule") == 0) {
		mol = atoi (value);
		set_atom_molecule (atom_number, mol);
	}
	else if (strcmp (field, "color") == 0) {
		if (isdigit(value[0]))
			col = atoi (value);
		else col = get_color_number(table, value);
		set_atom_color (atom_number, col);
	}
	else if (strcmp (field, "kind") == 0) {
		set_atom_kind (atom_number, value);
	}
	else {
		sprintf (message, "unable to assign to field %s of atom %d",
			field, atom_number);
		inform (message);
	}
}

double get_atom_radius (int atom_number)
{
	double radius;
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return (0.0);
	radius = atm_ptr -> radius;
	return (radius);
}

void set_atom_radius (int atom_number, double radius)
{
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return;
	atm_ptr -> radius = radius;
}

double get_atom_covalent (int atom_number)
{
	double covalent;
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return (0.0);
	covalent = atm_ptr -> covalent;
	return (covalent);
}

void set_atom_covalent (int atom_number, double covalent)
{
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return;
	atm_ptr -> covalent = covalent;
}

double get_atom_ball (int atom_number)
{
	double ball;
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return (0.0);
	ball = atm_ptr -> ball;
	return (ball);
}

void set_atom_ball (int atom_number, double ball)
{
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return;
	atm_ptr -> ball = ball;
}

double get_atom_angle (int atom_number)
{
	double angle;
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return (0.0);
	angle = atm_ptr -> angle;
	return (angle);
}

void set_atom_angle (int atom_number, double angle)
{
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return;
	atm_ptr -> angle = angle;
}

double get_atom_density (int atom_number)
{
	double density;
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return (0.0);
	density = atm_ptr -> density;
	return (density);
}

void set_atom_density (int atom_number, double density)
{
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return;
	atm_ptr -> density = density;
}

double get_atom_occupancy (int atom_number)
{
	double occupancy;
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return (0.0);
	occupancy = atm_ptr -> occupancy;
	return (occupancy);
}

void set_atom_occupancy (int atom_number, double occupancy)
{
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return;
	atm_ptr -> occupancy = occupancy;
}

double get_atom_opacity (int atom_number)
{
	double opacity;
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return (0.0);
	opacity = atm_ptr -> opacity;
	return (opacity);
}

void set_atom_opacity (int atom_number, double opacity)
{
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return;
	atm_ptr -> opacity = opacity;
}


double get_atom_tfactor (int atom_number)
{
	double tfactor;
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return (0.0);
	tfactor = atm_ptr -> tfactor;
	return (tfactor);
}

void set_atom_tfactor (int atom_number, double tfactor)
{
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return;
	atm_ptr -> tfactor = tfactor;
}

void get_atom_center (int atom_number, double center[3])
{
	int k;
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return;
	for (k = 0; k < 3; k++)
		center[k] = atm_ptr -> center[k];
}

void set_atom_center (int atom_number, double center[3])
{
	int k;
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return;
	for (k = 0; k < 3; k++)
		atm_ptr -> center[k] = center[k];
}

void set_atom_type (int atom_number, int type)
{
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return;
	atm_ptr -> type = (short) type;
}

int get_atom_type (int atom_number)
{
	int type;
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return (0);
	type = atm_ptr -> type;
	return (type);
}

void set_atom_element (int atom_number, int chemical_element)
{
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return;
	atm_ptr -> chemical_element = (short) chemical_element;
}

int get_atom_element (int atom_number)
{
	int chemical_element;
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return (0);
	chemical_element = atm_ptr -> chemical_element;
	return (chemical_element);
}

void set_atom_anumber (int atom_number, int anumber)
{
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return;
	atm_ptr -> anumber = (short) anumber;
}

int get_atom_anumber (int atom_number)
{
	int anumber;
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return (0);
	anumber = atm_ptr -> anumber;
	return (anumber);
}

void set_atom_srn (int atom_number, int srn)
{
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return;
	atm_ptr -> srn = (short) srn;
}

int get_atom_srn (int atom_number)
{
	int srn;
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return (0);
	srn = atm_ptr -> srn;
	return (srn);
}

void set_atom_molecule (int atom_number, int mol)
{
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return;
	atm_ptr -> mol = (short) mol;
}

int get_atom_molecule (int atom_number)
{
	int mol;
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return (0);
	mol = atm_ptr -> mol;
	return (mol);
}

void set_atom_color (int atom_number, int col)
{
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return;
	atm_ptr -> col = (short) col;
}

int get_atom_color (int atom_number)
{
	int col;
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return (0);
	col = atm_ptr -> col;
	return (col);
}


void get_atom_subunit (int atom_number, char subunit[MAX_ATNAME])
{
	int k;
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return;
	for (k = 0; k < MAX_ATNAME; k++)
		subunit[k] = atm_ptr -> subunit[k];
}

void get_atom_name (int atom_number, char name[MAX_ATNAME])
{
	int k;
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return;
	for (k = 0; k < MAX_ATNAME; k++)
		name[k] = atm_ptr -> name[k];
}

void get_atom_group (int atom_number, char group[MAX_ATNAME])
{
	int k;
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return;
	for (k = 0; k < MAX_ATNAME; k++)
		group[k] = atm_ptr -> group[k];
}

void get_atom_sequence (int atom_number, char sequence[MAX_ATNAME])
{
	int k;
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return;
	for (k = 0; k < MAX_ATNAME; k++)
		sequence[k] = atm_ptr -> sequence[k];
}

void get_atom_pdb (int atom_number, char pdb[MAX_ATNAME])
{
	int k;
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return;
	for (k = 0; k < MAX_ATNAME; k++)
		pdb[k] = atm_ptr -> pdb[k];
}

void get_atom_kind (int atom_number, char kind[MAX_ATNAME])
{
	int k;
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return;
	for (k = 0; k < MAX_ATNAME; k++)
		kind[k] = atm_ptr -> kind[k];
}


void set_atom_subunit (int atom_number, char subunit[MAX_ATNAME])
{
	int k;
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return;
	for (k = 0; k < MAX_ATNAME; k++)
		atm_ptr -> subunit[k] = subunit[k];
}

void set_atom_name (int atom_number, char name[MAX_ATNAME])
{
	int k;
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return;
	for (k = 0; k < MAX_ATNAME; k++)
		atm_ptr -> name[k] = name[k];
}

void set_atom_group (int atom_number, char group[MAX_ATNAME])
{
	int k;
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return;
	for (k = 0; k < MAX_ATNAME; k++)
		atm_ptr -> group[k] = group[k];
}

void set_atom_sequence (int atom_number, char sequence[MAX_ATNAME])
{
	int k;
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return;
	for (k = 0; k < MAX_ATNAME; k++)
		atm_ptr -> sequence[k] = sequence[k];
}

void set_atom_pdb (int atom_number, char pdb[MAX_ATNAME])
{
	int k;
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return;
	for (k = 0; k < MAX_ATNAME; k++)
		atm_ptr -> pdb[k] = pdb[k];
}

void set_atom_kind (int atom_number, char kind[MAX_ATNAME])
{
	int k;
	struct atom *atm_ptr;

	atm_ptr = atom_ptr (atom_number);
	if (atm_ptr == (struct atom *) NULL) return;
	for (k = 0; k < MAX_ATNAME; k++)
		atm_ptr -> kind[k] = kind[k];
}

