/* 
   Molecular Surface Package
   Protein Data Bank to MS format
   PDB2MS
   written by Michael L. Connolly
   January 27, 2006
   
   
   using:
   
   Protein Shape Library
   Copyright 1993 by Michael L. Connolly
   All Rights Reserved
   
*/

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

/* initialization of program */

void init_pdb ()
{
	long nri;
	char message[MAX_STRING];

	n_atom = 0; n_file = 0; n_molecule = 0;
	atoms = named_set ("atoms", ATOM);
	bonds = allocate_set (BOND);
	tubers = named_set ("tubers", ATOM);
	unknowns = named_set ("unknowns", ATOM);
	these = named_set ("these", ATOM);
	scratch = allocate_set (ATOM);
	
	init_patterns();
	init_radii();
	nri = create_residue_index ();
	if (error ()) return;
	sprintf (message, "%8ld residue indices for pattern table", nri);
	informd (message);
}


/* Input/Output */

int read_atom_file(struct array *atom_lines, char *format, char *setname)
{
	int result, format_type;
	int srn, mol;
	int atom_number, n_these;
	int set_number, chemical_element, vdw_type;
	int idx, jdx, best_match, found;
	int anumber;
	int ii;
	long l, n_line;
	long n_vdw_zero;
	double occupancy, tfactor, density;
	double radius, ball_radius, covalent, center[3];
	double min_r, max_r;
	char *in_line;
	char prev_subunit[6];
	char atom[6], res[6], seq[6], subunit[6];
	char pdb[8], ans[6];
	char xs[10], ys[10], zs[10];
	char os[8], tfs[8];
	char srns[6], ts[6], mols[6], rs[6], dens[6];
	char kind[MAX_ATNAME];
	char message[MAX_STRING];
	char shortened[MAX_STRING];
	char preface[MAX_STRING];
	struct radiusRecord *ar;
	struct patternRecord *pat;
	struct residueIndex *rin;
	
	if (format == NULL) format_type = PDB_FORMAT;
	else if (strlen (format) <= (unsigned) 0) format_type = PDB_FORMAT;
	else if (strcmp (format, "pdb") == 0) format_type = PDB_FORMAT;
	else if (strcmp (format, "ms") == 0) format_type = MS_FORMAT;
	else if (strcmp (format, "xyzr") == 0) format_type = XYZR_FORMAT;
	else if (strcmp (format, "ds") == 0) format_type = DS_FORMAT;
	else {
		sprintf (message, "read_atom_file: invalid input file format: %s", format);
		set_error1 (message);
		return (0);
	}
	set_number = lookup_set (setname);
	if (set_number == 0) {
		sprintf (message, "read_atom_file: invalid molecule atom set name: %s", setname);
		set_error1 (message);
		return (0);
	}
	clear_set (these);
	strcpy (prev_subunit, "@");	/* unlikely to be used as subunit */
	n_file++;
	n_line = atom_lines -> length;
	n_vdw_zero = 0;
	min_r =  10000000.0;
	max_r = -10000000.0;
	for (l = 0; l < n_line; l++) {
		in_line = fetch_string (atom_lines, l);
		if (in_line == NULL) break;
		/* initialize */
		strcpy (pdb, "");
		strcpy (ans, "");
		strcpy (atom, "");
		strcpy (res, "");
		strcpy (subunit, "");
		strcpy (seq, "");
		strcpy (xs, "");
		strcpy (ys, "");
		strcpy (zs, "");
		strcpy (os, "");
		strcpy (tfs, "");
		strcpy (ts, "");
		strcpy (srns, "");
		srn = DEFAULT_SRN;
		density = DEFAULT_DENSITY;
		strcpy (mols, "");
		strcpy (rs, "");
		strcpy (dens, "");
		
		switch (format_type) {
		case PDB_FORMAT:
			result = interpret_pdb (in_line,pdb,ans,atom,res,subunit,
				seq,xs,ys,zs,os,tfs);
			break;
		case MS_FORMAT:
			result = interpret_ms (in_line,xs, ys, zs, ts, srns, mols,
				res, seq, atom);
			break;
		case XYZR_FORMAT:
			result = interpret_xyzr (in_line, xs, ys, zs, rs,
				res, seq, atom);
			break;
		case DS_FORMAT:
			result = interpret_ds (in_line, xs, ys, zs, rs, srns, mols, dens,
				res, seq, atom);
			break;
		default:
			break;
		}
		if (result == 0) continue;
		if (result < 0) return (0);
		
		/* same for all: */
		center[0] = atof (xs);
		center[1] = atof (ys);
		center[2] = atof (zs);
		
		/* determine chemical element from first letter of atom name */
		if (atom[0] == 'H') chemical_element = 1;
		else if (atom[0] == 'C') chemical_element = 6;
		else if (atom[0] == 'N') chemical_element = 7;
		else if (atom[0] == 'O') chemical_element = 8;
		else if (atom[0] == 'P') chemical_element = 15;
		else if (atom[0] == 'S') chemical_element = 16;
		else if (atom[0] == 'F') chemical_element = 26;
		else chemical_element = 0;	/* not recognized */
	
		vdw_type = 0;
		radius = DEFAULT_ATOM_RADIUS;	/* default atomic vdw radius */
		covalent = DEFAULT_COVALENT;
		ball_radius = DEFAULT_BALL_RADIUS;
		strcpy (kind, "");			/* default kind */
		/* determine type from residue name and atom name */
		best_match = -1;
		for (jdx = 0; jdx < n_residue_index; jdx++) {
			rin = residue_table + jdx;
			if (pattern_match (res, rin -> residue)) {
				for (idx = rin -> first_pattern; idx <= rin -> last_pattern; idx++) {
					pat = pattern_table + idx;
					if (pattern_match (res,pat->residue) &&
						pattern_match(atom,pat->atom_name))
							best_match = idx;
				}
			}
		}
		if (best_match >= 0) {
			pat = pattern_table + best_match;
			vdw_type = pat -> type;
			strcpy (kind, pat -> kind);
			/* lookup radius */
			found = 0;
			for (idx = 0; idx < n_radii; idx++) {
				ar = radius_table + idx;
				if (ar -> type == pat -> type) {
					radius = ar -> radius;
					covalent = ar -> covalent;
					found = 1;
					break;
				}
			}
			if (!found) {
				sprintf(message,"radius not found for type %ld %s:",
					pat -> type, pat -> kind);
				inform (message);
				inform(in_line);
				sprintf(message,
					"using default radius: %8.3f", DEFAULT_ATOM_RADIUS);
				inform (message);
				radius = DEFAULT_ATOM_RADIUS;
				covalent = DEFAULT_COVALENT;
			}
		}
		if (format_type == PDB_FORMAT) {
			/* a change in the subunit character is interpreted as a new molecule */
			if (strcmp (subunit, prev_subunit) != 0) {
				strcpy (prev_subunit, subunit);
				n_molecule++;
			}
		}
		else {
			if (n_atom == 0) n_molecule = 1;
		}
	
		switch (format_type) {
		case PDB_FORMAT:
			anumber = atoi (ans);
			occupancy = atof (os);
			tfactor = atof (tfs);
			mol = n_molecule;
			break;
		case MS_FORMAT:
			anumber = n_atom + 1;
			vdw_type = atoi (ts);
			srn = atoi (srns);
			mol = atoi (mols);
			break;
		case XYZR_FORMAT:
			anumber = n_atom + 1;
			radius = atof (rs);
			mol = n_molecule;
			break;
		case DS_FORMAT:
			anumber = n_atom + 1;
			radius = atof (rs);
			srn = atoi (srns);
			mol = atoi (mols);
			density = atof (dens);
			break;
		default:
			break;
		}
		
		/* store in memory (using PSL) */
		atom_number = put_atom (center, radius, ball_radius, covalent, atom, res, seq,
			kind, pdb, subunit, chemical_element, vdw_type,
			anumber, srn, mol, occupancy, tfactor, density);
		if (atom_number == 0) break;	/* problem */
		
		n_atom++;
		/* add to standard sets */
		include_member (atoms, ATOM, atom_number);
		include_member (these, ATOM, atom_number);
		if (set_number > 0)
			include_member (set_number, ATOM, atom_number);
		if (vdw_type == 0) {
			include_member (unknowns, ATOM, atom_number);
			n_vdw_zero++;
			for (ii = 0; ii < 4; ii++)
				preface[ii] = *(in_line + 72 + ii);
			preface[4] = 0;
			for (ii = 0; ii < 26; ii++)
				shortened[ii] = *(in_line + ii);
			shortened[26] = 0;
			sprintf(message, "    %4s %26s unknown atom type given default radius %8.3f",
				preface, shortened, DEFAULT_ATOM_RADIUS);
			inform (message);
		}
		if (radius < min_r) min_r = radius;
		if (radius > max_r) max_r = radius;
	}

	/* finished with reading atomic coordinate file */
	n_these = count_set (these);
	sprintf (message,
		"%8ld atoms read from file, with radii varying between %6.3f and %6.3f",
		n_these, min_r, max_r);
	inform(message);
	if (n_vdw_zero > 0) {
		sprintf(message, "%8ld atoms of unknown type", n_vdw_zero);
		inform (message);
	}
	return (1);
}


int read_radii_file (struct array *rlines)
{
	long n_line, n_read, nw, l;
	struct radiusRecord *ar;
	char *in_line;
	char message[MAX_STRING];

	n_line = rlines -> length;
	if (n_line > MAX_RADII) {
		sprintf (message, "read_radius_file: %6ld > %6ld",
			n_line, (long) MAX_RADII);
		set_error1 (message);
		return (0);
	}
	n_radii = n_line;
	/* allocate memory */
	radius_table = (struct radiusRecord *) allocate_objects (RADIUSRECORD, n_radii);
	if (radius_table == NULL) {
		sprintf (message, "msroll: cannot allocate radius table");
		set_error1 (message);
		return (0);
	}
	n_read = 0;
	for (l = 0; l < n_radii; l++) {
		in_line = fetch_string (rlines, l);
		if (in_line == NULL) break;
		ar = radius_table + l;
		nw = sscanf (in_line, "%ld %lf %lf %s",
			&(ar -> type), &(ar -> radius), &(ar -> covalent), ar -> kind);
		if (nw < 2) {
			sprintf (message, "msroll: invalid radius format");
			set_error1 (message);
			return (0);
		}
		n_read++;
	}
	sprintf (message, "%8ld atom types and radii read", n_read);
	inform (message);
	return (1);
}

int read_pattern_file (struct array *plines)
{
	long n_line, n_read, nw, l;
	struct patternRecord *pat;
	char *in_line;
	char residue[MAX_STRING];
	char atom_name[MAX_STRING];
	long type;
	char kind[MAX_STRING];
	char message[MAX_STRING];

	n_line = plines -> length;
	if (n_line > MAX_PATTERNS) {
		sprintf (message, "read_pattern_file: %6ld > %6ld",
			n_line, (long) MAX_PATTERNS);
		set_error1 (message);
		return (0);
	}
	n_pattern = n_line;
	/* allocate memory */
	pattern_table = (struct patternRecord *) allocate_objects (PATTERNRECORD, n_pattern);
	if (pattern_table == NULL) {
		sprintf (message, "msroll: cannot allocate pattern table");
		set_error1 (message);
		return (0);
	}
	n_read = 0;
	for (l = 0; l < n_pattern; l++) {
		in_line = fetch_string (plines, l);
		if (in_line == NULL) break;
		strcpy (residue, "");
		strcpy (atom_name, "");
		strcpy (kind, "");
		pat = pattern_table + l;
		nw = sscanf (in_line, "%s %s %ld %s",
			residue, atom_name, &type, kind);
		if (nw < 3) {
			sprintf (message, "msroll: invalid pattern format");
			set_error1 (message);
			return (0);
		}
		if (strlen (residue) > (unsigned long) MAX_ATNAME) {
			sprintf (message, "msroll: residue name too long %s", residue);
			set_error1 (message);
			return (0);
		}
		if (strlen (atom_name) > (unsigned long) MAX_ATNAME) {
			sprintf (message, "msroll: atom name too long %s", atom_name);
			set_error1 (message);
			return (0);
		}
		if (strlen (kind) > (unsigned long) MAX_ATNAME) {
			sprintf (message, "msroll: atom kind too long %s", kind);
			set_error1 (message);
			return (0);
		}
		strcpy (pat -> residue, residue);
		strcpy (pat -> atom_name, atom_name);
		pat -> type = type;
		strcpy (pat -> kind, kind);
		n_read++;
	}
	sprintf (message, "%8ld patterns read", n_read);
	inform (message);
	return (1);
}

long create_residue_index ()
{
	long n, l;
	char previous_residue[MAX_ATNAME];
	struct patternRecord *pat;
	struct residueIndex *rin;

	strcpy (previous_residue, "");
	n_residue_index = 0;
	for (l = 0; l < n_pattern; l++) {
		pat = pattern_table + l;
		if (strcmp (pat -> residue, previous_residue) != 0) {
			n_residue_index++;
			strcpy (previous_residue, pat -> residue);
		}
	}
	if (n_residue_index <= 0) {
		set_error1 ("create_residue_index: underflow");
		return (0L);
	}
	residue_table = (struct residueIndex *)
		allocate_objects (RESIDUEINDEX, n_residue_index);
	if (residue_table == NULL) {
		set_error1 ("create_residue_index: memory failure");
		return (0L);
	}
	n = -1; rin = NULL;
	strcpy (previous_residue, "");
	for (l = 0; l < n_pattern; l++) {
		pat = pattern_table + l;
		if (strcmp (pat -> residue, previous_residue) != 0) {
			n++;
			rin = residue_table + n;
			rin -> first_pattern = l;
			strcpy (rin -> residue, pat -> residue);
			strcpy (previous_residue, pat -> residue);
		}
		rin -> last_pattern = l;
	}
	return (n_residue_index);
}

/* second string is from pattern table */
int pattern_match (char str1[6], char str2[6])
{
	int i, c1, c2;
	if (strcmp(str2,"*") == 0) return (1);
	for (i = 0; i < 6; i++) {
		c1 = str1[i]; c2 = str2[i];
		if (c1 == 0 || c2 == 0) break;
		if (c2 == '?') continue;
		/* always ignore case */
		c1 = tolower(c1);
		c2 = tolower(c2);
		if (c1 != c2) return (0);
	}
	return (c1 == 0 && c2 == 0);
}

int interpret_pdb (char in_line[100], char pdb[8], char ans[6], char atom[6], char res[6], char subunit[6], char seq[6], char xs[10], char ys[10], char zs[10], char os[8], char tfs[8])
{
	int i;

	/* first field is ATOM or HETATM for records we care about */
	for (i = 0; i < 6; i++)
		pdb[i] = in_line[i];
	pdb[6] = 0; remove_tb (pdb); remove_lb(pdb);
	if (strcmp (pdb, "ATOM") != 0 &&
		strcmp(pdb, "HETATM") != 0) return (0);

	/* atom number */
	for (i = 0; i < 5; i++)
		ans[i] = in_line[6+i];
	ans[5] = 0; remove_tb (ans); remove_lb(ans);

	/* atom name */
	for (i = 0; i < 4; i++)
		atom[i] = in_line[13+i];
	/* look for atom name shifted left */
	if (in_line[12] != ' ') {
		for (i = 0; i < 4; i++)
			atom[i] = in_line[12+i];
	}
	/* look for embedded blank in second position */
	if (atom[0] != ' ' && atom[1] == ' ' && atom[2] != ' ') {
		atom[1] = '_';
	}
	atom[4] = 0; remove_tb (atom); remove_lb(atom);
	
	/* residue name */
	for (i = 0; i < 4; i++)
		res[i] = in_line[17+i];
	res[4] = 0; remove_tb (res); remove_lb(res);
	
	/* subunit */
	subunit[0] = in_line[21];
	subunit[1] = 0; remove_tb (subunit); remove_lb(subunit);
	
	/* sequence number is stored as a character string */
	for (i = 0; i < 5; i++)
		seq[i] = in_line[22+i];
	seq[5] = 0; remove_tb (seq); remove_lb(seq);
	
	/* atom center */
	for (i = 0; i < 8; i++) {
		xs[i] = in_line[30+i];
		ys[i] = in_line[38+i];
		zs[i] = in_line[46+i];
	}
	xs[8] = 0;
	ys[8] = 0;
	zs[8] = 0;

	/* occupancy and temperature factor */
	for (i = 0; i < 6; i++) {
		os[i] = in_line[54+i];
		tfs[i] = in_line[60+i];
	}
	os[7] = 0;
	tfs[7] = 0;
	return(1);
}

int interpret_xyzr (char in_line[100], char xs[10], char ys[10], char zs[10], char rs[6], char res[6], char seq[6], char atom[6])
{
	sscanf (in_line, "%s %s %s %s %s %s %s",
		xs, ys, zs, rs, res, seq, atom);
	return (1);
}

int interpret_ms (char in_line[100], char xs[10], char ys[10], char zs[10], char ts[6], char srns[6], char mols[6], char res[6], char seq[6], char atom[6])
{
	sscanf (in_line, "%s %s %s %s %s %s %s %s %s",
		xs, ys, zs, ts, srns, mols, res, seq, atom);
	return (1);
}

int interpret_ds (char in_line[100], char xs[10], char ys[10], char zs[10], char rs[6], char srns[6], char mols[6], char dens[6], char res[6], char seq[6], char atom[6])
{
	sscanf (in_line, "%s %s %s %s %s %s %s %s %s %s",
		xs, ys, zs, rs, srns, mols, dens, res, seq, atom);
	return (1);
}

int select_bonds (int atom_set)
{
	int b, a1, a2, as[2];
	int bond_set;
	if (count_set (bonds) == 0) return (0);
	bond_set = allocate_set (BOND);
	/* go through bond set */
	for (b = init_for (bonds); b != 0; b = next_for (bonds)) {
		/* obtain atom object numbers */
		get_bond_atoms (b, as);
		if (error()) {
			return (0);
		}
		a1 = as[0];
		a2 = as[1];
		/* check for membership in set */
		if (!member_of (a1, atom_set)) continue;
		if (!member_of (a2, atom_set)) continue;
		include_member (bond_set, BOND, b);
	}
	return (bond_set);
}

long remove_duplicate_atoms (int atom_set)
{
	long n, i, j, n_duplicate;
	int o, k, ok;
	double d, radius1, radius2;
	double center1[3], center2[3];
	struct array *numbers, *coordinates;
	double   rect[3][2];
	atomnum *jnbr, *inrect;
	char message[MAX_STRING];



	if (atom_set == 0) return (0L);
	n = count_set (atom_set);
	if (n == 0L) return (0L);
	coordinates = new_array (SPHERE, n);
	if (coordinates == NULL) return (0L);
	/* kludge */
	natom = n;
	atmco = coordinates -> centers;
	atmrad = coordinates -> radii;
	atmatt = NULL;
	numbers = new_array (INTEGER, n);
	if (numbers == NULL) return (0L);

	i = 0;
	for (o = init_for (atom_set); o != 0; o = next_for (atom_set)) {
		get_atom_center (o, center1);
		radius1 = get_atom_radius (o);
		if (error()) return (0L);
		store_sphere (coordinates, i, center1, radius1);
		if (error ()) return (0L);
		i++;
	}
        /* determine number of entries in table */
        determine_ntab ();
	dsbtree ();     /* create box tree */
	if (errflg) {
		set_error1 (errstr);
		sprintf (message, "remove_duplicate_atoms: dsbtree returns error for atom %d", atom+1);
		set_error2 (message);
		return(0L);
	}

	n_duplicate = 0;
	for (i = 0, o = init_for (atom_set); o != 0; i++, o = next_for (atom_set)) {
		radius1 = fetch_sphere (coordinates, i, center1);
		if (error ()) return (0L);
        /* set up rectangular range for searching box tree */
        /* make it big enough to include any possible nbr of atom ia */
        for (k = 0; k < 3; k++) {
            rect[k][0] = center1[k] - EPSILON;
            rect[k][1] = center1[k] + EPSILON;
        }
        /* get atoms in rectangle */
        inrect = dsgair (rect, rootbox);
        if (errflg) {
        	set_error2 ("remove_duplicate_atoms: dsgair returns error");
            return(0L);
        }
        if (inrect == NULL) continue;
		/* check for duplicate */
		ok = 1;
        /* loop through atoms in rectangular range */
        for (jnbr = inrect; *jnbr != EOL; jnbr++) {
            j = *jnbr;         /* retrieve atom index from list */
			if (i <= j) continue;
			radius2 = fetch_sphere (coordinates, j, center2);
			if (error ()) return (0L);
			d = distance (center1, center2);
			if (d < EPSILON) {
				ok = 0;
				break;
			}
		}
		if (!ok) {
			store_integer (numbers, i, o);	/* flag as duplicate */
			if (error ()) return (0L);
			n_duplicate++;
		}
		/* free list from box tree search */
		free_atomnums (inrect);
	}
	dsfbox (rootbox);   /* free box (recursive) */
	free_cache(BOX);
	if (coordinates != NULL) free_array (coordinates);
	if (n_duplicate == 0) {
		if (numbers != NULL) free_array (numbers);
		return (n_duplicate);
	}

	atmco = NULL; atmrad = NULL;
	if (errflg) {
		set_error2 ("remove_duplicate_atoms: dsfbox returns error");
		return (0L);
	}
	if (n_duplicate >= n) {
		set_error1 ("remove_duplicate_atoms: none left");
		return (0L);
	}
	/* remove duplicate atoms */
	for (i = 0; i < n; i++) {
		o = fetch_integer (numbers, i);
		if (error ()) return (0L);
		if (o != 0) {
			exclude_member (atom_set, ATOM, o);
			if (error ()) return (0L);
		}
	}
	if (numbers != NULL) free_array (numbers);
	return (n_duplicate);
}

/* SPHERE AND PLANE CODE */

/* obtain object number of sphere or plane */

int lookup_region (char *vname)
{
	int region_number;
	
	region_number = object_number_of (vname);
	return (region_number);
}

/* define a sphere region */

void sphere_command (char *spherename, float x, float y, float z, float r)
{
	int region_number;
	double center[3], radius;
	
	define_region (spherename, SPHERE, 2);
	region_number = object_number_of (spherename);
	if (error ()) {
		return;
	}
	center[0] = x; center[1] = y; center[2] = z;
	radius = r;
	set_region_center (region_number, center);
	set_region_radius (region_number, radius);
}

/* define a plane region */

void plane_command (char *planename, float x, float y, float z, float xn, float yn, float zn)
{
	int region_number;
	double center[3], vector[3];
	
	define_region (planename, PLANE, 2);
	region_number = object_number_of (planename);
	if (error ()) {
		return;
	}
	center[0] = x; center[1] = y; center[2] = z;
	vector[0] = xn; vector[1] = yn; vector[2] = zn;
	set_region_center (region_number, center);
	normalize (vector);
	set_region_direction (region_number, vector);
}

void define_region (char *name, int type, int dim)
{
	int  object_number, symbol_number;
	struct region *vty_ptr;

	symbol_number = new_symbol (name, 1, REGION);
	if (error()) return;
	object_number = allocate_named (REGION, symbol_number);
	if (object_number == 0) return;
	vty_ptr = region_ptr (object_number);
	if (error()) return;
	vty_ptr -> type = (short) type;
	if (vty_ptr -> type == SPHERE || vty_ptr -> type == CIRCLE)
		vty_ptr -> radii[0] = 1.0;
	if (vty_ptr -> type == CIRCLE)
		vty_ptr -> axis[2] = 1.0;
}

void init_region ()
{
	memory_init (REGION);
	if (error ()) {
		return;
	}
}

/* second radius not currently accessible */
double get_region_radius (int variety_number)
{
	double radius;
	struct region *vty_ptr;

	vty_ptr = region_ptr (variety_number);
	if (vty_ptr == (struct region *) NULL) return (0.0);
	radius =  vty_ptr -> radii[0];
	return (radius);
}

void set_region_radius (int variety_number, double radius)
{
	struct region *vty_ptr;

	vty_ptr = region_ptr (variety_number);
	if (vty_ptr == (struct region *) NULL) return;
	vty_ptr -> radii[0] =  radius;
}

void get_region_center (int variety_number, double center[3])
{
	int k;
	struct region *vty_ptr;

	vty_ptr = region_ptr (variety_number);
	if (vty_ptr == (struct region *) NULL) return;
	for (k = 0; k < 3; k++)
		center[k] =  vty_ptr -> center[k];
}

void set_region_center (int variety_number, double center[3])
{
	int k;
	struct region *vty_ptr;

	vty_ptr = region_ptr (variety_number);
	if (vty_ptr == (struct region *) NULL) return;
	for (k = 0; k < 3; k++)
		vty_ptr -> center[k] =  center[k];
}

void get_region_direction (int variety_number, double direction[3])
{
	int k;
	struct region *vty_ptr;

	vty_ptr = region_ptr (variety_number);
	if (vty_ptr == (struct region *) NULL) return;
	for (k = 0; k < 3; k++)
		direction[k] =  vty_ptr -> axis[k];
}

void set_region_direction (int variety_number, double direction[3])
{
	int k;
	struct region *vty_ptr;

	vty_ptr = region_ptr (variety_number);
	if (vty_ptr == (struct region *) NULL) return;
	for (k = 0; k < 3; k++)
		vty_ptr -> axis[k] =  direction[k];
}

int region_size ()
{
	return (sizeof (struct region));
}

struct region *region_ptr (int number)
{
	return ((struct region *) generic_ptr (REGION, number));
}

void remove_tb (char *s)
{
	int l;
	char *p;
	
	l = strlen(s);
	if (l <= 0) return;
	p = s + l - 1;
	while (p >= s) {
		if (*p == ' ')
			*p-- = 0;
		else break;
	}
	return;
}

void remove_lb (char *s)
{
	int i, l;
	char *p;
	
	while (*s == ' ') {
		l = strlen(s);
		if (l <= 0) return;
		for (i = 0, p = s; i < l; i++, p++)
			*p = *(p+1);
	}
	return;
}

void help_scanner (char command_line[100])
{
	int i, j;
	char message[MAX_STRING];
	char temp_line[120];
	
	for (i = 0,j = 0; i < 100; i++,j++) {
		if (i > 0 && isalnum(command_line[i-1]) &&
			isoperator(command_line[i])) {
			temp_line[j] = ' '; j++;
		}
		temp_line[j] = command_line[i];
		if (temp_line[j] == (char) 0) break;
		if (isoperator(command_line[i]) &&
			isalnum(command_line[i+1])) {
			 j++; temp_line[j] = ' ';
		}
	}
	
	if (strlen (temp_line) > (unsigned) 99) {
		sprintf (message, "command line too long");
		set_error1 (message);
		return;
	}
	strcpy (command_line, temp_line);
}

int isoperator (char c)
{
	if (c == '+') return (1);	
	if (c == '-') return (1);	
	if (c == '*') return (1);	
	if (c == '=') return (1);
	if (c == '!') return (1);
	if (c == '<') return (1);
	if (c == '>') return (1);
	return (0);
}


int parse_sphere (struct msscene *ms, int n_words, char words[MAX_WORD][MAX_NAME])
{
	double x, y, z, r;
	/* sphere spherename x y z radius */
	if (n_words != 6) {
		set_error1 ("invalid sphere command");
		return (0);
	}
	x = atof (words[2]);
	y = atof (words[3]);
	z = atof (words[4]);
	r = atof (words[5]);
	
	sphere_command (words[1], x, y, z, r);
	return (1);
}

int parse_plane (struct msscene *ms, int n_words, char words[MAX_WORD][MAX_NAME])
{
	double x, y, z, xn, yn, zn;
	/* plane planename x y z xn yn zn */
	if (n_words != 8) {
		set_error1 ("invalid plane command");
		return (0);
	}
	x = atof (words[2]);
	y = atof (words[3]);
	z = atof (words[4]);
	xn = atof (words[5]);
	yn = atof (words[6]);
	zn = atof (words[7]);
	plane_command (words[1], x, y, z, xn, yn, zn);
	return (1);
}

int parse_clear (struct msscene *ms, int n_words, char words[MAX_WORD][MAX_NAME])
{
	/* clear setname */
	if (n_words != 1) {
		set_error1 ("user not allowed to modify protected set");
		return (0);
	}
	if (protected_set(words[1])) {
		set_error1 ("user not allowed to modify protected set");
		return (0);
	}
	clear_command (words[1]);
	return (1);
}

int pdb_command (char command_line[MAXLINE])
{
	int nw;
	char command__line[100];
	char setname[MAX_NAME], setname1[MAX_NAME];
	char setname2[MAX_NAME], setname3[MAX_NAME];
	char field[MAX_NAME], operator[MAX_NAME];
	char value[MAX_NAME];
	char message[MAX_STRING];

	
	/* purpose: to allow hyphens in filenames */
	strcpy (command__line, command_line);
	help_scanner (command__line);
	
	/* set *= field operator value */
	nw = sscanf (command__line, "%s *= %s %s %s",
		setname, field, operator, value);
	if (nw == 4) {
		times_equals (setname, field, operator, value);
		return (1);
	}
	
	/* set += field operator value */
	nw = sscanf (command__line, "%s += %s %s %s",
		setname, field, operator, value);
	if (nw == 4) {
		if (protected_set(setname)) {
			sprintf (message, "pdb_command: user not allowed to modify protected set");
			set_error1 (message);
			return (0);
		}
		plus_equals (setname, field, operator, value);
		return (1);
	}
	
	/* set -= field operator value */
	nw = sscanf (command__line, "%s -= %s %s %s",
		setname, field, operator, value);
	if (nw == 4) {
		if (protected_set(setname)) {
			sprintf (message, "pdb_command: user not allowed to modify protected set");
			set_error1 (message);
			return (0);
		}
		minus_equals (setname, field, operator, value);
		return (1);
	}
	
	/* set field = value */
	nw = sscanf (command__line, "%s %s = %s",
		setname, field, value);
	if (nw == 3) {
		field_assignment (setname, field, value);
		return (1);
	}

	/* set1 = set2 * set3 */
	nw = sscanf (command__line, "%s = %s * %s",
		setname1, setname2, setname3);
	if (nw == 3) {
		if (protected_set(setname1)) {
			sprintf (message,"pdb_command: user not allowed to modify protected set");
			set_error1 (message);
			return (0);
		}
		handle_intersection (setname1, setname2, setname3);
		return (1);
	}

	/* set1 = set2 + set3 */
	nw = sscanf (command__line, "%s = %s + %s",
		setname1, setname2, setname3);
	if (nw == 3) {
		if (protected_set(setname1)) {
			sprintf (message,"pdb_command: user not allowed to modify protected set");
			set_error1 (message);
			return (0);
		}
		handle_union (setname1, setname2, setname3);
		return (1);
	}

	/* set1 = set2 - set3 */
	nw = sscanf (command__line, "%s = %s - %s",
		setname1, setname2, setname3);
	if (nw == 3) {
		if (protected_set(setname1)) {
			sprintf (message,"pdb_command: user not allowed to modify protected set");
			set_error1 (message);
			return (0);
		}
		handle_subtraction (setname1, setname2, setname3);
		return (1);
	}

	/* set = field operator value */
	nw = sscanf (command__line, "%s = %s %s %s",
		setname, field, operator, value);
	if (nw == 4) {
		if (protected_set(setname)) {
			sprintf (message, "pdb_command: user not allowed to modify protected set");
			set_error1 (message);
			return (0);
		}
		assignment (setname, field, operator, value);
		return (1);
	}
	else if (nw == 2) {
		if (protected_set(setname1)) {
			sprintf (message, "pdb_command: user not allowed to modify protected set");
			set_error1 (message);
			return (0);
		}
		simple_set_copy (setname1, setname2);
		return (1);
	}

	return (0);
}


void read_atom_commands (struct msscene *ms, FILE *fp_command)
{
	int n_words, result, w;
	char command_line[MAXLINE+1];
	char words[MAX_WORD][MAX_NAME];
	char verb[MAXLINE];


	for (;;) {
		fgets (command_line, MAXLINE, fp_command);
		if (feof (fp_command)) break;
		informd (command_line);
		for (w = 0; w < MAX_WORD; w++)
			strcpy (words[w], " ");
		null_terminate (command_line, '#');
		n_words = sscanf (command_line, "%s %s %s %s %s %s %s %s %s %s",
			words[0], words[1], words[2], words[3], words[4],
			words[5], words[6], words[7], words[8], words[9]);
		if (n_words <= 0) continue;
		if (words[0][0] == '!') continue;
		if (words[0][0] == '#') continue;
		strcpy (verb, words[0]);
		if (strcmp (verb, "sphere") == 0) {
			result = parse_sphere (ms, n_words, words);
		}
		else if (strcmp (verb, "plane") == 0) {
			result = parse_plane (ms, n_words, words);
		}
		else if (strcmp (verb, "clear") == 0) {
			result = parse_clear (ms, n_words, words);
		}
		else {
			result = pdb_command (command_line);
		}
		if (result == 0 && !error ()) {
			set_error1("unrecognizable command");
			set_error2(command_line);
			return;
		}
		if (!result) break;
	}
}

int free_all_pdb ()
{
	free_objects (RESIDUEINDEX, (short *) residue_table);
	free_object_blocks (ATOM);
	free_object_blocks (BOND);
	return (1);
}

