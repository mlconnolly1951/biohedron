/*
   MSTran

   MS Translate to other formats

   Copyright 1995 by Michael L. Connolly
   All Rights Reserved

   December 20, 2001
*/

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

int mstran (struct msscene *ms, FILE *fpt, FILE *fpd, FILE *fpm, FILE *fpg, FILE *fpo, int fn, int which_value, double minval, double maxval, struct color_ramp *ramp, double ctrlev, int uniform_hue)
{
	int result;
	long norg;
	double tolerance = 0.3;
	char message[MAX_STRING];
	struct surface *phn, *tan, *btn, *orgphn;
	struct surface *den;
	struct oringe *org, *head, *tail;
	struct molecule *mol;

	if (fpo == NULL) return(0);
	if (fpt != NULL) {
		phn = read_polyhedron (fpt);
		fclose(fpt);
		if (phn == NULL) return (0);
		if (which_value < 0 || minval == maxval)
			phn -> scheme = define_scheme (UNIFORM_COLORING,
				uniform_hue, which_value, minval, maxval);
		else
			phn -> scheme = define_scheme (VERTEX_COLORING,
				0, which_value, minval, maxval);
		if (error()) return (0);
		if (phn -> scheme == NULL) {
			set_error2 ("mstran: define_scheme returns null");
			return (0);
		}
		if (fn == INVENTOR) {
			write_header_iv (ms, fpo);
			if (error()) return (0);
			write_tri_iv (ms, phn, fpo, ramp);
		}
		else if (fn == WHATIF)
			write_whatif (phn, fpo);
		else if (fn == O)
			write_Os (phn, fpo, "dot");
		if (error()) return (0);
	}
	if (fpd != NULL) {
		den = read_density (fpd);
		fclose (fpd);
		if (den == NULL) return(0);
		result = polygonize_density (den, ctrlev);
		if (result == 0) return (0);
		den -> scheme = define_scheme (UNIFORM_COLORING,
			uniform_hue, 0, 0.0, 1.0);
		if (den -> scheme == NULL) return (0);
		write_header_iv (ms, fpo);
		if (error()) return (0);
		write_den_iv (ms, den, ctrlev, fpo);
		if (error()) return (0);
		/* close output file */
		fclose (fpo);
	}

	if (fpg != NULL) {
		norg = 0;
		head = NULL;
		tail = NULL;
		for (;;) {
			org = read_oringe (fpg);
			if (error ()) break;
			if (org == NULL) break; /* End of File */
			norg++;
			if (head == NULL) head = org;
			else tail -> next = org; 
			tail = org;
		}
		fclose (fpg);
		if (error ()) return (0);
		sprintf (message, "%ld oringe records read", norg);
		inform (message);
		orgphn = polyhedron_oringes (head, 0.2, ORG_SURFACE);
		if (error()) return (0);
		orgphn -> scheme = define_scheme (UNIFORM_COLORING, 2, -1, 0.0, 0.0);
		if (error()) return (0);
		if (orgphn -> scheme == NULL) {
			set_error2 ("mstran: define_scheme returns null");
			return (0);
		}
		write_header_iv (ms, fpo);
		if (error()) return (0);
		write_edg_iv (ms, orgphn, fpo);
		if (error()) return (0);
	}

	if (fpm != NULL) {
		mol = read_molecule (ms, "", "", DEFAULT_ALPHA, fpm, NULL, NULL, NULL);
		fclose(fpm);
		if (mol == NULL) return (0);
		if (error()) return (0);
		connect_command ("atoms", tolerance);
		write_header_iv (ms, fpo);
		if (error()) return (0);
		write_mol_iv (ms, mol, fpo);
	}
	fclose (fpo);
	return(1);
}

void tran_mem ()
{
	unsigned long size;
	int type;
	char type_name[32];
	
	init_mem ((unsigned long) N_OBJECTS);
	if (error()) return;
	
	type = MSSCENE;
	size = sizeof (struct msscene);
	strcpy (type_name, "msscene");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = SPHERE;
	size = sizeof (struct sphere);
	strcpy (type_name, "sphere");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = TORUS;
	size = sizeof (struct torus);
	strcpy (type_name, "torus");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = PROBE;
	size = sizeof (struct probe);
	strcpy (type_name, "probe");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = VERTEX;
	size = sizeof (struct vertex);
	strcpy (type_name, "vertex");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = CIRCLE;
	size = sizeof (struct circle);
	strcpy (type_name, "circle");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = ARC;
	size = sizeof (struct arc);
	strcpy (type_name, "arc");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = FACE;
	size = sizeof (struct face);
	strcpy (type_name, "face");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = CYCLE;
	size = sizeof (struct cycle);
	strcpy (type_name, "cycle");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = EDGE;
	size = sizeof (struct edge);
	strcpy (type_name, "edge");
	define_type (type, size, type_name);
	if (error()) return;
		
	type = COMPONENT;
	size = sizeof (struct component);
	strcpy (type_name, "component");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = SURFACE;
	size = sizeof (struct surface);
	strcpy (type_name, "surface");
	define_type (type, size, type_name);
	if (error()) return;

	type = VARIETY;
	size = sizeof (struct variety);
	strcpy (type_name, "variety");
	define_type (type, size, type_name);
	if (error()) return;

	type = PHNVTX;
	size = phnvtx_size ();
	strcpy (type_name, "phnvtx");
	define_type (type, size, type_name);

	type = PHNEDG;
	size = phnedg_size ();
	strcpy (type_name, "phnedg");
	define_type (type, size, type_name);

	type = PHNTRI;
	size = phntri_size ();
	strcpy (type_name, "phntri");
	define_type (type, size, type_name);

	type = POLYGON;
	size = sizeof (struct polygon);
	strcpy (type_name, "polygon");
	define_type (type, size, type_name);

	type = HEDVTX;
	size = sizeof (struct hedvtx);
	strcpy (type_name, "hedvtx");
	define_type (type, size, type_name);
	if (error()) return;

	type = HEDEDG;
	size = sizeof (struct hededg);
	strcpy (type_name, "hededg");
	define_type (type, size, type_name);
	if (error()) return;

	type = HEDTRI;
	size = sizeof (struct hedtri);
	strcpy (type_name, "hedtri");
	define_type (type, size, type_name);
	if (error()) return;

	type = VTXGRP;
	size = sizeof (struct vtxgrp);
	strcpy (type_name, "vtxgrp");
	define_type (type, size, type_name);
	if (error()) return;

	type = HEDRON;
	size = sizeof (struct hedron);
	strcpy (type_name, "hedron");
	define_type (type, size, type_name);
	if (error()) return;

	type = MOLECULE;
	size = sizeof (struct molecule);
	strcpy (type_name, "molecule");
	define_type (type, size, type_name);

	type = OBJECT_SCHEME;
	size = sizeof (struct object_scheme);
	strcpy (type_name, "object_scheme");
	define_type (type, size, type_name);
	if (error()) return;

	type = COLOR_RAMP;
	size = sizeof (struct color_ramp);
	strcpy (type_name, "color_ramp");
	define_type (type, size, type_name);
	if (error()) return;

	type = ORINGE;
	size = sizeof (struct oringe);
	strcpy (type_name, "oringe");
	define_type (type, size, type_name);
	if (error()) return;

	type = RECORD;
	size = sizeof (struct record);
	strcpy (type_name, "record");
	define_type (type, size, type_name);
	if (error()) return;

	type = TOKEN;
	size = sizeof (struct token);
	strcpy (type_name, "token");
	define_type (type, size, type_name);
	if (error()) return;

	type = EDGER;
	size = sizeof (struct edger);
	strcpy (type_name, "edger");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = MATERIAL_TABLE;
	size = sizeof (struct material_table);
	strcpy (type_name, "material_table");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = PATTERNRECORD;
	size = sizeof (struct patternRecord);
	strcpy (type_name, "patternRecord");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = RADIUSRECORD;
	size = sizeof (struct radiusRecord);
	strcpy (type_name, "radiusRecord");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = RESIDUEINDEX;
	size = sizeof (struct residueIndex);
	strcpy (type_name, "residueIndex");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = ARRAY;
	size = sizeof (struct array);
	strcpy (type_name, "array");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = ATOM;
	size = sizeof (struct atom);
	strcpy (type_name, "atom");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = BOND;
	size = sizeof (struct bond);
	strcpy (type_name, "bond");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = OBJECT_HEADER;
	size = sizeof (struct object_header);
	strcpy (type_name, "object_header");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = OBJECT_BLOCK;
	size = sizeof (struct object_block);
	strcpy (type_name, "object_block");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = ELEMENT;
	size = sizeof (struct element);
	strcpy (type_name, "element");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = ELEMENT_BLOCK;
	size = sizeof (struct element_block);
	strcpy (type_name, "element_block");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = SET;
	size = sizeof (struct set);
	strcpy (type_name, "set");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = SYMBOL;
	size = sizeof (struct symbol);
	strcpy (type_name, "symbol");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = BOOLEAN;
	size = sizeof (struct boolean);
	strcpy (type_name, "boolean");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = INTEGER;
	size = sizeof (struct integer);
	strcpy (type_name, "integer");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = REAL;
	size = sizeof (struct real);
	strcpy (type_name, "real");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = STRING;
	size = sizeof (struct string);
	strcpy (type_name, "string");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = REGION;
	size = sizeof (struct region);
	strcpy (type_name, "region");
	define_type (type, size, type_name);
	if (error()) return;
	
	type = CEPT;
	size = sizeof (struct cept);
	strcpy (type_name, "cept");
	define_type (type, size, type_name);
	if (error()) return;
	
}


/*
   mstran

   Copyright 1998 by Michael L. Connolly
   All Rights Reserved

*/
