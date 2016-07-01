/*

	MSRoll

	Copyright 2006 by Michael L. Connolly
	All rights reserved

	Written by Michael L. Connolly.
	February 16, 2006

	The usage is:
	
	msroll -m molecule_file -r radii_file -y pattern_file 
			-p probe_radius
			-l coalesce_radius
			-n name
			-f fineness
			-x cusp_intersection_passes
			-a area_file [bac | bca] -v volume_file -q surface_file
			-e error_file
			-t polygon_file
			-c cavity_file
	
	The probe radius defaults to 1.5 (or whatever is in mspq.h).
	One should specify one of area or volume or surface,
	otherwise the run will not produce any useful output.
	The error_file contains error message and informative messages.
	It may be "stdout", "stderr" or a disk file.
	The area file and volume file can be "stdout" or a disk file.
	The optional second argument for -a concerns the sorting
	of the atomic area output file. The bca means by component and atom,
	the bac means by atom and component. Choose bca if you want
	each component (cavity) to be grouped together. If either bca or bac
	is chosen, the same atom may show up several times, once
	for each cavity (or external surface) it borders.
	The largest surface component is always component 1.
	It is always the outer surface. 
	Volumes for internal cavities are written as negative values,
	to distinguish them from solvent atoms separate from the
	main molecule.
*/

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"


struct argrec argarray[] = {
	{ 'a', 0, 0, 0, 0.0, ""},
	{ 'c', 0, 0, 0, 0.0, ""},
	{ 'd', 0, 0, 0, 0.0, ""},
	{ 'e', 0, 0, 0, 0.0, ""},
	{ 'f', 0, 0, 0, 0.0, ""},
	{ 'g', 0, 0, 0, 0.0, ""},
	{ 'j', 0, 0, 0, 0.0, ""},
	{ 'k', 0, 0, 0, 0.0, ""},
	{ 'l', 0, 0, 0, 0.0, ""},
	{ 'm', 0, 0, 0, 0.0, ""},
	{ 'n', 0, 0, 0, 0.0, ""},
	{ 'o', 0, 0, 0, 0.0, ""},
	{ 'p', 0, 0, 0, 0.0, ""},
	{ 'r', 0, 0, 0, 0.0, ""},
	{ 'q', 0, 0, 0, 0.0, ""},
	{ 't', 0, 0, 0, 0.0, ""},
	{ 'v', 0, 0, 0, 0.0, ""},
	{ 'w', 0, 0, 0, 0.0, ""},
	{ 'x', 0, 0, 0, 0.0, ""},
	{ 'y', 0, 0, 0, 0.0, ""},
};

char *intprop = "Copyright 2006 by Michael L. Connolly";
char *southam = "All rights reserved";

/* the main function does I/O and calls the functions that compute */

main (int argc, char *argv[])
{
	int nitems, ngotten, nw, atom_set;
	int arg, result, connected;
	int dox;
	unsigned long size1, size2;
	double pr,grid, coalesce;
	double alpha;
	char fname[64];
	char molname[64];
	char molfilename[MAXLINE];
	double surface_center[3];
	char order[MAXLINE];
	char message[MAXLINE];
	struct msscene *ms;
	struct molecule *mol;
	FILE *fpm;
	FILE *fpr;
	FILE *fpe;
	FILE *fpf;
	FILE *fpa;
	FILE *fpq;
	FILE *fpv;
	FILE *fpy;
	FILE *fpj;
	FILE *fpk;
	FILE *fpt;
	FILE *fpc;
	
	fperror = stderr;
	fpinform = stderr;
	fpdebug = stderr;

	sprintf (message,
	"msroll   version %d.%d.%d    Copyright 2006 by Michael L. Connolly",
		(int) VERSION, (int) SUBVERSION, (int) STUBVERSION);
	inform(message);

	/* no more interactive */	
	if (argc <= 2) {
		fprintf (stderr, "usage:\tmsroll -m pdb_file [-r radii_file] [-y pattern_file]\n");
		fprintf (stderr, "\t\t[-p probe_radius] [-a area_file [bac | bca]] [-v volume_file]\n");
		fprintf (stderr, "\t\t[-q surface_file] [-t polygon_file] [-c cavity_file]\n");
		fprintf (stderr, "\t\t[-j surface_point_file] [-k ds_area_file]\n");
		fprintf (stderr, "\t\t[-n name] [-f fineness_file | fineness_angle]\n");
		fprintf (stderr, "\t\t[-x cusp_intersection_passes]\n");
		fprintf (stderr, "\t\t[-e error_file] [-l coalesce_radius]\n");
		exit (0);
	}

	/* initialize */
	fpm = fpr = fpt = fpa = fpv = fpy = (FILE *) NULL;	/* no default */
	fpq = fpf = fpj = fpk = fpc = (FILE *) NULL;	/* no default */
	fpe = stderr;			/* default */
	pr = -1.0; 				/* if not specified, msroll will use default */
	grid = 0.0; 			/* no neighbor or collision grid */
	connected = 0;
	coalesce = DEFAULT_COALESCE;
	strcpy(order,"");		/* empty will cause msroll to use default (bca-bac) */
	strcpy(molfilename, "");			
	strcpy(molname, "");					
	alpha = 0.0; 			/* if not specified, trb will use default */
	dox = DEFAULT_PASS;
	
	
	/* streams and flags */
	/* get flags from command line */
	size1 = sizeof(argarray);
	size2 = sizeof(struct argrec);
	nitems = size1/size2;
	ngotten = get_arguments (argc, argv, argarray, nitems);
	if (ngotten <= 0) {
		if (ngotten < -1)
			fprintf (stderr, "msroll: invalid flag: %c\n", -ngotten);
		else fprintf (stderr, "msroll: problem reading command-line arguments\n");
		exit (1);
	}
	for (arg = 0; arg < nitems; arg++) {
		if (!argarray[arg].on) continue;
		switch (argarray[arg].flag) {
		case 'm':
			if (strcmp(argarray[arg].stringarg,"stdin") == 0) fpm = stdin;
			else {
				strcpy (molfilename, argarray[arg].stringarg);
				fpm = fopen (molfilename, "r");
				if (fpm == NULL) {
					fprintf (stderr,
						"msroll: cannot open molecule file: %s\n",
							argarray[arg].stringarg);
					exit (1);
				}
			}
			break;
		case 'r':
			if (strcmp(argarray[arg].stringarg,"stdin") == 0) fpr = stdin;
			else {
				fpr = fopen (argarray[arg].stringarg, "r");
				if (fpr == NULL) {
					fprintf (stderr,
						"msroll: cannot open radii file: %s\n",
							argarray[arg].stringarg);
					exit (1);
				}
			}
			break;
		case 'e':
			if (strcmp(argarray[arg].stringarg,"stderr") == 0) fpe = stderr;
			else if (strcmp(argarray[arg].stringarg,"stdout") == 0) fpe = stdout;
			else {
				fpe = fopen (argarray[arg].stringarg, "w");
				if (fpe == NULL) {
					fprintf (stderr,
						"msroll: cannot open error file: %s\n",
						argarray[arg].stringarg);
					exit (1);
				}
			}
			break;
		case 'q':
			fpq = fopen (argarray[arg].stringarg, "wb");	/* binary */
			if (fpq == NULL) {
				fprintf (stderr,
					"msroll: cannot open surface file: %s\n",
					argarray[arg].stringarg);
				exit (1);
			}
			break;
		case 'a':
			nw = sscanf(argarray[arg].stringarg, "%s %s", fname, order);
			if (nw < 1) continue;
			if (strcmp(fname,"stderr") == 0) fpa = stderr;
			else if (strcmp(fname,"stdout") == 0) fpa = stdout;
			else {
				fpa = fopen (fname, "w");
				if (fpa == NULL) {
					fprintf (stderr, "msroll: cannot open area file: %s\n", fname);
					exit (1);
				}
			}
			break;
		case 'v':
			if (strcmp(argarray[arg].stringarg,"stderr") == 0) fpv = stderr;
			else if (strcmp(argarray[arg].stringarg,"stdout") == 0) fpv = stdout;
			else {
				fpv = fopen (argarray[arg].stringarg, "w");
				if (fpv == NULL) {
					fprintf (stderr,"msroll: cannot open volume file: %s\n",
						argarray[arg].stringarg);
					exit (1);
				}
			}
			break;
		case 'j':
			if (strcmp(argarray[arg].stringarg,"stderr") == 0) fpj = stderr;
			else if (strcmp(argarray[arg].stringarg,"stdout") == 0) fpj = stdout;
			else {
				fpj = fopen (argarray[arg].stringarg, "w");
				if (fpj == NULL) {
					fprintf (stderr,"msroll: cannot open surface point file: %s\n",
						argarray[arg].stringarg);
					exit (1);
				}
			}
			break;
		case 'k':
			if (strcmp(argarray[arg].stringarg,"stderr") == 0) fpk = stderr;
			else if (strcmp(argarray[arg].stringarg,"stdout") == 0) fpk = stdout;
			else {
				fpk = fopen (argarray[arg].stringarg, "w");
				if (fpk == NULL) {
					fprintf (stderr,"msroll: cannot open ds area file: %s\n",
						argarray[arg].stringarg);
					exit (1);
				}
			}
			break;
		case 'p':
			pr = argarray[arg].doublearg;
			break;
		case 'l':
			coalesce = argarray[arg].doublearg;
			if (coalesce <= 0.0) coalesce = 0.0;
			break;
		case 'g':
			grid = argarray[arg].doublearg;
			break;
		case 'y':
			if (strcmp(argarray[arg].stringarg,"stdin") == 0) fpy = stdin;
			else {
				fpy = fopen (argarray[arg].stringarg, "r");
				if (fpy == NULL) {
					fprintf (stderr,
						"msroll: cannot open pattern file: %s\n",
							argarray[arg].stringarg);
					exit (1);
				}
			}
			break;
		case 't':
			nw = sscanf(argarray[arg].stringarg, "%s", fname);
			if (nw < 1) continue;
			if (strcmp(fname,"stdout") == 0) fpt = stdout;
			else {
				fpt = fopen (fname, "w");
				if (fpt == NULL) {
					fprintf (stderr,
						"msroll: cannot open polyhedron file: %s\n", fname);
					exit (1);
				}
			}
			break;
		case 'c':
			nw = sscanf(argarray[arg].stringarg, "%s", fname);
			if (nw < 1) continue;
			if (strcmp(fname,"stdout") == 0) fpc = stdout;
			else {
				fpc = fopen (fname, "w");
				if (fpc == NULL) {
					fprintf (stderr,
						"msroll: cannot open cavity file: %s\n", fname);
					exit (1);
				}
			}
			break;
		case 'f':
			/* check whether real, stdin or file name */
			if (argarray[arg].doublearg != 0.0)
				alpha = argarray[arg].doublearg;
			else if (strcmp(argarray[arg].stringarg,"stdin") == 0) fpf = stdin;
			else {
				fpf = fopen (argarray[arg].stringarg, "r");
				if (fpf == NULL) {
					fprintf (stderr, "msroll: cannot open fineness file: %s\n",
						argarray[arg].stringarg);
					exit (1);
				}
			}
			break;
		case 'n':
			strcpy (molname, argarray[arg].stringarg);
			break;
		case 'x':
			dox = argarray[arg].longarg;
			break;
		case 'd':
			debug = argarray[arg].longarg;
			break;
		default:
			fprintf (stderr,
				"msroll: invalid flag: %s\n", argarray[arg].stringarg);
			exit (1);
		}
	}
	if (fpm == NULL) {
		fprintf (stderr, "msroll: missing molecule file name\n");
		exit (1);
	}

	init_routing (fpe, fpe, fpe);

	roll_mem ();
	if (error()) {
		print_error ();
		exit (1);
	}
	ms = new_msscene ();
	if (ms == NULL || error ()) {
		set_error2 ("msroll: new_msscene fails");
		print_error ();
		exit (1);
	}
	ms -> this_srf = (struct surface *) NULL;
	ms -> debug = 0;
	ms -> coalesce = coalesce;

	/* initialization of Protein Shape Library */
	psl_init ();
	if (error()) {
		print_error ();
		exit (1);
	}
	init_mol();
	if (error()) {
		print_error ();
		exit (1);
	}
	init_region();
	if (error()) {
		print_error ();
		exit (1);
	}
	init_pdb ();
	if (error()) {
		print_error ();
		exit (1);
	}

	ms -> table = NULL;

	mol = read_molecule (ms, molname, molfilename, alpha,fpm, fpr, fpy, fpf);
	if (error()) {
		print_error ();
		exit (1);
	}
	atom_set = mol -> atom_set;
	
	if (fpj != NULL || fpk != NULL) {
		result = ds_stream (ms, atom_set, fpe, fpj, fpk, pr, connected);
		if (error() || ! result) {
			print_error ();
			exit (1);
		}
	}

	if (fpq != NULL || fpa != NULL || fpv != NULL || fpt != NULL || fpc != NULL) {
		result = pqms_stream (ms,mol,fpe,fpq,fpa,fpv,pr,grid, dox, order, surface_center);
		if (error() || ! result) {
			print_error ();
			exit (1);
		}
	}

	if (fpt != NULL || fpc != NULL) {
		result = trb_stream (ms, fpe, fpt, fpc, mol);
		if (error() || ! result) {
			print_error ();
			exit (1);
		}
	}
	free_scene (ms);
	if (error()) exit (1);
	free_all_pdb ();
	if (error()) exit (1);
	free_all_psl ();
	if (error()) exit (1);
	print_counts();
	return(1);
}

