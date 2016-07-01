/* Molecular Surface Package Copyright 1989 by Michael L. Connolly */
/* Written by Michael L. Connolly */
/* Last revised: January 27, 2006 */

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

 /* Developed originally on MASSCOMP 5500 */
 /* Further developed on Macintosh SE during 1987 and 1988 */
 /* Ported to Alliant & Raster -- June, July & November 1988 */
 /* Further Modified - November 1988 at Alliant */
 /* Further debugged on Macintosh SE during 1989-90 */
 /* In September 1990 kludge was added to skip over problem concave faces */
 /* Further debugged on Macintosh IIsi during 1991 */
 /* Note: remember that the Mac uses 16-bit integers */
 /* April 1991: some ints changed to longs */
 /* July 1991: some shorts changed to longs */
 /* July 1991: sunraster output formatted added and made default */
 /* October 1991: common msp.h file for binary surface file format */
 /* May 1992: common mspmath.c file for low-level math routines */
 /* January 1993: transluscent surfaces */
 /* February 1993: tubes */
 /* December 1994: AVS output formats, multiple writes of same image */
 /* January 1995: alpha transparency */
 /* October 1995: modularity: rearranging source code into source files */
 /* January 1996: merged srf and msplot to form msdraw */


int shorty;

struct material_table *table;

/* GLOBAL VARIABLES */
int current_object_type;
int current_object_number;
int current_symbol;
int n_ranges;
unsigned long scratch1, scratch2;
double pi_constant;
struct object_header *object_table;
struct element_block *first_element_block;
char function_name[MAX_NAME];

int error_flag;
int fatal_flag;
long debug;
char error_string[MAXLINE];
char error_message[MAXLINE];
FILE *fperror;
FILE *fpinform;
FILE *fpdebug;
struct cept *head_cept;
struct cept *tail_cept;



long max_type;
long bytes_per_short;
long lword_size;
long mask_size;
long block_size;
long extension_size;
unsigned long current_memory;
unsigned long maximum_memory;
struct class *classes;
long memory_count[MEMORY_TYPES];
long variable_count[VARIABLE_TYPES];

int n_atom;				/* total number of atoms */
int n_file;				/* number of files read */
int n_molecule;			/* number of molecules read */
int atoms;				/* all atoms - set number */
int these;				/* last file read */
int tubers;				/* set of atoms using tubes */
int unknowns;			/* set of atoms with unknown vdw type */
int scratch; 			/* temp set */
int bonds;				/* set of bonds */
int n_radii;
int n_pattern;
long n_residue_index;
struct residueIndex *residue_table;
struct patternRecord *pattern_table;
struct radiusRecord *radius_table;

/* global variables corresponding to dot surface descriptor record: */

/* input fields */
int     maxatom;                /* maximum number of atoms (for future use) */
int     natom;                  /* (current) number of atoms for molecule */
double  *atmco;                  /* pointer to atomic coordinate array */
double  *atmrad;                 /* pointer to atomic radii array */
double  *atmden;                 /* pointer to array of surface densities */
short  *atmatt;                 /* pointer to array of attention numbers */
double   pradius;                /* probe radius */
int     connected;              /* 1 = connected surface; 0 = complete surf */
/* output fields */
int     errflg;                 /* non-zero value means fatal error */
char   *errstr;                 /* brief description of error */
char    stage[80];              /* current stage in calculation */
int     atom;                   /* probable atom involved in error */
int     verbose;                /* true: print out stages as they start */
int     nsrfpnt;                /* total number of surface points for mol */
double   cvxarea;                /* total convex area for molecule */
double   renarea;                /* total reentrant area for molecule */
struct cluster **cvxsp;         /* array of pointers to contact clusters */
struct cluster **rensp;         /* array of pointers to reentrant clusters */
struct hidden  *intern;         /* pointer to record hidden from user */

/* global variables corresponding to hidden record: */

/* counters: */
int     ntori;
int     nprobe;
int     nlow;
int     ndead;
/* pointers to linked lists */
struct torus   *hedtor;
struct torus   *taltor;
struct probe   *hedprb;
struct probe   *talprb;
/* pointers to arrays of pointers */
struct probe  **lowp;
struct run **hedrun;
struct cluster **hedcav;
/* for connected rolling */
struct plist  **hedpl;
int     nplist;
short *problem;
/* for octree algorithm */
double   minwid;
struct box *rootbox;

/* other global variables: */

double    density;       /* current density, for communication with dsdivid() */
int ntab;

struct errmsg errtab[] = {
	500, "memory allocation fails",
	555, "dsd contains null pointer",
	570, "invalid attention number",
	571, "invalid atom radius",
	572, "invalid probe radius",
	573, "invalid surface point density",
	574, "invalid complete/connected parameter",
	575, "invalid maximum number of atoms",
	576, "invalid number of atoms",
	577, "atomic coordinate too large",
	585, "atom centers too close",
	610, "first endpoint not arc begin",
	611, "arc end angle <= arc begin angle",
	612, "null first endpoint",
	613, "null second endpoint",
	622, "confused probe links",
	624, "unexpected null probe link",
	625, "cluster nmem-points inconsistency",
	627, "cluster normals inconsistency",
	630, "cluster unexpectedly without points",
	658, "bad case in switch",
	665, "probe hash key outside limits",
	672, "probe link already used",
	685, "atom list not in ascending order",
	702, "minimum not found",
	712, "non-torus neighbor count inconsistency",
	713, "atom cycle formation unexpected zero count",
	714, "free atom cycle count inconsistency",
	716, "atom cycle formation count inconsistency",
	750, "null pointer returned unexpectedly",
	760, "null pointer as function argument",
	765, "structure contents inconsistent with request",
	770, "null pointer to global list",
	775, "null pointer to input array",
	780, "surface points mistakenly calculated",
	782, "unexpected concave surface points",
	800, "cannot normalize 0 vector",
	810, "arbitrary perpendicular creation fails",
	850, "unit vector norm not 1.0",
	852, "atoms with coincident centers",
	853, "atom unexpectedly outside box",
	855, "initial probe requires > 1 tori",
	859, "negative or zero radius",
	860, "zero angle arc",
	861, "imaginary radius",
	862, "dot product outside expected range",
	863, "arc dots in wrong place",
	864, "zero or negative current density",
	865, "unexpected absence of latitudes",
	866, "probe placed on low attention atoms",
	868, "clockwise concave face",
	869, "distant low probe midplane",
	870, "no mutual neighbors to step onto",
	875, "step-off atom not mutual neighbor",
	876, "step torus calculation error",
	877, "arc does not start at step-off atom",
	880, "inaccessible circle for step",
	882, "arc does not start at step-off probe",
	920, "surface point inside atom",
	925, "probe placed by normal not tangent"
};

/*
The pattern tables was last revised: October 31, 1996.
Last match found is used,
so less-specific early matches will be overriden by
more specific later matches.
A question mark matches any one character.
A string consisting of only an asterisk matches any string.
There is no feature for, C*, for example,
to match any string beginning with C.
You would have to use as series of lines: C, C?, C??, C???, C????.
The residue and atom names can be at most 5 characters.
*/

struct radiusRecord default_radius_table[] = {
    1,  1.60,      0.57,   "O=C",                    
    2,  1.70,      0.66,   "OH",
    3,  1.60,      0.57,   "OCO",
    4,  1.65,      0.70,   "NH",
    5,  1.70,      0.70,   "NH2",
    6,  1.75,      0.70,   "NH3",
    7,  1.85,      0.77,   "CH",
    8,  1.90,      0.77,   "CH2",
    9,  1.95,      0.77,   "CH3",
   10,  1.80,      0.67,   "CAr",
   11,  1.90,      0.70,   "CHAr",
   12,  1.90,      1.04,   "S",
   21,  1.50,      1.50,   "Metal",                   
   31,  2.00,      0.77,   "CAlNu",
   32,  1.77,      0.67,   "CArNu",
   33,  1.40,      0.66,   "OSuNu",
   34,  1.64,      0.57,   "O=CNu",
   35,  1.64,      0.57,   "OPONu",
   36,  1.55,      0.65,   "NArNu",
   37,  1.86,      0.70,   "NAlNu",
   38,  1.80,      0.95,   "P",
   99,  1.00,      0.50,   "H"
};


struct patternRecord default_pattern_table[] ={
"*",    "H",   99,  "H",
"*",    "H?",  99,  "H",
"*",    "H??", 99,  "H",
"*",    "H???", 99, "H",
"*",    "1H", 99,  "H",
"*",    "1H?", 99,  "H",
"*",    "1H??", 99, "H",
"*",    "2H", 99,  "H",
"*",    "2H?", 99,  "H",
"*",    "2H??", 99, "H",
"*",    "3H", 99,  "H",
"*",    "3H?", 99,  "H",
"*",    "3H??", 99, "H",
"*",    "LP",   99,  "H",
"*",    "LP?",  99,  "H",
"*",    "LP??",  99,  "H",
"*",	"C?",	7,	"CH",
"*",	"C??",	7,	"CH",
"*",	"C???",	7,	"CH",
"*",	"Q?",	7,	"CH",
"*",	"Q??",	7,	"CH",
"*",	"N?",	4,	"NH",
"*",	"N??",	4,	"NH",
"*",	"N???",	4,	"NH",
"*",	"O?",	1,	"O=C",
"*",	"O??",  1,	"O=C",
"*",	"O???",  1,	"O=C",
"*",	"S?",	12,	"S",
"*",	"S??",	12,	"S",
"*",	"P?",	38,	"P",
"*",	"P??",	38,	"P",
"*",	"N",	4,	"NH",
"*",	"CA",	7,	"CH",
"*",	"C",	10,	"CAr",
"*",	"O",	1,	"O=C",
"*",	"CB",	8,	"CH2",
"*",	"OXT",	3,	"OCO",
"*",	"OXO",	2,	"OH",
"*",	"OT",	3,	"OCO",
"*",	"CU",	21,	"Metal",
"*",	"ZN",	21,	"Metal",
"*",	"ZN1",	21,	"Metal",
"*",	"ZN2",	21,	"Metal",
"*",	"FE",	21,	"Metal",
"*",	"FE1",	21,	"Metal",
"*",	"FE2",	21,	"Metal",
"*",	"FE3",	21,	"Metal",
"*",	"FE4",	21,	"Metal",
"*",	"MN",	21,	"Metal",
"*",	"NA",	21,	"Metal",
"*",	"CD",	21,	"Metal",
"*",	"MG",	21,	"Metal",
"*",	"S1",	12,	"S",
"*",	"S2",	12,	"S",
"*",	"S3",	12,	"S",
"*",	"S4",	12,	"S",
"*",	"P",	38,	"P",
"*",	"O1P",	35,	"OPONu",
"*",	"O2P",	35,	"OPONu",
"*",	"C1*",	31,	"CAlNu",
"*",	"C2*",	31,	"CAlNu",
"*",	"O2*",	33,	"OSuNu",
"*",	"C3*",	31,	"CAlNu",
"*",	"O3*",	33,	"OSuNu",
"*",	"C4*",	31,	"CAlNu",
"*",	"O4*",	33,	"OSuNu",
"*",	"C5*",	31,	"CAlNu",
"*",	"O5*",	33,	"OSuNu",
"CA",	"CA",	21,	"Metal",
"ACE",	"C",	10,	"CAr",
"ACE",	"O",	1,	"O=C",
"ACE",	"CH3",	9,	"CH3",
"ALA",	"CB",	9,	"CH3",
"ARG",	"CG",	8,	"CH2",
"ARG",	"CD",	8,	"CH2",
"ARG",	"NE",	4,	"NH",
"ARG",	"CZ",	10,	"CAr",
"ARG",	"NH1",	5,	"NH2",
"ARG",	"NH2",	5,	"NH2",
"ASN",	"CG",	7,	"CH",
"ASN",	"OD1",	3,	"OCO",
"ASN",	"ND2",	5,	"NH2",
"ASN",	"AD1",	3,	"OCO",
"ASN",	"AD2",	3,	"OCO",
"ASP",	"CG",	7,	"CH",
"ASP",	"OD1",	3,	"OCO",
"ASP",	"OD2",	3,	"OCO",
"CYS",	"SG",	12,	"S",
"GLU",	"CG",	8,	"CH2",
"GLU",	"CD",	7,	"CH",
"GLU",	"OE1",	3,	"OCO",
"GLU",	"OE2",	3,	"OCO",
"GLN",	"CG",	8,	"CH2",
"GLN",	"CD",	7,	"CH",
"GLN",	"OE1",	3,	"OCO",
"GLN",	"NE2",	5,	"NH2",
"GLN",	"AE1",	3,	"OCO",
"GLN",	"AE2",	3,	"OCO",
"GLN",	"CG_A",	8,	"CH2",
"GLN",	"CD_A",	7,	"CH",
"GLN",	"OE1A",	3,	"OCO",
"GLN",	"NE2A",	5,	"NH2",
"GLN",	"AE1A",	3,	"OCO",
"GLN",	"AE2A",	3,	"OCO",
"GLN",	"CG_B",	8,	"CH2",
"GLN",	"CD_B",	7,	"CH",
"GLN",	"OE1B",	3,	"OCO",
"GLN",	"NE2B",	5,	"NH2",
"GLN",	"AE1B",	3,	"OCO",
"GLN",	"AE2B",	3,	"OCO",
"GLX",	"CG",	8,	"CH2",
"GLX",	"CD",	7,	"CH",
"GLX",	"OE1",	3,	"OCO",
"GLX",	"NE2",	5,	"NH2",
"GLX",	"AE1",	3,	"OCO",
"GLX",	"AE2",	3,	"OCO",
"HIS",	"CG",	7,	"CH",
"HIS",	"ND1",	4,	"NH",
"HIS",	"CD2",	7,	"CH",
"HIS",	"CE1",	7,	"CH",
"HIS",	"NE2",	4,	"NH",
"HIS",	"AD1",	4,	"NH",
"HIS",	"AD2",	4,	"NH",
"HIS",	"AE1",	4,	"NH",
"HIS",	"AE2",	4,	"NH",
"ILE",	"CB",	7,	"CH",
"ILE",	"CG1",	8,	"CH2",
"ILE",	"CG2",	9,	"CH3",
"ILE",	"CD1",	9,	"CH3",
"LEU",	"CG",	7,	"CH",
"LEU",	"CD1",	9,	"CH3",
"LEU",	"CD2",	9,	"CH3",
"LYS",	"CG",	8,	"CH2",
"LYS",	"CD",	8,	"CH2",
"LYS",	"CE",	8,	"CH2",
"LYS",	"NZ",	6,	"NH3",
"MET",	"CG",	8,	"CH2",
"MET",	"SD",	12,	"S",
"MET",	"CE",	9,	"CH3",
"PHE",	"CG",	10,	"CAr",
"PHE",	"CD1",	11,	"CHAr",
"PHE",	"CD2",	11,	"CHAr",
"PHE",	"CE1",	11,	"CHAr",
"PHE",	"CE2",	11,	"CHAr",
"PHE",	"CZ",	11,	"CHAr",
"PRO",	"CG",	8,	"CH2",
"PRO",	"CD",	8,	"CH2",
"SER",	"OG",	2,	"OH",
"THR",	"CB",	7,	"CH",
"THR",	"OG1",	2,	"OH",
"THR",	"CG2",	9,	"CH3",
"TRP",	"CG",	10,	"CAr",
"TRP",	"CD1",	11,	"CHAr",
"TRP",	"CD2",	10,	"CAr",
"TRP",	"CE2",	10,	"CAr",
"TRP",	"NE1",	4,	"NH",
"TRP",	"CE3",	11,	"CHAr",
"TRP",	"CZ2",	11,	"CHAr",
"TRP",	"CZ3",	11,	"CHAr",
"TRP",	"CH2",	11,	"CHAr",
"TYR",	"CG",	10,	"CAr",
"TYR",	"CD1",	11,	"CHAr",
"TYR",	"CD2",	11,	"CHAr",
"TYR",	"CE1",	11,	"CHAr",
"TYR",	"CE2",	11,	"CHAr",
"TYR",	"CZ",	10,	"CAr",
"TYR",	"OH",	2,	"OH",
"VAL",	"CG1",	9,	"CH3",
"VAL",	"CG2",	9,	"CH3",
"VAL",	"CB",	7,	"CH",
"HOH",	"O",	2,	"OH",
"HOH",	"O1",	2,	"OH",
"HOH",	"O2",	2,	"OH",
"HOH",	"O3",	2,	"OH",
"HOH",	"O4",	2,	"OH",
"HOH",	"O7",	2,	"OH",
"HEM",	"FE",	21,	"Metal",
"HEM",	"CHA",	11,	"CHAr",
"HEM",	"CHB",	11,	"CHAr",
"HEM",	"CHC",	11,	"CHAr",
"HEM",	"CHD",	11,	"CHAr",
"HEM",	"N A",	4,	"NH",
"HEM",	"N_A",	4,	"NH",
"HEM",	"NA",	4,	"NH",
"HEM",	"C1A",	11,	"CHAr",
"HEM",	"C2A",	11,	"CHAr",
"HEM",	"C3A",	11,	"CHAr",
"HEM",	"C4A",	11,	"CHAr",
"HEM",	"CMA",	9,	"CH3",
"HEM",	"CAA",	8,	"CH2",
"HEM",	"CBA",	8,	"CH2",
"HEM",	"CGA",	10,	"CAr",
"HEM",	"O1A",	3,	"OCO",
"HEM",	"O2A",	3,	"OCO",
"HEM",	"N B",	4,	"NH",
"HEM",	"N_B",	4,	"NH",
"HEM",	"NB",	4,	"NH",
"HEM",	"C1B",	11,	"CHAr",
"HEM",	"C2B",	11,	"CHAr",
"HEM",	"C3B",	11,	"CHAr",
"HEM",	"C4B",	11,	"CHAr",
"HEM",	"CMB",	9,	"CH3",
"HEM",	"CAB",	8,	"CH2",
"HEM",	"CBB",	8,	"CH2",
"HEM",	"N C",	4,	"NH",
"HEM",	"N_C",	4,	"NH",
"HEM",	"NC",	4,	"NH",
"HEM",	"C1C",	11,	"CHAr",
"HEM",	"C2C",	11,	"CHAr",
"HEM",	"C3C",	11,	"CHAr",
"HEM",	"C4C",	11,	"CHAr",
"HEM",	"CMC",	9,	"CH3",
"HEM",	"CAC",	8,	"CH2",
"HEM",	"CBC",	8,	"CH2",
"HEM",	"N D",	4,	"NH",
"HEM",	"N_D",	4,	"NH",
"HEM",	"ND",	4,	"NH",
"HEM",	"C1D",	11,	"CHAr",
"HEM",	"C2D",	11,	"CHAr",
"HEM",	"C3D",	11,	"CHAr",
"HEM",	"C4D",	11,	"CHAr",
"HEM",	"CMD",	9,	"CH3",
"HEM",	"CAD",	8,	"CH2",
"HEM",	"CBD",	8,	"CH2",
"HEM",	"CGD",	10,	"CAr",
"HEM",	"O1D",	3,	"OCO",
"HEM",	"O2D",	3,	"OCO",
"HEM",	"OH2",	2,	"OH",
"A",	"N1",	36,	"NArNu",
"A",	"C2",	32,	"CArNu",
"A",	"N3",	36,	"NArNu",
"A",	"C4",	32,	"CArNu",
"A",	"C5",	32,	"CArNu",
"A",	"C6",	32,	"CArNu",
"A",	"N6",	37,	"NAlNu",
"A",	"N7",	36,	"NArNu",
"A",	"C8",	32,	"CArNu",
"A",	"N9",	36,	"NArNu",
"1MA",	"N1",	36,	"NArNu",
"1MA",	"C1",	31,	"CAlNu",
"1MA",	"C2",	32,	"CArNu",
"1MA",	"NE",	36,	"NArNu",
"1MA",	"C4",	32,	"CArNu",
"1MA",	"C5",	32,	"CArNu",
"1MA",	"C6",	32,	"CArNu",
"1MA",	"N6",	37,	"NAlNu",
"1MA",	"N7",	36,	"NArNu",
"1MA",	"C8",	32,	"CArNu",
"1MA",	"N9",	36,	"NArNu",
"G",	"N1",	36,	"NArNu",
"G",	"C2",	32,	"CArNu",
"G",	"N2",	37,	"NAlNu",
"G",	"N3",	36,	"NArNu",
"G",	"C4",	32,	"CArNu",
"G",	"C5",	32,	"CArNu",
"G",	"C6",	32,	"CArNu",
"G",	"O6",	34,	"O=CNu",
"G",	"N7",	36,	"NArNu",
"G",	"C8",	32,	"CArNu",
"G",	"N9",	36,	"NArNu",
"2MG",	"N1",	36,	"NArNu",
"2MG",	"C2",	32,	"CArNu",
"2MG",	"N2",	37,	"NAlNu",
"2MG",	"C2A",	31,	"CAlNu",
"2MG",	"N3",	36,	"NArNu",
"2MG",	"C4",	32,	"CArNu",
"2MG",	"C5",	32,	"CArNu",
"2MG",	"C6",	32,	"CArNu",
"2MG",	"O6",	34,	"O=CNu",
"2MG",	"N7",	36,	"NArNu",
"2MG",	"C8",	32,	"CArNu",
"2MG",	"N9",	36,	"NArNu",
"M2G",	"N1",	36,	"NArNu",
"M2G",	"C2",	32,	"CArNu",
"M2G",	"N2",	37,	"NAlNu",
"M2G",	"C2A",	31,	"CAlNu",
"M2G",	"C2B",	31,	"CAlNu",
"M2G",	"N3",	36,	"NArNu",
"M2G",	"C4",	32,	"CArNu",
"M2G",	"C5",	32,	"CArNu",
"M2G",	"C6",	32,	"CArNu",
"M2G",	"O6",	34,	"O=CNu",
"M2G",	"N7",	36,	"NArNu",
"M2G",	"C8",	32,	"CArNu",
"M2G",	"N9",	36,	"NArNu",
"7MG",	"N1",	36,	"NArNu",
"7MG",	"C2",	32,	"CArNu",
"7MG",	"N2",	37,	"NAlNu",
"7MG",	"N3",	36,	"NArNu",
"7MG",	"C4",	32,	"CArNu",
"7MG",	"C5",	32,	"CArNu",
"7MG",	"C6",	32,	"CArNu",
"7MG",	"O6",	34,	"O=CNu",
"7MG",	"N7",	36,	"NArNu",
"7MG",	"C7",	31,	"CAlNu",
"7MG",	"C8",	32,	"CArNu",
"7MG",	"N9",	36,	"NArNu",
"OMG",	"C2A",	31,	"CAlNu",
"OMG",	"N1",	36,	"NArNu",
"OMG",	"C2",	32,	"CArNu",
"OMG",	"N2",	37,	"NAlNu",
"OMG",	"N3",	36,	"NArNu",
"OMG",	"C4",	32,	"CArNu",
"OMG",	"C5",	32,	"CArNu",
"OMG",	"C6",	32,	"CArNu",
"OMG",	"O6",	34,	"O=CNu",
"OMG",	"N7",	36,	"NArNu",
"OMG",	"C8",	32,	"CArNu",
"OMG",	"N9",	36,	"NArNu",
"C",	"N1",	36,	"NArNu",
"C",	"C2",	32,	"CArNu",
"C",	"O2",	34,	"O=CNu",
"C",	"N3",	36,	"NArNu",
"C",	"C4",	32,	"CArNu",
"C",	"N4",	37,	"NAlNu",
"C",	"C5",	32,	"CArNu",
"C",	"C6",	32,	"CArNu",
"OMC",	"C2A",	31,	"CAlNu",
"OMC",	"N1",	36,	"NArNu",
"OMC",	"C2",	32,	"CArNu",
"OMC",	"O2",	34,	"O=CNu",
"OMC",	"N3",	36,	"NArNu",
"OMC",	"C4",	32,	"CArNu",
"OMC",	"N4",	37,	"NAlNu",
"OMC",	"C5",	32,	"CArNu",
"OMC",	"C6",	32,	"CArNu",
"5MC",	"N1",	36,	"NArNu",
"5MC",	"C2",	32,	"CArNu",
"5MC",	"O2",	34,	"O=CNu",
"5MC",	"N3",	36,	"NArNu",
"5MC",	"C4",	32,	"CArNu",
"5MC",	"N4",	37,	"NAlNu",
"5MC",	"C5",	32,	"CArNu",
"5MC",	"C5A",	31,	"CAlNu",
"5MC",	"C6",	32,	"CArNu",
"U",	"N1",	36,	"NArNu",
"U",	"C2",	32,	"CArNu",
"U",	"O2",	34,	"O=CNu",
"U",	"N3",	36,	"NArNu",
"U",	"C4",	32,	"CArNu",
"U",	"OR",	34,	"O=CNu",
"U",	"C5",	32,	"CArNu",
"U",	"C6",	32,	"CArNu",
"T",	"N1",	36,	"NArNu",
"T",	"C2",	32,	"CArNu",
"T",	"O2",	34,	"O=CNu",
"T",	"N3",	36,	"NArNu",
"T",	"C4",	32,	"CArNu",
"T",	"O4",	34,	"O=CNu",
"T",	"C5",	32,	"CArNu",
"T",	"C5M",	31,	"CAlNu",
"T",	"C6",	32,	"CArNu",
"H2U",	"N1",	37,	"NAlNu",
"H2U",	"C2",	31,	"CAlNu",
"H2U",	"O2",	34,	"O=CNu",
"H2U",	"N3",	37,	"NAlNu",
"H2U",	"C4",	31,	"CAlNu",
"H2U",	"O4",	34,	"O=CNu",
"H2U",	"C5",	31,	"CAlNu",
"H2U",	"C6",	31,	"CAlNu",
"5MU",	"N1",	36,	"NArNu",
"5MU",	"C2",	32,	"CArNu",
"5MU",	"O2",	34,	"O=CNu",
"5MU",	"N3",	36,	"NArNu",
"5MU",	"C4",	32,	"CArNu",
"5MU",	"O4",	34,	"O=CNu",
"5MU",	"C5",	32,	"CArNu",
"5MU",	"C5A",	31,	"CAlNu",
"5MU",	"C6",	32,	"CArNu",
"PSU",	"N1",	36,	"NArNu",
"PSU",	"C2",	32,	"CArNu",
"PSU",	"O2",	34,	"O=CNu",
"PSU",	"N3",	36,	"NArNu",
"PSU",	"C4",	32,	"CArNu",
"PSU",	"O4",	34,	"O=CNu",
"PSU",	"C5",	32,	"CArNu",
"PSU",	"C6",	32,	"CArNu",
"YG",	"N1",	36,	"NArNu",
"YG",	"C2",	32,	"CArNu",
"YG",	"N2",	36,	"NArNu",
"YG",	"C3",	31,	"CAlNu",
"YG",	"N3",	36,	"NArNu",
"YG",	"C4",	32,	"CArNu",
"YG",	"C5",	32,	"CArNu",
"YG",	"C6",	32,	"CArNu",
"YG",	"O6",	33,	"OSuNu",
"YG",	"N7",	36,	"NArNu",
"YG",	"C8",	32,	"CArNu",
"YG",	"N9",	36,	"NArNu",
"YG",	"C10",	31,	"CAlNu",
"YG",	"C11",	32,	"CArNu",
"YG",	"C12",	32,	"CArNu",
"YG",	"C13",	31,	"CAlNu",
"YG",	"C14",	31,	"CAlNu",
"YG",	"C15",	31,	"CAlNu",
"YG",	"C16",	31,	"CAlNu",
"YG",	"O17",	33,	"OSuNu",
"YG",	"O18",	33,	"OSuNu",
"YG",	"C19",	31,	"CAlNu",
"YG",	"N20",	37,	"NAlNu",
"YG",	"C21",	31,	"CAlNu",
"YG",	"O22",	33,	"OSuNu",
"YG",	"O23",	33,	"OSuNu",
"YG",	"C24",	31,	"CAlNu",
"PC",	"O4",	35,	"OPONu",
"PC",	"P1",	38,	"P",
"PC",	"O1",	35,	"OPONu",
"PC",	"O3",	35,	"OPONu",
"PC",	"O2",	35,	"OPONu",
"PC",	"C1",	8,	"CH2",
"PC",	"C2",	8,	"CH2",
"PC",	"N1",	4,	"NH",
"PC",	"C3",	9,	"CH3",
"PC",	"C5",	9,	"CH3",
"PC",	"C4",	9,	"CH3",
"AMP",  "O1A",  35, "OPONu",
"AMP",  "O2A",  35, "OPONu",
"AMP",  "O3A",  35, "OPONu",
"ATP",  "O1A",  35, "OPONu",
"ATP",  "O1B",  35, "OPONu",
"ATP",  "O1G",  35, "OPONu",
"ATP",  "O2A",  35, "OPONu",
"ATP",  "O2B",  35, "OPONu",
"ATP",  "O2G",  35, "OPONu",
"ATP",  "O3A",  35, "OPONu",
"ATP",  "O3B",  35, "OPONu",
"ATP",  "O3G",  35, "OPONu",
"DMS",  "S",    12, "S",
"SAD",  "SE1",  12, "S",
"SO4",  "S",    12, "S",
};

void init_patterns()
{
	pattern_table = default_pattern_table;
	n_pattern = sizeof(default_pattern_table) / sizeof(struct patternRecord);
}

void init_radii ()
{
	radius_table = default_radius_table;
	n_radii = sizeof(default_radius_table) / sizeof(struct radiusRecord);
}

int determine_ntab ()
{
    ntab = sizeof (errtab) / sizeof (struct errmsg);
    return (ntab);
}

int init_routing (FILE *e, FILE *i, FILE *d)
{
	if (e != NULL) fperror = e; else fpinform = stderr;
	if (i != NULL) fpinform = i; else fpinform = stderr;
	if (d != NULL) fpdebug = d; else fpdebug = stderr;
	return(1);
}

int get_arguments (int argc, char *argv[], struct argrec array[], int nitems)
{
	char *str;
	double doublearg;
	short flag, ch;
	long longarg, len, slen;
	int idx, a, found, nflag, c, nalpha, ndigit, npunct;
	int isalphanumeric, isinteger, isreal;
	
	if (argc <= 1) return (0);
	nflag = 0; idx = -1;  flag = ' ';
	/* loop through argv command-line arguments */
	for (a = 1; a < argc; a++) {
		str = argv[a]; len = strlen (str);
		if (len > 63) return(-1);
		if (*str == '-' && isalpha(*(str+1))) {
			flag = *(str+1);
			found = 0;
			for (idx = 0; idx < nitems; idx++) {
				if (array[idx].flag == flag) {
					found = 1;
					break;
				}
			}
			if (!found) return (-flag);
			array[idx].on = 1; nflag++;
			continue;
		}
		if (idx < 0) continue;	/* arg not preceded by flag */
		/* determine whether integer, real or string */
		isalphanumeric = 0; isinteger = 0; isreal = 0;
		nalpha = 0; ndigit = 0; npunct = 0;
		for (c = 0; c < len; c++) {
			ch = *(str+c);
			if (isalpha(ch)) nalpha++;
			if (isdigit(ch)) ndigit++;
			if (ispunct(ch)) npunct++;
		}
		/* assuming for now that all numbers are >= 0 */
		if (nalpha > 0)
			isalphanumeric++;
		else if (npunct == 0)
			isinteger++;
		else if (npunct > 0)
			isreal++;
		else return (-1);
		if (isinteger) {
			longarg = atoi(str);
			array[idx].longarg = longarg;
		}
		else if (isreal) {
			doublearg = atof(str);
			array[idx].doublearg = doublearg;
		}
		else if (isalphanumeric) {
			slen = strlen(array[idx].stringarg);
			if (slen == 0)
				strcpy(array[idx].stringarg, str);
			else {
				if (slen + 1 + len > 63) return(-1);
				strcat(array[idx].stringarg, " ");
				strcat(array[idx].stringarg, str);
			}
		}
		else return(-1);
	}
	if (nflag > 0) return(nflag);
	/* look for positional parameters */
	/* loop through argv command-line arguments */
	for (a = 1; a < argc; a++) {
		str = argv[a]; len = strlen (str);
		if (len > 63) return(-1);
		found = 0; idx = -1;
		for (idx = 0; idx < nitems; idx++) {
			if (array[idx].position == a) {
				found = 1;
				break;
			}
		}
		if (!found) return (-1);
		if (idx < 0) return(-1);
		array[idx].on = 1; nflag++;
		/* determine whether integer, real or string */
		isalphanumeric = 0; isinteger = 0; isreal = 0;
		nalpha = 0; ndigit = 0; npunct = 0;
		for (c = 0; c < len; c++) {
			ch = *(str+c);
			if (isalpha(ch)) nalpha++;
			if (isdigit(ch)) ndigit++;
			if (ispunct(ch)) npunct++;
		}
		/* assuming for now that all integers are >= 0 */
		if (nalpha > 0)
			isalphanumeric++;
		else if (npunct == 0)
			isinteger++;
		else if (npunct > 0)
			isreal++;
		else return (-1);
		if (isinteger) {
			longarg = atoi(str);
			array[idx].longarg = longarg;
		}
		else if (isreal) {
			doublearg = atof(str);
			array[idx].doublearg = doublearg;
		}
		else if (isalphanumeric) {
			slen = strlen(array[idx].stringarg);
			if (slen == 0)
				strcpy(array[idx].stringarg, str);
			else {
				if (slen + 1 + len > 63) return(-1);
				strcat(array[idx].stringarg, " ");
				strcat(array[idx].stringarg, str);
			}
		}
		else return(-1);
	}
	return(nflag);
}

void null_terminate (char *line, char chr)
{
	char *c;
	

	for (c = line; *c != (char) 0; c++)
		if (*c == chr) {
			*c = (char) 0;
			break;
		}
}

void setup_default_name (char *atom_file, char name[])
{
	int i, l;
	char *p1, *p2;
	char word[MAXLINE];

	if (atom_file == NULL) {
		strcpy(name,"unknown");
		return;
	}
	i = sscanf(atom_file, "%s", word);
	if (i < 1) {
		strcpy(name,"unknown");
		return;
	}
	p1 = atom_file;
	l = strlen (atom_file);
	p1 += l;
	/* move p1 backward to first / */
	while (p1 > atom_file) {
		if (*p1 == '/' || *p1 == ':' || *p1 == ']') {
			p1++;
			break;
		}
		p1--;
	}
	/* move p1 forward past non-alpha */
	while (p1 < atom_file + l) {
		if (*p1 == 0) break;
		if (isalpha(*p1)) break;
		p1++;
	}
	p2 = p1;
	/* move p2 forward to first non-alphanumeric */
	while (p2 <= atom_file + l) {
		if (*p2 == 0) break;
		if (!isalnum(*p2)) break;
		p2++;
	}
	if ((p2-p1) >= 64) {
		strcpy(name,"unknown");
		return;
	}
	if ((p2-p1) <= 0) {
		strcpy(name,"unknown");
		return;
	}
	/* copy to new array */
	for (i=0; p1 < p2; p1++, i++)
		name[i] = *p1;
	name[i] = (char) 0;
}


struct msscene *new_msscene ()
{
	int j, k;
	struct msscene *ms;

	ms = (struct msscene *) allocate_object (MSSCENE);
	if (ms == NULL) return (ms);
	ms -> fperror = stderr;
	ms -> fpinform = stderr;
	ms -> fpdebug = stderr;
	for (j = 0; j < 3; j++)
		for (k = 0; k < 3; k++)
			ms -> rotation[j][k] = ((j == k) ? 1.0 : 0.0);
	ms -> overlap_hue = 1;
	return (ms);
}


struct component *get_component_ptr (struct surface *srf, int component_number)
{
	struct component *cmp_ptr;

	if (srf -> component_handles == (struct component **) NULL)
		cmp_ptr = (struct component *) NULL;
	else
		cmp_ptr = *(srf -> component_handles + (component_number - 1));
	return (cmp_ptr);
}


int cavity_component (struct surface *srf, int component_number) {
	struct component *cmp_ptr;
	cmp_ptr = get_component_ptr (srf, component_number);
	return (cmp_ptr -> volume < 0.0);
}



struct surface *new_surface ()
{
	struct surface *srf_ptr;
    struct cept *ex;
	
	srf_ptr = (struct surface *) allocate_object (SURFACE);
	if (srf_ptr == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
		add_object (ex,  SURFACE, "srf_ptr");
		add_function (ex, "new_surface");
		return (NULL);
	}
	srf_ptr -> surface_thickness = DEFAULT_THICKNESS;
	return (srf_ptr);
}



void get_format_name (int output_format, char format_name[MAX_NAME])
{
	strcpy (format_name, " ");
	if (output_format == MS_FORMAT)
		strcpy (format_name, "ms");
	else if (output_format == DS_FORMAT)
		strcpy (format_name, "ds");
	else if (output_format == PDB_FORMAT)
		strcpy (format_name, "pdb");
	else if (output_format == SUNRASTER)
		strcpy (format_name, "sun");
	else if (output_format == SGIIMAGE)
		strcpy (format_name, "sgi");
	else if (output_format == WHATIF)
		strcpy (format_name, "whatif");
	else if (output_format == O)
		strcpy (format_name, "O");
	else if (output_format == AVSIMAGE)
		strcpy (format_name, "avsimage");
	else if (output_format == AVSFIELD)
		strcpy (format_name, "avsfield");
	else if (output_format == AVSDENSITY)
		strcpy (format_name, "avsdensity");
	else if (output_format == INVENTOR)
		strcpy (format_name, "inventor");
	else if (output_format == HPGL)
		strcpy (format_name, "hpgl");
	else if (output_format == POSTSCRIPT)
		strcpy (format_name, "ps");
	else if (output_format == VET)
		strcpy (format_name, "vet");
}

int get_format_number(char *format_name)
{
	int output_format;
	if (strcmp (format_name, "xyzr") == 0)
		output_format = XYZR_FORMAT;
	else if (strcmp (format_name, "ms") == 0)
		output_format = MS_FORMAT;
	else if (strcmp (format_name, "pdb") == 0)
		output_format = PDB_FORMAT;
	else if (strcmp (format_name, "ds") == 0)
		output_format = DS_FORMAT;
	else if (strcmp (format_name, "hpgl") == 0)
		output_format = HPGL;
	else if (strcmp (format_name, "ps") == 0)
		output_format = POSTSCRIPT;
	else if (strcmp (format_name, "sun") == 0)
		output_format = SUNRASTER;
	else if (strcmp (format_name, "sgi") == 0)
		output_format = SGIIMAGE;
	else if (strcmp (format_name, "avsimage") == 0)
		output_format = AVSIMAGE;
	else if (strcmp (format_name, "avsfield") == 0)
		output_format = AVSFIELD;
	else if (strcmp (format_name, "avsdensity") == 0)
		output_format = AVSDENSITY;
	else if (strcmp (format_name, "inventor") == 0)
		output_format = INVENTOR;
	else if (strcmp (format_name, "iv") == 0)
		output_format = INVENTOR;
	else if (strcmp (format_name, "inventor") == 0)
		output_format = INVENTOR;
	else if (strcmp (format_name, "whatif") == 0)
		output_format = WHATIF;
	else if (strcmp (format_name, "o") == 0)
		output_format = O;
	else if (strcmp (format_name, "vet") == 0)
		output_format = VET;
	else if (strcmp (format_name, "bmp") == 0)
		output_format = BMP;
	else {
		output_format = 0;
	}
	return(output_format);
}

void free_scene (struct msscene *ms) {
	struct square *squares;
	struct molecule *head_molecule;
	struct surface *head_srf;
	struct material_table *table;
	struct color_ramp *head_ramp;
	struct object_scheme *scheme;
	struct surface *srf, *next_srf;
	struct color_ramp *rmp, *next_rmp;
	struct molecule *mol, *next_mol;

	head_molecule = ms -> head_molecule;
	head_ramp = ms -> head_ramp;
	table = ms -> table;
	squares = ms -> squares;

	for (mol = head_molecule; mol != NULL; mol = next_mol) {
			next_mol = mol -> next;
			head_srf = mol -> head_surface;
			for (srf = head_srf; srf != NULL; srf = next_srf) {
				next_srf = srf -> next;
				scheme = srf -> scheme;
				if (scheme != NULL) {
					if (scheme -> atom_opacities != NULL)
						free_doubles (scheme -> atom_opacities, 0, ATOM_OPACITIES);
					if (scheme -> atom_colors != NULL)
						free_longs (scheme -> atom_colors, 0, ATOM_COLORS);
					free_object (OBJECT_SCHEME, (short *) (srf -> scheme));
					if (error()) return;
				}
				clear_pqms (srf);
				if (error()) return;
				free_surface (srf);
				if (error()) return;
				free_phn (srf);
				if (error()) return;
				free_object (SURFACE, (short *) srf);
				if (error()) return;
			}
			if (mol -> atom_alphas != NULL)
				free_doubles (mol -> atom_alphas, 0, ATOM_ALPHAS);
			free_object (MOLECULE, (short *) mol);
			if (error()) return;
	}
	for (srf = ms -> this_srf; srf != NULL; srf = next_srf) {
		next_srf = srf -> next;
		clear_pqms (srf);
		if (error()) return;
		free_surface (srf);
		if (error()) return;
		free_phn (srf);
		if (error()) return;
		free_object (SURFACE, (short *) srf);
		if (error()) return;
	}
	for (rmp = head_ramp; rmp != NULL; rmp = next_rmp) {
		next_rmp = rmp -> next;
		free_object (COLOR_RAMP, (short *) rmp);
		if (error()) return;
	}
	free_object (MATERIAL_TABLE, (short *) (ms -> table));
	free_object (MSSCENE, (short *) ms);
}



/* Molecular Surface Package Copyright 1989 by Michael L. Connolly */
/* Written by Michael L. Connolly */
