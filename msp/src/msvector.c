#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

/* Molecular Surface Package Copyright 1993 by Michael L. Connolly */
/* September 4, 1999 */


int write_whatif (struct surface *phn, FILE *fp_tri)
{
	int hue;
	long e, v0, v1, comp, n_component;
	long n_edge_written;
	double hue360, red, green, blue;
	struct phnvtx *pv0, *pv1;
	char message[MAXLINE];
	struct phnedg *pe;

	if (phn -> phnvtx_handles == NULL) return (0);
	if (phn -> phnedg_handles == NULL) return (0);
	if (phn -> phntri_handles == NULL) return (0);
	if (phn -> n_component <= 0) {
		set_error1 ("write_whatif: no components");
		return (0);
	}
	n_component = phn -> n_component;

	/* write edges */
	n_edge_written = 0;
	for (e = 0; e < phn -> n_phnedg; e++) {
		pe = num2phnedg (phn, e + 1);
		if (pe == NULL) return (0);
		hue = 0;
		comp = pe -> comp;
		hue = pe -> hue;
		v0 = pe -> vtxnum[0];
		v1 = pe -> vtxnum[1];
		pv0 = num2phnvtx (phn, v0);
		if (pv0 == NULL) return (0);
		pv1 = num2phnvtx (phn, v1);
		if (pv1 == NULL) return (0);
		red = table -> red[hue];
		green = table -> green[hue];
		blue = table -> blue[hue];
		hue360 = rgb_to_wheel (red, green, blue);
		fprintf (fp_tri,
		" %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",
		pv0 -> center[0], pv0 -> center[1], pv0 -> center[2],
		pv1 -> center[0], pv1 -> center[1], pv1 -> center[2], hue360);
		n_edge_written++;
	}
	sprintf (message,"%8ld edges written to disk in whatif format", n_edge_written);
	inform(message);
	return (1);
}

int write_Os (struct surface *phn, FILE *fpo, char *name)
{
	int icomp;
	char component_name[MAXLINE];
	char message[MAXLINE];
	char shortened[4];

	if (phn -> n_component <= 0) {
		set_error1 ("write_Os: no components");
		return (0);
	}
	/* write files of dots */
	for (icomp = 1; icomp <= phn -> n_component; icomp++) {
		strncpy(shortened,name,3);
		shortened[3] = 0;
		sprintf(component_name, "%s%d", shortened, icomp);
		sprintf(message, "write_Os: %s", component_name);
		inform(message);
		write_O (phn, fpo, icomp, component_name);
	}
	return(1);
}

void write_O (struct surface *phn, FILE *fpo, int icomp, char *component_name)
{
	int k;
	int comp, ccomp, hue;
	long n_dot_written;
	long e;
	double dot_co[3];
	char message[MAXLINE];
	char colour_name[MAXLINE];
	struct phnedg *pe;
	struct phnvtx *pv0, *pv1;
	struct object_scheme *cs;
	int ind;

	/* O format:
	begin s_2
	colour oringe                         
	dot 0.119 -5.280 31.350                                                      
	dot 0.119 -5.280 31.776                                                       
	....                                                                        
	end   
	*/
	cs = phn -> scheme;
	if (cs == NULL) {
		set_error1 ("write_O: null object scheme");
		return;
	}
	hue = icomp;
	fprintf(fpo,"begin %s\n", component_name);

	if (table == NULL) {
		set_error1("write_o: null material table");
		return;
	}
	if (table -> name[hue] == NULL) {
		set_error1("write_o: null material table hue name");
		return;
	}
	strcpy(colour_name, table -> name[hue]);
	fprintf(fpo,"colour %s\n", colour_name);

	/* write midpoints of edges as dots */
	n_dot_written = 0;
	for (e = 0; e < phn -> n_phnedg; e++) {
		pe = num2phnedg (phn, e + 1);
		if (pe == NULL) return;
		comp = pe -> comp;
		if (comp != icomp) continue;
		pv0 = num2phnvtx (phn, pe -> vtxnum[0]);
		if (pv0 == NULL) return;
		pv1 = num2phnvtx (phn, pe -> vtxnum[1]);
		if (pv1 == NULL) return;
		for (k = 0; k < 3; k++)
			dot_co[k] = (pv0 -> center[k] + pv1 -> center[k])/2;
		fprintf (fpo, "dot %8.3f %8.3f %8.3f\n",
			dot_co[0], dot_co[1], dot_co[2]);
		n_dot_written++;
	}
	fprintf(fpo,"end\n");
	sprintf (message,"%8ld dots written for %s",
		n_dot_written, component_name);
	inform(message);
}

