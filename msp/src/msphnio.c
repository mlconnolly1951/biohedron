/* Molecular Surface Package Copyright 1995 Michael L. Connolly */
/* January 8, 2002 */

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"


struct surface *read_polyhedron (FILE *fp_polyhedron)
{
	long i, max_comp;
	struct surface *phn;
	struct phnvtx *pv;
	phn = read_vet (fp_polyhedron);
	/* count components */
	max_comp = 0;
	for (i = 0; i < phn -> n_phnvtx; i++) {
		pv = num2phnvtx (phn, i + 1);
		if (pv == NULL) return (NULL);
		if (pv -> comp > max_comp) max_comp = pv -> comp;
	}
	phn -> n_component = max_comp;
	return (phn);
}

struct surface *read_vet (FILE *fp_polyhedron)
{
	int j, k, nscan, c, nperiod, len, orn0, orn1, orn2;
	short comp, atm, hue;
	long i, nv, ne, nt, en0, en1, en2;
	long minatm, maxatm;
	unsigned long uen0, uen1, uen2, vn0, vn1, vn2;
	double vtxco0[3], vtxco1[3], vtxco2[3];
	double evect0[3], evect1[3], evect2[3];
	double outward[3];
	double elength, tlength, radius;
	char message[MAXLINE];
	char input_line[MAXLINE+1];
	struct surface *phn;
	struct phnvtx *pv, *pv0, *pv1, *pv2;
	struct phnedg *pe, *pe0, *pe1, *pe2;
	struct phntri *tri;
	
	minatm = 1000000;
	maxatm = 0;
	/* read counts of vertices, edges and triangles */

	fgets (input_line, MAXLINE, fp_polyhedron);
	if (feof (fp_polyhedron)) {
		set_error1("(read_vet): premature EOF");
		return(NULL);
	}
	radius = 0.0;
	nscan = sscanf (input_line, "%ld %ld %ld", &nv, &ne, &nt);
	if (nscan < 3) {
		set_error1("(read_vet): bad format for first line");
		return (NULL);
	}
	phn = init_phn (nv, ne, nt);
	if (phn == NULL) return (NULL);

	/* read vertices, one by one */
	for (i = 0; i < phn -> n_phnvtx; i++) {
		fgets (input_line, MAXLINE, fp_polyhedron);
		if (feof (fp_polyhedron)) {
			set_error1 ("(read_vet): premature end of file");
			return(NULL);
		}
		pv = num2phnvtx (phn, i + 1);
		if (pv == NULL) return (NULL);
		len = strlen (input_line);
		nperiod = 0;
		for (c = 0; c < len; c++)
			if (input_line[c] == '.') nperiod++;
		if (nperiod == 7) {
			nscan = sscanf (input_line, "%lf %lf %lf %lf %lf %lf %lf %hd %hd %hd",
				&(pv -> center[0]), &(pv -> center[1]), &(pv -> center[2]),
				&(pv -> outward[0]), &(pv -> outward[1]), &(pv -> outward[2]),
				&(pv -> values[0]), &(pv -> comp), &(pv -> atm), &(pv -> hue));
			if (nscan < 10) {
				set_error1 ("(read_vet): bad vertex format (10)");
				return(NULL);
			}
			if (pv -> values[0] < phn -> minvals[0])
				phn -> minvals[0] = pv -> values[0];
			if (pv -> values[0] > phn -> maxvals[0])
				phn -> maxvals[0] = pv -> values[0];
		}
		else if (nperiod == 9) {
			nscan = sscanf (input_line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %hd %hd %hd",
				&(pv -> center[0]), &(pv -> center[1]), &(pv -> center[2]),
				&(pv -> outward[0]), &(pv -> outward[1]), &(pv -> outward[2]),
				&(pv -> values[0]), &(pv -> values[1]), &(pv -> values[2]), &(pv -> comp), &(pv -> atm), &(pv -> hue));
			if (nscan < 12) {
				set_error1 ("(read_vet): bad vertex format (12)");
				return(NULL);
			}
			for (j = 0; j < 3; j++) {
			if (pv -> values[j] < phn -> minvals[j])
				phn -> minvals[j] = pv -> values[j];
			if (pv -> values[j] > phn -> maxvals[j])
				phn -> maxvals[j] = pv -> values[j];
			}
		}
		else {
			set_error1 ("(read_vet): bad vertex format");
			return(NULL);
		}
		for (k = 0; k < 3; k++) {
			if (pv -> center[k] < phn -> bounds[0][k])
				phn -> bounds[0][k] = pv -> center[k];
			if (pv -> center[k] > phn -> bounds[1][k])
				phn -> bounds[1][k] = pv -> center[k];
		}
		pv -> number = i + 1;
	}
	if (phn -> minvals[0] < phn -> maxvals[0]) {
		sprintf (message, "%8.3f minimum u value", phn -> minvals[0]);
		inform(message);
		sprintf (message, "%8.3f maximum u value", phn -> maxvals[0]);
		inform(message);
	}
	if (phn -> minvals[1] < phn -> maxvals[1]) {
		sprintf (message, "%8.3f minimum v value", phn -> minvals[1]);
		inform(message);
		sprintf (message, "%8.3f maximum v value", phn -> maxvals[1]);
		inform(message);
	}
	if (phn -> minvals[2] < phn -> maxvals[2]) {
		sprintf (message, "%8.3f minimum w value", phn -> minvals[2]);
		inform(message);
		sprintf (message, "%8.3f maximum w value", phn -> maxvals[2]);
		inform(message);
	}
	for (k = 0; k < 3; k++)
		phn -> center[k] = (phn -> bounds[0][k] + phn -> bounds[1][k]) / 2.0;

	/* read edges, one by one */
	for (i = 0; i < phn -> n_phnedg; i++) {
		fgets (input_line, MAXLINE, fp_polyhedron);
		if (feof (fp_polyhedron)) {
			set_error1 ("(read_vet): premature end of file");
			return(NULL);
		}
		pe = num2phnedg(phn, i + 1);
		if (pe == NULL) return (NULL);
		nscan = sscanf (input_line, "%ld %ld %hd %hd %hd",
			&vn0, &vn1, &comp, &atm, &hue);
		if (nscan < 5) {
			set_error1 ("(read_vet): bad edge format");
			return(NULL);
		}
		if (vn0 < 1 || vn0 > phn -> n_phnvtx) {
			sprintf (message, "(read_vet): bad vertex number: %6ld", vn0);
			set_error1(message);
			return(NULL);
		}
		if (vn1 < 1 || vn1 > phn -> n_phnvtx) {
			sprintf (message, "(read_vet): bad vertex number: %6ld", vn1);
			set_error1(message);
			return(NULL);
		}
		pv0 = num2phnvtx (phn, vn0);
		if (pv0 == NULL) return (NULL);
		pv1 = num2phnvtx (phn, vn1);
		if (pv1 == NULL) return (NULL);
		/* check for degeneracy */
		for (k = 0; k < 3; k++) {
			vtxco0[k] = pv0 -> center[k];
			vtxco1[k] = pv1 -> center[k];
		}
		elength = distance (vtxco0, vtxco1);
		if (elength <= 0.0) {
			sprintf (message,
					"(read_vet): warning, degenerate edge, vertices: %6ld %6ld", vn0+1, vn1+1);
			inform(message);
			/* kludge: add a bit to one cooridinate of one vertex */
			pv0 -> center[vn0 % 3] += 0.000001 * (1 + vn0 % 7);
		}
		pe -> comp = comp;
		pe -> atm = atm;
		pe -> hue = hue;
		pe -> vtxnum[0] = vn0;
		pe -> vtxnum[1] = vn1;
		pe -> pvt[0] = pv0;
		pe -> pvt[1] = pv1;
		pe -> number = i + 1;
	}

	/* read triangles, one by one */
	for (i = 0; i < phn -> n_phntri; i++) {
		fgets (input_line, MAXLINE, fp_polyhedron);
		if (feof (fp_polyhedron)) {
			set_error1 ("(read_vet): premature end of file");
			return (NULL);
		}
		tri = num2phntri (phn, i + 1);
		if (tri == NULL) return (NULL);
		nscan = sscanf (input_line, "%ld %ld %ld %ld %ld %ld %hd %hd %hd",
			&en0, &en1, &en2, &vn0, &vn1, &vn2, &comp, &atm, &hue);
		if (nscan < 9) {
			set_error1 ("(read_vet): bad triangle format");
			return (NULL);
		}
		uen0 = abs (en0);
		uen1 = abs (en1);
		uen2 = abs (en2);
		orn0 = (en0 < 0);
		orn1 = (en1 < 0);
		orn2 = (en2 < 0);
		if (uen0 < 1 || uen0 > phn -> n_phnedg) {
			set_error1 ("(read_vet): bad edge number");
			return (NULL);
		}
		if (uen1 < 1 || uen1 > phn -> n_phnedg) {
			set_error1 ("(read_vet): bad edge number");
			return (NULL);
		}
		if (uen2 < 1 || uen2 > phn -> n_phnedg) {
			set_error1 ("(read_vet): bad edge number");
			return (NULL);
		}
		if (vn0 < 1 || vn0 > phn -> n_phnvtx) {
			sprintf (message, "(read_vet): bad vertex number: %6ld", vn0);
			set_error1(message);
			return (NULL);
		}
		if (vn1 < 1 || vn1 > phn -> n_phnvtx) {
			sprintf (message, "(read_vet): bad vertex number: %6ld", vn1);
			set_error1(message);
			return (NULL);
		}
		if (vn2 < 1 || vn2 > phn -> n_phnvtx) {
			sprintf (message, "(read_vet): bad vertex number: %6ld", vn2);
			set_error1(message);
			return (NULL);
		}
		pe0 = num2phnedg(phn, uen0);
		if (pe0 == NULL) return (NULL);
		pe1 = num2phnedg(phn, uen1);
		if (pe1 == NULL) return (NULL);
		pe2 = num2phnedg(phn, uen2);
		if (pe2 == NULL) return (NULL);
		tri -> edg[0] = pe0;
		tri -> edg[1] = pe1;
		tri -> edg[2] = pe2;
		tri -> orn[0] = (short) orn0;
		tri -> orn[1] = (short) orn1;
		tri -> orn[2] = (short) orn2;
		tri -> comp = comp;
		tri -> atm = atm;
		tri -> hue = hue;
		tri -> shape = 0;
		tri -> edgnum[0] = en0;
		tri -> edgnum[1] = en1;
		tri -> edgnum[2] = en2;
		tri -> vtxnum[0] = vn0;
		tri -> vtxnum[1] = vn1;
		tri -> vtxnum[2] = vn2;
		if (atm < minatm) minatm = atm;
		if (atm > maxatm) maxatm = atm;

		pv0 = num2phnvtx (phn, vn0);
		if (pv0 == NULL) return (NULL);
		pv1 = num2phnvtx (phn, vn1);
		if (pv1 == NULL) return (NULL);
		pv2 = num2phnvtx (phn, vn2);
		if (pv2 == NULL) return (NULL);
		/* compute triangle normal vector */
		for (k = 0; k < 3; k++) {
			vtxco0[k] = pv0 -> center[k];
			vtxco1[k] = pv1 -> center[k];
			vtxco2[k] = pv2 -> center[k];
		}
		for (k = 0; k < 3; k++) {
			evect0[k] = vtxco1[k] - vtxco0[k];
			evect1[k] = vtxco2[k] - vtxco1[k];
			evect2[k] = vtxco0[k] - vtxco2[k];
		}
		cross (evect1, evect2, outward);
		tlength = norm (outward);
		if (tlength <= 0.0) {
			sprintf (message,
					"(read_vet): warning, degenerate triangle, vertices: %6ld %6ld %6ld", 
					vn0+1, vn1+1, vn2+1);
			inform(message);
		}
		else normalize (outward);
		for (k = 0; k < 3; k++)
			tri -> axis[k] = outward[k];
	}
	if (maxatm > minatm) phn -> n_atom = maxatm;

	link_polyhedron (phn);
	if (error()) return (NULL);

	return (phn);
}



struct surface *read_outer (FILE *fp_polyhedron)
{
	struct surface *oldphn;
	struct surface *newphn;

	oldphn = read_vet (fp_polyhedron);
	if (oldphn -> n_component >= 2) {
		newphn = filter_phn (oldphn);
		free_phn (oldphn);
		if (newphn -> next != NULL) {
			free_phn (newphn -> next);
			newphn -> next = NULL;
		}
	}
	else newphn = oldphn;
	return (newphn);
}

int write_vet (struct surface *phn, FILE *fp_tri)
{
	int comp, atm, hue;
	long v, e, t;
	char message[MAXLINE];
	struct phnvtx *pv;
	struct phnedg *pe;
	struct phntri *pt;

	if (phn -> phnvtx_handles == NULL) return (0);
	if (phn -> phnedg_handles == NULL) return (0);
	if (phn -> phntri_handles == NULL) return (0);
	fprintf (fp_tri, "%ld %ld %ld\n", phn -> n_phnvtx, phn -> n_phnedg, phn -> n_phntri);

	/* write vertex coordinates */
	for (v = 0; v < phn -> n_phnvtx; v++) {
		pv = num2phnvtx (phn, v + 1);
		if (pv == NULL) return (0);
		comp = pv -> comp;
		atm = pv -> atm;
		hue = pv -> hue;
		fprintf (fp_tri, "%12.6f %12.6f %12.6f %7.4f %7.4f %7.4f %10.6f %10.6f %10.6f %3d %5d %3d\n",
			pv -> center[0], pv -> center[1], pv -> center[2],
			pv -> outward[0], pv -> outward[1], pv -> outward[2],
			pv -> values[0], pv -> values[1], pv -> values[2], comp, atm, hue);
		/*
                fprintf (fp_tri,
			"%6d %3d %10.6f %10.6f %10.6f %12.6f %12.6f %12.6f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
			atm, comp, pv -> values[0], pv -> values[1], pv -> values[2],
			pv -> center[0], pv -> center[1], pv -> center[2],
			pv -> outward[0], pv -> outward[1], pv -> outward[2],
			pv -> base[0], pv -> base[1], pv -> base[2],
			pv -> zenith[0], pv -> zenith[1], pv -> zenith[2]);
                */
	}
	sprintf (message,"%8ld vertices written to disk", phn -> n_phnvtx);
	inform(message);

	/* write edges */
	for (e = 0; e < phn -> n_phnedg; e++) {
		pe = num2phnedg (phn, e + 1);
		if (pe == NULL) return (0);
		comp = pe -> comp;
		atm = pe -> atm;
		hue = pe -> hue;
		fprintf (fp_tri, "%6ld %6ld %3d %5d %3d\n",
			pe -> vtxnum[0], pe -> vtxnum[1], comp, atm, hue);
	}
	sprintf (message,"%8ld edges written to disk", phn -> n_phnedg);
	inform(message);

	for (t = 0; t < phn -> n_phntri; t++) {
		pt = num2phntri (phn, t + 1);
		if (pt == NULL) return (0);
		comp = pt -> comp;
		atm = pt -> atm;
		hue = pt -> hue;
		fprintf (fp_tri, "%7ld %7ld %7ld %6ld %6ld %6ld %3d %5d %3d\n",
			pt -> edgnum[0], pt -> edgnum[1], pt -> edgnum[2],
			pt -> vtxnum[0], pt -> vtxnum[1], pt -> vtxnum[2],
			comp, atm, hue);
	}
	sprintf (message,"%8ld triangles written to disk", phn -> n_phntri);
	inform(message);
	return (1);
}


struct surface *filter_phn (struct surface *oldphn)
{
	int j, k;
	long v, e, t;
	long vold, eold, told;
	long vnew0, enew0, tnew0;
	long vnew1, enew1, tnew1;
	long oldnv, oldne, oldnt;
	long newnv0, newne0, newnt0;
	long newnv1, newne1, newnt1;
	long *vnumbers0, *enumbers0;
	long *vnumbers1, *enumbers1;
	struct surface *newphn0, *newphn1;
	struct phnvtx *opv, *npv0, *npv1;
	struct phnedg *ope, *npe0, *npe1;
	struct phntri *opt, *npt0, *npt1;
	char message[MAXLINE];

	if (oldphn -> n_component <= 1) {
		set_error1 ("filter_phn: too few components");
		return (NULL);
	}

	oldnv = oldphn -> n_phnvtx;
	oldne = oldphn -> n_phnedg;
	oldnt = oldphn -> n_phntri;
	newnv0 = 0; newne0 = 0; newnt0 = 0;
	newnv1 = 0; newne1 = 0; newnt1 = 0;
	for (v = 0; v < oldnv; v++) {
		opv = num2phnvtx (oldphn, v + 1);
		if (opv == NULL) {
			set_error1 ("filter_phn: invalid old vertex number");
			return (NULL);
		}
		if (!cavity_component (oldphn, (int) opv -> comp)) newnv0++; else newnv1++;
	}
	for (e = 0; e < oldne; e++) {
		ope = num2phnedg (oldphn, e + 1);
		if (ope == NULL) {
			set_error1 ("filter_phn: invalid old edge number");
			return (NULL);
		}
		if (!cavity_component (oldphn, (int) ope -> comp)) newne0++; else newne1++;
	}
	for (t = 0; t < oldnt; t++) {
		opt = num2phntri (oldphn, t + 1);
		if (opt == NULL) {
			set_error1 ("filter_phn: invalid old triangle number");
			return (NULL);
		}
		if (!cavity_component (oldphn, (int) opt -> comp)) newnt0++; else newnt1++;
	}
	if (newnv0 <= 0 || newne0 <= 0 || newnt0 <= 0) {
		sprintf (message, "filter_phn: missing outer component");
		set_error1 (message);
		return (NULL);
	}
	if (newnv1 <= 0 || newne1 <= 0 || newnt1 <= 0) {
		sprintf (message, "filter_phn: missing cavity components");
		set_error1 (message);
		return (NULL);
	}
	newphn0 = init_phn (newnv0, newne0, newnt0);
	if (newphn0 == NULL) return (NULL);
	newphn1 = init_phn (newnv1, newne1, newnt1);
	if (newphn1 == NULL) return (NULL);

	/* set up number translation tables */
	vnumbers0 = allocate_longs (oldnv, 0, VNUMBERS0);
	if (vnumbers0 == NULL) {
		sprintf (message, "filter_phn: memory allocation failure for new vertices 0");
		return (NULL);
	}
	vnumbers1 = allocate_longs (oldnv, 0, VNUMBERS1);
	if (vnumbers1 == NULL) {
		sprintf (message, "filter_phn: memory allocation failure for new vertices 1");
		return (NULL);
	}
	enumbers0 = allocate_longs (oldne, 0, ENUMBERS0);
	if (enumbers0 == NULL) {
		sprintf (message, "filter_phn: memory allocation failure for new edges 0");
		return (NULL);
	}
	enumbers1 = allocate_longs (oldne, 0, ENUMBERS1);
	if (enumbers1 == NULL) {
		sprintf (message, "filter_phn: memory allocation failure for new edges 1");
		return (NULL);
	}

	vnew0 = 1; vnew1 = 1;
	for (v = 0; v < oldnv; v++) {
		opv = num2phnvtx (oldphn, v + 1);
		if (!cavity_component (oldphn, (int) opv -> comp)) {
			*(vnumbers0 + v) = vnew0;
			vnew0++;
		}
		else {
			*(vnumbers1 + v) = vnew1;
			vnew1++;
		}
	}
	vnew0--; vnew1--;
	enew0 = 1; enew1 = 1;
	for (e = 0; e < oldne; e++) {
		ope = num2phnedg (oldphn, e +1);
		if (ope == NULL) {
			set_error1 ("filter_phn: invalid old edge number");
			return (NULL);
		}
		if (!cavity_component (oldphn, (int) ope -> comp)) {
			*(enumbers0 + e) = enew0;
			enew0++;
		}
		else {
			*(enumbers1 + e) = enew1;
			enew1++;
		}
	}
	enew0--; enew1--;

	/* transfer information, but omitting all but the selected component */
	for (k = 0; k < 3; k++) {
		newphn0 -> bounds[0][k] = oldphn -> bounds[0][k];
		newphn0 -> bounds[1][k] = oldphn -> bounds[1][k];
		newphn0 -> center[k] = oldphn -> center[k];
	}
	for (j = 0; j < 3; j++) {
		newphn0 -> minvals[k] = 1000000.0;
		newphn0 -> maxvals[k] = -1000000.0;
	}
	for (k = 0; k < 3; k++) {
		newphn1 -> bounds[0][k] = oldphn -> bounds[0][k];
		newphn1 -> bounds[1][k] = oldphn -> bounds[1][k];
		newphn1 -> center[k] = oldphn -> center[k];
	}
	for (j = 0; j < 3; j++) {
		newphn1 -> minvals[k] = 1000000.0;
		newphn1 -> maxvals[k] = -1000000.0;
	}
	vnew0 = 0; vnew1 = 0;
	for (v = 0; v < oldnv; v++) {
		opv = num2phnvtx (oldphn, v + 1);
		if (opv == NULL) {
			set_error1 ("filter_phn: invalid vertex number");
			return (NULL);
		}
		if (!cavity_component (oldphn, (int) opv -> comp)) {
			npv0 = num2phnvtx (newphn0, vnew0 + 1);
			if (npv0 == NULL) {
				set_error1 ("filter_phn: invalid vertex number");
				return (NULL);
			}
			for (k = 0; k < 3; k++) {
				npv0 -> center[k] = opv -> center[k];
				npv0 -> outward[k] = opv -> outward[k];
				npv0 -> normal[k] = opv -> normal[k];
				npv0 -> base[k] = opv -> base[k];
				npv0 -> zenith[k] = opv -> zenith[k];
			}
			for (j = 0; j < 3; j++) {
				npv0 -> values[j] = opv -> values[j];
				if (npv0 -> values[j] < newphn0 -> minvals[j])
					newphn0 -> minvals[j] = npv0 -> values[j];
				if (npv0 -> values[j] > newphn0 -> maxvals[j])
					newphn0 -> maxvals[j] = npv0 -> values[j];
			}
			npv0 -> comp = opv -> comp;
			npv0 -> atm = opv -> atm;
			npv0 -> hue = opv -> hue;
			vnew0++;
			npv0 -> number = vnew0;
		}
		else {
			npv1 = num2phnvtx (newphn1, vnew1 + 1);
			if (npv1 == NULL) {
				set_error1 ("filter_phn: invalid vertex number");
				return (NULL);
			}
			for (k = 0; k < 3; k++) {
				npv1 -> center[k] = opv -> center[k];
				npv1 -> outward[k] = opv -> outward[k];
				npv1 -> normal[k] = opv -> normal[k];
				npv1 -> base[k] = opv -> base[k];
				npv1 -> zenith[k] = opv -> zenith[k];
			}
			for (j = 0; j < 3; j++) {
				npv1 -> values[j] = opv -> values[j];
				if (npv1 -> values[j] < newphn1 -> minvals[j])
					newphn1 -> minvals[j] = npv1 -> values[j];
				if (npv1 -> values[j] > newphn1 -> maxvals[j])
					newphn1 -> maxvals[j] = npv1 -> values[j];
			}
			npv1 -> comp = opv -> comp;
			npv1 -> atm = opv -> atm;
			npv1 -> hue = opv -> hue;
			vnew1++;
			npv1 -> number = vnew1;
		}
	}
	enew0 = 0; enew1 = 0;
	for (e = 0; e < oldne; e++) {
		ope = num2phnedg (oldphn, e + 1);
		if (ope == NULL) {
			set_error1 ("filter_phn: invalid edge number");
			return (NULL);
		}

		if (!cavity_component (oldphn, (int) ope -> comp)) {
			npe0 = num2phnedg (newphn0, enew0 + 1);
			if (npe0 == NULL) {
				set_error1 ("filter_phn: invalid edge number");
				return (NULL);
			}
			for (j = 0; j < 2; j++) {
				vold = ope -> vtxnum[j];
				vnew0 = *(vnumbers0 + (vold - 1));
				npe0 -> vtxnum[j] = vnew0;
			}
			npe0 -> comp = ope -> comp;
			npe0 -> atm = ope -> atm;
			npe0 -> hue = ope -> hue;
			enew0++;
		}
		else {
			npe1 = num2phnedg (newphn1, enew1 + 1);
			if (npe1 == NULL) {
				set_error1 ("filter_phn: invalid edge number");
				return (NULL);
			}
			for (j = 0; j < 2; j++) {
				vold = ope -> vtxnum[j];
				vnew1 = *(vnumbers1 + (vold - 1));
				npe1 -> vtxnum[j] = vnew1;
			}
			npe1 -> comp = ope -> comp;
			npe1 -> atm = ope -> atm;
			npe1 -> hue = ope -> hue;
			enew1++;
		}
	}
	tnew0 = 0; tnew1 = 0;
	for (t = 0; t < oldnt; t++) {
		opt = num2phntri (oldphn, t + 1);
		if (opt == NULL) {
			set_error1 ("filter_phn: invalid triangle number");
			return (NULL);
		}
		if (!cavity_component (oldphn, (int) opt -> comp)) {
			npt0 = num2phntri (newphn0, tnew0 + 1);
			if (npt0 == NULL) {
				set_error1 ("filter_phn: invalid triangle number");
				return (NULL);
			}
			for (j = 0; j < 3; j++) {
				eold = opt -> edgnum[j];
				if (eold > 0) enew0 = *(enumbers0 + (eold - 1));
				else {
					eold = (-eold);
					enew0 = *(enumbers0 + (eold - 1));
					enew0 = (-enew0);
				}
				npt0 -> edgnum[j] = enew0;
			}
			for (j = 0; j < 3; j++) {
				vold = opt -> vtxnum[j];
				vnew0 = *(vnumbers0 + (vold - 1));
				npt0 -> vtxnum[j] = vnew0;
			}
			npt0 -> comp = opt -> comp;
			npt0 -> atm = opt -> atm;
			npt0 -> hue = opt -> hue;
			npt0 -> shape = opt -> shape;
			tnew0++;
		}
		else {
			npt1 = num2phntri (newphn1, tnew1 + 1);
			if (npt1 == NULL) {
				set_error1 ("filter_phn: invalid triangle number");
				return (NULL);
			}
			for (j = 0; j < 3; j++) {
				eold = opt -> edgnum[j];
				if (eold > 0) enew1 = *(enumbers1 + (eold - 1));
				else {
					eold = (-eold);
					enew1 = *(enumbers1 + (eold - 1));
					enew1 = (-enew1);
				}
				npt1 -> edgnum[j] = enew1;
			}
			for (j = 0; j < 3; j++) {
				vold = opt -> vtxnum[j];
				vnew0 = *(vnumbers1 + (vold - 1));
				npt1 -> vtxnum[j] = vnew1;
			}
			npt1 -> comp = opt -> comp;
			npt1 -> atm = opt -> atm;
			npt1 -> hue = opt -> hue;
			npt1 -> shape = opt -> shape;
			tnew1++;
		}
	}
	free_longs (vnumbers0, 0, VNUMBERS0);
	free_longs (vnumbers1, 0, VNUMBERS1);
	free_longs (enumbers0, 0, ENUMBERS0);
	free_longs (enumbers1, 0, ENUMBERS1);
	sprintf (message, "%8ld vertices, %8ld edges, %8ld triangles first part",
		newnv0, newne0, newnt0);
	inform (message);
	sprintf (message, "%8ld vertices, %8ld edges, %8ld triangles second part",
		newnv1, newne1, newnt1);
	inform (message);
	newphn0 -> next = newphn1;
	return (newphn0);
}



