/*
 * Molecular Surface Package
 * Copyright 1986 by Michael L. Connolly
 * All Rights Reserved
 * March 15, 2000
 */

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"


struct msdata *newmsdata ()
{
	struct msdata *msd;
	
	msd = (struct msdata *) allocate_object (MSDATA);
	return (msd);
}

int freemsdata (struct msdata *msd)
{
	int i;
	for (i = 0; i < msd -> narray; i++) {
		if (msd -> floats[i] != NULL)
			free_floats (msd -> floats[i]);
		if (msd -> longs[i] != NULL)
			free_longs (msd -> longs[i], 0, msd -> variable);
		if (msd -> bytes[i] != NULL)
			free_bytes (msd -> bytes[i]);
	}
	free_object (MSDATA, (short *) msd);
	free_cache (MSDATA);
	return (1);
}

int addfloats (struct msdata *msd, long width, long count, float *floats)
{
	int n;
	if (msd -> narray >= MaxMSData) return (0);
	n = msd -> narray;
	msd -> narray++;
	msd -> types[n] = FLOATDATA;
	msd -> widths[n] = width;
	msd -> counts[n] = count;
	msd -> floats[n] = floats;
	return (1);
}

int addlongs (struct msdata *msd, long width, long count, long *longs, int variable)
{
	int n;
	if (msd -> narray >= MaxMSData) return (0);
	n = msd -> narray;
	msd -> narray++;
	msd -> types[n] = LONGDATA;
	msd -> widths[n] = width;
	msd -> counts[n] = count;
	msd -> longs[n] = longs;
	msd -> variable = variable;
	return (1);
}

int addbytes (struct msdata *msd, long width, long count, unsigned char *bytes)
{
	int n;
	if (msd -> narray >= MaxMSData) return (0);
	n = msd -> narray;
	msd -> narray++;
	msd -> types[n] = BYTEDATA;
	msd -> widths[n] = width;
	msd -> counts[n] = count;
	msd -> bytes[n] = bytes;
	return (1);
}



/*
 * Molecular Surface Package
 * Copyright 1986 by Michael L. Connolly
 * All Rights Reserved
 * June 6, 1998
 */

