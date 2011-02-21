/*
 Authors 
 Martin Schlather, martin.schlather@math.uni-goettingen.de

 Definition of correlation functions and derivatives (spectral measures, 
 tbm operators)

Note:
 * Never use the below functions directly, but only by the functions indicated 
   in RFsimu.h, since there is no error check (e.g. initialization of RANDOM)
 * VARIANCE, SCALE are not used here 
 * definitions for the random coin method can be found in MPPFcts.cc
 * definitions for genuinely anisotropic or nonstationary models are in
   SophisticatedModel.cc; hyper models also in Hypermodel.cc


 Copyright (C) 2001 -- 2003 Martin Schlather
 Copyright (C) 2004 -- 2004 Yindeng Jiang & Martin Schlather
 Copyright (C) 2005 -- 2011 Martin Schlather

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  
*/


#include <math.h>
#include <assert.h>
#include "RF.h"

#ifdef DoublePrecision
#define boolean int
#define uchar int
#define uword int
#else
#define boolean bool
#define uchar unsigned char
#define uword short unsigned int
#endif
//#define gz unsigned short int


void threshold(int n, res_type *in, double percent, boolean* out) {
    int i;
    double threshold;
    
    threshold = qnorm(percent, 0.0, 1.0, 1, 0);
    for (i=0; i<n; i++) {
	out[i] = in[i] < threshold;
    }
}

void coarsening(int nr, int nc, boolean *in, int res, uchar *out) {
    if (res > 15 && sizeof(uchar)==1) error("max coarse is 15");
    
    int threshold = (res * res) / 2,
	ncDres = nc / res,
	nrres = nr * res,
	nrDres = nr / res;
    boolean *pin, *pc,  *pinbase, *lastres,
	*lastin = in + ncDres * nrres;
    uchar *proutende, *pr,
	*pout = out, 
	*lastout = out + nrDres * ncDres;

    
    for(pinbase = in; pinbase<lastin; pout += nrDres) {
	pc = pinbase + nrres;

	for (; pinbase<pc;  pinbase+=nr) {
	    pin = pinbase;
	    proutende = pout + nrDres;
	    for (pr = pout; pr<proutende; pr++) {	    
		for (lastres = pin + res; pin<lastres; pin++) {
		    *pr += *pin;  // *pr = *pr + *pin;
		}
	    }
	}
    }

    for (pout = out; pout<lastout; pout++) {
	*pout = *pout > threshold;	
    }
    
}


void reference_area(int nr, int nc, uchar *in, int edge, uword *out) {
    int i, j, old, laststart,
	n = nc * nr,
	edgeM1 = edge - 1,
	edgeM1nr = nr * edgeM1,
	last = n - edge,
	edgenr = edge * nr,
	nrMedge1 = nr - edge +1;
    uword dummy;
    
    out[0] = 0.0;
    for(i=0; i<edge; i++) out[0] += in[i];
    for(i=1; i<=last; i++) {
	out[i] = out[i-1] - in[i-1] + in[i + edgeM1];
    }


    last = n - 1 - edgeM1nr;
    for(i=0; i<nrMedge1; i++) {
	old = out[i];
	for(j=i + nr, laststart = i + edgenr; j<laststart; j+=nr) {
	    out[i] += out[j];
	}
	for(j=i + nr; j<=last; j+=nr) {
	    dummy = out[j];
	    out[j] = out[j - nr] - old + out[j + edgeM1nr];
	    old = dummy;
	}
    }	   
}


void analyse_internal(int* keynr, 
		      double *gausspercent, int *ngaussp,
		      int *coarse, int *ncrs, 
		      int *edges, int *nedg,
		      double *percent, int  *nperc,
		      res_type *image,
		      boolean *binary,
		      uchar * decreased,
		      int *Nrdecr,
		      int *Ncdecr,
		      uword *refarea,
		      int *Areathreshold,
		      double *result
    ) {
    int k, Err, nettoedge, igaussp, ires, iedg, iperc, 
	lastnc, forest, land, j, i, lastnr,
	nrdecr=-1, 
	ncdecr=-1, 
        areathreshold=-1, 
	n = KEY[*keynr].loc.totalpoints,
	nr = KEY[*keynr].loc.length[0],
	nc = KEY[*keynr].loc.length[1],	
	one = 1,
	zero = 0;
    double percent_covered;

    bool debug = false;

    DoSimulateRF(keynr, &one, &zero, image, &Err);
    if (Err != 0) {
	error("DoSimulateRF does not work");
    }
    if (debug) PRINTF("simulation done with %dx%d pixels\n", nr, nc);

    k = 0;
    for (igaussp = 0; igaussp < *ngaussp; igaussp++) {
	threshold(n, image, gausspercent[igaussp], binary);
	for (ires = 0; ires<*ncrs;  ires++) {
	    nrdecr = nr / coarse[ires];
	    ncdecr = nc / coarse[ires];
	    if (debug) PRINTF("thresh=%2.2f coarse=%d new size:%dx%d pixels\n", 
		   gausspercent[igaussp], coarse[ires],
		   nrdecr, ncdecr);
	    if (coarse[ires] == 1) {
	        for (i=0; i<n; i++) decreased[i] = binary[i];
	    } else {
	        for (i=0; i<n; i++) decreased[i] = 0;
	        coarsening(nr, nc, binary, coarse[ires], decreased);
	    }
	    for (iedg = 0; iedg<*nedg; iedg++) {
		nettoedge = (int) round((double) edges[iedg] / 
					(double) coarse[ires]);
		if (nettoedge > nrdecr) nettoedge = nrdecr;
		if (nettoedge > ncdecr) nettoedge = ncdecr;
		if (debug) PRINTF("   netto.referenz=%d\n", nettoedge);
		for (iperc=0; iperc<*nperc; iperc++) {
		    areathreshold = 
			(int) ceil(nettoedge * nettoedge * percent[iperc]);
		    if (debug) PRINTF("       min=%2.0f%(=%dpx) ",
			   percent[iperc] * 100, areathreshold 
			);
		    reference_area(nrdecr, ncdecr, decreased, nettoedge, 
				   refarea);
		    lastnc = (ncdecr - nettoedge +1) * nrdecr;	    
		    forest = 0;
		    land = (nrdecr - nettoedge + 1) * (ncdecr - nettoedge + 1);
		    for (j=0; j<lastnc; j+=nrdecr) {
			lastnr = j + nrdecr - nettoedge;		    
			for (i=j; i<=lastnr; i++) {
			    forest += (refarea[i] >= areathreshold);
			}
		    }

		    percent_covered = (double) forest / (double) land;
		    if (debug) PRINTF("-> forest=%d land=(%d-%d+1)^2=%d coverd=%2.2f%\n",
			   forest, nrdecr, nettoedge, land, percent_covered);
		    result[k++] = percent_covered;
 		}
	    }
	}
    }
    *Nrdecr = nrdecr;
    *Ncdecr = ncdecr;
    *Areathreshold = areathreshold;
}



void analyseForst(int* keynr, 
		  double *gausspercent, int *ngaussp,
		  int *coarse, int *ncrs, 
		  int *edges, int *nedg,
		  double *percent, int *nperc,
		  double *result) {

    int nrdecr, ncdecr, areathreshold,
	n = KEY[*keynr].loc.totalpoints;
    res_type *image;
    boolean *binary;
    uchar *decreased;
    uword *refarea;

    PRINTF("memory needed %e + %e + %e + %e\n", 
	   (double) n * sizeof(res_type), 
	   (double) n * sizeof(boolean),
	   (double) n * sizeof(uchar),
	   (double) n * sizeof(uword)
	);
    if ((image = (res_type*) malloc(n * sizeof(res_type))) == NULL) 
	error("mem err 1");
    if ((binary = (boolean*) malloc(n * sizeof(boolean))) == NULL) 
	error("mem err 2");
    if ((decreased = (uchar*) calloc(n, sizeof(uchar))) == NULL) 
	error("mem err 3");
    if ((refarea   = (uword*) malloc(n * sizeof(uword))) == NULL)
	error("mem err 4");

    analyse_internal(keynr, gausspercent, ngaussp, coarse, ncrs, 
		     edges, nedg, percent, nperc, 
		     image, binary, decreased, 
		     &nrdecr, &ncdecr, refarea, 
		     &areathreshold,
		     result);
  
    free(refarea);
    free(decreased);
    free(binary);
    free(image);
}
    

void analyseForstImages(int* keynr, 
			double *gausspercent, int *ngaussp,
			int *coarse, int *ncrs, 
			int *edges, int *nedg,
			double *percent, int *nperc,
			//
			double *image,
			int *binary,
			int *decreased,
			int *nrdecr,
			int *ncdecr,
			int *refarea,
			//
			int *areathreshold,
			double *result
			
    ) {
    

    if (sizeof(res_type) != sizeof(double)) {
	error("images cannot be returned (res_type)\n");	
    }
    if (sizeof(boolean) != sizeof(int)) {
	error("images cannot be returned (bool)\n");	
    }
    if (sizeof(uchar) != sizeof(int)) {
	error("images cannot be returned (uchar)\n");	
    }
    if (sizeof(uword) != sizeof(int)) {
	error("images cannot be returned (uword)\n");	
    }
    analyse_internal(keynr, gausspercent, ngaussp, coarse, ncrs, 
		     edges, nedg, percent, nperc,
		     (res_type *)image, 
		     (boolean*) binary, 
		     (uchar*) decreased, 
		     nrdecr, ncdecr,  
		     (uword*) refarea, areathreshold,
		     result);
    // assert(*areathreshold < 10);
}

