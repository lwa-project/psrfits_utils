// This code average data from two LWA stations together.
// J. Dowell  Dec 2021

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <fitsio.h>
#include <libgen.h>
#include "psrfits.h"
#include "average_lwa_cmd.h"

static void print_percent_complete(int current, int number, int reset)
{
   static int newper = 0, oldper = -1;

   if (reset) {
      oldper = -1;
      newper = 0;
   } else {
      newper = (int) (current / (float) (number) * 100.0);
      if (newper < 0)
         newper = 0;
      if (newper > 100)
         newper = 100;
      if (newper > oldper) {
         printf("\r%3d%% ", newper);
         fflush(stdout);
         oldper = newper;
      }
   }
}

void weighted_sum(unsigned char *x1, float* s1, float* z1, float* w1,
                  unsigned char *x2, float* s2, float* z2, float* w2,
                  unsigned char *x,  float* s,  float* z,  float* w,
                  int nwin, int nchan, int npol, unsigned char clip) {
  long i, j;
  float va, vb;
  float sa, sb;
  float za, zb;
  float wa, wb;
  
  float* temp;
  temp = (float*) malloc(nwin * sizeof(float));
  memset(temp, 0, nwin*sizeof(float));

  for(i=0; i<nchan*npol; i++) {
    // Load in the scale factors, zero offsets, and weights
    sa = *(s1 + i);
    sb = *(s2 + i);
    za = *(z1 + i);
    zb = *(z2 + i);
    if(i % npol == 0) {
      wa = *(w1 + i/npol);
      wb = *(w2 + i/npol);
      *(w + i/npol) = (wa + wb) / 2.0;
    }
    
    // Compute the weighted average values and save the min/max value encountered
    float vmin, vmax;
    vmin =  1e20;
    vmax = -1e20;
    for(j=0; j<nwin; j++) {
      // Apply the scale and zero offsets
      va = *(x1 + j*nchan*npol + i) * sa + za;
      vb = *(x2 + j*nchan*npol + i) * sb + zb;
      
      // Average
      *(temp + j) = (va*wa + vb*wb) / (wa + wb);
      
      // Min/max comparisions
      if( *(temp + j) < vmin ) {
        vmin = *(temp + j);
      } else if( *(temp + j) > vmax ) {
        vmax = *(temp + j);
      }
      
      // Clean up NaNs - we need to do this after the min/max step
      if( *(temp + j) != *(temp + j) ) {
        *(temp + j) = 0.0;
      }
    }
    
    // Compute the scale factor and zero offsets for the averaged data
    if( vmin > vmax ) {
      // But first, deal with the case of no valid data populating vmin and vmax
      vmin = 0.0;
      vmax = 1.0;
    }
    *(s + i) = (vmax - vmin) / clip;
    *(z + i) = vmin;
    
    // Save the avearged data
    for(j=0; j<nwin; j++) {
      *(temp + j) -= *(z + i);
      *(temp + j) /= *(s + i);
      *(temp + j) = round(*(temp + j));
      *(x + j*nchan*npol + i) = ((unsigned char) *(temp + j)) & clip;
    }
  }
  
  free(temp);
}

int main(int argc, char *argv[])
{
   Cmdline *cmd;
   struct psrfits pfs1, pfs2, pfo;
   fitsfile *infits, *outfits;
   char *pc1, *pc2;
   char outfilename[200];       //Name of outfile if not specified on command line
   int stat = 0, padding = 0, userN = 0, status;

   // Call usage() if we have no command line arguments
   if (argc == 1) {
      Program = argv[0];
      usage();
      exit(0);
   }
   // Parse the command line using the excellent program Clig
   cmd = parseCmdline(argc, argv);
   pfs1.tot_rows = pfs1.N = pfs1.T = pfs1.status = 0;       //Initialize first station
   pfs2.tot_rows = pfs2.N = pfs2.T = pfs2.status = 0;       //Initialize second station
   pfs1.filenum = pfs2.filenum = 0;
   pfo.tot_rows = pfo.N = pfo.T = pfo.status = pfo.multifile = 0;       //Initialize output
   pfs1.filenames = (char **)malloc(sizeof(char *));
   pfs2.filenames = (char **)malloc(sizeof(char *));
   char tmpfilename[200],tmpfilename2[200];
   sprintf(tmpfilename, cmd->argv[0]);     //Copy filename specified on command line to
   pfs1.filenames[0]=tmpfilename;
   pfs1.numfiles=1;
   sprintf(tmpfilename2, cmd->argv[1]);     //both inputs, will correct filenames shortly
   pfs2.filenames[0]=tmpfilename2;
   pfs2.numfiles=1;
   int rv = psrfits_open(&pfs1);   //Open upper band
   if (rv) {
      fits_report_error(stderr, rv);
      exit(1);
   }
   rv = psrfits_open(&pfs2);       //Open lower band
   if (rv) {
      fits_report_error(stderr, rv);
      exit(1);
   }
   pfo = pfs2;               //Copy all variables from the second station into the output struct
   if (!cmd->outputbasenameP)
      sprintf(pfo.basefilename, "averaged");
   else
      sprintf(pfo.basefilename, cmd->outputbasename);
   pfo.filenum = 0;
   sprintf(pfo.filename, "\0"); //Set filename to null so psrfits_open will create the filename for me
   pfo.rownum = 1;
   pfo.tot_rows = 0;
   pfo.N = 0;
   printf("first rows_per_file=%d\n",pfs2.rows_per_file);
   printf("second rows_per_file=%d\n",pfs1.rows_per_file);
   if (pfs1.rows_per_file != pfs2.rows_per_file) {        //Sanity check for the two input frequency bands
      fprintf(stderr, "rows_per_file in input files do not match!\n");
      exit(1);
   }

   //Variables used to make code cleaner
   double df = pfs2.hdr.df;
   int nchan = pfs2.hdr.nchan;
   int npol = pfs2.hdr.npol;
   int nbits = pfs2.hdr.nbits;
   int nsblk = pfs2.hdr.nsblk;
   //Allocate memory for the data from both stations
   pfs2.sub.dat_freqs = (double *) malloc(sizeof(double) * nchan);
   pfs2.sub.dat_weights = (float *) malloc(sizeof(float) * nchan);
   pfs2.sub.dat_offsets = (float *) malloc(sizeof(float) * nchan * npol);
   pfs2.sub.dat_scales = (float *) malloc(sizeof(float) * nchan * npol);
   pfs2.sub.rawdata = (unsigned char *) malloc(pfs2.sub.bytes_per_subint);
   pfs2.sub.data = (unsigned char *) malloc(pfs2.sub.bytes_per_subint*2);

   pfs1.sub.dat_freqs = (double *) malloc(sizeof(double) * nchan);
   pfs1.sub.dat_weights = (float *) malloc(sizeof(float) * nchan);
   pfs1.sub.dat_offsets = (float *) malloc(sizeof(float) * nchan * npol);
   pfs1.sub.dat_scales = (float *) malloc(sizeof(float) * nchan * npol);
   pfs1.sub.rawdata = (unsigned char *) malloc(pfs1.sub.bytes_per_subint);
   pfs1.sub.data = (unsigned char *) malloc(pfs1.sub.bytes_per_subint*2);

   int firsttime = 1;           //First time through do while loop
   do {
      print_percent_complete(pfs2.rownum, pfs2.rows_per_file,
                             pfs2.rownum == 1 ? 1 : 0);
      psrfits_read_subint(&pfs2);
      psrfits_read_subint(&pfs1);
      if(abs(pfs2.hdr.MJD_epoch-pfs1.hdr.MJD_epoch)>1e-10)
      {
        printf("Error: MJDs do not match...Exiting.\n");
        exit(1);
      }
      if(abs(pfs2.sub.dat_freqs[0]-pfs1.sub.dat_freqs[0])>1e-10)
      {
        printf("Error: Starting frequencies do not match...Exiting.\n");
        exit(1);
      }
      if(abs(pfs2.sub.dat_freqs[1]-pfs1.sub.dat_freqs[1])>1e-10)
      {
        printf("Error: Channel widths do not match...Exiting.\n");
        exit(1);
      }
      if(abs(pfs2.hdr.BW-pfs1.hdr.BW)>1e-10)
      {
        printf("Error: BWs do not match...Exiting.\n");
        exit(1);
      }
      if (firsttime) {
         firsttime = 0;
         pfo.hdr.nchan = nchan;
         printf("outnchan=%d\n",nchan);
         //if(overlap)
         //  exit(1);
         pfo.hdr.BW = pfs2.hdr.BW;
         pfo.hdr.fctr = pfs2.hdr.fctr;
         pfo.sub.bytes_per_subint = pfo.hdr.nchan * nsblk * nbits / 8 * npol;
         pfo.sub.dat_freqs = (double *) malloc(sizeof(double) * nchan);
         pfo.sub.dat_weights = (float *) malloc(sizeof(float) * nchan);
         pfo.sub.dat_offsets = (float *) malloc(sizeof(float) * nchan * npol);
         pfo.sub.dat_scales = (float *) malloc(sizeof(float) * nchan * npol);
         pfo.sub.rawdata = (unsigned char *) malloc(pfo.sub.bytes_per_subint);
         pfo.sub.data = (unsigned char *) malloc(pfo.sub.bytes_per_subint*2);
      }
      if (pfs2.status == 0 && pfs1.status == 0) {
         pfo.sub.tsubint = pfs2.sub.tsubint;
         pfo.sub.offs = pfs2.sub.offs;
         pfo.sub.lst = pfs2.sub.lst;
         pfo.sub.ra = pfs2.sub.ra;
         pfo.sub.dec = pfs2.sub.dec;
         pfo.sub.glon = pfs2.sub.glon;
         pfo.sub.glat = pfs2.sub.glat;
         pfo.sub.feed_ang = pfs2.sub.feed_ang;
         pfo.sub.pos_ang = pfs2.sub.pos_ang;
         pfo.sub.par_ang = pfs2.sub.par_ang;
         pfo.sub.tel_az = pfs2.sub.tel_az;
         pfo.sub.tel_zen = pfs2.sub.tel_zen;
         pfo.sub.FITS_typecode = pfs2.sub.FITS_typecode;

         //Create variables to reduce column width of lines below
         double *dat_freqs = pfo.sub.dat_freqs;
         double *dat_freqs2 = pfs2.sub.dat_freqs;
         float *dat_weights = pfo.sub.dat_weights;
         float *dat_weights1 = pfs1.sub.dat_weights;
         float *dat_weights2 = pfs2.sub.dat_weights;
         float *dat_offsets = pfo.sub.dat_offsets;
         float *dat_offsets1 = pfs1.sub.dat_offsets;
         float *dat_offsets2 = pfs2.sub.dat_offsets;
         float *dat_scales = pfo.sub.dat_scales;
         float *dat_scales1 = pfs1.sub.dat_scales;
         float *dat_scales2 = pfs2.sub.dat_scales;
         unsigned char *data = pfo.sub.data;
         unsigned char *rawdata = pfo.sub.rawdata;
         unsigned char *data1 = pfs1.sub.data;
         unsigned char *rawdata1 = pfs1.sub.rawdata;
         unsigned char *data2 = pfs2.sub.data;
         unsigned char *rawdata2 = pfs2.sub.rawdata;
         memcpy(dat_freqs, dat_freqs2, sizeof(double)*nchan);
         
         if(nbits==4)
         {
           weighted_sum(data1, dat_scales1, dat_offsets1, dat_weights1,\
                        data2, dat_scales2, dat_offsets2, dat_weights2,\
                        data,  dat_scales,  dat_offsets,  dat_weights,\
                        nsblk, nchan, npol, 0x0F);
         }
         else if(nbits==8)
         {
           weighted_sum(rawdata1, dat_scales1, dat_offsets1, dat_weights1,\
                        rawdata2, dat_scales2, dat_offsets2, dat_weights2,\
                        rawdata,  dat_scales,  dat_offsets,  dat_weights,\
                        nsblk, nchan, npol, 0xFF);
         }
         psrfits_write_subint(&pfo);
      }
   } while (pfo.rownum <= pfo.rows_per_file && pfs1.status==0 && pfs2.status==0);
   printf("Closing file '%s'\n", pfs2.filename);
   fits_close_file(pfs1.fptr, &status);
   printf("Closing file '%s'\n", pfs1.filename);
   fits_close_file(pfs2.fptr, &status);
   if(pfs2.status!=0||pfs1.status!=0)
   {
     fprintf(stderr,"An error occurred when averaging the two LWA files!\n");
     if(pfs2.status==108||pfs1.status==108)
       fprintf(stderr,"One or both of the files is incomplete.\n");
     exit(1);
   }
   exit(0);
}
