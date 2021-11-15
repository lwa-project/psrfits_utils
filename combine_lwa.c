//This code combines the two frequency bands from the Mock 
//spectrometers at Arecibo. K. Stovall  Oct 2010

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <fitsio.h>
#include <libgen.h>
#include "psrfits.h"
#include "combine_mocks_cmd.h"

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

//Routine taken from PRESTO
void avg_var(float *x, int n, double *mean, double *var)
/* For a float vector, *x, of length n, this routine  */
/* returns the mean and variance of *x.               */
{
   long i;
   double an = 0.0, an1 = 0.0, dx;

   /*  Modified (29 June 98) C version of the following:        */
   /*  ALGORITHM AS 52  APPL. STATIST. (1972) VOL.21, P.226     */
   /*  Returned values were checked with Mathematica 3.01       */

   if (n < 1) {
      printf("\vVector length must be > 0 in avg_var().  Exiting\n");
      exit(1);
   } else {
      *mean = (double) x[0];
      *var = 0.0;
   }

   for (i = 1; i < n; i++) {
      an = (double) (i + 1);
      an1 = (double) (i);
      dx = (x[i] - *mean) / an;
      *var += an * an1 * dx * dx;
      *mean += dx;
   }

   if (n > 1)
      *var /= an1;

   return;
}

//End of routine taken from PRESTO

int main(int argc, char *argv[])
{
   Cmdline *cmd;
   struct psrfits pfupper, pflower, pfo;
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
   pfupper.tot_rows = pfupper.N = pfupper.T = pfupper.status = 0;       //Initialize upper band
   pflower.tot_rows = pflower.N = pflower.T = pflower.status = 0;       //Initialize lower band
   pfupper.filenum = pflower.filenum = 1;
   pfo.tot_rows = pfo.N = pfo.T = pfo.status = pfo.multifile = 0;       //Initialize output
   sprintf(pfupper.filename, cmd->argv[0]);     //Copy filename specified on command line to
   sprintf(pflower.filename, cmd->argv[1]);     //upper and lower bands, will correct filenames shortly
//   if ((pc2 = strstr(pfupper.filename, "s1")) != NULL)  //Upper contains s1, change to s0
//      strncpy(pc2, "s0", 2);
//   else if ((pc2 = strstr(pflower.filename, "s0")) != NULL)     //Lower contains s0, change to s1
//      strncpy(pc2, "s1", 2);
//   else {
//      printf("Unable to determine which sideband is which\n");
//      exit(EXIT_FAILURE);
//   }
   //Setting the name of the output file, setting as same name as input file, but removing s0/s1. 
//   pc1 = strstr(pflower.filename, "s1");
//   pc2 = strrchr(pflower.filename, '.');        //At '.fits'
//   pc2--;
//   while ((pc2 >= pflower.filename) && isdigit(*pc2))   //Move through the digits to the separation char.
//      pc2--;
//   strncpy(outfilename, pflower.filename, pc1 - pflower.filename);      //Copy everything up to s1 into outfilename
//   strncpy(outfilename + (pc1 - pflower.filename), pc1 + 2, pc2 - pc1 - 2);     //Concatenate from after s1 to char before the separation char.
//   pc1 = outfilename + (pc2 - pflower.filename - 2);
//   *pc1 = 0;
   int rv = psrfits_open(&pfupper);   //Open upper band
   if (rv) {
      fits_report_error(stderr, rv);
      exit(1);
   }
   rv = psrfits_open(&pflower);       //Open lower band
   if (rv) {
      fits_report_error(stderr, rv);
      exit(1);
   }
   pfo = pflower;               //Copy all lower band variables into the output struct
   if (!cmd->outputbasenameP)
      sprintf(pfo.basefilename, "test");
   else
      sprintf(pfo.basefilename, cmd->outputbasename);
   pfo.filenum = 0;
   sprintf(pfo.filename, "\0"); //Set filename to null so psrfits_open will create the filename for me
   pfo.rownum = 1;
   pfo.tot_rows = 0;
   pfo.N = 0;
   printf("lower rows_per_file=%d\n",pflower.rows_per_file);
   printf("upper rows_per_file=%d\n",pfupper.rows_per_file);
   if (pfupper.rows_per_file != pflower.rows_per_file) {        //Sanity check for the two input frequency bands
      fprintf(stderr, "rows_per_file in input files do not match!\n");
      exit(1);
   }

   double upperfreqoflower, nextfromlower, lowerfreqofupper, numchandiff;       //Used to find which frequencies to take from each band
   double offsetfactor, scalefactor;    //Factors which will be applied to offsets and scales
   int upchanskip, lowchanskip; //Number of channels to skip in each banda

   //Variables used to make code cleaner
   int extrachanoffset, outoffset, upperoffset, numtocopyupper, loweroffset_skip,
       loweroffset, numtocopylower, newuppernchan, newlowernchan;
   double df = pflower.hdr.df;
   int nchan = pflower.hdr.nchan;
   int outnchan;
   int npol = pflower.hdr.npol;
   int nbits = pflower.hdr.nbits;
   int nsblk = pflower.hdr.nsblk;
   //Allocate memory for all upper and lower data
   pflower.sub.dat_freqs = (double *) malloc(sizeof(double) * nchan);
   pflower.sub.dat_weights = (float *) malloc(sizeof(float) * nchan);
   pflower.sub.dat_offsets = (float *) malloc(sizeof(float) * nchan * npol);
   pflower.sub.dat_scales = (float *) malloc(sizeof(float) * nchan * npol);
   pflower.sub.rawdata = (unsigned char *) malloc(pflower.sub.bytes_per_subint);
   pflower.sub.data = (unsigned char *) malloc(pflower.sub.bytes_per_subint*2);

   pfupper.sub.dat_freqs = (double *) malloc(sizeof(double) * nchan);
   pfupper.sub.dat_weights = (float *) malloc(sizeof(float) * nchan);
   pfupper.sub.dat_offsets = (float *) malloc(sizeof(float) * nchan * npol);
   pfupper.sub.dat_scales = (float *) malloc(sizeof(float) * nchan * npol);
   pfupper.sub.rawdata = (unsigned char *) malloc(pfupper.sub.bytes_per_subint);
   pfupper.sub.data = (unsigned char *) malloc(pfupper.sub.bytes_per_subint*2);

   int firsttime = 1;           //First time through do while loop
   do {
      print_percent_complete(pflower.rownum, pflower.rows_per_file,
                             pflower.rownum == 1 ? 1 : 0);
      psrfits_read_subint(&pflower);
      //printf("pflower[0]=%f\n",pflower.sub.dat_freqs[0]);
      //printf("pflower[nchan-1]=%f\n",pflower.sub.dat_freqs[nchan - 1]);
      psrfits_read_subint(&pfupper);
      //printf("pfupper[0]=%f\n",pfupper.sub.dat_freqs[0]);
      //printf("pfupper[nchan-1]=%f\n",pfupper.sub.dat_freqs[nchan - 1]);
      double lowerchantoremove=(pfupper.sub.dat_freqs[0]-pflower.sub.dat_freqs[nchan-1])/(pflower.hdr.BW/nchan)-1;
      //printf("lowerchantoremove=%f\n",lowerchantoremove);
      if(abs(pflower.hdr.MJD_epoch-pfupper.hdr.MJD_epoch)>1e-10)
      {
        printf("Error: MJDs do not match...Exiting.\n");
        exit(1);
      }
      if(abs(pflower.hdr.BW-pfupper.hdr.BW)>1e-10)
      {
        printf("Error: BWs do not match...Exiting.\n");
        exit(1);
      }
      int overlap=0;
      int chaninsert=0;
      if(lowerchantoremove<0)
        lowerchantoremove=lowerchantoremove*-1;
      overlap=(int)(lowerchantoremove);
      float remainder=lowerchantoremove-(float)(overlap);
      if(remainder>0.5)
        overlap++;
      if(abs(lowerchantoremove-(float)(overlap))>0.001)
      {
        printf("Error: Frequency channels do not line up well, %f ...Exiting.\n",lowerchantoremove);
        exit(1);
      }
      if (firsttime) {
         firsttime = 0;
         outnchan = nchan+nchan-overlap;
         pfo.hdr.nchan = outnchan;
         printf("outnchan=%d\n",outnchan);
         if(overlap)
           exit(1);
         pfo.hdr.BW = (double) pflower.hdr.BW*2.0;      //New bandwidth
         pfo.hdr.fctr =         //New center frequency
             pflower.sub.dat_freqs[nchan-1]+pflower.hdr.BW/nchan;
         pfo.sub.bytes_per_subint = pfo.hdr.nchan * nsblk * nbits / 8 * npol;
         pfo.sub.dat_freqs = (double *) malloc(sizeof(double) * outnchan);
         pfo.sub.dat_weights = (float *) malloc(sizeof(float) * outnchan);
         pfo.sub.dat_offsets = (float *) malloc(sizeof(float) * outnchan * npol);
         pfo.sub.dat_scales = (float *) malloc(sizeof(float) * outnchan * npol);
         pfo.sub.rawdata = (unsigned char *) malloc(pfo.sub.bytes_per_subint);
         pfo.sub.data = (unsigned char *) malloc(pfo.sub.bytes_per_subint*2);
      }
      if (pflower.status == 0 && pfupper.status == 0) {
         pfo.sub.tsubint = pflower.sub.tsubint;
         pfo.sub.offs = pflower.sub.offs;
         pfo.sub.lst = pflower.sub.lst;
         pfo.sub.ra = pflower.sub.ra;
         pfo.sub.dec = pflower.sub.dec;
         pfo.sub.glon = pflower.sub.glon;
         pfo.sub.glat = pflower.sub.glat;
         pfo.sub.feed_ang = pflower.sub.feed_ang;
         pfo.sub.pos_ang = pflower.sub.pos_ang;
         pfo.sub.par_ang = pflower.sub.par_ang;
         pfo.sub.tel_az = pflower.sub.tel_az;
         pfo.sub.tel_zen = pflower.sub.tel_zen;
         pfo.sub.FITS_typecode = pflower.sub.FITS_typecode;

         //Create variables to reduce column width of lines below
         double *dat_freqs = pfo.sub.dat_freqs;
         double *udat_freqs = pfupper.sub.dat_freqs;
         double *ldat_freqs = pflower.sub.dat_freqs;
         float *dat_weights = pfo.sub.dat_weights;
         float *udat_weights = pfupper.sub.dat_weights;
         float *ldat_weights = pflower.sub.dat_weights;
         float *dat_offsets = pfo.sub.dat_offsets;
         float *udat_offsets = pfupper.sub.dat_offsets;
         float *ldat_offsets = pflower.sub.dat_offsets;
         float *dat_scales = pfo.sub.dat_scales;
         float *udat_scales = pfupper.sub.dat_scales;
         float *ldat_scales = pflower.sub.dat_scales;
         unsigned char *data = pfo.sub.data;
         unsigned char *rawdata = pfo.sub.rawdata;
         unsigned char *udata = pfupper.sub.data;
         unsigned char *urawdata = pfupper.sub.rawdata;
         unsigned char *ldata = pflower.sub.data;
         unsigned char *lrawdata = pflower.sub.rawdata;
         memcpy(dat_freqs,ldat_freqs, sizeof(double) *nchan);
         memcpy(dat_freqs+nchan,udat_freqs, sizeof(double) *nchan);
         memcpy(dat_weights,ldat_weights, sizeof(float) *nchan);
         memcpy(dat_weights+nchan,udat_weights, sizeof(float) *nchan);
         memcpy(dat_offsets,ldat_offsets, sizeof(float) *nchan);
         memcpy(dat_offsets+nchan,udat_offsets, sizeof(float) *nchan);
         memcpy(dat_scales,ldat_scales, sizeof(float)*nchan);
         memcpy(dat_scales+nchan,udat_scales, sizeof(float)*nchan);
         int ii;
         for (ii = 0; ii < nsblk; ++ii)
         {
           if(nbits==4)
           {
             memcpy(data + ii * outnchan, ldata + ii * nchan, nchan*npol);
             memcpy(data + ii * outnchan + nchan, udata + ii * nchan, nchan*npol);
           }
           else if(nbits==8)
           {
             memcpy(rawdata + ii * outnchan, lrawdata + ii * nchan, nchan*npol);
             memcpy(rawdata + ii * outnchan + nchan, urawdata + ii * nchan, nchan*npol);
           }
         }
         psrfits_write_subint(&pfo);
      }
   } while (pfo.rownum <= pfo.rows_per_file && pfupper.status==0 && pflower.status==0);
   printf("Closing file '%s'\n", pflower.filename);
   fits_close_file(pfupper.fptr, &status);
   printf("Closing file '%s'\n", pfupper.filename);
   fits_close_file(pflower.fptr, &status);
   if(pflower.status!=0||pfupper.status!=0)
   {
     fprintf(stderr,"An error occurred when combining the two Mock files!\n");
     if(pflower.status==108||pfupper.status==108)
       fprintf(stderr,"One or both of the files is incomplete.\n");
     exit(1);
   }
   exit(0);
}
