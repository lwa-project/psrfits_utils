//This code creates beams for palfalfa drift scan
//K. Stovall  Feb 2011

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fitsio.h>
#include <libgen.h>
#include "psrfits.h"
#include "psrfits_drifttool_cmd.h"

int main(int argc, char *argv[])
{
   Cmdline *cmd;
   struct psrfits pfin, pfout;
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
   pfin.tot_rows = pfin.N = pfin.T = pfin.status = 0;       //Initialize input file
   pfin.filenum = 0;
   pfout.tot_rows = pfout.N = pfout.T = pfout.status = 0;       //Initialize output
   sprintf(pfin.filename, cmd->argv[0]);     //Copy filename specified on command line to
   //Get the file basename and number from command-line argument
   //(code taken from psrfits2fil)
   strcpy(outfilename, pfin.filename);       //Using outfilename as temp. variable
   pc2 = strrchr(outfilename, '.');     // at .fits
   *pc2 = 0;                    // terminate string
   pc1 = pc2 - 1;
   while ((pc1 >= outfilename) && isdigit(*pc1))
      pc1--;
   if (pc1 <= outfilename) {    // need at least 1 char before filenum
      puts("Illegal input filename. must have chars before the filenumber");
      exit(1);
   }
   strcpy(pfin.basefilename, outfilename);
   //(end of code taken from psrfits2fil)
   //Setting the name of the output file, setting as same name as input file, but removing s0/s1. 
   int rv = psrfits_open(&pfin);   //Open input file
   if (rv) {
      fits_report_error(stderr, rv);
      exit(1);
   }
   pfout = pfin;               //Copy all variables into the output struct
   if (!cmd->outputbasenameP) {   
     fprintf(stderr,"You must specify the output basename using -o\n");
     exit(1);
   }
   else
      sprintf(pfout.basefilename, cmd->outputbasename);
   int subintperfile=17,subintadvance=15;
   if(cmd->nsubintP)
     subintperfile=cmd->nsubint;
   if(cmd->advP)
     subintadvance=cmd->adv;
   pfout.filenum = 0;
   sprintf(pfout.filename, "\0"); //Set filename to null so psrfits_open will create the filename for me
   pfout.rownum = 1;
   pfout.tot_rows = 0;
   pfout.N = 0;
   //Get the beam number from the file name
   char *ibeam;
   ibeam = strrchr(pfin.basefilename, 'b');
   ibeam = ibeam+1;
   *(ibeam+1) = 0;  //terminate string
   printf("subintsperfile=%d,subintadvance=%d,rows_per_file=%d\n",subintperfile,subintadvance,pfin.rows_per_file);
   if(pfin.rows_per_file<subintperfile)
   {
      fprintf(stderr,"File does not contain enough subints\n");
      exit(1);
   }
   if(subintadvance>subintperfile)
   {
      fprintf(stderr,"Number of subints to advance is larger than the number of subints per file?\n");
      exit(1);
   }
   int numfiles=(pfin.rows_per_file-subintperfile)/subintadvance+1;
   if(pfin.rows_per_file%subintadvance!=0)
     ++numfiles;
   printf("Creating %d files\n",numfiles);
   char templatename[200];
   sprintf(templatename,"%s.template.fits",pfout.basefilename);
   printf("%s\n",templatename);
   status = 0; //fits_create_file fails if this is not set to zero
//   fits_create_file(&outfits, templatename, &status);

   //Instead of copying HDUs one by one, move to the SUBINT HDU
   //and copy all the HDUs preceding it
//   infits=pfin.fptr;
//   fits_movnam_hdu(infits, BINARY_TBL, "SUBINT", 0, &status);
//   fits_copy_file(infits, outfits, 1, 0, 0, &status);

   //Copy the SUBINT table header
//   fits_copy_header(infits, outfits, &status);
   //fprintf(stderr,"After fits_copy_header, status: %d\n", status);
//   fits_flush_buffer(outfits, 0, &status);

   //Set NAXIS2 in the output SUBINT table to 0 b/c we haven't
   //written any rows yet
   int dummy = 0;
//   fits_update_key(outfits, TINT, "NAXIS2", &dummy, NULL, &status);
//   fits_movabs_hdu(outfits, 1, NULL, &status);
//   fits_update_key(outfits, TSTRING, "IBEAM", ibeam, "Beam number for multibeam systems", &status);

//   fits_close_file(outfits, &status);
printf("1\n");   
   int nchan=pfin.hdr.nchan;
   int npol=pfin.hdr.npol;

   pfin.sub.dat_freqs = (float *) malloc(sizeof(float) * nchan);
   pfin.sub.dat_weights = (float *) malloc(sizeof(float) * nchan);
   pfin.sub.dat_offsets = (float *) malloc(sizeof(float) * nchan * npol);
   pfin.sub.dat_scales = (float *) malloc(sizeof(float) * nchan * npol);
   pfin.sub.data = (unsigned char *) malloc(pfin.sub.bytes_per_subint);
   pfin.sub.rawdata = (unsigned char *) malloc(pfin.sub.bytes_per_subint);
   int count=0;
   char ralist[numfiles][16];
   double radeglist[numfiles];
   char declist[numfiles][16];
   double decdeglist[numfiles];
   char basefilenames[numfiles][100];
   int ii;
printf("2\n");   
   while (count<pfin.rows_per_file)
   {
//printf("3\n");   
     psrfits_read_subint(&pfin);
//printf("4\n");   
     for(ii=0;ii<numfiles;++ii)
     {
       if(count==subintadvance*ii+subintperfile/2)
       {
          printf("%f %f\n",pfin.sub.ra,pfin.sub.dec);
          radeglist[ii]=pfin.sub.ra;
          decdeglist[ii]=pfin.sub.dec;
          int rah=pfin.sub.ra/15;
          int ram=((pfin.sub.ra/15)-rah)*60;
          double ras=((((pfin.sub.ra/15)-rah)*60)-ram)*60;
          int ded,dem;
          double des;
          if(decdeglist[ii]<0)
          {
            ded=abs(pfin.sub.dec);
            dem=abs((pfin.sub.dec+ded)*60);
            des=fabs((((pfin.sub.dec+ded)*60)+dem)*60);
          }
          else
          {
            ded=pfin.sub.dec;
            dem=(pfin.sub.dec-ded)*60;
            des=(((pfin.sub.dec-ded)*60)-dem)*60;
          }
          if(ras>10)
            sprintf(ralist[ii],"%02d:%02d:%.4f",rah,ram,ras);
          else
            sprintf(ralist[ii],"%02d:%02d:0%.4f",rah,ram,ras);
          if(des>10)
            if(decdeglist[ii]<0)
              sprintf(declist[ii],"-%02d:%02d:%.4f",ded,dem,des);
            else
              sprintf(declist[ii],"%02d:%02d:%.4f",ded,dem,des);
          else
            if(decdeglist[ii]<0)
              sprintf(declist[ii],"-%02d:%02d:0%.4f",ded,dem,des);
            else
              sprintf(declist[ii],"%02d:%02d:0%.4f",ded,dem,des);
          if(decdeglist[ii]<0)
            sprintf(basefilenames[ii],"%s.D%02d%02d%02d-%02d%02d",cmd->outputbasename,rah,ram,(int)ras,ded,dem);
          else
            sprintf(basefilenames[ii],"%s.D%02d%02d%02d+%02d%02d",cmd->outputbasename,rah,ram,(int)ras,ded,dem);
       }
     }
     count++;
   }
   struct psrfits pfo[numfiles];
   int filesactive[numfiles];
   for(ii=0;ii<numfiles;++ii)
   {
     printf("%s\n",basefilenames[ii]);
     pfo[ii].tot_rows = pfo[ii].N = pfo[ii].T = pfo[ii].status = 0;       //Initialize output
     pfo[ii] = pfin;               //Copy all variables into the output struct
     sprintf(pfo[ii].basefilename, basefilenames[ii]);
     pfo[ii].filenum = 0;
     sprintf(pfo[ii].filename, "\0"); //Set filename to null so psrfits_open will create the filename for me
     pfo[ii].rownum = 1;
     pfo[ii].tot_rows = 0;
     pfo[ii].N = 0;
     sprintf(pfo[ii].hdr.ra_str,ralist[ii]);
     sprintf(pfo[ii].hdr.dec_str,declist[ii]);
     pfo[ii].hdr.ra2000=radeglist[ii];
     pfo[ii].hdr.dec2000=decdeglist[ii];
     int nchan=pfo[ii].hdr.nchan;
     int npol=pfo[ii].hdr.npol;
     printf("bytes_per_subint=%d/%d\n",pfin.sub.bytes_per_subint,(nchan*npol*pfo[ii].hdr.nsblk));
     pfo[ii].sub.dat_freqs = (float *) malloc(sizeof(float) * nchan);
     pfo[ii].sub.dat_weights = (float *) malloc(sizeof(float) * nchan);
     pfo[ii].sub.dat_offsets = (float *) malloc(sizeof(float) * nchan * npol);
     pfo[ii].sub.dat_scales = (float *) malloc(sizeof(float) * nchan * npol);
     pfo[ii].sub.data = (unsigned char *) malloc(pfin.sub.bytes_per_subint);
     pfo[ii].sub.rawdata = (unsigned char *) malloc(pfin.sub.bytes_per_subint);
     pfo[ii].sub.offs=0;
     filesactive[ii]=0;
   }
   psrfits_close(&pfin);
   rv = psrfits_open(&pfin);   //Open input file again
   count=0;
   while(count<pfin.rows_per_file)
   {
     printf("Count:%d/%d\n",count,pfin.rows_per_file);
     psrfits_read_subint(&pfin);
     printf("%d\n",pfin.sub.rawdata[0]);
     if(count%subintadvance==0)
     {
       if(count>=subintperfile)
       {
         psrfits_close(&pfo[count/subintadvance-2]);
         filesactive[count/subintadvance-2]=0;
         printf("Deactivating %d\n",(count/subintadvance-2));
       }
       if(count<(pfin.rows_per_file-subintadvance))
       {
         psrfits_create(&pfo[count/subintadvance]);
         filesactive[count/subintadvance]=1;
         printf("Activating %d\n",(count/subintadvance));
       }
     }
     for(ii=0;ii<numfiles;++ii)
     {
       if(filesactive[ii]==1)
       {
         pfo[ii].sub.tsubint=pfin.sub.tsubint;
         pfo[ii].sub.offs=pfin.sub.offs;
         pfo[ii].sub.lst=pfin.sub.lst;
         pfo[ii].sub.ra=pfin.sub.ra;
         pfo[ii].sub.dec=pfin.sub.dec;
         pfo[ii].sub.glon=pfin.sub.glon;
         pfo[ii].sub.glat=pfin.sub.glat;
         pfo[ii].sub.feed_ang=pfin.sub.feed_ang;
         pfo[ii].sub.pos_ang=pfin.sub.pos_ang;
         pfo[ii].sub.par_ang=pfin.sub.par_ang;
         pfo[ii].sub.tel_az=pfin.sub.tel_az;
         pfo[ii].sub.tel_zen=pfin.sub.tel_zen;
         pfo[ii].sub.FITS_typecode=pfin.sub.FITS_typecode;
         memcpy(pfo[ii].sub.dat_freqs,pfin.sub.dat_freqs,sizeof(float)*nchan);
         memcpy(pfo[ii].sub.dat_weights,pfin.sub.dat_weights,sizeof(float)*nchan*npol);
         memcpy(pfo[ii].sub.dat_offsets,pfin.sub.dat_offsets,sizeof(float)*nchan*npol);
         memcpy(pfo[ii].sub.dat_scales,pfin.sub.dat_scales,sizeof(float)*nchan*npol);
         memcpy(pfo[ii].sub.data,pfin.sub.data,pfin.sub.bytes_per_subint);
         printf("%d:%d ",ii,pfin.sub.rawdata[0]);
         memcpy(pfo[ii].sub.rawdata,pfin.sub.rawdata,pfin.sub.bytes_per_subint);
         printf("%d\n",pfo[ii].sub.rawdata[0]);
         psrfits_write_subint(&pfo[ii]);
       }
     }
     count++;
   }
   if(pfin.rows_per_file%subintadvance!=0)
     psrfits_close(&pfo[numfiles-1]);
   exit(0);
}
