#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "io_serial.h"
#include "../libutility/utility.h"


/*================================================================================
 * write particles to a file (either dumpfile or output file)
 *================================================================================*/
void output(char *outfile_name, int dumpflag)
{
  partptr cur_part;                  /* current particle               */
  int     i,j;                       /* loop index                     */
  int     slen;                      /* string length                  */
  char    file_ncpy[MAXSTRING];      /* copy of file name              */
  char    file_no[10];               /* file number                    */
  char   *f_name;                    /* file name pointer              */
  FILE   *outstream;                 /* file pointer                   */
  int     machine_sizeof_long;       /* store sizeof(long) in outfile  */
  flouble pos, mom, fweight;
   
  
  /* copy timestep dependent quantities over to io.header */
  io.header.no_timestep = global.no_timestep;
  io.header.a_current   = global.a;
  io.header.K_current   = energy.K_current;
  io.header.U_current   = energy.U_current;
  io.header.Eintegral   = energy.integral;
  io.header.Econst      = energy.econst;
  
  
  if(dumpflag == 0)   /* normal output */
    {
     f_name = strcpy(file_ncpy, outfile_name);
	  
#ifdef REDSHIFTNAME
     /* OUTPUT CONVENTION: 3 digits*/
     sprintf(file_no,"z%.3f",fabs(global.z));
#else
     sprintf(file_no,"%05ld",global.no_timestep);
#endif
     strcat(f_name, file_no);
     
       /* open file */
       if ((outstream = fopen(f_name,"wb")) == NULL) 
       {
          fprintf(io.logfile,"output: could not open file %s\n", f_name);
          exit(1);
       }
    }
  else
    {
      if((outstream = fopen(outfile_name,"wb")) == NULL)
        {
         fprintf(io.logfile,"output: could not open dump file %s\n",
                 outfile_name);
         exit(1);
        }
    }
  
  
  /* write a simple "1" to file (for BYTESWAP testing) */
  machine_sizeof_long = sizeof(long);
  fwrite(&machine_sizeof_long, sizeof(int), 1, outstream);

  /* write IO header */
  if(fwrite(&(io.header), sizeof(io.header), 1, outstream) != 1)
    {
      fprintf(io.logfile,"\n\noutput: could not write io.header\n");
      fflush(io.logfile);
      fclose(io.logfile);
      exit(1);
    }

  cur_part = global.fst_part;
  for(i = 0; i < global.no_part; i++)
   {
     for(j = X; j <= Z; j++)
       {
        cur_part->pos[j] = fmod((double)cur_part->pos[j] + 1.0, 1.0);
        
        pos = cur_part->pos[j];
        mom = cur_part->mom[j];
        
        fwrite((void*)&pos, sizeof(flouble), 1, outstream);
        fwrite((void*)&mom, sizeof(flouble), 1, outstream);
       }
#ifdef MULTIMASS
     fweight  = cur_part->weight;
     fwrite((void*)&fweight, sizeof(flouble), 1, outstream);
#endif  // MULTIMASS

     /* move to next particle */
     cur_part++;
   }
  fflush(outstream);
  fclose(outstream);
}


/*==============================================================================
 * here we decide whether to write an output, a dumpfile or nothing...
 *==============================================================================*/
void manage_outputs(double timecounter, double timestep)
{
  int     iout;
  double  a_out;             /* check for output file               */
  double  a_left, a_right;
  
  if(global.terminate == FALSE)
    {
#ifdef LIGHTCONE 
      /* write the lightcone */
      if(global.z<=simu.z_lightcone)
        output_lc();
#endif /* LIGHTCONE */
      
#ifndef OUTDUMPS
      /* write output ? */
      if(global.output_count < io.no_outputs)
        {
          /* bracket current expansion factor [a_left, a_right] */
          a_left  = calc_super_a(timecounter-timestep/2.0);
          a_right = calc_super_a(timecounter+timestep/2.0);
          
          /* check whether there is an a_out within [a_left, a_right] */
          for(iout=0; iout<io.no_outputs; iout++)
            {
              a_out   = io.a_out[iout];
              
              if(a_left <= a_out && a_out < a_right)
                {
                  fprintf(io.logfile,"\nstarting to write output files...");
                  fflush(io.logfile);
                  output(io.outfile_prefix, 0);
#ifdef HYDRO
                  output_grid(global.dom_grid, 0);
#if (NO_DM || DM_GAS || HYDRO_TEST==7)
                  write_hydro(global.dom_grid);
#endif
#endif
                  global.output_count++;
                  fprintf(io.logfile, "done\n");
                  fflush(io.logfile);
                  
                  /* we wrote an output file... */
                  global.ioflag      = TRUE;
                  PkSpectrum.dump_Pk = TRUE;
                }
            }
        }
#endif /* OUTDUMPS */
      
      /* write dumpfile ? */
      if((global.no_timestep % io.out_dumps) == 0 && global.no_timestep > 0)
        {
          fprintf(io.logfile,"\nstarting to dump output files...");
          fflush(io.logfile);
          output(io.dumpfile_name, 1);
#ifdef HYDRO
          output_grid(global.dom_grid, 1);
#if (NO_DM || DM_GAS || HYDRO_TEST==7)
          write_hydro(global.dom_grid);
#endif
#endif
          fprintf(io.logfile, "done\n");
          fflush(io.logfile);
          
#ifdef OUTDUMPS
          fprintf(io.logfile,"\nstarting to write output files...");
          fflush(io.logfile);
          output(io.outfile_prefix, 0);
#ifdef HYDRO
          output_grid(global.dom_grid, 0);
#endif
          fprintf(io.logfile, "done\n");
          fflush(io.logfile);
          
          /* we wrote an output file... */
          global.ioflag      = TRUE;
          PkSpectrum.dump_Pk = TRUE;
#endif /* OUTDUMPS */
          
#ifdef LIGHTCONE 
          /* synchronize the dumpfile with the partial lightcone files
           * (if the lightcone file exists) */
          if(io.lightcone!=NULL) 
            {
              fclose(io.lightcone);
              io.lightcone=NULL;
            }
#endif /* LIGHTCONE */
        }      
    }
}
