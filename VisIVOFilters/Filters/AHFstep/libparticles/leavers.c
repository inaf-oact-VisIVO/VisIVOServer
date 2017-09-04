#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "particles.h"
#include "../libutility/utility.h"
#include "../libgrids/grids.h"

/*
 * NOTE:  if you want the code to terminate on any suspicious events
 *        related to undrifting and unkicking, please compile the code
 *        with: 
 *              -DTERMINATE2
 */

/*===============================================================================
* store cur_part within cur_grid leavers list
*===============================================================================*/
void store_leaver(partptr cur_part, gridls *cur_grid)
{
   /* 1.-2. send back that particle to have it moved on next coarser level */
   if(cur_grid->multistep < 3)
     {
      cur_grid->leavers.no_sendback++;
      cur_grid->leavers.send_back = 
         (partptr *) realloc(cur_grid->leavers.send_back, 
                             cur_grid->leavers.no_sendback*sizeof(partptr));
      
      *(cur_grid->leavers.send_back+cur_grid->leavers.no_sendback-1) = cur_part;
      
      /* 
         * do *NOT* transfer cur_part to next coarser grid:
       *
       *    it's better to undrift it beforehand !
       */
     }
   
   /* 3. coa_grid was already kicked and cur_part needs synchronisation */
   else if(cur_grid->multistep == 3)
     {
      cur_grid->leavers.no_keepmoving++;
      cur_grid->leavers.keep_moving = 
         (partptr *) realloc(cur_grid->leavers.keep_moving, 
                             cur_grid->leavers.no_keepmoving*sizeof(partptr));
      
      *(cur_grid->leavers.keep_moving+cur_grid->leavers.no_keepmoving-1) = cur_part;
      
      /* 
         * transfer particle to next coarser grid is needed:
       *
       *   otherwise it's density contribution is lost !
       */
#ifdef VERBOSELOG2
      fprintf(io.logfile,"store_leaver: trying to link particle to grid %ld\n",cur_grid->l1dim/2);
#endif
      coagrid_llsearch(cur_part, cur_grid-1);
     }
   else
     {
#ifdef VERBOSELOG2
      fprintf(io.logfile,"store_leaver: trying to link particle to grid %ld\n",cur_grid->l1dim/2);
#endif
      /* end of multistep: just transfer the particle to the next coarser grid */
      coagrid_llsearch(cur_part, cur_grid-1);
     }
}

/*-------------------------------------------------------------------------------
* undo latest DRIFT for supported particle cur_part
*-------------------------------------------------------------------------------*/
void undrift_leaver(partptr cur_part, double DRIFT)
{
   int    i;
   double distance;
   
   for(i = X; i <= Z; i++)
     {
      distance = (double)cur_part->mom[i] * DRIFT;

      cur_part->pos[i] = ((double)cur_part->pos[i] - distance);
      cur_part->pos[i] = f1mod((double)cur_part->pos[i]+1.0, 1.0);
     }
}

/*-------------------------------------------------------------------------------
* do DRIFT for supported particle cur_part
*-------------------------------------------------------------------------------*/
void drift_leaver(partptr cur_part, double DRIFT, gridls *cur_grid)
{
   int    i;
   double distance;
   
   for(i = X; i <= Z; i++)
     {
      distance          = (double)cur_part->mom[i] * DRIFT;

      cur_part->pos[i]  = ((double)cur_part->pos[i] + distance);
      cur_part->pos[i]  = f1mod((double)cur_part->pos[i]+1.0, 1.0);
       
#ifdef TERMINATE2
      /* check distance travelled by particle */
      if(distance > CELLFRAC_MAX * cur_grid->spacing)
        {
         fprintf(io.logfile,"drift_leaver: %f (%ld = %d)\n", 
                 distance/cur_grid->spacing, 
                 cur_grid->l1dim, cur_grid->multistep);
         fflush(io.logfile);
#ifndef ATS
         global.terminate = TRUE;
#endif
        }
#endif
     }
}

/*-------------------------------------------------------------------------------
* undo latest KICK for supported particle cur_part
*-------------------------------------------------------------------------------*/
void unkick_leaver(partptr cur_part, double KICK, gridls *cur_grid)
{
   nptr   tsc_nodes[3][3][3];  /* nodes to assign to                 */
   gridls *for_grid;
   int    idim;                /* coord changing index               */
   int    i,j,k;               /* indices for 3D arrays              */
   double pnarg_a;             /* pn sep temp arg                    */
   double pnarg_b;             /* pn sep temp arg                    */
   double tpnarg;              /* temp to calc pn args               */
   dvect  temp_coords;         /* double vector - calc un_mod coords */
   dvect  xyz_coords;          /* actual double coords of node       */
   dvect  pn_sep;              /* particle-node separation / spacing */
   dvect  weights[3];          /* weights in each dimension          */
   double force_comp[NDIM];    /* array for force components         */
   double momupdate;
   
   long   n[NDIM];             /* particle - node coords        */
   pqptr  cur_pquad;           /* current pquad                 */
   cqptr  cur_cquad;           /* current cquad                 */
   nqptr  cur_nquad;           /* current nquad                 */
   int    exp;                 /* exponent for frexp            */
   double edge_shift;          /* coords of edges of first cell */
   double cur_shift;
   
   /* shift of cell centre as compared to edge of box [grid units] */
   cur_shift = 0.5/(double)cur_grid->l1dim;
   
   /* find node orginally carrying cur_part */
   edge_shift = 0.0;
   
   /* coordintes of that node */
   for(i = X; i <= Z; i++)
      n[i] = (long)((double)cur_grid->l1dim * 
                    f1mod((double)cur_part->pos[i]-edge_shift+1.0, 1.0));
   
   
   /*---------------------------------------------------------
    * find node within QUAD's by running down linked lists... 
    *---------------------------------------------------------*/
   for(cur_pquad = cur_grid->pquad;
       (n[Z] > cur_pquad->z + cur_pquad->length) && (cur_pquad->next != NULL); 
       cur_pquad = cur_pquad->next)
      ;
   
   if((n[Z] >= cur_pquad->z) && (n[Z] < cur_pquad->z + cur_pquad->length))
     {
      for(cur_cquad = cur_pquad->loc + (n[Z] - cur_pquad->z);
          (n[Y] > cur_cquad->y + cur_cquad->length) && (cur_cquad->next != NULL); 
          cur_cquad = cur_cquad->next)
         ;
      
      if((n[Y] >= cur_cquad->y) && (n[Y] < cur_cquad->y + cur_cquad->length))
        {
         for(cur_nquad = cur_cquad->loc + (n[Y] - cur_cquad->y);
             (n[X] > cur_nquad->x + cur_nquad->length)&&(cur_nquad->next != NULL);
             cur_nquad = cur_nquad->next)
            ;
         
         if((n[X] >= cur_nquad->x) && (n[X] < cur_nquad->x + cur_nquad->length))
           {
            /*------------------------------------
             * eventually reached correct node...
             *------------------------------------*/
            tsc_nodes[1][1][1] = cur_nquad->loc + (n[X] - cur_nquad->x);
            
            get_TSCnodes(cur_grid, cur_pquad, cur_cquad, cur_nquad, tsc_nodes,
                         &n[Z], &n[Y], &n[X]);
            
            if(test_tsc(tsc_nodes) == FALSE)
              {
#ifdef TERMINATE2
               fprintf(io.logfile,
                       "\nunkick_leaver: particle %ld can't be un-kicked on %ld grid\n",
                       cur_part,cur_grid->l1dim);
               fflush(io.logfile);
               global.terminate = TRUE;
#endif
#ifdef KICK_ANYWAYS
               /* forces on next coarser level have NOT been calculated yet 
                  *  ...and hence kicking it on cur_grid-1 is sort of dangerous */
               unkick_leaver(cur_part, KICK, cur_grid-1);
#endif
               return;
              }
            
            /* calculate realspace coords of cur_node */
            temp_coords[X] = (((double)n[X])/(double)cur_grid->l1dim)+cur_shift;
            temp_coords[Y] = (((double)n[Y])/(double)cur_grid->l1dim)+cur_shift;
            temp_coords[Z] = (((double)n[Z])/(double)cur_grid->l1dim)+cur_shift;
            xyz_coords[X]  = f1mod(temp_coords[X]+1.0, 1.0);
            xyz_coords[Y]  = f1mod(temp_coords[Y]+1.0, 1.0);
            xyz_coords[Z]  = f1mod(temp_coords[Z]+1.0, 1.0);
            
            /* calc fraction to be un-assigned from each tsc node */
            for(idim = 0; idim < NDIM; idim++)
              {
               pn_sep[idim] = (double)cur_grid->l1dim *
               ((double)cur_part->pos[idim] - xyz_coords[idim]);
               
               if(fabs(pn_sep[idim]) > 0.5*(double)cur_grid->l1dim)
                 {
                  tpnarg          = (double)cur_part->pos[idim] + 0.5;
                  pnarg_a         = f1mod(tpnarg+1.0, 1.0);
                  tpnarg          = xyz_coords[idim] + 0.5;
                  pnarg_b         = f1mod(tpnarg+1.0, 1.0);
                  pn_sep[idim]    = (pnarg_a-pnarg_b) * (double)cur_grid->l1dim;
                 }
#ifdef TSC
               weights[1][idim] = 0.75 - pow2(pn_sep[idim]);
               weights[0][idim] = pow2((0.5 - pn_sep[idim]))/2;
               weights[2][idim] = pow2((0.5 + pn_sep[idim]))/2;
#endif
#ifdef CIC
               if(pn_sep[idim] > 0.)
                 {
                  weights[0][idim] = 0.0;
                  weights[1][idim] = 1.0 - pn_sep[idim];
                  weights[2][idim] =       pn_sep[idim];
                 }
               else
                 {
                  weights[0][idim] =     - pn_sep[idim];
                  weights[1][idim] = 1.0 + pn_sep[idim];
                  weights[2][idim] = 0.0;
                 }
#endif
#ifdef NGP
               weights[0][idim] = 0.0;
               weights[1][idim] = 1.0;
               weights[2][idim] = 0.0;
#endif
               force_comp[idim] = 0.0;     /* initialize comp. array */
              }
            
            for(k = 0; k < 3; k++)
               for(j = 0; j < 3; j++)
                  for(i = 0; i < 3; i++)
                     for(idim = X; idim <= Z; idim++)
                       {
                        force_comp[idim] +=
                        weights[k][Z]*weights[j][Y]*weights[i][X]
                        * (double)tsc_nodes[k][j][i]->force.forces[idim];
                       }  
                        
            /* un-kick particle */
            for(idim = X; idim <= Z; idim++)
              {
               momupdate   = (double)force_comp[idim] * KICK;

               /* eventually un-kick momentum */
               cur_part->mom[idim] -= momupdate;
              }
           }
         else
           {
#ifdef TERMINATE2
            fprintf(io.logfile,
                    "\nunkick_leaver[X]: particle %ld not found in %ld grid ?!\n",
                    cur_part,cur_grid->l1dim);
            fflush(io.logfile);
            global.terminate = TRUE;
#endif
           }
        }
      else
        {
#ifdef TERMINATE2
         fprintf(io.logfile,
                 "\nunkick_leaver[Y]: particle %ld not found in %ld grid ?!\n",
                 cur_part,cur_grid->l1dim);
         fflush(io.logfile);
         global.terminate = TRUE;
#endif
        }
     }
   else
     {
#ifdef TERMINATE2
      fprintf(io.logfile,"\nunkick_leaver[Z]: %ld particle not found in %ld grid ?!\n",
              cur_part,cur_grid->l1dim);
      fflush(io.logfile);
      global.terminate = TRUE;
#endif
     }
}

/*-------------------------------------------------------------------------------
* kick particle linked to cur_grid with KICK
*-------------------------------------------------------------------------------*/
void kick_leaver(partptr cur_part, double KICK, gridls *cur_grid)
{
   nptr    tsc_nodes[3][3][3];  /* nodes to assign to                 */
   gridls *for_grid;
   int     idim;                /* coord changing index               */
   int     i,j,k;               /* indices for 3D arrays              */
   double  pnarg_a;             /* pn sep temp arg                    */
   double  pnarg_b;             /* pn sep temp arg                    */
   double  tpnarg;              /* temp to calc pn args               */
   dvect   temp_coords;         /* double vector - calc un_mod coords */
   dvect   xyz_coords;          /* actual double coords of node       */
   dvect   pn_sep;              /* particle-node separation / spacing */
   dvect   weights[3];          /* weights in each dimension          */
   double  force_comp[NDIM];    /* array for force components         */
   double  momupdate;
   
   long   n[NDIM];             /* particle - node coords        */
   pqptr  cur_pquad;           /* current pquad                 */
   cqptr  cur_cquad;           /* current cquad                 */
   nqptr  cur_nquad;           /* current nquad                 */
   int    exp;                 /* exponent for frexp            */
   double edge_shift;          /* coords of edges of first cell */
   double cur_shift;
   
   /* shift of cell centre as compared to edge of box [grid units] */
   cur_shift = 0.5/(double)cur_grid->l1dim;
   
   /* find node orginally carrying cur_part */
   edge_shift = 0.0;
   
   /* coordintes of that node */
   for(i = X; i <= Z; i++)
      n[i] = (long)((double)cur_grid->l1dim * f1mod((double)cur_part->pos[i]-edge_shift+1.0, 1.0));
   
   
   /*---------------------------------------------------------
      * find node within QUAD's by running down linked lists... 
      *---------------------------------------------------------*/
   for(cur_pquad = cur_grid->pquad;
       (n[Z] > cur_pquad->z + cur_pquad->length) && (cur_pquad->next != NULL); 
       cur_pquad = cur_pquad->next)
      ;
   
   if((n[Z] >= cur_pquad->z) && (n[Z] < cur_pquad->z + cur_pquad->length))
     {
      for(cur_cquad = cur_pquad->loc + (n[Z] - cur_pquad->z);
          (n[Y] > cur_cquad->y + cur_cquad->length) && (cur_cquad->next != NULL); 
          cur_cquad = cur_cquad->next)
         ;
      
      if((n[Y] >= cur_cquad->y) && (n[Y] < cur_cquad->y + cur_cquad->length))
        {
         for(cur_nquad = cur_cquad->loc + (n[Y] - cur_cquad->y);
             (n[X] > cur_nquad->x + cur_nquad->length)&&(cur_nquad->next != NULL);
             cur_nquad = cur_nquad->next)
            ;
         
         if((n[X] >= cur_nquad->x) && (n[X] < cur_nquad->x + cur_nquad->length))
           {
            /*------------------------------------
            * eventually reached correct node...
            *------------------------------------*/
            tsc_nodes[1][1][1] = cur_nquad->loc + (n[X] - cur_nquad->x);
            
            get_TSCnodes(cur_grid, cur_pquad, cur_cquad, cur_nquad, tsc_nodes,
                         &n[Z], &n[Y], &n[X]);
            
            if(test_tsc(tsc_nodes) == FALSE)
              {
#ifdef TERMINATE2
               fprintf(io.logfile,
                       "\nkick_leaver: particle %ld can't be kicked on %ld grid\n",
                       cur_part,cur_grid->l1dim);
               fflush(io.logfile);
               /*global.terminate = TRUE;*/
#endif
#ifdef KICK_ANYWAYS
               /* forces on next coarser level have NOT been calculated yet 
                  *  ...and hence kicking it on cur_grid-1 is sort of dangerous */
               kick_leaver(cur_part, KICK, cur_grid-1);
#endif
               return;
              }
            
            /* calculate realspace coords of cur_node */
            temp_coords[X] = (((double)n[X])/(double)cur_grid->l1dim)+cur_shift;
            temp_coords[Y] = (((double)n[Y])/(double)cur_grid->l1dim)+cur_shift;
            temp_coords[Z] = (((double)n[Z])/(double)cur_grid->l1dim)+cur_shift;
            xyz_coords[X]  = f1mod(temp_coords[X]+1.0, 1.0);
            xyz_coords[Y]  = f1mod(temp_coords[Y]+1.0, 1.0);
            xyz_coords[Z]  = f1mod(temp_coords[Z]+1.0, 1.0);
            
            /* calc fraction to be un-assigned from each tsc node */
            for(idim = 0; idim < NDIM; idim++)
              {
               pn_sep[idim] = (double)cur_grid->l1dim *
               ((double)cur_part->pos[idim] - xyz_coords[idim]);
               
               if(fabs(pn_sep[idim]) > 0.5*(double)cur_grid->l1dim)
                 {
                  tpnarg          = (double)cur_part->pos[idim] + 0.5;
                  pnarg_a         = f1mod(tpnarg+1.0, 1.0);
                  tpnarg          = xyz_coords[idim] + 0.5;
                  pnarg_b         = f1mod(tpnarg+1.0, 1.0);
                  pn_sep[idim]    = (pnarg_a-pnarg_b) * (double)cur_grid->l1dim;
                 }
#ifdef TSC
               weights[1][idim] = 0.75 - pow2(pn_sep[idim]);
               weights[0][idim] = pow2((0.5 - pn_sep[idim]))/2;
               weights[2][idim] = pow2((0.5 + pn_sep[idim]))/2;
#endif
#ifdef CIC
               if(pn_sep[idim] > 0.)
                 {
                  weights[0][idim] = 0.0;
                  weights[1][idim] = 1.0 - pn_sep[idim];
                  weights[2][idim] =       pn_sep[idim];
                 }
               else
                 {
                  weights[0][idim] =     - pn_sep[idim];
                  weights[1][idim] = 1.0 + pn_sep[idim];
                  weights[2][idim] = 0.0;
                 }
#endif
#ifdef NGP
               weights[0][idim] = 0.0;
               weights[1][idim] = 1.0;
               weights[2][idim] = 0.0;
#endif
               force_comp[idim] = 0.0;     /* initialize comp. array */
              }
            
            for(k = 0; k < 3; k++)
               for(j = 0; j < 3; j++)
                  for(i = 0; i < 3; i++)
                     for(idim = X; idim <= Z; idim++)
                       {
                        force_comp[idim] +=
                        weights[k][Z]*weights[j][Y]*weights[i][X]
                        * (double)tsc_nodes[k][j][i]->force.forces[idim];
                       }
                        /* kick particle */
                        for(idim = X; idim <= Z; idim++)
                          {
                           momupdate   = (double)force_comp[idim] * KICK;

                           /* eventually kick momentum */
                           cur_part->mom[idim] +=  momupdate;
                          }
           }
                  else
                    {
#ifdef TERMINATE2
                     fprintf(io.logfile,
                             "\nkick_leaver: particle %ld not found in %ld grid ?!\n",
                             cur_part,cur_grid->l1dim);
                     fflush(io.logfile);
                     global.terminate = TRUE;
#endif
                    }
        }
            else
              {
#ifdef TERMINATE2
               fprintf(io.logfile,
                       "\nkick_leaver: particle %ld not found in %ld grid ?!\n",
                       cur_part,cur_grid->l1dim);
               fflush(io.logfile);
               global.terminate = TRUE;
#endif
              }
     }
         else
           {
#ifdef TERMINATE2
            fprintf(io.logfile,
                    "\nkick_leaver: particle %ld not found in %ld grid ?!\n",
                    cur_part,cur_grid->l1dim);
            fflush(io.logfile);
            global.terminate = TRUE;
#endif
           }
}


/*=============================================================================
* unstep DRIFT's and KICK's for leavers
*=============================================================================*/
void unstep_leavers(gridls *cur_grid, double DRIFT1, double DRIFT2, double KICK)
{
   partptr cur_part;
   int     ileaver;
   
   /* keep track of percentage of particles crossing grid boundaries */
   global.no_leavers += cur_grid->leavers.no_sendback;
   global.no_lsteps++;
   
   /* send back early leavers */
   if(cur_grid->leavers.no_sendback > 0)
     {
      for(ileaver = 0; ileaver < cur_grid->leavers.no_sendback; ileaver++)
        {
         cur_part = *(cur_grid->leavers.send_back+ileaver);
         
         if(cur_grid->multistep == 1)
           {
            undrift_leaver(cur_part, DRIFT1);
           }
         else if(cur_grid->multistep == 2)
           {
            undrift_leaver(cur_part, DRIFT2);
            unkick_leaver (cur_part, KICK, cur_grid);
            undrift_leaver(cur_part, DRIFT1);
           }
         else
           {
#ifdef TERMINATE2
            fprintf(io.logfile,
                    "\nunstep_leavers:  send-back for multistep > 2 ?!\n");
            fflush(io.logfile);
            global.terminate = TRUE;
#endif
           }
         
#ifdef VERBOSELOG2
         fprintf(io.logfile,"unstep_leavers: trying to link particle to grid %ld\n",cur_grid->l1dim/2);
#endif
         /* now it's safe to transfer the particle to the next coarser grid */
         coagrid_llsearch(cur_part, cur_grid-1);
        }
      
      /* remove 'sent-back' particles from list */
      free(cur_grid->leavers.send_back);
      cur_grid->leavers.send_back   = NULL;
      cur_grid->leavers.no_sendback = 0;
     }
}

/*===============================================================================
* keep on moving particles according to fine grid stepper
*===============================================================================*/
void advance_leavers(gridls *cur_grid, double DRIFT, double KICK)
{
   partptr cur_part;
   int     ileaver;
   
   /* keep track percentage of particles crossing grid boundaries */
   global.no_leavers += cur_grid->leavers.no_keepmoving;
   global.no_lsteps++;
   
   /* keep on moving late leavers */
   if(cur_grid->leavers.no_keepmoving > 0)
     {
      /* advance all leavers */
      for(ileaver = 0; ileaver < cur_grid->leavers.no_keepmoving; ileaver++)
        {
         cur_part = *(cur_grid->leavers.keep_moving+ileaver);
         
         if(cur_grid->multistep == 3)
           {
            kick_leaver (cur_part, KICK,  cur_grid-1);
            drift_leaver(cur_part, DRIFT, cur_grid-1);
           }
        }
      
      /* make sure all particles are at correct nodes on next coarser level */
      NULL_newll(cur_grid-1);
      move_part(cur_grid-1);
      
#ifdef TERMINATE2
      /* check, if everything is still right on (cur_grid-1) */
      if((cur_grid-1)->leavers.no_sendback   != 0 || 
         (cur_grid-1)->leavers.no_keepmoving != 0)
        {
         fprintf(io.logfile,
                 "advance_leavers: send-back after move_part(%ld) = %d\n",
                 (cur_grid-1)->l1dim,((cur_grid-1)->leavers.no_sendback));
         fprintf(io.logfile,
                 "advance_leavers: keep-moving after move_part(%ld) = %d\n",
                 (cur_grid-1)->l1dim,((cur_grid-1)->leavers.no_keepmoving));
         fflush(io.logfile);
         global.terminate = TRUE;
        }
#endif
      
      /* remove 'keep-moving' particles from list */
      free(cur_grid->leavers.keep_moving);
      cur_grid->leavers.keep_moving   = NULL;
      cur_grid->leavers.no_keepmoving = 0;
     }
}




