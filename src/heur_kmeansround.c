/**@file    heur_kmeansround.c
 * @ingroup DEFPLUGINS_HEUR
 * @brief   kmeans rounding primal heuristic
 * @author  Carina Moreira Costa
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <stdio.h>

#include "heur_kmeansround.h"
#include "datapoints.h"
#include "getProbdata.h"
#include "typedefs.h"

#define HEUR_NAME             "kmeansround"
#define HEUR_DESC             "LP rounding heuristic for kmeans clustering problem"
#define HEUR_DISPCHAR         'Y'
#define HEUR_PRIORITY         0
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERLPNODE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */


/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   /* data of the problem */
   SCIP_Real**           points;             /**< (ndatapoints x dimension)-array of data points */
   int                   nclusters;          /**< number of clusters */
   int                   ndatapoints;        /**< number of data points */
   int                   dimension;          /**< dimension of data points */

   /* data of the model */
   SCIP_VAR*             objvar;             /**< variable to linearize objective */
   SCIP_VAR***           evars;              /**< (ndatapoints x nclusters)-array of objective variables */
   SCIP_VAR***           xvars;              /**< (ndatapoints x nclusters)-array of assignment variables */
   SCIP_VAR***           cvars;              /**< (nclusters x dimension)-array of centroid variables */
   SCIP_VAR***           quadcvars;          /**< (nclusters x dimension)-array of squared centroid variables */
   SCIP_VAR**            kappavars;          /**< ncluster-array to encode cardinality of clusters */
};


/*
 * Local methods
 */


/*
 * Callback methods of primal heuristic
 */


 /** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeKmeansround)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( scip != NULL );
   assert( heur != NULL );

   heurdata = SCIPheurGetData(heur);

   if ( heurdata != NULL )
   {
      SCIPfreeBlockMemory(scip, &heurdata);
   }

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitKmeansround)
{  /*lint --e{715}*/
   SCIP_PROBDATA* probdata;
   SCIP_HEURDATA* heurdata;
   Datapoints* datapoints;
   int model;

   assert( scip != NULL );
   assert( heur != NULL );

   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   SCIP_CALL( SCIPgetIntParam(scip, "clustering/model", &model) );

   datapoints = SCIPprobdataGetDatapoints(probdata, model);
   heurdata->nclusters = SCIPprobdataGetNClusters(probdata, model);
   heurdata->ndatapoints = datapoints->ndatapoints;
   heurdata->dimension = datapoints->dimension;
   heurdata->points = datapoints->points;
   heurdata->objvar = SCIPprobdataGetObjvar(probdata, model);
   heurdata->evars = SCIPprobdataGetEvars(probdata, model);
   heurdata->xvars = SCIPprobdataGetXvars(probdata, model);
   heurdata->cvars = SCIPprobdataGetCvars(probdata, model);
   heurdata->quadcvars = SCIPprobdataGetQuadCvars(probdata, model);
   heurdata->kappavars = SCIPprobdataGetKappaVars(probdata, model);

   return SCIP_OKAY;
}


/** round LP solution to the closest feasible binary solution */
static
SCIP_RETCODE roundSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic */
   SCIP_HEURDATA*        heurdata,           /**< data of heuristic */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   SCIP_Real** points;
   SCIP_Real** barycenter;
   SCIP_VAR* objvar;
   SCIP_VAR*** evars;
   SCIP_VAR*** xvars;
   SCIP_VAR*** cvars;
   SCIP_VAR*** quadcvars;
   SCIP_VAR** kappavars;
   SCIP_SOL* sol;
   int nclusters;
   int ndatapoints;
   int dimension;
   int j;
   int i;
   int l;
   int model;
   int maxidx;
   int* npointsincluster;
   SCIP_Real maxval;
   SCIP_Real curval;
   SCIP_Real eval;
   SCIP_Real objval = 0.0;
   SCIP_Bool* ispointassigned;
   SCIP_Bool stored;
   SCIP_Bool uselocalizedcut;

   assert( scip != NULL );
   assert( heur != NULL );
   assert( heurdata != NULL );

   points = heurdata->points;
   objvar = heurdata->objvar;
   evars = heurdata->evars;
   xvars = heurdata->xvars;
   cvars = heurdata->cvars;
   quadcvars = heurdata->quadcvars;
   kappavars = heurdata->kappavars;
   nclusters = heurdata->nclusters;
   ndatapoints = heurdata->ndatapoints;
   dimension = heurdata->dimension;

   assert( points != NULL );
   assert( objvar != NULL || evars != NULL );
   assert( xvars != NULL );
   assert( cvars != NULL );
   assert( kappavars != NULL );
   assert( nclusters > 0 );
   assert( ndatapoints > 0 );
   assert( dimension > 0 );

   /* copy the current LP solution to the working solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );
   SCIP_CALL( SCIPlinkLPSol(scip, sol) );

   SCIP_CALL( SCIPallocClearBufferArray(scip, &ispointassigned, ndatapoints) );

   SCIP_CALL( SCIPallocClearBufferArray(scip, &npointsincluster, nclusters) );

   SCIP_CALL( SCIPallocBufferArray(scip, &barycenter, nclusters) );
   for (j = 0; j < nclusters; ++j)
      SCIP_CALL( SCIPallocClearBufferArray(scip, &barycenter[j], dimension) );

   /* ensure that all clusters will have at least one point each */
   for (j = 0; j < nclusters; ++j)
   {
      maxidx = -1;
      maxval = -1.0;

      /* find value and index of maximum x-entry */
      for (i = 0; i < ndatapoints; ++i)
      {
         if ( ispointassigned[i] )
            continue;

         curval = SCIPgetSolVal(scip, sol, xvars[i][j]);
         if ( SCIPisGT(scip, curval, maxval) )
         {
            maxval = curval;
            maxidx = i;
         }
      }

      assert( 0 <= maxidx && maxidx < ndatapoints );

      /* stop if point cannot be assigned to cluster*/
      if ( SCIPvarGetUbLocal(xvars[maxidx][j]) < 0.5 )
      {
         *result = SCIP_DIDNOTFIND;
         SCIP_CALL( SCIPfreeSol(scip, &sol) );
         goto FREEMEMEMORY;
      }

      /* assign point with largest fractional entry to cluster j and ensure this point is only in one cluster */
      SCIP_CALL( SCIPsetSolVal(scip, sol, xvars[maxidx][j], 1) );
      ispointassigned[maxidx] = TRUE;

      for (int c = 0; c < nclusters; ++c)
      {
         if ( c == j )
            continue;

         SCIP_CALL( SCIPsetSolVal(scip, sol, xvars[maxidx][c], 0) );
      }

      /* update barycenter, which will be exactly the point chosen above */
      for (l = 0; l < dimension; ++l)
         barycenter[j][l] = points[maxidx][l];
      ++npointsincluster[j];
   }

   /* round LP solution to the nearest binary one */
   for (i = 0; i < ndatapoints; ++i)
   {
      if ( ispointassigned[i] )
         continue;

      maxidx = -1;
      maxval = -1.0;

      for (j = 0; j < nclusters; ++j)
      {
         curval = SCIPgetSolVal(scip, sol, xvars[i][j]);
         if ( SCIPisGT(scip, curval, maxval) )
         {
            maxval = curval;
            maxidx = j;
         }
      }

      assert( 0 <= maxidx && maxidx < nclusters );

      /* stop if point cannot be assigned to cluster */
      if ( SCIPvarGetUbLocal(xvars[i][maxidx]) < 0.5 )
      {
         *result = SCIP_DIDNOTFIND;
         SCIP_CALL( SCIPfreeSol(scip, &sol) );
         goto FREEMEMEMORY;
      }

      SCIP_CALL( SCIPsetSolVal(scip, sol, xvars[i][maxidx], 1) );
      ispointassigned[i] = TRUE;

      for (j = 0; j < nclusters; ++j)
      {
         if ( j == maxidx )
            continue;

         SCIP_CALL( SCIPsetSolVal(scip, sol, xvars[i][j], 0) );
      }

      /* update barycenter */
      for (l = 0; l < dimension; ++l)
         barycenter[maxidx][l] = (points[i][l] + npointsincluster[maxidx] * barycenter[maxidx][l]) / (npointsincluster[maxidx] + 1);
      ++npointsincluster[maxidx];
   }

   /* check which model is used to update variables correctly */
   SCIP_CALL( SCIPgetIntParam(scip, "clustering/model", &model) );

   if ( model == MODEL_QUADRATIC )
   {
      assert( quadcvars != NULL );
      assert( objvar != NULL );
   }
   else
      assert( evars != NULL );

   /* update cvals */
   for (j = 0; j < nclusters; ++j)
   {
      for (l = 0; l < dimension; ++l)
      {
         SCIP_CALL( SCIPsetSolVal(scip, sol, cvars[j][l], barycenter[j][l]) );

         if ( model == MODEL_QUADRATIC )
         {
            SCIP_CALL( SCIPsetSolVal(scip, sol, quadcvars[j][l], barycenter[j][l] * barycenter[j][l]) );
         }
      }
   }

   /* compute obj val or evals given the rounded solution */
   for (int c = 0; c < nclusters; ++c)
   {
      for (i = 0; i < ndatapoints; ++i)
      {
         if ( SCIPgetSolVal(scip, sol, xvars[i][c]) == 0 )
            continue;

         eval = 0.0;

         for (l = 0; l < dimension; ++l)
            eval += points[i][l] * points[i][l] - 2.0 * points[i][l] * barycenter[c][l] + barycenter[c][l] * barycenter[c][l];

         if ( model == MODEL_QUADRATIC )
            objval += eval;

         else
         {
            SCIP_CALL( SCIPsetSolVal(scip, sol, evars[i][c], eval) );

            for (j = 0; j < nclusters; ++j)
            {
               if ( j == c )
                  continue;
               else
               {
                  SCIP_CALL( SCIPsetSolVal(scip, sol, evars[i][j], 0.0) );
               }
            }
         }
      }
   }

   /* set kappa vars */
   SCIP_CALL( SCIPgetBoolParam(scip, "clustering/uselocalizedcardcut", &uselocalizedcut) );

   // set kappavars
   if ( uselocalizedcut )
   {
      int ncomplement;
      for (j = 0; j < nclusters; ++j)
      {
         ncomplement = 0;

         for (int c = 0; c < nclusters; ++c)
         {
            if ( c == j )
               continue;

            for (i = 0; i < ndatapoints; ++i)
            {
               if ( SCIPgetSolVal(scip, sol, xvars[i][c]) == 1 )
                  ++ncomplement;
            }
         }
         SCIP_CALL( SCIPsetSolVal(scip, sol, kappavars[j], ncomplement) );
      }
   }
   else
   {
      for (j = 0; j < nclusters; ++j)
      {
         SCIP_CALL( SCIPsetSolVal(scip, sol, kappavars[j], 0.0) );
      }
   }

   if ( model == MODEL_QUADRATIC )
      SCIP_CALL( SCIPsetSolVal(scip, sol, objvar, objval) );

   SCIP_CALL( SCIPtrySol(scip, sol, FALSE, FALSE, TRUE, TRUE, TRUE, &stored) );

   if ( stored )
      *result = SCIP_FOUNDSOL;
   else
      *result = SCIP_DIDNOTFIND;

   SCIP_CALL( SCIPfreeSol(scip, &sol) );

 FREEMEMEMORY:
   for (j = nclusters - 1; j >= 0; --j)
   {
      SCIPfreeBufferArray(scip, &barycenter[j]);
   }
   SCIPfreeBufferArray(scip, &barycenter);

   SCIPfreeBufferArray(scip, &npointsincluster);
   SCIPfreeBufferArray(scip, &ispointassigned);

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecKmeansround)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandssol;
   int nfrac;
   int nlpcands;

   assert( scip != NULL );
   assert( heur != NULL );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   *result = SCIP_DIDNOTFIND;

   /* Check that everything is there */
   assert( scip );
   assert( result != NULL );
   assert( SCIPhasCurrentNodeLP(scip) );

   *result = SCIP_DIDNOTRUN;

   /* only call heuristic if an optimal LP solution is at hand */
   if ( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only call heuristic if the LP objective value is smaller than the cutoff bound */
   if ( SCIPisGE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)) )
      return SCIP_OKAY;

   /* get fractional variables, that should be binary */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, NULL, &nlpcands, NULL, NULL) );
   nfrac = nlpcands;

   /* only call heuristic if LP solution is fractional */
   if ( nfrac == 0 )
      return SCIP_OKAY;

   *result= SCIP_DIDNOTFIND;

   /* construct feasible solution by rounding LP solution */
   SCIP_CALL( roundSolution(scip, heur, heurdata, result) );

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the improvement primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurKmeansround(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create heuristic data */
   heurdata = NULL;
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   heurdata->nclusters = 0;
   heurdata->ndatapoints = 0;
   heurdata->dimension = 0;
   heurdata->xvars = NULL;
   heurdata->cvars = NULL;
   heurdata->objvar = NULL;
   heurdata->evars = NULL;
   heurdata->quadcvars = NULL;
   heurdata->kappavars = NULL;

   heur = NULL;

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecKmeansround, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeKmeansround) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitKmeansround) );

   return SCIP_OKAY;
}
