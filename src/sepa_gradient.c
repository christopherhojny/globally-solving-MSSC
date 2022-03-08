/**@file   sepa_gradient.c
 * @brief  gradient cut separator
 * @author Carina Costa
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <stdio.h>

#include "datapoints.h"
#include "getProbdata.h"
#include "sepa_gradient.h"


#define SEPA_NAME              "gradient"
#define SEPA_DESC              "gradient cut separator"
#define SEPA_PRIORITY              1000
#define SEPA_FREQ                    -1
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

#define DEFAULT_MAXNGEN              10 /**< maximum number of gradient cuts to be generated per round (-1: unlimited) */


/*
 * Data structures
 */

/** separator data */
struct SCIP_SepaData
{
   Datapoints*           datapoints;         /**< pointer to underlying data points */
   int                   nclusters;          /**< number of clusters */
   int                   ndatapoints;        /**< number of data points */
   int                   dimension;          /**< dimension of data points */
   SCIP_VAR***           xvars;              /**< (ndatapoints x nclusters)-array of assignment var */
   SCIP_VAR***           evars;              /**< (ndatapoints x nclusters)-array of epigraph var */
   SCIP_VAR***           cvars;              /**< (nclusters x dimension)-array of centroid variables */
   int                   maxngen;            /**< maximum number of gradient cuts to be generated per round */
};


/*
 * Local methods
 */

/** frees data of separator */
static
SCIP_RETCODE sepadataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA**       sepadata            /**< pointer to sepadata */
   )
{
   assert( scip != NULL );
   assert( sepadata != NULL );

   SCIPfreeBlockMemory(scip, sepadata);

   return SCIP_OKAY;
}


/** creates data structure of separator */
static
SCIP_RETCODE sepadataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA**       sepadata            /**< pointer to store separator data */
   )
{
   assert( scip != NULL );
   assert( sepadata != NULL );

   SCIP_CALL( SCIPallocBlockMemory(scip, sepadata) );

   return SCIP_OKAY;
}


/** set data structure of separator */
static
SCIP_RETCODE sepadataSet(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< pointer to store separator data */
   Datapoints*           datapoints,         /**< array of data points */
   int                   nclusters,          /**< number of clusters */
   int                   ndatapoints,        /**< number of data points in datapoints */
   int                   dimension,          /**< dimension of points in data points */
   SCIP_VAR***           xvars,              /**< array of point-cluster-assignment-variables */
   SCIP_VAR***           evars,              /**< (ndatapoints x nclusters)-array of epigraph var */
   SCIP_VAR***           cvars               /**< (nclusters x dimension)-array of centroid variables */
   )
{
   assert( scip != NULL );
   assert( sepadata != NULL );
   assert( datapoints != NULL );
   assert( nclusters > 0 );
   assert( ndatapoints > 0 );
   assert( dimension > 0 );
   assert( xvars != NULL );
   assert( cvars != NULL );

   sepadata->datapoints = datapoints;
   sepadata->nclusters = nclusters;
   sepadata->ndatapoints = ndatapoints;
   sepadata->dimension = dimension;
   sepadata->xvars = xvars;
   sepadata->evars = evars;
   sepadata->cvars = cvars;

   return SCIP_OKAY;
}


/** checks whether a given primal or LP solution is integral in xvars */
static
SCIP_RETCODE xvarsIntegral(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_SEPADATA*        sepadata,           /**< data of separator */
   SCIP_SOL*             sol,                /**< primal solution that should be separated, or NULL for LP solution */
   SCIP_Bool*            allintegral         /**< pointer to store whether all xvars are integral */
   )
{
   SCIP_VAR*** xvars;
   SCIP_Real solval;
   int nclusters;
   int ndatapoints;
   int i;
   int j;

   assert( scip != NULL );
   assert( sepadata != NULL );
   assert( allintegral != NULL );

   *allintegral = TRUE;

   xvars = sepadata->xvars;
   nclusters = sepadata->nclusters;
   ndatapoints = sepadata->ndatapoints;

   assert( xvars != NULL );
   assert( nclusters > 0 );
   assert( ndatapoints > 0 );

   // check if x is integral
   for (i = 0; i < ndatapoints; ++i)
   {
      for (j = 0; j < nclusters; ++j)
      {
         solval = SCIPgetSolVal(scip, sol, xvars[i][j]);

         if ( !SCIPisFeasIntegral(scip, solval) )
         {
            *allintegral = FALSE;
            return SCIP_OKAY;
         }
      }
   }

   return SCIP_OKAY;
}

/** separate given solution */
static
SCIP_RETCODE separate(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPADATA*        sepadata,           /**< data of separator */
   SCIP_SOL*             sol,                /**< primal solution that should be separated, or NULL for LP solution */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   Datapoints* datapoints;
   SCIP_VAR*** xvars;
   SCIP_VAR*** evars;
   SCIP_VAR*** cvars;
   SCIP_Real** points;
   SCIP_Real** createdxsol;
   SCIP_Real** barycenter;
   SCIP_Real** boxubs;
   SCIP_Real** boxlbs;
   SCIP_Real coordub;
   SCIP_Real coordlb;
   SCIP_Real** bigM;
   int nclusters;
   int ndatapoints;
   int dimension;
   int npointsincluster;
   int i;
   int j;
   int k;
   int l;
   int p;
   int maxngen;
   int ngen = 0;
   SCIP_Bool allintegral;

   assert( scip != NULL );
   assert( sepadata != NULL );
   assert( *result == SCIP_DIDNOTRUN );

   datapoints = sepadata->datapoints;
   xvars = sepadata->xvars;
   evars = sepadata->evars;
   cvars = sepadata->cvars;
   nclusters = sepadata->nclusters;
   ndatapoints = sepadata->ndatapoints;
   dimension = sepadata->dimension;
   points = datapoints->points;
   maxngen = sepadata->maxngen;

   assert( datapoints != NULL );
   assert( xvars != NULL );
   assert( evars != NULL );
   assert( cvars != NULL );
   assert( nclusters > 0 );
   assert( ndatapoints > 0 );
   assert( dimension > 0 );
   assert( points != NULL );

   SCIP_CALL( SCIPallocBufferArray(scip, &barycenter, nclusters) );
   for (j = 0; j < nclusters; ++j)
      SCIP_CALL( SCIPallocClearBufferArray(scip, &barycenter[j], dimension) );

   SCIP_CALL( SCIPallocBufferArray(scip, &bigM, ndatapoints) );
   for (i = 0; i < ndatapoints; ++i)
      SCIP_CALL( SCIPallocClearBufferArray(scip, &bigM[i], nclusters) );

   /* check if x is integral */
   SCIP_CALL( xvarsIntegral(scip, sepadata, sol, &allintegral) );

   if ( allintegral )
   {
      /* compute barycenters */
      for (j = 0; j < nclusters; ++j)
      {
         npointsincluster = 0;

         for (i = 0; i < ndatapoints; ++i)
         {
            /* check whether data point is in cluster */
            if ( SCIPgetSolVal(scip, sol, xvars[i][j]) < 0.5 )
               continue;

            /* update barycenter for each coordinate */
            for (l = 0; l < dimension; ++l)
               barycenter[j][l] = (points[i][l] + npointsincluster * barycenter[j][l]) / (npointsincluster + 1);
            ++npointsincluster;
         }
      }
   }
   else
   {
      /* guess an integral x-solution from the LP solution */
      SCIP_CALL( SCIPallocBufferArray(scip, &createdxsol, ndatapoints) );
      for (i = 0; i < ndatapoints; ++i)
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &createdxsol[i], nclusters) );
      }

      for (i = 0; i < ndatapoints; ++i)
      {
         int idxmax = -1;
         SCIP_Real maxval = -1.0;
         SCIP_Real curval;

         /* find cluster that gives the largest fractional entry */
         for (j = 0; j < nclusters; ++j)
         {
            curval = SCIPgetSolVal(scip, sol, xvars[i][j]);

            if ( SCIPisGT(scip, curval, maxval) )
            {
               maxval = curval;
               idxmax = j;
            }
         }

         assert( 0 <= idxmax && idxmax < nclusters );
         assert( SCIPisGE(scip, maxval, 0.0) && SCIPisLE(scip, maxval, 1.0) );

         /* assign point to the cluster and ensure this point is only in one cluster */
         createdxsol[i][idxmax] = 1;

         for (k = 0; k < nclusters; ++k)
         {
            if ( k == idxmax )
               continue;
            createdxsol[i][k] = 0;
         }
      }

      /* compute barycenters */
      for (j = 0; j < nclusters; ++j)
      {
         npointsincluster = 0;
         for (i = 0; i < ndatapoints; ++i)
         {
            if ( createdxsol[i][j] == 1 )
            {
               /* update barycenter for each coordinate */
               for (l = 0; l < dimension; ++l)
                  barycenter[j][l] = (points[i][l] + npointsincluster * barycenter[j][l]) / (npointsincluster + 1);
               ++npointsincluster;
            }
            else
               continue;
         }
      }

      for (i = ndatapoints - 1; i >= 0; --i)
      {
         SCIPfreeBufferArray(scip, &createdxsol[i]);
      }
      SCIPfreeBufferArray(scip, &createdxsol);
   }

   /* collect maximum coordinate for a cluster k in each dimension */
   SCIP_CALL( SCIPallocBufferArray(scip, &boxubs, nclusters) );
   for (k = 0; k < nclusters; ++k)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &boxubs[k], dimension) );
      for (i = 0; i < dimension; ++i)
      {
         coordub = -SCIPinfinity(scip);

         for (p = 0; p < ndatapoints; ++p)
         {
            /* skip variables that cannot be part of cluster k */
            if ( SCIPvarGetUbLocal(xvars[p][k]) < 0.5 )
               continue;

            if (points[p][i] > coordub)
               coordub = points[p][i];
         }

         if ( !SCIPisEQ(scip, coordub, -SCIPinfinity(scip)) && SCIPisLT(scip, coordub, SCIPvarGetUbLocal(cvars[k][i])) )
         {
            boxubs[k][i] = coordub;
         }
         else
            boxubs[k][i] = SCIPvarGetUbLocal(cvars[k][i]);
      }
   }

   /* collect minimum coordinate for a cluster k in each dimension */
   SCIP_CALL( SCIPallocBufferArray(scip, &boxlbs, nclusters) );
   for (k = 0; k < nclusters; ++k)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &boxlbs[k], dimension) );
      for (i = 0; i < dimension; ++i)
      {
         coordlb = SCIPinfinity(scip);

         for (p = 0; p < ndatapoints; ++p)
         {
            /* skip variables that cannot be part of cluster k */
            if ( SCIPvarGetUbLocal(xvars[p][k]) < 0.5 )
               continue;

            if (points[p][i] < coordlb)
               coordlb = points[p][i];
         }

         if ( !SCIPisEQ(scip, coordlb, SCIPinfinity(scip)) && SCIPisGT(scip, coordlb, SCIPvarGetLbLocal(cvars[k][i])) )
         {
            boxlbs[k][i] = coordlb;
         }
         else
            boxlbs[k][i] = SCIPvarGetLbLocal(cvars[k][i]);
      }
   }

   /* compute maximum squared distance of point i to bounding box of cluster j */
   for (i = 0; i < ndatapoints; ++i)
   {
      for (j = 0; j < nclusters; ++j)
      {
         for (l = 0; l < dimension; ++l)
         {
            SCIP_Real lowerdist;
            SCIP_Real upperdist;

            lowerdist = (points[i][l] - boxlbs[j][l]) * (points[i][l] - boxlbs[j][l]);
            upperdist = (boxubs[j][l] - points[i][l]) * (boxubs[j][l] - points[i][l]);

            if ( SCIPisLT(scip, lowerdist, upperdist) )
               bigM[i][j] += upperdist;
            else
               bigM[i][j] += lowerdist;
         }
      }
   }

   for (k = nclusters - 1; k >= 0; --k)
   {
      SCIPfreeBufferArray(scip, &boxlbs[k]);
   }
   SCIPfreeBufferArray(scip, &boxlbs);
   for (k = nclusters - 1; k >= 0; --k)
   {
      SCIPfreeBufferArray(scip, &boxubs[k]);
   }
   SCIPfreeBufferArray(scip, &boxubs);


   /* generate gradient cuts */
   for (i = 0; i < ndatapoints; ++i)
   {
      for (j = 0; j < nclusters; ++j)
      {
         /* pointer to store cut */
         SCIP_ROW* cut;

         SCIP_Real cutlhs = 0.0;

         cutlhs -= bigM[i][j];

         for (l = 0; l < dimension; ++l)
            cutlhs += points[i][l] * points[i][l] - barycenter[j][l] * barycenter[j][l];

         SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut, sepa, "gradient", cutlhs, SCIPinfinity(scip), FALSE, FALSE, TRUE) );
         SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );

         SCIP_CALL( SCIPaddVarToRow(scip, cut, xvars[i][j], - bigM[i][j]) );
         SCIP_CALL( SCIPaddVarToRow(scip, cut, evars[i][j], 1.0) );

         for (l = 0; l < dimension; ++l)
         {
            SCIP_CALL( SCIPaddVarToRow(scip, cut, cvars[j][l], 2 * points[i][l] - 2 * barycenter[j][l]) );
         }

         SCIP_CALL( SCIPflushRowExtensions(scip, cut) );

         /* checks whether cut is sufficiently violated with respect to the given solution */
         if ( SCIPisCutEfficacious(scip, sol, cut) )
         {
            SCIP_Bool infeasible;

            SCIP_CALL( SCIPaddRow(scip, cut, FALSE, &infeasible) );

            if ( infeasible )
            {
               *result = SCIP_CUTOFF;
                SCIP_CALL( SCIPreleaseRow(scip, &cut) );
                goto FREEDATA;
            }

            *result = SCIP_SEPARATED;
            ++ngen;

            /* if not already existing, adds row to global cut pool */
            SCIP_CALL( SCIPaddPoolCut(scip, cut) );
         }

         SCIP_CALL( SCIPreleaseRow(scip, &cut) );

         /* possibly terminate early */
         if ( maxngen != -1 && ngen >= maxngen )
            goto FREEDATA;
      }
   }

 FREEDATA:
   for (i = ndatapoints -1; i >= 0; --i)
   {
      SCIPfreeBufferArray(scip, &bigM[i]);
   }
   SCIPfreeBufferArray(scip, &bigM);

   for (j = nclusters - 1; j >= 0; --j)
   {
      SCIPfreeBufferArray(scip, &barycenter[j]);
   }
   SCIPfreeBufferArray(scip, &barycenter);

   return SCIP_OKAY;
}

/*
 * Callback methods of separator
 */

/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeGradient)
{
   SCIP_SEPADATA* sepadata;

   assert( scip != NULL );
   assert( sepa != NULL );

   sepadata = SCIPsepaGetData(sepa);
   assert( sepadata != NULL );

   SCIP_CALL( sepadataFree(scip, &sepadata) );

   return SCIP_OKAY;
}


/** initialization method of separator (called after problem was transformed) */
static
SCIP_DECL_SEPAINIT(sepaInitGradient)
{  /*lint --e{715}*/
   SCIP_PROBDATA* probdata;
   SCIP_SEPADATA* sepadata;
   Datapoints* datapoints;
   int nclusters;
   int ndatapoints;
   int dimension;
   SCIP_VAR*** xvars;
   SCIP_VAR*** evars;
   SCIP_VAR*** cvars;
   int model;

   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   sepadata = SCIPsepaGetData(sepa);
   assert( sepadata != NULL );

   SCIP_CALL( SCIPgetIntParam(scip, "clustering/model", &model) );

   datapoints = SCIPprobdataGetDatapoints(probdata, model);
   nclusters = SCIPprobdataGetNClusters(probdata, model);
   ndatapoints = SCIPprobdataGetNDatapoints(probdata, model);
   dimension = SCIPprobdataGetDimension(probdata, model);
   xvars = SCIPprobdataGetXvars(probdata, model);
   cvars = SCIPprobdataGetCvars(probdata, model);
   evars = SCIPprobdataGetEvars(probdata, model);

   /* create convexity propagator data */
   SCIP_CALL( sepadataSet(scip, sepadata, datapoints, nclusters, ndatapoints, dimension, xvars, evars, cvars) );

   return SCIP_OKAY;
}


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpGradient)
{
   SCIP_SEPADATA* sepadata;

   assert( scip != NULL );
   assert( result != NULL );
   assert( sepa != NULL );

   *result = SCIP_DIDNOTRUN;

   sepadata = SCIPsepaGetData(sepa);
   assert( sepadata != NULL );

   /* separate cuts on the LP solution */
   SCIP_CALL( separate(scip, sepa, sepadata, NULL, result) );

   return SCIP_OKAY;
}


/** arbitrary primal solution separation method of separator */
static
SCIP_DECL_SEPAEXECSOL(sepaExecsolGradient)
{
   SCIP_SEPADATA* sepadata;

   assert( scip != NULL );
   assert( result != NULL );
   assert( sepa != NULL );

   sepadata = SCIPsepaGetData(sepa);
   assert( sepadata != NULL );

   *result = SCIP_DIDNOTRUN;

   /* separate cuts on the given primal solution */
   SCIP_CALL( separate(scip, sepa, sepadata, sol, result) );

   return SCIP_OKAY;
}


/*
 * separator specific interface methods
 */

/** creates the gradient cut separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaGradient(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_SEPA* sepa;

   /* create gradient separator data */
   SCIP_CALL( sepadataCreate(scip, &sepadata) );

   /* include separator */
   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY,
         sepaExeclpGradient, sepaExecsolGradient,
         sepadata) );

   assert(sepa != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeGradient) );
   SCIP_CALL( SCIPsetSepaInit(scip, sepa, sepaInitGradient) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/maxngen",
         "maximum number of cuts to be generated per separation round",
         &sepadata->maxngen, TRUE, DEFAULT_MAXNGEN, -1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
