/**@file   heur_improve.c
 * @ingroup DEFPLUGINS_HEUR
 * @brief  improvement primal heuristic
 * @author Carina Moreira Costa
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "datapoints.h"
#include "getProbdata.h"
#include "heur_improve.h"
#include "typedefs.h"
#include "math.h"


#define HEUR_NAME             "improve"
#define HEUR_DESC             "improvement heuristic"
#define HEUR_DISPCHAR         'X'
#define HEUR_PRIORITY         0
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */


/*
 * Data structures
 */


/** primal heuristic data */
struct SCIP_HeurData
{
   Datapoints*           datapoints;         /**< pointer to underlying data points */
   SCIP_Real**           points;             /**< (ndatapoints x dimension)-array of data points */
   int                   nclusters;          /**< number of clusters */
   int                   ndatapoints;        /**< number of data points */
   int                   dimension;          /**< dimension of data points */
   SCIP_VAR*             objvar;             /**< variable to linearize objective */
   SCIP_VAR***           evars;              /**< (ndatapoints x nclusters)-array of objective variables */
   SCIP_VAR***           xvars;              /**< (ndatapoints x nclusters)-array of assignment variables */
   SCIP_VAR***           cvars;              /**< (nclusters x dimension)-array of centroid variables */
   SCIP_VAR***           quadcvars;          /**< (nclusters x dimension)-array of squared centroid variables */
   SCIP_VAR**            kappavars;          /**< ncluster-array to encode cardinality of clusters */
   int                   bestsolidx;         /**< best solution during the previous run */
};

/*
 * Local methods
 */

/** computes candidate clusters to be changed */
static
SCIP_Bool computeCandidates(
   SCIP*                 scip,               /**< SCIP data structure */
   Datapoints*           datapoints,         /**< array of data points */
   int                   nclusters,          /**< number of clusters */
   int*                  pointincluster,     /**< ndatapoints-array to inform cluster in which the point is assigned to */
   int*                  npointsincluster,   /**< nclusters-array to inform number of points per cluster */
   SCIP_Real*            weightedvariances,  /**< nclusters-array of weighted variances */
   SCIP_Real**           barycenter          /**< (nclusters x dim)-array of centroids */
   )
{
   SCIP_Real** newbarycenter;
   SCIP_Real* tmpbarycenter;
   SCIP_Real* ratios;
   SCIP_Real tmpweightedvariance = 0.0;
   SCIP_Real currentratio;
   SCIP_Real maxdistance;
   SCIP_Real dist1;
   SCIP_Real dist2;
   SCIP_Real tolupdatecentroids = 1e-2;
   SCIP_Real infinitynormdif = 1;
   int maxNrIterations = 10;
   int iterationscounter;
   SCIP_Bool* isclusterupdated;
   SCIP_Bool foundcandidatesol = FALSE;
   int* clustertojoin;
   int* clustertosplit;
   int i;
   int j1;
   int l;
   int j;
   int ncandidatechanges = 0;
   int point1;
   int point2;
   int npointsinclusterc1;

   assert( scip != NULL );
   assert( datapoints != NULL );
   assert( nclusters > 0 );
   assert( pointincluster != NULL );
   assert( npointsincluster != NULL );
   assert( weightedvariances != NULL );
   assert( barycenter != NULL );

   SCIP_CALL( SCIPallocClearBufferArray(scip, &tmpbarycenter, datapoints->dimension) );
   SCIP_CALL( SCIPallocBufferArray(scip, &newbarycenter, nclusters) );
   for (j = 0; j < nclusters; ++j)
   {
      SCIP_CALL( SCIPallocClearBufferArray(scip, &newbarycenter[j], datapoints->dimension) );
   }
   SCIP_CALL( SCIPallocBufferArray(scip, &ratios, nclusters - 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &clustertojoin, nclusters - 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &clustertosplit, nclusters - 1) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &isclusterupdated, nclusters) );

   /* join each two clusters, compute its new barycenter and new weighted variance */
   for (j1 = 0; j1 < nclusters - 1; ++j1)
   {
      int j2;

      /* initialize ratio (between weighted joint variance and the weighted variance from a third cluster) */
      ratios[j1] = 1.0;
      for (j2 = j1 + 1; j2 < nclusters; ++j2)
      {
         int j3;

         /* compute temporary barycenter by joining clusters j1 and j2 */
         for (i = 0; i < datapoints->ndatapoints; ++i)
         {
            /* check if point i is in cluster j1 or j2 */
            if ( pointincluster[i] == j1 || pointincluster[i] == j2 )
            {
               for (l = 0; l < datapoints->dimension; ++l)
                  tmpbarycenter[l] += datapoints->points[i][l];
            }
         }
         for (l = 0; l < datapoints->dimension; ++l)
            tmpbarycenter[l] = tmpbarycenter[l] / (npointsincluster[j1] + npointsincluster[j2]);

         /* compute temporary joint weighted variance */
         for (i = 0; i < datapoints->ndatapoints; ++i)
         {
            if ( pointincluster[i] == j1 || pointincluster[i] == j2 )
            {
               SCIP_Real distpointclusterval = 0.0;
               for (l = 0; l < datapoints->dimension; ++l)
                  distpointclusterval += datapoints->points[i][l] * datapoints->points[i][l] - 2.0 * datapoints->points[i][l] * tmpbarycenter[l] + tmpbarycenter[l] * tmpbarycenter[l];
               tmpweightedvariance += distpointclusterval;
            }
         }
         tmpweightedvariance = tmpweightedvariance /(npointsincluster[j1] + npointsincluster[j2]);

         for (j3 = 0; j3 < nclusters; ++j3)
         {
            if ( j3 == j1 || j3 == j2 )
               continue;

            /* check if joint weighted variance is smaller than the weighted variance of a third cluster */
            currentratio = tmpweightedvariance / weightedvariances[j3];
            if ( SCIPisLT(scip, currentratio, ratios[j1]) )
            {
               ratios[j1] = currentratio;
               clustertosplit[j1] = j3;
               clustertojoin[j1] = j2;

               /* store temporarily */
               for (l = 0; l < datapoints->dimension; ++l)
                  newbarycenter[j1][l] = tmpbarycenter[l];
            }
         }

         tmpweightedvariance = 0.0;
         for (l = 0; l < datapoints->dimension; ++l)
            tmpbarycenter[l] = 0.0;
      }

      /* check if there were canditate clusters to join j1 and to be splitted */
      if ( ratios[j1] != 1.0 )
         ++ncandidatechanges;
   }

   if ( ncandidatechanges == 0 )
   {
      goto FREEMEMORYHERE;
   }

   foundcandidatesol = TRUE;

   /* choose smallest ratio and update barycenters and assignments */
   while ( ncandidatechanges != 0 )
   {
      SCIP_Real minratio = 1.0;
      int c1 = -1;
      int c2 = -1;
      int c3 = -1;
      for (j = 0; j < nclusters - 1; ++j)
      {
         if ( SCIPisLT(scip, ratios[j], minratio) && !isclusterupdated[j] )
         {
            if ( isclusterupdated[clustertojoin[j]] || isclusterupdated[clustertosplit[j]] )
            {
               --ncandidatechanges;
               continue;
            }
            minratio = ratios[j];
            c1 = j;
            c2 = clustertojoin[j];
            c3 = clustertosplit[j];
         }
         else
            continue;
      }

      /* check if while loop continues */
      if ( minratio == 1.0 )
      {
         goto FREEMEMORYHERE;
      }

      /* change will be performed */
      --ncandidatechanges;

      assert( 0 <= c1 && c1 < nclusters && !isclusterupdated[c1] );
      isclusterupdated[c1] = TRUE;

      assert( 0 <= c2 && c2 < nclusters && !isclusterupdated[c2] );
      isclusterupdated[c2] = TRUE;

      assert( 0 <= c3 && c3 < nclusters && !isclusterupdated[c3] );
      isclusterupdated[c3] = TRUE;

     /* update barycenter cluster c1 */
     for (l = 0; l < datapoints->dimension; ++l)
        barycenter[c1][l] = newbarycenter[c1][l];

     /* update points in cluster c1 */
     npointsinclusterc1 = 0;
     for (i = 0; i < datapoints->ndatapoints; ++i)
     {
        if ( pointincluster[i] == c1 || pointincluster[i] == c2 )
        {
           pointincluster[i] = c1;
           ++npointsinclusterc1;
        }
     }
     assert( npointsinclusterc1 == npointsincluster[c1] + npointsincluster[c2] );
     npointsincluster[c1] = npointsinclusterc1;

     /* compute new barycenters of c2 and c3 */
     /* start with the two furthest points in the third cluster to be the new barycenters */
      maxdistance = -1.0;
      point1 = -1;
      point2 = -1;

      for (i = 0; i < datapoints->ndatapoints - 1; ++i)
      {
         int p;
         if ( pointincluster[i] != c3 )
            continue;
         for (p = i + 1; p < datapoints->ndatapoints; ++p)
         {
            if ( pointincluster[p] == c3 )
            {
               SCIP_Real currentdist = 0;
               for (l = 0; l < datapoints->dimension; ++l)
                  currentdist += datapoints->points[i][l] * datapoints->points[i][l] - 2 * datapoints->points[i][l] *
                                 datapoints->points[p][l] + datapoints->points[p][l] * datapoints->points[p][l];
               if ( SCIPisGT(scip, currentdist, maxdistance) )
               {
                  point1 = i;
                  point2 = p;
                  maxdistance = currentdist;
               }
            }
         }
      }

      assert( maxdistance > 0 && point1 >= 0 && point2 >= 0 );

      /* if cluster to be splitted is cluster 0, ensure that point 0 remains in cluster 0 due to symmetry breaking */
      if ( c3 == 0 )
      {
         dist1 = 0;
         dist2 = 0;
         for (l = 0; l < datapoints->dimension; ++l)
         {
            dist1 += datapoints->points[0][l] * datapoints->points[0][l] - 2 * datapoints->points[0][l] *
                     datapoints->points[point1][l] + datapoints->points[point1][l] * datapoints->points[point1][l];
            dist2 += datapoints->points[0][l] * datapoints->points[0][l] - 2 * datapoints->points[0][l] *
                     datapoints->points[point2][l] + datapoints->points[point2][l] * datapoints->points[point2][l];
         }
         if ( SCIPisLT(scip, dist1, dist2) )
         {
            for (l = 0; l < datapoints->dimension; ++l)
            {
               barycenter[c3][l] = 0.5 * (datapoints->points[point1][l] + datapoints->points[0][l]);
               barycenter[c2][l] = datapoints->points[point2][l];
            }
         }
         else
         {
            for (l = 0; l < datapoints->dimension; ++l)
            {
               barycenter[c3][l] = 0.5 * (datapoints->points[point2][l] + datapoints->points[0][l]);
               barycenter[c2][l] = datapoints->points[point1][l];
            }
         }
      }
      else
      {
         for (l = 0; l < datapoints->dimension; ++l)
         {
            barycenter[c2][l] = datapoints->points[point1][l];
            barycenter[c3][l] = datapoints->points[point2][l];
         }
      }

      /* iterate over points in third cluster to see which barycenter is closest, and update barycenter */
      iterationscounter = 0;
      while ( infinitynormdif > tolupdatecentroids && maxNrIterations > iterationscounter )
      {
         infinitynormdif = -1;
         npointsincluster[c2] = 0;
         npointsincluster[c3] = 0;
         for (i = 0; i < datapoints->ndatapoints; ++i)
         {
            if ( pointincluster[i] != c3 && pointincluster[i] != c2 )
               continue;

            dist1 = 0;
            dist2 = 0;
            for (l = 0; l < datapoints->dimension; ++l)
            {
               dist1 += datapoints->points[i][l] * datapoints->points[i][l] - 2 * datapoints->points[i][l] *
                  barycenter[c2][l] + barycenter[c2][l] * barycenter[c2][l];
               dist2 += datapoints->points[i][l] * datapoints->points[i][l] - 2 * datapoints->points[i][l] *
                  barycenter[c3][l] + barycenter[c3][l] * barycenter[c3][l];
            }

            if ( SCIPisLT(scip, dist1, dist2) )
            {
               pointincluster[i] = c2;
               for (l = 0; l < datapoints->dimension; ++l)
                  newbarycenter[c2][l] = (datapoints->points[i][l] + npointsincluster[c2] * newbarycenter[c2][l]) / (npointsincluster[c2] + 1);
               ++npointsincluster[c2];
            }
            else if ( SCIPisLE(scip, dist2, dist1) )
            {
               pointincluster[i] = c3;
               for (l = 0; l < datapoints->dimension; ++l)
                  newbarycenter[c3][l] = (datapoints->points[i][l] + npointsincluster[c3] * newbarycenter[c3][l]) / (npointsincluster[c3] + 1);
               ++npointsincluster[c3];
            }
         }
         /* compute difference */
         for (l = 0; l < datapoints->dimension; ++l)
         {
            dist1 = fabs(newbarycenter[c2][l] - barycenter[c2][l]);
            dist2 = fabs(newbarycenter[c3][l] - barycenter[c3][l]);

            if ( SCIPisGT(scip, dist1, dist2) && SCIPisGT(scip, dist1, infinitynormdif) )
               infinitynormdif = dist1;
            else if ( SCIPisGE(scip, dist2, dist1) && SCIPisGT(scip, dist2, infinitynormdif) )
               infinitynormdif = dist2;
         }

         /* update barycenter */
         for (l = 0; l < datapoints->dimension; ++l)
         {
            barycenter[c2][l] = newbarycenter[c2][l];
            barycenter[c3][l] = newbarycenter[c3][l];
         }
         ++iterationscounter;
      }
   }

 FREEMEMORYHERE:
   SCIPfreeBufferArray(scip, &isclusterupdated);
   SCIPfreeBufferArray(scip, &clustertosplit);
   SCIPfreeBufferArray(scip, &clustertojoin);
   SCIPfreeBufferArray(scip, &ratios);
   for (j = nclusters - 1; j >= 0; --j)
   {
      SCIPfreeBufferArray(scip, &newbarycenter[j]);
   }
   SCIPfreeBufferArray(scip, &newbarycenter);
   SCIPfreeBufferArray(scip, &tmpbarycenter);

   return foundcandidatesol;
}


/*
 * Callback methods of primal heuristic
 */

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeImprove)
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
SCIP_DECL_HEURINIT(heurInitImprove)
{  /*lint --e{715}*/
   SCIP_PROBDATA* probdata;
   SCIP_HEURDATA* heurdata;
   int model;

   assert( scip != NULL );
   assert( heur != NULL );

   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   SCIP_CALL( SCIPgetIntParam(scip, "clustering/model", &model) );

   heurdata->datapoints = SCIPprobdataGetDatapoints(probdata, model);
   heurdata->nclusters = SCIPprobdataGetNClusters(probdata, model);
   heurdata->ndatapoints = SCIPprobdataGetNDatapoints(probdata, model);
   heurdata->dimension = SCIPprobdataGetDimension(probdata, model);
   heurdata->points = heurdata->datapoints->points;
   heurdata->objvar = SCIPprobdataGetObjvar(probdata, model);
   heurdata->evars = SCIPprobdataGetEvars(probdata, model);
   heurdata->xvars = SCIPprobdataGetXvars(probdata, model);
   heurdata->cvars = SCIPprobdataGetCvars(probdata, model);
   heurdata->quadcvars = SCIPprobdataGetQuadCvars(probdata, model);
   heurdata->kappavars = SCIPprobdataGetKappaVars(probdata, model);

   return SCIP_OKAY;
}


/** try to construct better solution */
static
SCIP_RETCODE constructSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic */
   SCIP_HEURDATA*        heurdata,           /**< data of heuristic */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   Datapoints* datapoints;
   SCIP_Real** points;
   SCIP_VAR* objvar;
   SCIP_VAR*** evars;
   SCIP_VAR*** xvars;
   SCIP_VAR*** cvars;
   SCIP_VAR*** quadcvars;
   SCIP_VAR** kappavars;
   int nclusters;
   int ndatapoints;
   int dimension;
   int k;
   int j;
   int i;
   int l;
   int model;
   int* npointsincluster;
   int* pointincluster;
   SCIP_SOL* bestsol;
   SCIP_SOL* newsol;
   SCIP_Real** barycenter;
   SCIP_Real* weightedvariances;
   SCIP_Real* disttocluster;
   SCIP_Real bestub;
   SCIP_Real newub = 0;
   SCIP_Bool foundcandidatesol;
   SCIP_Bool uselocalizedcut;
   SCIP_Bool success;

   assert( scip != NULL );
   assert( heur != NULL );
   assert( heurdata != NULL );

   datapoints = heurdata->datapoints;
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

   SCIP_CALL( SCIPgetIntParam(scip, "clustering/model", &model) );

   assert( datapoints != NULL );
   assert( points != NULL );
   assert( objvar != NULL || evars != NULL );
   assert( xvars != NULL );
   assert( cvars != NULL );
   assert( kappavars != NULL );
   assert( nclusters > 0 );
   assert( ndatapoints > 0 );
   assert( dimension > 0 );

   if ( model == MODEL_QUADRATIC )
      assert( quadcvars != NULL );

   /* get best feasible primal solution found so far */
   bestsol = SCIPgetBestSol(scip);

   /* only call heuristic, if a solution is at hand and this solution is new */
   if ( bestsol == NULL || SCIPsolGetIndex(bestsol) == heurdata->bestsolidx )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &npointsincluster, nclusters) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pointincluster, ndatapoints) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &weightedvariances, nclusters) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &disttocluster, ndatapoints) );

   SCIP_CALL( SCIPallocBufferArray(scip, &barycenter, nclusters) );
   for (j = 0; j < nclusters; ++j)
      SCIP_CALL( SCIPallocClearBufferArray(scip, &barycenter[j], dimension) );

   /* compute barycenter for each cluster and coordinate */
   for (j = 0; j < nclusters; ++j)
   {
      npointsincluster[j] = 0;

      for (i = 0; i < ndatapoints; ++i)
      {
         /* check if point is in cluster j */
         if ( SCIPgetSolVal(scip, bestsol, xvars[i][j]) > 0.5 )
         {
            pointincluster[i] = j;

            /* update barycenter of cluster */
            for (l = 0; l < dimension; ++l)
               barycenter[j][l] = (points[i][l] + npointsincluster[j] * barycenter[j][l]) / (npointsincluster[j] + 1);
            ++npointsincluster[j];
         }
      }

      /* there cannot be empty clusters */
      if ( SCIPisLE(scip, npointsincluster[j], 0.0) )
      {
         *result = SCIP_DIDNOTFIND;
         goto freemememory;
      }
   }

   /* compute weighted variance within each cluster */
   for (j = 0; j < nclusters; ++j)
   {
      for (i = 0; i < ndatapoints; ++i)
      {
         /* check if point is in cluster j */
         if ( SCIPgetSolVal(scip, bestsol, xvars[i][j]) > 0.5 )
         {
            /* update weighted variance within cluster j */
            for (l = 0; l < dimension; ++l)
            {
               weightedvariances[j] += points[i][l] * points[i][l] - 2.0 * points[i][l] * barycenter[j][l] + barycenter[j][l] * barycenter[j][l];
            }
         }
      }
      weightedvariances[j] = weightedvariances[j] / npointsincluster[j];

      /* if variance is zero the heuristic is not applicable */
      if ( SCIPisLE(scip, weightedvariances[j], 0.0) )
      {
         *result = SCIP_DIDNOTFIND;
         goto freemememory;
      }
   }

   foundcandidatesol = computeCandidates(scip, datapoints, nclusters, pointincluster, npointsincluster, weightedvariances, barycenter);

   /* check if there is a candidate solution to be constructed */
   if ( !foundcandidatesol )
   {
      *result = SCIP_DIDNOTFIND;
      goto freemememory;
   }

   /* try to create new feasible solution by improving best sol */
   SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );

   /* compute new objective upper bound */
   for (i = 0; i < ndatapoints; ++i)
   {
      j = pointincluster[i];
      for (l = 0; l < dimension; ++l)
         disttocluster[i] += points[i][l] * points[i][l] - 2 * points[i][l] * barycenter[j][l] + barycenter[j][l] * barycenter[j][l];
      newub += disttocluster[i];
   }

   /* check if new upper bound is smaller than best known in original space */
   bestub = SCIPgetSolOrigObj(scip, bestsol);

   if ( SCIPisLT(scip, newub, bestub) )
   {
      if ( model == MODEL_QUADRATIC )
      {
         /* set objective function value */
         SCIP_CALL( SCIPsetSolVal(scip, newsol, objvar, newub) );
      }
   }
   else
   {
      *result = SCIP_DIDNOTFIND;
      SCIP_CALL( SCIPfreeSol(scip, &newsol) );
      goto freemememory;
   }

   /* set xvars */
   for (i = 0; i < ndatapoints; ++i)
   {
      j = pointincluster[i];

      /* if point is already fixed, we do not continue */
      if ( SCIPisEQ(scip, SCIPvarGetUbLocal(xvars[i][j]), 0.0) && SCIPisEQ(scip, SCIPvarGetLbLocal(xvars[i][j]), 0.0) )
      {
         *result = SCIP_DIDNOTFIND;
          SCIP_CALL( SCIPfreeSol(scip, &newsol) );
          goto freemememory;
      }
      SCIP_CALL( SCIPsetSolVal(scip, newsol, xvars[i][j], 1) );
      for (k = 0; k < nclusters; ++k)
      {
         if ( k == j )
            continue;
         SCIP_CALL( SCIPsetSolVal(scip, newsol, xvars[i][k], 0) );
      }
   }

   /* set cvars and quadcvars */
   for (k = 0; k < nclusters; ++k)
   {
      for (l = 0; l < dimension; ++l)
      {
         SCIP_CALL( SCIPsetSolVal(scip, newsol, cvars[k][l], barycenter[k][l]) );

         if ( model == MODEL_QUADRATIC )
         {
            SCIP_CALL( SCIPsetSolVal(scip, newsol, quadcvars[k][l], barycenter[k][l] * barycenter[k][l]) );
         }
      }
   }

   /* set evars in epigraph model */
   if ( model == MODEL_QUADSOC )
   {
      for (i = 0; i < ndatapoints; ++i)
      {
         j = pointincluster[i];
         SCIP_CALL( SCIPsetSolVal(scip, newsol, evars[i][j], disttocluster[i]) );
         for (k = 0; k < nclusters; ++k)
         {
            if ( k == j )
               continue;
            SCIP_CALL( SCIPsetSolVal(scip, newsol, evars[i][k], 0.0) );
         }
      }
   }

   /* set kappa vars */
   SCIP_CALL( SCIPgetBoolParam(scip, "clustering/uselocalizedcardcut", &uselocalizedcut) );

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
               if ( SCIPgetSolVal(scip, newsol, xvars[i][c]) == 1 )
                  ++ncomplement;
            }
         }
         SCIP_CALL( SCIPsetSolVal(scip, newsol, kappavars[j], ncomplement) );
      }
   }
   else
   {
      for (j = 0; j < nclusters; ++j)
      {
         SCIP_CALL( SCIPsetSolVal(scip, newsol, kappavars[j], 0.0) );
      }
   }

   /* try and set solution */
   SCIP_CALL( SCIPtrySolFree(scip, &newsol,
      TRUE, /* printreason */
      TRUE, /* completely */
      TRUE, /* checkbounds */
      TRUE, /* checkintegrality */
      TRUE, /* checklprows */
      &success) );

   if ( success )
   {
      *result = SCIP_FOUNDSOL;
      goto freemememory;
   }
   else
   {
      *result = SCIP_DIDNOTFIND;
      goto freemememory;
   }

freemememory:
   for (j = nclusters - 1; j >= 0; --j)
   {
      SCIPfreeBufferArray(scip, &barycenter[j]);
   }
   SCIPfreeBufferArray(scip, &barycenter);
   SCIPfreeBufferArray(scip, &disttocluster);
   SCIPfreeBufferArray(scip, &weightedvariances);
   SCIPfreeBufferArray(scip, &pointincluster);
   SCIPfreeBufferArray(scip, &npointsincluster);

   heurdata->bestsolidx = SCIPsolGetIndex(bestsol);

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecImprove)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( scip != NULL );
   assert( heur != NULL );

   *result = SCIP_DIDNOTRUN;

   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( constructSolution(scip, heur, heurdata, result) );

   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the improvement primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurImprove(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create heuristic data */
   heurdata = NULL;
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   heurdata->datapoints = NULL;
   heurdata->nclusters = 0;
   heurdata->ndatapoints = 0;
   heurdata->dimension = 0;
   heurdata->objvar = NULL;
   heurdata->evars = NULL;
   heurdata->xvars = NULL;
   heurdata->cvars = NULL;
   heurdata->quadcvars = NULL;
   heurdata->kappavars = NULL;
   heurdata->bestsolidx = -1;

   heur = NULL;

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecImprove, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeImprove) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitImprove) );

   return SCIP_OKAY;
}
