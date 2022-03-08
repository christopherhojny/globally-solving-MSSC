/**@file   warmStart.cpp
 * @brief  compute a warm start solution for the k-means problem (using the usual k-means method)
 * @author Carina Costa
 */

#include <iostream>
#include "warmStart.h"
#include "datapoints.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include "getProbdata.h"
#include "typedefs.h"


/** initializes centroids */
static
std::vector<std::vector<SCIP_Real>> guessInitCentroids(
   SCIP*                 scip,               //!< SCIP data structure
   Datapoints*           datapoints,         //!< pointer to data point structure
   int                   nclusters           //!< number of clusters
   )
{
   assert( scip != NULL );
   assert( datapoints != NULL );
   assert( nclusters > 0 );

   int ndatapoints;
   int dimension;
   int i;
   int l;
   int c;
   int indexMaxDist;
   int nCentroidsAlreadyFound = 0;
   SCIP_Real** points;

   ndatapoints = datapoints->ndatapoints;
   dimension = datapoints->dimension;
   points = datapoints->points;

   assert( ndatapoints > 0 );
   assert( dimension > 0 );
   assert( points != NULL );

   std::vector<std::vector<SCIP_Real>> centroids(nclusters, std::vector<SCIP_Real>(dimension));
   std::vector<SCIP_Real> distPointToClosestCentroid(ndatapoints, 0.0);

   // start with an initial guess for the centroids, take the first point as first centroid
   for (l = 0; l < dimension; ++l)
      centroids[0][l] = points[0][l];
   ++nCentroidsAlreadyFound;

   // compute distance between points and first centroid
   distPointToClosestCentroid[0] = 0.0; // first point is exactly the first centroid

   for (i = 1; i < ndatapoints; ++i)
   {
      for (l = 0; l < dimension; ++l)
      {
         distPointToClosestCentroid[i] += points[i][l] * points[i][l] - 2.0 * points[i][l] * centroids[0][l] + centroids[0][l] * centroids[0][l];
      }
   }

   // second centroid is the furthest point from the first centroid
   auto findMaxIndex = std::max_element(distPointToClosestCentroid.begin(), distPointToClosestCentroid.end());
   indexMaxDist = std::distance(distPointToClosestCentroid.begin(), findMaxIndex);

   for (l = 0; l < dimension; ++l)
      centroids[1][l] = points[indexMaxDist][l];
   ++nCentroidsAlreadyFound;

   // if k > 2, the next centroid is the furthest point from its nearest centroid
   while ( nCentroidsAlreadyFound < nclusters )
   {
      // update distance between point and closest centroid
      c = nCentroidsAlreadyFound - 1;
      for (i = 1; i < ndatapoints; ++i)
      {
         SCIP_Real tmpDist = 0.0;
         for (l = 0; l < dimension; ++l)
         {
            tmpDist += points[i][l] * points[i][l] - 2.0 * points[i][l] * centroids[c][l] + centroids[c][l] * centroids[c][l];
         }

         // check if point is closest to the last found cluster than to the other clusters
         if ( SCIPisLT(scip, tmpDist, distPointToClosestCentroid[i]) )
            distPointToClosestCentroid[i] = tmpDist;
      }

      // obtain new centroid
      auto itMaxIndex = std::max_element(distPointToClosestCentroid.begin(), distPointToClosestCentroid.end());
      indexMaxDist = std::distance(distPointToClosestCentroid.begin(), itMaxIndex);
      for (l = 0; l < dimension; ++l)
         centroids[c + 1][l] = points[indexMaxDist][l];
      ++nCentroidsAlreadyFound;
   }

   return centroids;
}


/** create and add warm start solution for clustering problem */
SCIP_RETCODE addClusteringWarmStart(
   SCIP*                 scip,               /**< SCIP data structure */
   Datapoints*           datapoints,         /**< pointer to data point structure */
   int                   nclusters           /**< number of clusters */
   )
{
   SCIP_PROBDATA* probdata;

   assert( scip != NULL );
   assert( datapoints != NULL );
   assert( nclusters > 0 );

   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   SCIP_VAR*** xvars;
   SCIP_VAR*** evars;
   SCIP_VAR* objvar;
   SCIP_VAR** kappavars;
   SCIP_VAR*** cvars;
   SCIP_VAR*** quadcvars;
   SCIP_SOL* warmstartsol;
   int model;
   int ndatapoints;
   int dimension;
   int i;
   int l;
   int j;
   int c;
   int indexClosestCentroid = 0;
   int npointsincluster;
   int iterationscounter;
   int ncomplement;
   SCIP_Real** points;
   SCIP_Real minDist;
   SCIP_Real toleranceTermination = 1e-4;
   int maxNrIterations = 100;
   SCIP_Real infinityNormVal = SCIPinfinity(scip);
   SCIP_Real eval;
   SCIP_Real objval = 0.0;
   SCIP_Bool stored;
   SCIP_Bool alternativesolstored;
   SCIP_Bool uselocalizedcut;

   SCIP_CALL( SCIPgetIntParam(scip, "clustering/model", &model) );

   minDist = SCIPinfinity(scip);

   ndatapoints = datapoints->ndatapoints;
   dimension = datapoints->dimension;
   points = datapoints->points;
   objvar = SCIPprobdataGetObjvar(probdata, model);
   evars = SCIPprobdataGetEvars(probdata, model);
   xvars = SCIPprobdataGetXvars(probdata, model);
   kappavars = SCIPprobdataGetKappaVars(probdata, model);
   cvars = SCIPprobdataGetCvars(probdata, model);
   quadcvars = SCIPprobdataGetQuadCvars(probdata, model);

   assert( points != NULL );
   assert( objvar != NULL || evars != NULL );
   assert( xvars != NULL );
   assert( cvars != NULL );
   assert( kappavars != NULL );
   assert( ndatapoints > 0 );
   assert( dimension > 0 );

   std::vector<std::vector<SCIP_Real>> centroids(nclusters, std::vector<SCIP_Real>(dimension));
   std::vector<std::vector<int>> assignments(ndatapoints, std::vector<int>(nclusters, 0));

   std::vector<std::vector<SCIP_Real>> centroidsNew(nclusters, std::vector<SCIP_Real>(dimension, 0));
   std::vector<std::vector<int>> assignmentsNew(ndatapoints, std::vector<int>(nclusters, 0));
   std::vector<SCIP_Real> normValues(nclusters * dimension);

   // create warm start solution
   SCIP_CALL( SCIPcreateOrigSol(scip, &warmstartsol, NULL) );

   // guess initial centroids
   centroids = guessInitCentroids(scip, datapoints, nclusters);

   // compute assignments
   assignments[0][0] = 1; // first point was set to be the first cluster

   for (i = 1; i < ndatapoints; ++i)
   {
      for (j = 0; j < nclusters; ++j)
      {
         SCIP_Real currentDist = 0.0;

         for (l = 0; l < dimension; ++l)
            currentDist += points[i][l] * points[i][l] - 2.0 * points[i][l] * centroids[j][l] + centroids[j][l] * centroids[j][l];

         if ( j == 0 )
         {
            minDist = currentDist;
            indexClosestCentroid = j;
         }

         if ( SCIPisLT(scip, currentDist, minDist) )
         {
            minDist = currentDist;
            indexClosestCentroid = j;
         }
      }

      assert( 0 <= minDist && indexClosestCentroid < nclusters && 0 <= indexClosestCentroid );

      assignments[i][indexClosestCentroid] = 1;
   }

   iterationscounter = 1;

   // update centroids and assignments alternatingly
   while ( SCIPisGT(scip, infinityNormVal, toleranceTermination) && SCIPisLT(scip, iterationscounter, maxNrIterations) )
   {
      // update centroids
      for (j = 0; j < nclusters; ++j)
      {
         npointsincluster = 0;

         for (i = 0; i < ndatapoints; ++i)
         {
            if ( assignments[i][j] != 1 )
               continue;

            for (l = 0; l < dimension; ++l)
               centroidsNew[j][l] = (points[i][l] + npointsincluster * centroidsNew[j][l]) / (npointsincluster + 1);
            ++npointsincluster;
         }
      }

      // compute new assignment
      assignmentsNew[0][0] = 1;

      for (i = 1; i < ndatapoints; ++i)
      {
         for (j = 0; j < nclusters; ++j)
         {
            SCIP_Real currentDist = 0.0;

            for (l = 0; l < dimension; ++l)
            {
               currentDist += points[i][l] * points[i][l] - 2.0 * points[i][l] * centroidsNew[j][l] + centroidsNew[j][l] * centroidsNew[j][l];
            }

            if ( j == 0 )
            {
               minDist = currentDist;
               indexClosestCentroid = j;
            }

            if ( SCIPisLT(scip, currentDist, minDist) )
            {
               minDist = currentDist;
               indexClosestCentroid = j;
            }
         }

         assert( 0 <= minDist && indexClosestCentroid < nclusters && 0 <= indexClosestCentroid );

         assignmentsNew[i][indexClosestCentroid] = 1;
      }

      // check if centroids changed in the last two iterations
      for (j = 0; j < nclusters; ++j)
      {
         for (l = 0; l < dimension; ++l)
         {
            normValues[l + (dimension * j)] = std::abs(centroidsNew[j][l] - centroids[j][l]);
         }
      }

      infinityNormVal = *std::max_element(normValues.begin(), normValues.end());

      // store new centroids and assignments
      centroids = centroidsNew;
      assignments = assignmentsNew;
      ++iterationscounter;

      // reset
      for (i = 0; i < ndatapoints; ++i)
      {
         for (j = 0; j < nclusters; ++j)
         {
             assignmentsNew[i][j] = 0;
         }
      }
   }

   // build warm start solution
   if ( model == MODEL_QUADRATIC )
   {
      assert( quadcvars != NULL );
      assert( objvar != NULL );
   }
   else
      assert( evars != NULL );

   // set cvals
   for (j = 0; j < nclusters; ++j)
   {
      for (l = 0; l < dimension; ++l)
      {
         SCIP_CALL( SCIPsetSolVal(scip, warmstartsol, cvars[j][l], centroids[j][l]) );

         if ( model == MODEL_QUADRATIC )
         {
            SCIP_CALL( SCIPsetSolVal(scip, warmstartsol, quadcvars[j][l], centroids[j][l] * centroids[j][l]) );
         }
      }
   }

   // set xvars
   for (i = 0; i < ndatapoints; ++i)
   {
      for (j = 0; j < nclusters; ++j)
      {
         SCIP_CALL( SCIPsetSolVal(scip, warmstartsol, xvars[i][j], assignments[i][j]) );
      }
   }

   SCIP_CALL( SCIPgetBoolParam(scip, "clustering/uselocalizedcardcut", &uselocalizedcut) );

   // set kappavars
   if ( uselocalizedcut )
   {
      for (j = 0; j < nclusters; ++j)
      {
         ncomplement = 0;

         for (c = 0; c < nclusters; ++c)
         {
            if ( c == j )
               continue;

            for (i = 0; i < ndatapoints; ++i)
            {
               if ( assignments[i][c] == 1 )
                  ++ncomplement;
            }
         }
         SCIP_CALL( SCIPsetSolVal(scip, warmstartsol, kappavars[j], ncomplement) );
      }
   }
   else
   {
      for (j = 0; j < nclusters; ++j)
      {
         SCIP_CALL( SCIPsetSolVal(scip, warmstartsol, kappavars[j], 0.0) );
      }
   }

   // compute obj val or evals
   for (j = 0; j < nclusters; ++j)
   {
      for (i = 0; i < ndatapoints; ++i)
      {
         if ( SCIPgetSolVal(scip, warmstartsol, xvars[i][j]) == 0 )
            continue;

         eval = 0.0;

         for (l = 0; l < dimension; ++l)
            eval += points[i][l] * points[i][l] - 2.0 * points[i][l] * centroids[j][l] + centroids[j][l] * centroids[j][l];

         if ( model == MODEL_QUADRATIC )
            objval += eval;
         else
         {
            SCIP_CALL( SCIPsetSolVal(scip, warmstartsol, evars[i][j], eval) );

            for (c = 0; c < nclusters; ++c)
            {
               if ( c == j )
                  continue;
               else
               {
                  SCIP_CALL( SCIPsetSolVal(scip, warmstartsol, evars[i][c], 0.0) );
               }
            }
         }
      }
   }

   if ( model == MODEL_QUADRATIC )
   {
      SCIP_CALL( SCIPsetSolVal(scip, warmstartsol, objvar, objval) );
   }

   SCIP_CALL( SCIPaddSol(scip, warmstartsol, &stored) );

   if ( stored )
      std::cout << "Warmstart solution was added." << std::endl;

   // create alternative solution to avoid warmstartsol being discarded due to numerical tolerance issues
   SCIP_SOL* alternativewarmstartsol;
   SCIP_CALL( SCIPcreateSolCopyOrig(scip, &alternativewarmstartsol, warmstartsol) );

   if ( model == MODEL_QUADRATIC )
   {
      SCIP_CALL( SCIPsetSolVal(scip, alternativewarmstartsol, objvar, objval + 130) );
   }
   else
   {
      for (j = 0; j < nclusters; ++j)
      {
         for (i = 0; i < ndatapoints; ++i)
         {
            SCIP_Real currenteval = SCIPgetSolVal(scip, warmstartsol, evars[i][j]);
            SCIP_CALL( SCIPsetSolVal(scip, alternativewarmstartsol, evars[i][j], currenteval + 1) );
         }
      }
   }

   SCIP_CALL( SCIPaddSol(scip, alternativewarmstartsol, &alternativesolstored) );
   SCIP_CALL( SCIPfreeSol(scip, &alternativewarmstartsol) );

   SCIP_CALL( SCIPfreeSol(scip, &warmstartsol) );

   return SCIP_OKAY;
}
