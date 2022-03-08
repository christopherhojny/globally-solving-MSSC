/**@file   branch_pairs_midmost.c
 * @brief  branching rule based on pairs of points (midmost) and clusters
 * @author Carina Costa
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <math.h>

#include "scip/cons_linear.h"

#include "branch_pairs_midmost.h"
#include "datapoints.h"
#include "getProbdata.h"

#define BRANCHRULE_NAME            "pairs_midmost"
#define BRANCHRULE_DESC            "branching rule based on pairs of points (midmost) and clusters"
#define BRANCHRULE_PRIORITY        -500000
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0


/** computes a penalty distance between a point and two clusters */
static
SCIP_Real computePenaltyDistance(
   SCIP*                 scip,               /**< SCIP data structure */
   Datapoints*           datapoints,         /**< array of data points */
   SCIP_VAR***           cvars,              /**< (nclusters x dimension)-array of centroid variables */
   int                   curpointidx,        /**< index of data point */
   int                   clusteridx1,        /**< index of first cluster */
   int                   clusteridx2         /**< index of second cluster */
   )
{
   SCIP_Real penaltydist;
   SCIP_Real squareddist1 = 0.0;
   SCIP_Real squareddist2 = 0.0;
   SCIP_Real partdist;
   int i;

   assert( scip != NULL );
   assert( datapoints != NULL );
   assert( cvars != NULL );

   for (i = 0; i < datapoints->dimension; ++i)
   {
      partdist = datapoints->points[curpointidx][i] - SCIPgetSolVal(scip, NULL, cvars[clusteridx1][i]);
      squareddist1 += partdist * partdist;
      partdist = datapoints->points[curpointidx][i] - SCIPgetSolVal(scip, NULL, cvars[clusteridx2][i]);
      squareddist2 += partdist * partdist;
   }
   penaltydist = squareddist1 + squareddist2;

   return penaltydist;
}


/** branching execution method for fractional LP solutions */
static SCIP_DECL_BRANCHEXECLP(branchExeclpPairsMidmost)
{
   SCIP_PROBDATA* probdata;
   Datapoints* datapoints;
   SCIP_VAR*** cvars;
   SCIP_VAR*** xvars;
   int nclusters;
   int ndatapoints;
   int dimension;
   int model;
   int nmidpoints;

   SCIP_VARDATA* vardata;
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandsfrac;
   int nlpcands;

   SCIP_Real penaltydist;
   SCIP_Real minval = SCIPinfinity(scip);
   SCIP_Real sum1;
   SCIP_Real sum2;
   int clusteridx1 = -1;
   int clusteridx2 = -1;
   int clusteridx = -1;
   int pointidx1 = -1;
   int pointidx2 = -1;
   int curclusteridx;
   int curpointidx;
   SCIP_Real downbound;

   int** candsbycluster;
   int* ncandsincluster;

   int* pointsnearmiddle;
   SCIP_Bool* pointisinset;

   int v;
   int k;
   int c;
   int i;
   int p;

   assert( scip != NULL );
   assert( branchrule != NULL );
   assert( strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   SCIP_CALL( SCIPgetIntParam(scip, "clustering/model", &model) );

   datapoints = SCIPprobdataGetDatapoints(probdata, model);
   cvars = SCIPprobdataGetCvars(probdata, model);
   xvars = SCIPprobdataGetXvars(probdata, model);
   nclusters = SCIPprobdataGetNClusters(probdata, model);
   dimension = SCIPprobdataGetDimension(probdata, model);
   ndatapoints = SCIPprobdataGetNDatapoints(probdata, model);

   assert( datapoints != NULL );
   assert( cvars != NULL );
   assert( xvars != NULL );
   assert( nclusters > 0 );
   assert( ndatapoints > 0 );
   assert( dimension > 0 );

   /* get fractional LP candidates */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, NULL, &lpcandsfrac, NULL, &nlpcands, NULL) );
   assert( nlpcands > 0 );

   /* store branching candidates by cluster */
   SCIP_CALL( SCIPallocBufferArray(scip, &candsbycluster, nclusters) );
   for (v = 0; v < nclusters; ++v)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &candsbycluster[v], ndatapoints) );
   }
   SCIP_CALL( SCIPallocClearBufferArray(scip, &ncandsincluster, nclusters) );

   for (v = 0; v < nlpcands; ++v)
   {
      assert( lpcands[v] != NULL );
      vardata = SCIPvarGetData(lpcands[v]);

      /* skip variables which are not of x-type */
      if ( vardata == NULL )
         continue;

      curclusteridx = SCIPvardataGetCluster(vardata, model);
      curpointidx = SCIPvardataGetPoint(vardata, model);
      candsbycluster[curclusteridx][ncandsincluster[curclusteridx]] = curpointidx;
      ncandsincluster[curclusteridx] += 1;
   }

   /* find the pair of clusters and the candidate point p such that p is the midmost point between the two clusters */
   for (k = 0; k < nclusters - 1; ++k)
   {
      for (c = k + 1; c < nclusters; ++c)
      {
         for (v = 0; v < ncandsincluster[k]; ++v)
         {
            curpointidx = candsbycluster[k][v];
            penaltydist = computePenaltyDistance(scip, datapoints, cvars, curpointidx, k, c);

            if ( penaltydist < minval )
            {
               minval = penaltydist;
               pointidx1 = curpointidx;
               clusteridx1 = k;
               clusteridx2 = c;
            }
         }

         for (v = 0; v < ncandsincluster[c]; ++v)
         {
            curpointidx = candsbycluster[c][v];
            penaltydist = computePenaltyDistance(scip, datapoints, cvars, curpointidx, k, c);

            if ( penaltydist < minval )
            {
               minval = penaltydist;
               pointidx1 = curpointidx;
               clusteridx1 = k;
               clusteridx2 = c;
            }
         }
      }
   }

   if ( pointidx1 == -1 || SCIPisEQ(scip, minval, SCIPinfinity(scip)) )
   {
      *result = SCIP_DIDNOTFIND;
       goto FREEBRANCHINGPOINTS;
   }

   /* for the chosen pair of clusters, find the other L midmost points between the two clusters */
   SCIP_CALL( SCIPgetIntParam(scip, "clustering/branchpairsnmidpoints", &nmidpoints) );
   nmidpoints = MIN(nmidpoints, ncandsincluster[clusteridx1] + ncandsincluster[clusteridx2] - 1);

   SCIP_CALL( SCIPallocBufferArray(scip, &pointsnearmiddle, nmidpoints + 1) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &pointisinset, ndatapoints) );
   pointsnearmiddle[0] = pointidx1;
   pointisinset[pointidx1] = TRUE;

   for (int m = 1; m < nmidpoints + 1; ++m)
   {
      minval = SCIPinfinity(scip);
      pointidx2 = -1;

      for (v = 0; v < ncandsincluster[clusteridx1]; ++v)
      {
         curpointidx = candsbycluster[clusteridx1][v];

         if ( pointisinset[curpointidx] )
            continue;
         penaltydist = computePenaltyDistance(scip, datapoints, cvars, curpointidx, clusteridx1, clusteridx2);

         if ( penaltydist < minval )
         {
            minval = penaltydist;
            pointidx2 = curpointidx;
         }
      }

      for (v = 0; v < ncandsincluster[clusteridx2]; ++v)
      {
         curpointidx = candsbycluster[clusteridx2][v];

         if ( pointisinset[curpointidx] )
            continue;
         penaltydist = computePenaltyDistance(scip, datapoints, cvars, curpointidx, clusteridx1, clusteridx2);

         if ( penaltydist < minval )
         {
            minval = penaltydist;
            pointidx2 = curpointidx;
         }
      }

      if ( pointidx2 == -1 || SCIPisEQ(scip, minval, SCIPinfinity(scip)) )
      {
         *result = SCIP_DIDNOTFIND;
         SCIPfreeBufferArray(scip, &pointisinset);
         SCIPfreeBufferArray(scip, &pointsnearmiddle);
         goto FREEBRANCHINGPOINTS;
      }

      pointsnearmiddle[m] = pointidx2;
      pointisinset[pointidx2] = TRUE;
   }

   /* find the two nearest points among the middle points */
   if ( nmidpoints > 1 )
   {
      minval = SCIPinfinity(scip);
      pointidx1 = -1;
      pointidx2 = -2;

      for (p = 0; p < nmidpoints; ++p)
      {
         for (v = p + 1; v < nmidpoints + 1; ++v)
         { 
            SCIP_Real dist = 0.0;
            sum1 = SCIPgetSolVal(scip, NULL, xvars[pointsnearmiddle[p]][clusteridx1]) + SCIPgetSolVal(scip, NULL, xvars[pointsnearmiddle[v]][clusteridx1]);
            sum2 = SCIPgetSolVal(scip, NULL, xvars[pointsnearmiddle[p]][clusteridx2]) + SCIPgetSolVal(scip, NULL, xvars[pointsnearmiddle[v]][clusteridx2]);
            if ( SCIPisEQ(scip, sum1, 0.0) || SCIPisEQ(scip, sum1, 1.0) || SCIPisEQ(scip, sum2, 0.0) || SCIPisEQ(scip, sum2, 1.0) )
               continue;

            for (i = 0; i < dimension; ++i)
            {
               SCIP_Real partdist = datapoints->points[pointsnearmiddle[p]][i] - datapoints->points[pointsnearmiddle[v]][i];
               dist += partdist * partdist;
            }
            if ( dist < minval )
            {
               minval = dist;
               pointidx1 = pointsnearmiddle[p];
               pointidx2 = pointsnearmiddle[v];
            }
         }
      }
   }
   SCIPfreeBufferArray(scip, &pointisinset);
   SCIPfreeBufferArray(scip, &pointsnearmiddle);

   if ( pointidx1 == -1 || pointidx2 == -1 || SCIPisEQ(scip, minval, SCIPinfinity(scip)) )
      *result = SCIP_DIDNOTFIND;
   else
   {
      SCIP_NODE* downchild;
      SCIP_NODE* upchild;

      SCIP_CONS* consdownchild;
      SCIP_CONS* consupchild;

      SCIP_VAR** vars;
      SCIP_Real* vals;
      char name[SCIP_MAXSTRLEN];

      SCIP_Real fracval1;
      SCIP_Real fracval2;
      SCIP_Bool branchmostfractional;

      /* choose cluster based on most fractional or least fractional part */
      sum1 = SCIPgetSolVal(scip, NULL, xvars[pointidx1][clusteridx1]) + SCIPgetSolVal(scip, NULL, xvars[pointidx2][clusteridx1]);
      sum2 = SCIPgetSolVal(scip, NULL, xvars[pointidx1][clusteridx2]) + SCIPgetSolVal(scip, NULL, xvars[pointidx2][clusteridx2]);
      fracval1 = SCIPisLT(scip, sum1, 1.0) ? MIN(sum1, 1.0 - sum1) : MIN(sum1 - 1.0, 2.0 - sum1);
      fracval2 = SCIPisLT(scip, sum2, 1.0) ? MIN(sum2, 1.0 - sum2) : MIN(sum2 - 1.0, 2.0 - sum2);

      if ( SCIPisEQ(scip, fracval1, 0.0) && SCIPisEQ(scip, fracval2, 0.0) )
      {
         *result = SCIP_DIDNOTFIND;
         goto FREEBRANCHINGPOINTS;
      }

      SCIP_CALL( SCIPgetBoolParam(scip, "clustering/branchpairsmostfractional", &branchmostfractional) );

      if ( branchmostfractional )
      {
         if ( SCIPisGE(scip, fracval1, fracval2) )
         {
            clusteridx = clusteridx1;
            downbound = SCIPisLT(scip, sum1, 1.0) ? 0.0 : 1.0;
         }
         else
         {
            clusteridx = clusteridx2;
            downbound = SCIPisLT(scip, sum2, 1.0) ? 0.0 : 1.0;
         }
      }
      else
      {
         if ( SCIPisLE(scip, fracval1, fracval2) )
         {
            clusteridx = clusteridx1;
            downbound = SCIPisLT(scip, sum1, 1.0) ? 0.0 : 1.0;
         }
         else
         {
            clusteridx = clusteridx2;
            downbound = SCIPisLT(scip, sum2, 1.0) ? 0.0 : 1.0;
         }
      }

      assert( clusteridx != -1 );

      /* create two children: x[idx1][k] + x[idx2][k] <= downbound or x[idx1][k] + x[idx2][k] >= downbound + 1 */
      SCIP_CALL( SCIPallocBufferArray(scip, &vars, 2) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vals, 2) );

      vars[0] = xvars[pointidx1][clusteridx];
      vars[1] = xvars[pointidx2][clusteridx];

      vals[0] = 1.0;
      vals[1] = 1.0;

      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "down_%d_%d_%d", pointidx1, pointidx2, (int)downbound);
      SCIP_CALL( SCIPcreateConsLinear(scip, &consdownchild, name, 2, vars, vals, 0.0, downbound,
            TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "up_%d_%d_%d", pointidx1, pointidx2, (int)downbound);
      SCIP_CALL( SCIPcreateConsLinear(scip, &consupchild, name, 2, vars, vals, downbound + 1.0, 2.0,
            TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

      /* add constraints to nodes */
      SCIP_CALL( SCIPcreateChild(scip, &downchild, 0.0, SCIPgetLocalTransEstimate(scip)) );
      SCIP_CALL( SCIPcreateChild(scip, &upchild, 0.0, SCIPgetLocalTransEstimate(scip)) );

      SCIP_CALL( SCIPaddConsNode(scip, downchild, consdownchild, NULL) );
      SCIP_CALL( SCIPaddConsNode(scip, upchild, consupchild, NULL) );

      /* release constraints */
      SCIP_CALL( SCIPreleaseCons(scip, &consdownchild) );
      SCIP_CALL( SCIPreleaseCons(scip, &consupchild) );

      *result = SCIP_BRANCHED;

      SCIPfreeBufferArray(scip, &vals);
      SCIPfreeBufferArray(scip, &vars);
   }


 FREEBRANCHINGPOINTS:
   SCIPfreeBufferArray(scip, &ncandsincluster);
   for (v = nclusters - 1; v >= 0; --v)
   {
      SCIPfreeBufferArray(scip, &candsbycluster[v]);
   }
   SCIPfreeBufferArray(scip, &candsbycluster);

   return SCIP_OKAY;
}

/** creates the pairsmidmost branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchrulePairsMidmost(
   SCIP *scip /**< SCIP data structure */
)
{
   SCIP_BRANCHRULEDATA *branchruledata;
   SCIP_BRANCHRULE *branchrule;

   /* create branching rule data */
   branchruledata = NULL;
   branchrule = NULL;

   /* include branching rule */
   SCIP_CALL(SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata));
   assert(branchrule != NULL);

   SCIP_CALL(SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpPairsMidmost));

   return SCIP_OKAY;
}
