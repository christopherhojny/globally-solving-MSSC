/**@file   branch_maxdist.c
 * @brief  branching rule based on maximum distance points
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <math.h>

#include "branch_maxdist.h"
#include "datapoints.h"
#include "getProbdata.h"

#define BRANCHRULE_NAME            "maxdist"
#define BRANCHRULE_DESC            "branching rule based on maximum distance points"
#define BRANCHRULE_PRIORITY        -450000
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpMaxdist)
{
   SCIP_PROBDATA* probdata;
   Datapoints* datapoints;
   int dimension;
   int model;

   SCIP_VAR*** xvars;
   SCIP_VAR*** cvars;

   SCIP_VARDATA* vardata;
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandsfrac;
   int nlpcands;

   SCIP_Real curdist;
   SCIP_Real maxdist = 0.0;
   int clusteridx = -1;
   int pointidx = -1;

   SCIP_NODE* upchild;
   SCIP_NODE* downchild;
   SCIP_VAR* branchvar;
   SCIP_Real branchval;

   SCIP_Real partdist;
   int v;
   int i;
   int branchcand;
   int branchpoint;
   int branchcluster;

   assert( scip != NULL );
   assert( branchrule != NULL );
   assert( strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   SCIP_CALL( SCIPgetIntParam(scip, "clustering/model", &model) );

   xvars = SCIPprobdataGetXvars(probdata, model);
   cvars = SCIPprobdataGetCvars(probdata, model);
   dimension = SCIPprobdataGetDimension(probdata, model);
   datapoints = SCIPprobdataGetDatapoints(probdata, model);

   assert( xvars != NULL );
   assert( cvars != NULL );
   assert( dimension > 0 );
   assert( datapoints != NULL );

   /* get fractional LP candidates */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, NULL, &lpcandsfrac, NULL, &nlpcands, NULL) );
   assert(nlpcands > 0);

   for (v = 0; v < nlpcands; ++v)
   {
      assert(lpcands[v] != NULL);
      vardata = SCIPvarGetData(lpcands[v]);

      /* skip variables which are not of x-type */
      if ( vardata == NULL )
         continue;

      /* get point index of branch candidate */
      branchpoint = SCIPvardataGetPoint(vardata, model);
      branchcluster = SCIPvardataGetCluster(vardata, model);

      /* compute distance to cluster */
      curdist = 0.0;
      for (i = 0; i < dimension; ++i)
      {
         partdist = datapoints->points[branchpoint][i] - SCIPgetSolVal(scip, NULL, cvars[branchcluster][i]);
         curdist += partdist * partdist;
      }

      if ( SCIPisGT(scip, curdist, maxdist) )
      {
         maxdist = curdist;
         clusteridx = branchcluster;
         pointidx = branchpoint;
         branchcand = v;
      }
   }

   if ( clusteridx == -1 )
      *result = SCIP_DIDNOTFIND;
   else
   {
      branchvar = xvars[pointidx][clusteridx];
      branchval = lpcandsfrac[branchcand];

      SCIP_CALL( SCIPbranchVarVal(scip, branchvar, branchval, &downchild, NULL, &upchild) );
      assert(downchild != NULL);
      assert(upchild != NULL);

      *result = SCIP_BRANCHED;
   }

   return SCIP_OKAY;
}

/** creates the entropy branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleMaxdist(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   /* create branching rule data */
   branchruledata = NULL;
   branchrule = NULL;

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata) );
   assert( branchrule != NULL );

   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpMaxdist) );

   return SCIP_OKAY;
}
