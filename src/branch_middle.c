/**@file   branch_middle.c
 * @brief  branching rule based on centrality of data point
 * @author
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "branch_middle.h"
#include "datapoints.h"
#include "getProbdata.h"

#define BRANCHRULE_NAME            "middle"
#define BRANCHRULE_DESC            "branching rule based on centrality of data point"
#define BRANCHRULE_PRIORITY        -50000
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpMiddle)
{
   SCIP_PROBDATA* probdata;
   Datapoints* datapoints;
   int dimension;
   int nclusters;
   int model;

   SCIP_VAR*** cvars;

   SCIP_VARDATA* vardata;
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandsfrac;
   int nlpcands;

   SCIP_VAR* branchvar;
   SCIP_Real branchval;
   SCIP_NODE* upchild;
   SCIP_NODE* downchild;

   SCIP_Real minvalue = SCIPinfinity(scip);
   SCIP_Real centralityval;
   int v;
   int k;
   int i;
   int branchpoint;
   SCIP_Real partdist;

   assert( scip != NULL );
   assert( branchrule != NULL );
   assert( strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   SCIP_CALL( SCIPgetIntParam(scip, "clustering/model", &model) );

   cvars = SCIPprobdataGetCvars(probdata, model);
   dimension = SCIPprobdataGetDimension(probdata, model);
   datapoints = SCIPprobdataGetDatapoints(probdata, model);
   nclusters = SCIPprobdataGetNClusters(probdata, model);

   assert( cvars != NULL );
   assert( dimension > 0 );
   assert( datapoints != NULL );
   assert( nclusters > 0 );

   /* get fractional LP candidates */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, NULL, &lpcandsfrac, NULL, &nlpcands, NULL) );
   assert( nlpcands > 0 );

   for (v = 0; v < nlpcands; ++v)
   {
      assert( lpcands[v] != NULL );
      vardata = SCIPvarGetData(lpcands[v]);

      /* skip variables which are not of x-type */
      if ( vardata == NULL )
         continue;

      /* get point index of branch candidate */
      branchpoint = SCIPvardataGetPoint(vardata, model);

      /* iterate over all clusters to compute centrality of the branchpoint candidate */
      centralityval = 0.0;
      for (k = 0; k < nclusters; ++k)
      {
         /* compute distance to cluster */
         for (i = 0; i < dimension; ++i)
         {
            partdist = datapoints->points[branchpoint][i] - SCIPgetSolVal(scip, NULL, cvars[k][i]);
            centralityval += partdist * partdist;
         }
      }
      centralityval = sqrt(centralityval);

      if ( SCIPisLT(scip, centralityval, minvalue) )
      {
         minvalue = centralityval;
         branchvar = lpcands[v]; /* the variable to branch on */
         branchval = lpcandsfrac[v]; /* the fractional value of the variable */
      }
   }

   assert( SCIPisLT(scip, minvalue, SCIPinfinity(scip)) );

   /* create the branch-and-bound tree child nodes of the current node */
   SCIP_CALL( SCIPbranchVarVal(scip, branchvar, branchval, &downchild, NULL, &upchild) );
   assert(downchild != NULL);
   assert(upchild != NULL);

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}

/** creates the middle branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleMiddle(
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

   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpMiddle) );

   return SCIP_OKAY;
}
