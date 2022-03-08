/**@file   branch_entropy.c
 * @brief  branching rule based on entropy
 * @author
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <math.h>

#include "branch_entropy.h"
#include "datapoints.h"
#include "getProbdata.h"

#define BRANCHRULE_NAME            "entropy"
#define BRANCHRULE_DESC            "branching rule based on entropy"
#define BRANCHRULE_PRIORITY        -50000
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpEntropy)
{
   SCIP_PROBDATA* probdata;
   SCIP_VARDATA* vardata;
   SCIP_VAR*** xvars;
   SCIP_VAR** lpcands;
   SCIP_VAR* branchvar;
   SCIP_NODE* upchild;
   SCIP_NODE* downchild;
   SCIP_Real* lpcandsfrac;
   SCIP_Real bestvalue;
   SCIP_Real entropyvalue;
   SCIP_Real branchval;
   int nclusters;
   int model;
   int nlpcands;
   int v;
   int k;
   int branchpoint;
   SCIP_Bool useminentropy;

   assert( scip != NULL );
   assert( branchrule != NULL );
   assert( strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   SCIP_CALL( SCIPgetIntParam(scip, "clustering/model", &model) );

   xvars = SCIPprobdataGetXvars(probdata, model);
   nclusters = SCIPprobdataGetNClusters(probdata, model);

   assert( xvars != NULL );
   assert( nclusters > 0 );

   /* get fractional LP candidates */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, NULL, &lpcandsfrac, NULL, &nlpcands, NULL) );
   assert( nlpcands > 0 );

   SCIP_CALL( SCIPgetBoolParam(scip, "clustering/useminentropy", &useminentropy) );

   if ( useminentropy )
      bestvalue = SCIPinfinity(scip);
   else
      bestvalue = -1.0;

   for (v = 0; v < nlpcands; ++v)
   {
      entropyvalue = 0.0;

      assert( lpcands[v] != NULL );
      vardata = SCIPvarGetData(lpcands[v]);
      if ( vardata == NULL )
         continue;

      /* get point index of branch candidate */
      branchpoint = SCIPvardataGetPoint(vardata, model);

      /* iterate over all clusters to compute entropy of the branchpoint candidate */
      for (k = 0; k < nclusters; ++k)
      {
         SCIP_Real solval;
         solval = SCIPgetSolVal(scip, NULL, xvars[branchpoint][k]);

         if ( SCIPisEQ(scip, solval, 0.0) || SCIPisLT(scip, solval, 0.0) ) /* the log2 (0) is defined as 0 */
            continue;

         entropyvalue -= solval * log2(solval);
      }

      /* whether we select least uncertain point or most uncertain point based on entropy */
      if ( useminentropy )
      {
         if ( SCIPisLT(scip, entropyvalue, bestvalue) )
         {
            bestvalue = entropyvalue;
            branchvar = lpcands[v]; /* the variable to branch on */
            branchval = lpcandsfrac[v]; /* the fractional value of the variable */
         }
         assert( SCIPisGT(scip, entropyvalue, 0.0) );
      }
      else
      {
         if ( SCIPisGT(scip, entropyvalue, bestvalue) )
         {
            bestvalue = entropyvalue;
            branchvar = lpcands[v]; /* the variable to branch on */
            branchval = lpcandsfrac[v]; /* the fractional value of the variable */
         }
      }
   }

   /* if entropy is less or equal than zero, then we have no fractional values */
   if ( SCIPisLE(scip, bestvalue, 0.0) )
   {
      *result = SCIP_DIDNOTRUN;
   }
   else
   {
      /* create the branch-and-bound tree child nodes of the current node */
      SCIP_CALL( SCIPbranchVarVal(scip, branchvar, branchval, &downchild, NULL, &upchild) );
      assert(downchild != NULL);
      assert(upchild != NULL);

      *result = SCIP_BRANCHED;
   }

   return SCIP_OKAY;
}

/** creates the entropy branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleEntropy(
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

   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpEntropy) );

   return SCIP_OKAY;
}
