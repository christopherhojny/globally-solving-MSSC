/**@file   branch_pairs.c
 * @brief  branching rule based on pairs of points
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <math.h>

#include "scip/cons_linear.h"

#include "branch_pairs.h"
#include "datapoints.h"
#include "getProbdata.h"

#define BRANCHRULE_NAME "pairs"
#define BRANCHRULE_DESC "branching rule based on pairs of points"
#define BRANCHRULE_PRIORITY -500000
#define BRANCHRULE_MAXDEPTH -1
#define BRANCHRULE_MAXBOUNDDIST 1.0

/** branching execution method for fractional LP solutions */
static SCIP_DECL_BRANCHEXECLP(branchExeclpPairs)
{
	SCIP_PROBDATA *probdata;
	SCIP_VAR ***xvars;
	int nclusters;
	int ndatapoints;
	int model;

	SCIP_VARDATA *vardata;
	SCIP_VAR **lpcands;
	SCIP_Real *lpcandsfrac;
	int nlpcands;

	SCIP_Real curval;
	SCIP_Real fracval;
	SCIP_Real bestval = 0.0;
	int clusteridx = -1;
	int pointidx1 = -1;
	int pointidx2 = -1;
	SCIP_Real downbound;
	int curclusteridx;
	int curpointidx1;
	int curpointidx2;

	int **candsbycluster;
	int *ncandsincluster;

	int v;
	int w;
	int k;

	assert(scip != NULL);
	assert(branchrule != NULL);
	assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
	assert(result != NULL);

	*result = SCIP_DIDNOTRUN;

	probdata = SCIPgetProbData(scip);
	assert(probdata != NULL);

	SCIP_CALL(SCIPgetIntParam(scip, "clustering/model", &model));

	xvars = SCIPprobdataGetXvars(probdata, model);
	nclusters = SCIPprobdataGetNClusters(probdata, model);
	ndatapoints = SCIPprobdataGetNDatapoints(probdata, model);

	assert(xvars != NULL);
	assert(nclusters > 0);
	assert(ndatapoints > 0);

	/* get fractional LP candidates */
	SCIP_CALL(SCIPgetLPBranchCands(scip, &lpcands, NULL, &lpcandsfrac, NULL, &nlpcands, NULL));
	assert(nlpcands > 0);

	/* store branching candidates by cluster */
	SCIP_CALL(SCIPallocBufferArray(scip, &candsbycluster, nclusters));
	for (v = 0; v < nclusters; ++v)
	{
		SCIP_CALL(SCIPallocBufferArray(scip, &candsbycluster[v], ndatapoints));
	}
	SCIP_CALL(SCIPallocClearBufferArray(scip, &ncandsincluster, nclusters));

	for (v = 0; v < nlpcands; ++v)
	{
		assert(lpcands[v] != NULL);
		vardata = SCIPvarGetData(lpcands[v]);

		/* skip variables which are not of x-type */
		if (vardata == NULL)
			continue;

		curclusteridx = SCIPvardataGetCluster(vardata, model);
		curpointidx1 = SCIPvardataGetPoint(vardata, model);
		candsbycluster[curclusteridx][ncandsincluster[curclusteridx]] = curpointidx1;
		ncandsincluster[curclusteridx] += 1;
	}

	/* per cluster k, iterate over all pairs (p,q) of candidates and find one such that
    * xvars[p][k] + xvars[q][k] is as fractional as possible
    */

	for (k = 0; k < nclusters; ++k)
	{
		for (v = 0; v < ncandsincluster[k]; ++v)
		{
			curpointidx1 = candsbycluster[k][v];

			for (w = v + 1; w < ncandsincluster[k]; ++w)
			{
				curpointidx2 = candsbycluster[k][w];

				curval = SCIPgetSolVal(scip, NULL, xvars[curpointidx1][k]) + SCIPgetSolVal(scip, NULL, xvars[curpointidx2][k]);

				if (SCIPisLT(scip, curval, 1.0))
					fracval = MIN(curval, 1.0 - curval);
				else
					fracval = MIN(curval - 1.0, 2.0 - curval);

				if (fracval > bestval)
				{
					bestval = fracval;
					pointidx1 = curpointidx1;
					pointidx2 = curpointidx2;
					clusteridx = k;
					downbound = SCIPisLT(scip, curval, 1.0) ? 0.0 : 1.0;
				}
			}
		}
	}

	if (pointidx1 == -1 || SCIPisEQ(scip, bestval, 0.0))
		*result = SCIP_DIDNOTFIND;
	else
	{
		SCIP_NODE *downchild;
		SCIP_NODE *upchild;
		SCIP_CONS *consdownchild;
		SCIP_CONS *consupchild;

		SCIP_VAR **vars;
		SCIP_Real *vals;
		char name[SCIP_MAXSTRLEN];

		/* create two children enforcing x[idx1][k] + x[idx2][k] <= downbound or x[idx1][k] + x[idx2][k] >= downbound + 1 */
		SCIP_CALL(SCIPallocBufferArray(scip, &vars, 2));
		SCIP_CALL(SCIPallocBufferArray(scip, &vals, 2));

		vars[0] = xvars[pointidx1][clusteridx];
		vars[1] = xvars[pointidx2][clusteridx];

		vals[0] = 1.0;
		vals[1] = 1.0;

		(void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "down_%d_%d_%d", pointidx1, pointidx2, (int)downbound);
		SCIP_CALL(SCIPcreateConsLinear(scip, &consdownchild, name, 2, vars, vals, 0.0, downbound,
												 TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE));
		(void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "up_%d_%d_%d", pointidx1, pointidx2, (int)downbound);
		SCIP_CALL(SCIPcreateConsLinear(scip, &consupchild, name, 2, vars, vals, downbound + 1.0, 2.0,
												 TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE));

		/* add constraints to nodes */
		SCIP_CALL(SCIPcreateChild(scip, &downchild, 0.0, SCIPgetLocalTransEstimate(scip)));
		SCIP_CALL(SCIPcreateChild(scip, &upchild, 0.0, SCIPgetLocalTransEstimate(scip)));

		SCIP_CALL(SCIPaddConsNode(scip, downchild, consdownchild, NULL));
		SCIP_CALL(SCIPaddConsNode(scip, upchild, consupchild, NULL));

		/* release constraints */
		SCIP_CALL(SCIPreleaseCons(scip, &consdownchild));
		SCIP_CALL(SCIPreleaseCons(scip, &consupchild));

		*result = SCIP_BRANCHED;

		SCIPfreeBufferArray(scip, &vals);
		SCIPfreeBufferArray(scip, &vars);
	}

	SCIPfreeBufferArray(scip, &ncandsincluster);
	for (v = nclusters - 1; v >= 0; --v)
	{
		SCIPfreeBufferArray(scip, &candsbycluster[v]);
	}
	SCIPfreeBufferArray(scip, &candsbycluster);

	return SCIP_OKAY;
}

/** creates the entropy branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchrulePairs(
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

	SCIP_CALL(SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpPairs));

	return SCIP_OKAY;
}
