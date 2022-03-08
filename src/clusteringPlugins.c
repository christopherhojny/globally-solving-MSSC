/**@file   clusterPlugins.c
 * @brief  load SCIP plugins for the k-means clustering problem
 * @author Christopher Hojny
 */

#include "clusteringPlugins.h"
#include "heur_improve.h"
#include "heur_kmeansround.h"
#include "prop_barycenter.h"
#include "prop_convexity.h"
#include "prop_convexity_qhull.h"
#include "prop_distance.h"
#include "branch_entropy.h"
#include "sepa_gradient.h"
#include "branch_maxdist.h"
#include "branch_pairs.h"
#include "branch_middle.h"
#include "branch_pairs_midmost.h"
#include "typedefs.h"

#include "scip/scipdefplugins.h"


/** Include basic plugins needed for the k-means clustering problem */
SCIP_RETCODE includeClusteringPlugins(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nclusters           /**< number of clusters */
   )
{
   int model;
   int libconvexhull;

   assert( scip != NULL );

   SCIP_CALL( SCIPgetIntParam(scip, "clustering/model", &model) );
   SCIP_CALL( SCIPgetIntParam(scip, "clustering/libconvexhull", &libconvexhull) );

   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   if ( libconvexhull == QHULL_LIB )
   {
      SCIP_CALL( SCIPincludePropConvexityQhull(scip) );
   }
   else
   {
      SCIP_CALL( SCIPincludePropConvexity(scip) );
   }

   SCIP_CALL( SCIPincludePropBarycenter(scip) );
   SCIP_CALL( SCIPincludePropDistance(scip) );

   SCIP_CALL( SCIPincludeHeurKmeansround(scip) );

   if ( nclusters > 2 )
   {
      SCIP_CALL( SCIPincludeHeurImprove(scip) );
   }

   SCIP_CALL( SCIPincludeBranchruleEntropy(scip) );
   SCIP_CALL( SCIPincludeBranchruleMaxdist(scip) );
   SCIP_CALL( SCIPincludeBranchrulePairs(scip) );
   SCIP_CALL( SCIPincludeBranchruleMiddle(scip) );
   SCIP_CALL( SCIPincludeBranchrulePairsMidmost(scip) );

   if ( model == MODEL_QUADSOC )
   {
      SCIP_CALL( SCIPincludeSepaGradient(scip) );
   }

   return SCIP_OKAY;
}
