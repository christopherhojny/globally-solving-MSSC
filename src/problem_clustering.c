#include <scip/scip.h>

#include "problem_clustering.h"
#include "probdata_clustering.h"
#include "probdata_clustering_quadsoc.h"
#include "typedefs.h"

/* creates initial model of the kidney exchange problem */
SCIP_RETCODE SCIPcreateModel(
   SCIP*                 scip,               /**< SCIP data structure */
   Datapoints*           datapoints,         /**< pointer to data point structure */
   int                   nclusters           /**< number of clusters */
   )
{
   int model;

   assert( scip != NULL );
   assert( datapoints != NULL );
   assert( nclusters > 0 );

   SCIP_CALL( SCIPgetIntParam(scip, "clustering/model", &model) );

   if ( model == MODEL_QUADRATIC )
   {
      SCIP_CALL( SCIPprobdataCreateQuadratic(scip, "name", datapoints, nclusters) );
   }
   else if ( model == MODEL_QUADSOC )
   {
      SCIP_CALL( SCIPprobdataCreateQuadSOC(scip, "name", datapoints, nclusters) );
   }
   else
   {
      // later will be SOC
   }

   return SCIP_OKAY;
}

/* free kidney exchange problem data */
SCIP_RETCODE SCIPfreeModel(
   SCIP*                 scip,               /**< SCIP data structure */
   Datapoints**          datapoints,         /**< pointer to data point structure */
   int                   nclusters           /**< number of clusters */
   )
{
   int i;
   int ndatapoints;

   assert( scip != NULL );
   assert( *datapoints != NULL );

   ndatapoints = (*datapoints)->ndatapoints;

   SCIPfreeBlockMemoryArrayNull(scip, &(*datapoints)->maxvals, (*datapoints)->dimension);
   SCIPfreeBlockMemoryArrayNull(scip, &(*datapoints)->minvals, (*datapoints)->dimension);
   for (i = ndatapoints - 1; i >= 0; --i)
   {
      SCIPfreeBlockMemoryArrayNull(scip, &(*datapoints)->points[i], (*datapoints)->dimension);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*datapoints)->points, ndatapoints);

   SCIPfreeBlockMemory(scip, datapoints);

   return SCIP_OKAY;
}