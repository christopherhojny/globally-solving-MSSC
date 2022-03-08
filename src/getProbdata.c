/**@file   getProbdata.c
 * @brief  Problem data for clustering problem
 * @author Carina Costa
 *
 * This file checks which model is used and then calls the interface methods of the problem.
 *
 * The list of all interface methods can be found in probdata_clustering.h. and probdata_clustering_quadsoc.h.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "datapoints.h"
#include "getProbdata.h"
#include "probdata_clustering.h"
#include "probdata_clustering_quadsoc.h"
#include "typedefs.h"

#include "scip/scip.h"


/** returns array pointer to data points */
Datapoints* SCIPprobdataGetDatapoints(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   model               /**< model used */
   )
{
   assert( probdata != NULL );

   if ( model == MODEL_QUADSOC )
   {
      return SCIPprobdataGetDatapointsQuadSOC(probdata);
   }
   else
   {
      return SCIPprobdataGetDatapointsQuadratic(probdata);
   }
}


/** returns number of clusters */
int SCIPprobdataGetNClusters(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   model               /**< model used */
   )
{
   assert( probdata != NULL );

   if ( model == MODEL_QUADSOC )
   {
      return SCIPprobdataGetNClustersQuadSOC(probdata);
   }
   else
   {
      return SCIPprobdataGetNClustersQuadratic(probdata);
   }
}

/** returns dimension */
int SCIPprobdataGetDimension(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   model               /**< model used */
   )
{
   assert( probdata != NULL );

   if ( model == MODEL_QUADSOC )
   {
      return SCIPprobdataGetDimensionQuadSOC(probdata);
   }
   else
   {
      return SCIPprobdataGetDimensionQuadratic(probdata);
   }
}

/** returns number of data points */
int SCIPprobdataGetNDatapoints(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   model               /**< model used */
   )
{
   assert( probdata != NULL );

   if ( model == MODEL_QUADSOC )
   {
      return SCIPprobdataGetNDatapointsQuadSOC(probdata);
   }
   else
   {
      return SCIPprobdataGetNDatapointsQuadratic(probdata);
   }
}

/** returns objective variable */
SCIP_VAR* SCIPprobdataGetObjvar(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   model               /**< model used */
   )
{
   assert( probdata != NULL );

   if ( model == MODEL_QUADSOC )
   {
      return NULL;
   }
   else
   {
      return SCIPprobdataGetObjvarQuadratic(probdata);
   }
}

/** returns evars variables */
SCIP_VAR*** SCIPprobdataGetEvars(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   model               /**< model used */
   )
{
   assert( probdata != NULL );

   if ( model == MODEL_QUADSOC )
   {
      return SCIPprobdataGetEvarsQuadSOC(probdata);
   }
   else
   {
      return NULL;
   }
}

/** returns xvars variables */
SCIP_VAR*** SCIPprobdataGetXvars(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   model               /**< model used */
   )
{
   assert( probdata != NULL );

   if ( model == MODEL_QUADSOC )
   {
      return SCIPprobdataGetXvarsQuadSOC(probdata);
   }
   else
   {
      return SCIPprobdataGetXvarsQuadratic(probdata);
   }
}

/** returns cvars variables */
SCIP_VAR*** SCIPprobdataGetCvars(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   model               /**< model used */
   )
{
   assert( probdata != NULL );

   if ( model == MODEL_QUADSOC )
   {
      return SCIPprobdataGetCvarsQuadSOC(probdata);
   }
   else
   {
      return SCIPprobdataGetCvarsQuadratic(probdata);
   }
}

/** returns quadcvars variables */
SCIP_VAR*** SCIPprobdataGetQuadCvars(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   model               /**< model used */
   )
{
   assert( probdata != NULL );

   if ( model == MODEL_QUADSOC )
   {
      return NULL;
   }
   else
   {
      return SCIPprobdataGetQuadCvarsQuadratic(probdata);
   }
}

/** returns kappavars variables */
SCIP_VAR** SCIPprobdataGetKappaVars(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   model               /**< model used */
   )
{
   assert( probdata != NULL );

   if ( model == MODEL_QUADSOC )
   {
      return SCIPprobdataGetKappaVarsQuadSOC(probdata);
   }
   else
   {
      return SCIPprobdataGetKappaVarsQuadratic(probdata);
   }
}

/** returns index of data point associated with a x-variable */
int SCIPvardataGetPoint(
   SCIP_VARDATA*         vardata,            /**< variable data */
   int                   model               /**< model used */
   )
{
   assert( vardata != NULL );

   if ( model == MODEL_QUADSOC )
   {
      return SCIPvardataGetPointQuadSOC(vardata);
   }
   else
   {
      return SCIPvardataGetPointQuadratic(vardata);
   }
}

/** returns index of cluster associated with a x-variable */
int SCIPvardataGetCluster(
   SCIP_VARDATA*         vardata,            /**< variable data */
   int                   model               /**< model used */
   )
{
   assert( vardata != NULL );

   if ( model == MODEL_QUADSOC )
   {
      return SCIPvardataGetClusterQuadSOC(vardata);
   }
   else
   {
      return SCIPvardataGetClusterQuadratic(vardata);
   }
}
