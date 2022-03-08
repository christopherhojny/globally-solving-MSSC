/**@file   probdata_clustering.h
 * @brief  Problem data for clustering problem
 * @author Christopher Hojny
 *
 * This file handles the main problem data used in that project.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef PROBDATA_CLUSTERING_H
#define PROBDATA_CLUSTERING_H

#include "datapoints.h"

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** sets up the problem data */
SCIP_RETCODE SCIPprobdataCreateQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           probname,           /**< problem name */
   Datapoints*           datapoints,         /**< pointer to underlying data points */
   int                   nclusters           /**< number of clusters */
   );

/** returns array pointer to data points */
Datapoints* SCIPprobdataGetDatapointsQuadratic(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns number of clusters */
int SCIPprobdataGetNClustersQuadratic(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns kappavars variables */
SCIP_VAR** SCIPprobdataGetKappaVarsQuadratic(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns dimension */
int SCIPprobdataGetDimensionQuadratic(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns number of data points */
int SCIPprobdataGetNDatapointsQuadratic(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns objective variable */
SCIP_VAR* SCIPprobdataGetObjvarQuadratic(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns index of data point associated with a x-variable */
int SCIPvardataGetPointQuadratic(
   SCIP_VARDATA*         vardata             /**< variable data */
   );

/** returns index of cluster associated with a x-variable */
int SCIPvardataGetClusterQuadratic(
   SCIP_VARDATA*         vardata             /**< variable data */
   );

/** returns xvars variables */
SCIP_VAR*** SCIPprobdataGetXvarsQuadratic(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns cvars variables */
SCIP_VAR*** SCIPprobdataGetCvarsQuadratic(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns quadcvars variables */
SCIP_VAR*** SCIPprobdataGetQuadCvarsQuadratic(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

#ifdef __cplusplus
}
#endif

#endif
