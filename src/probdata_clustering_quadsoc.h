/**@file   probdata_clustering_quadsoc.h
 * @brief  Problem data for clustering problem
 * @author Carina Costa
 *
 * This file handles the main problem data used in that project.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef PROBDATA_CLUSTERING_QUADSOC_H
#define PROBDATA_CLUSTERING_QUADSOC_H

#include "datapoints.h"

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** sets up the problem data */
SCIP_RETCODE SCIPprobdataCreateQuadSOC(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           probname,           /**< problem name */
   Datapoints*           datapoints,         /**< pointer to underlying data points */
   int                   nclusters           /**< number of clusters */
   );

/** returns array pointer to data points */
Datapoints* SCIPprobdataGetDatapointsQuadSOC(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns number of clusters */
int SCIPprobdataGetNClustersQuadSOC(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns dimension */
int SCIPprobdataGetDimensionQuadSOC(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns number of data points */
int SCIPprobdataGetNDatapointsQuadSOC(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns evars variables */
SCIP_VAR*** SCIPprobdataGetEvarsQuadSOC(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns kappa vars variables */
SCIP_VAR** SCIPprobdataGetKappaVarsQuadSOC(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns index of data point associated with a x-variable */
int SCIPvardataGetPointQuadSOC(
   SCIP_VARDATA*         vardata             /**< variable data */
   );

/** returns index of cluster associated with a x-variable */
int SCIPvardataGetClusterQuadSOC(
   SCIP_VARDATA*         vardata             /**< variable data */
   );

/** returns xvars variables */
SCIP_VAR*** SCIPprobdataGetXvarsQuadSOC(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );


/** returns cvars variables */
SCIP_VAR*** SCIPprobdataGetCvarsQuadSOC(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );


#ifdef __cplusplus
}
#endif

#endif