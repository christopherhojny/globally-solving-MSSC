/**@file   getProbdata.h
 * @brief  Problem data for clustering problem
 * @author Carina Costa
 *
 * This file handles the main problem data used in that project.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GETPROBDATA_H
#define GETPROBDATA_H

#include "datapoints.h"

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** returns array pointer to data points */
Datapoints* SCIPprobdataGetDatapoints(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   model               /**< model used */
   );

/** returns number of clusters */
int SCIPprobdataGetNClusters(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   model               /**< model used */
   );

/** returns dimension */
int SCIPprobdataGetDimension(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   model               /**< model used */
   );

/** returns number of data points */
int SCIPprobdataGetNDatapoints(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   model               /**< model used */
   );

/** returns objective variable */
SCIP_VAR* SCIPprobdataGetObjvar(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   model               /**< model used */
   );

/** returns index of data point associated with a x-variable */
int SCIPvardataGetPoint(
   SCIP_VARDATA*         vardata,            /**< variable data */
   int                   model               /**< model used */
   );

/** returns index of cluster associated with a x-variable */
int SCIPvardataGetCluster(
   SCIP_VARDATA*         vardata,            /**< variable data */
   int                   model               /**< model used */
   );

/** returns evars variables */
SCIP_VAR*** SCIPprobdataGetEvars(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   model               /**< model used */
   );

/** returns xvars variables */
SCIP_VAR*** SCIPprobdataGetXvars(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   model               /**< model used */
   );

/** returns cvars variables */
SCIP_VAR*** SCIPprobdataGetCvars(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   model               /**< model used */
   );

/** returns quadcvars variables */
SCIP_VAR*** SCIPprobdataGetQuadCvars(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   model               /**< model used */
   );

/** returns kappavars variables */
SCIP_VAR** SCIPprobdataGetKappaVars(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   model               /**< model used */
   );

#ifdef __cplusplus
}
#endif

#endif
