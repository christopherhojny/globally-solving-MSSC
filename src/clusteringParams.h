/**@file   clusteringParams.h
 * @brief  set up parameters for the clustering problem
 * @author Christopher Hojny
 */

#ifndef CLUSTERINGPARAMS_H
#define CLUSTERINGPARAMS_H

// SCIP include
#include <scip/scip.h>

#include "datapoints.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Set basic SCIP parameters that are relevant for the clusterin problem */
extern
SCIP_RETCODE setSCIPParameters(
   SCIP*                 scip                /**< SCIP data structure */
   );


/** Set parameters based on problem data */
extern
SCIP_RETCODE setDatabasedParameters(
   SCIP*                 scip,               /**< SCIP data structure */
   Datapoints*           datapoints          /**< data points */
   );

/** Introduce parameters that are relevant for the clustering problem */
extern
SCIP_RETCODE addClusteringParameters(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
