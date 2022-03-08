/**@file   problem_clustering.h
 * @brief  Problem data for clustering problem
 * @author Christopher Hojny
 *
 * This file handles the main problem data used in that project.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef PROBLEM_CLUSTERING_H
#define PROBLEM_CLUSTERING_H

#include "datapoints.h"
#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/* creates initial model of the clustering problem */
extern
SCIP_RETCODE SCIPcreateModel(
   SCIP*                 scip,               /**< SCIP data structure */
   Datapoints*           datapoints,         /**< pointer to data point structure */
   int                   nclusters           /**< number of clusters */
   );

/* free clustering problem data */
extern
SCIP_RETCODE SCIPfreeModel(
   SCIP*                 scip,               /**< SCIP data structure */
   Datapoints**          datapoints,         /**< pointer to data point structure */
   int                   nclusters           /**< number of clusters */
   );

#ifdef __cplusplus
}
#endif

#endif

