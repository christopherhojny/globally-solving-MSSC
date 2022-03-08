/**@file   clusteringPlugins.h
 * @brief  load SCIP plugins for the k-means clustering problem
 * @author Christopher Hojny
 */

#ifndef CLUSTERINGPLUGINS_H
#define CLUSTERINGPLUGINS_H

// SCIP include
#include <scip/scip.h>


#ifdef __cplusplus
extern "C" {
#endif

/** Include basic plugins needed for the kidney exchange problem */
extern
SCIP_RETCODE includeClusteringPlugins(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   ncluster            /**< number of clusters */
   );

#ifdef __cplusplus
}
#endif

#endif
