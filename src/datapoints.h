/**@file   datapoints.h
 * @brief Declaration of data points
 * @author Christopher Hojny
 *
 * This file contains the definitions of a struct to store data points
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef DATAPOINTS_H
#define DATAPOINTS_H

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** data of a node */
typedef struct Datapoints
{
   int                   dimension;          /**< dimension of set of data points */
   int                   ndatapoints;        /**< number of data points */
   SCIP_Real**           points;             /**< (ndatapoints x dimension)-array of data points */
   SCIP_Real*            minvals;            /**< for each dimension, the minimum value of a data point */
   SCIP_Real*            maxvals;            /**< for each dimension, the maximum value of a data point */
} Datapoints;

#ifdef __cplusplus
}
#endif

#endif

