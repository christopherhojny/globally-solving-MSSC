// -*- C++ -*-

#ifndef WARM_START_H
#define WARM_START_H

#include "datapoints.h"


/** create and add warmstart solution for the clustering problem */
extern
SCIP_RETCODE addClusteringWarmStart(
   SCIP*                 scip,               //!< SCIP data structure
   Datapoints*           datapoints,         //!< pointer to data point structure
   int                   nclusters           //!< number of clusters
   );

#endif