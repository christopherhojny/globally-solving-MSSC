/**@file   prop_distance.c
 * @brief  distance propagator
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "datapoints.h"
#include "getProbdata.h"
#include "prop_distance.h"

/* fundamental propagator properties */
#define PROP_NAME              "distance"
#define PROP_DESC              "excludes points from a cluster if they are too far away from a bounding box of the centroid"
#define PROP_PRIORITY                 0 /**< propagator priority */
#define PROP_FREQ                     -1 /**< propagator frequency */
#define PROP_DELAY                FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define PROP_TIMING             SCIP_PROPTIMING_BEFORELP/**< propagation timing mask */


/*
 * Data structures
 */

/** propagator data */
struct SCIP_PropData
{
   Datapoints*           datapoints;         /**< pointer to underlying data points */
   int                   nclusters;          /**< number of clusters */
   int                   ndatapoints;        /**< number of data points */
   int                   dimension;          /**< dimension of data points */
   SCIP_VAR***           xvars;              /**< (ndatapoints x dimension)-array of assignment var */
   SCIP_VAR***           cvars;              /**< (nclusters x dimension)-array of centroid variables */
};


/*
 * Local methods
 */

/** frees data of propagator */
static
SCIP_RETCODE propdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA**       propdata            /**< pointer to propdata */
   )
{
   assert( scip != NULL );
   assert( propdata != NULL );

   SCIPfreeBlockMemory(scip, propdata);

   return SCIP_OKAY;
}


/** set data structure of propagator */
static
SCIP_RETCODE propdataSet(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< pointer to store propagator data */
   Datapoints*           datapoints,         /**< array of data points */
   int                   nclusters,          /**< number of clusters */
   int                   ndatapoints,        /**< number of data points in datapoints */
   int                   dimension,          /**< dimension of points in data points */
   SCIP_VAR***           xvars,              /**< array of point-cluster-assignment-variables */
   SCIP_VAR***           cvars               /**< array of cluster-centroid-variables */
   )
{
   assert( scip != NULL );
   assert( propdata != NULL );
   assert( datapoints != NULL );
   assert( nclusters > 0 );
   assert( ndatapoints > 0 );
   assert( dimension > 0 );
   assert( xvars != NULL );
   assert( cvars != NULL );

   propdata->datapoints = datapoints;
   propdata->nclusters = nclusters;
   propdata->ndatapoints = ndatapoints;
   propdata->dimension = dimension;
   propdata->xvars = xvars;
   propdata->cvars = cvars;

   return SCIP_OKAY;
}


/** fixes x-variables for a cluster to 0 if the corresponding data point is too far away from
 *    a bounding box of the cluster's centroid
 */
static
SCIP_RETCODE propagate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of propagator */
   int*                  nprop,              /**< pointer to store number of bound tightenings */
   SCIP_Bool*            infeasible          /**< pointer to store if infeasibility has been detected */
   )
{
   Datapoints* datapoints;
   SCIP_VAR*** xvars;
   SCIP_VAR*** cvars;
   SCIP_Real** boxubs;
   SCIP_Real** boxlbs;
   SCIP_Real** points;
   SCIP_Real* mindists;
   SCIP_Real minmaxdist;
   SCIP_Real tmpmaxdist;
   int nclusters;
   int ndatapoints;
   int dimension;
   int k;
   int p;
   int i;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( nprop != NULL );
   assert( infeasible != NULL );

   *nprop = 0;
   *infeasible = FALSE;

   datapoints = propdata->datapoints;
   xvars = propdata->xvars;
   cvars = propdata->cvars;
   nclusters = propdata->nclusters;
   ndatapoints = propdata->ndatapoints;
   dimension = propdata->dimension;
   points = datapoints->points;

   assert( datapoints != NULL );
   assert( xvars != NULL );
   assert( cvars != NULL );
   assert( nclusters > 0 );
   assert( ndatapoints > 0 );
   assert( dimension > 0 );
   assert( points != NULL );

   /* collect upper and lower bounds on centroid variables in each dimension */
   SCIP_CALL( SCIPallocBufferArray(scip, &boxubs, nclusters) );
   for (k = 0; k < nclusters; ++k)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &boxubs[k], dimension) );
      for (i = 0; i < dimension; ++i)
         boxubs[k][i] = SCIPvarGetUbLocal(cvars[k][i]);
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &boxlbs, nclusters) );
   for (k = 0; k < nclusters; ++k)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &boxlbs[k], dimension) );
      for (i = 0; i < dimension; ++i)
         boxlbs[k][i] = SCIPvarGetLbLocal(cvars[k][i]);
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &mindists, nclusters) );

   /* iterate over all data points and compute the minimum and maximum distance to a cluster */
   for (p = 0; p < ndatapoints; ++p)
   {
      minmaxdist = SCIPinfinity(scip);

      for (k = 0; k < nclusters; ++k)
      {
         /* compute min/max squared distance of point p to bounding box of centroid for cluster k */
         mindists[k] = 0.0;
         tmpmaxdist = 0.0;
         for (i = 0; i < dimension; ++i)
         {
            /* if the point's i-th coordinate is outside the box... */
            if ( SCIPisGE(scip, points[p][i], boxubs[k][i]) )
            {
               mindists[k] += (points[p][i] - boxubs[k][i])*(points[p][i] - boxubs[k][i]);
               tmpmaxdist += (points[p][i] - boxlbs[k][i])*(points[p][i] - boxlbs[k][i]);
            }
            else if ( SCIPisLE(scip, points[p][i], boxlbs[k][i]) )
            {
               mindists[k] += (points[p][i] - boxlbs[k][i])*(points[p][i] - boxlbs[k][i]);
               tmpmaxdist += (points[p][i] - boxubs[k][i])*(points[p][i] - boxubs[k][i]);
            }
            else
            {
               /* ...or inside the box */
               SCIP_Real lowerdist;
               SCIP_Real upperdist;

               lowerdist = points[p][i] - boxlbs[k][i];
               upperdist = boxubs[k][i] - points[p][i];

               if ( SCIPisLT(scip, lowerdist, upperdist) )
                  tmpmaxdist += upperdist * upperdist;
               else
                  tmpmaxdist += lowerdist * lowerdist;
            }
         }

         /* update the minimum maximum distance */
         if ( SCIPisLT(scip, tmpmaxdist, minmaxdist) )
            minmaxdist = tmpmaxdist;
      }

      /* fix xvars[p][k] to 0 if the minimum distance to cluster is larger than the minimum maximum distance */
      for (k = 0; k < nclusters; ++k)
      {
         if ( SCIPisGT(scip, mindists[k], minmaxdist) )
         {
            if ( SCIPvarGetLbLocal(xvars[p][k]) > 0.5 )
            {
               *infeasible = TRUE;
               goto FREEDATA;
            }
            else if ( SCIPvarGetUbLocal(xvars[p][k]) > 0.5 )
            {
               SCIP_CALL( SCIPchgVarUb(scip, xvars[p][k], 0.0) );
               *nprop += 1;
            }
         }
      }
   }

 FREEDATA:
   SCIPfreeBufferArray(scip, &mindists);
   for (k = nclusters - 1; k >= 0; --k)
   {
      SCIPfreeBufferArray(scip, &boxlbs[k]);
   }
   SCIPfreeBufferArray(scip, &boxlbs);
   for (k = nclusters - 1; k >= 0; --k)
   {
      SCIPfreeBufferArray(scip, &boxubs[k]);
   }
   SCIPfreeBufferArray(scip, &boxubs);

   return SCIP_OKAY;
}


/*
 * Callback methods of propagator
 */


/** initialization method of propagator (called after problem was transformed) */
static
SCIP_DECL_PROPINIT(propInitDistance)
{  /*lint --e{715}*/
   SCIP_PROBDATA* probdata;
   SCIP_PROPDATA* propdata;
   Datapoints* datapoints;
   int nclusters;
   int ndatapoints;
   int dimension;
   SCIP_VAR*** xvars;
   SCIP_VAR*** cvars;
   int model;

   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   SCIP_CALL( SCIPgetIntParam(scip, "clustering/model", &model) );

   datapoints = SCIPprobdataGetDatapoints(probdata, model);
   nclusters = SCIPprobdataGetNClusters(probdata, model);
   ndatapoints = SCIPprobdataGetNDatapoints(probdata, model);
   dimension = SCIPprobdataGetDimension(probdata, model);
   xvars = SCIPprobdataGetXvars(probdata, model);
   cvars = SCIPprobdataGetCvars(probdata, model);

   /* create barycenter propagator data */
   SCIP_CALL( propdataSet(scip, propdata, datapoints, nclusters, ndatapoints, dimension, xvars, cvars) );

   return SCIP_OKAY;
}


/** destructor of propagator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PROPFREE(propFreeDistance)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   assert( scip != NULL );
   assert( prop != NULL );

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   SCIP_CALL( propdataFree(scip, &propdata) );

   return SCIP_OKAY;
}


/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecDistance)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_Bool infeasible = FALSE;
   int nfixings = 0;

   assert( scip != NULL );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   /* do not run if we are in the root or not yet solving */
   if ( SCIPgetDepth(scip) <= 0 || SCIPgetStage(scip) < SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   /* get data */
   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( propagate(scip, propdata, &nfixings, &infeasible) );

   if ( infeasible )
   {
      *result = SCIP_CUTOFF;

      return SCIP_OKAY;
   }
   else if ( nfixings > 0 )
      *result = SCIP_REDUCEDDOM;

   return SCIP_OKAY;
}


/** propagation conflict resolving method of propagator */
static
SCIP_DECL_PROPRESPROP(propRespropDistance)
{  /*lint --e{715}*/
   assert( result != NULL );

   *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/*
 * propagator specific interface methods
 */

/** creates the barycenter propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropDistance(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROP* prop;
   SCIP_PROPDATA* propdata;

   SCIP_CALL( SCIPallocBlockMemory(scip, &propdata) );

   /* include propagator */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecDistance, propdata) );

   assert( prop != NULL );

   /* set optional callbacks via setter functions */
   SCIP_CALL( SCIPsetPropInit(scip, prop, propInitDistance) );
   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeDistance) );
   SCIP_CALL( SCIPsetPropResprop(scip, prop, propRespropDistance) );

   return SCIP_OKAY;
}
