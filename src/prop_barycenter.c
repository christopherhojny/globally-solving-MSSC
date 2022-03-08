/**@file   prop_barycenter.c
 * @brief  barycenter propagator
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "datapoints.h"
#include "getProbdata.h"
#include "prop_barycenter.h"
#include "typedefs.h"


/* fundamental propagator properties */
#define PROP_NAME              "barycenter"
#define PROP_DESC              "finds bounds on the centroids and the objective"
#define PROP_PRIORITY                 0 /**< propagator priority */
#define PROP_FREQ                     -1 /**< propagator frequency */
#define PROP_DELAY                FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define PROP_TIMING             SCIP_PROPTIMING_BEFORELP/**< propagation timing mask */

#define DEFAULT_USESTRONGCPROP     FALSE /**< whether strong centroid propagation shall be used */

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
   SCIP_VAR***           xvars;              /**< (ndatapoints x nclusters)-array of assignment var */
   SCIP_VAR***           cvars;              /**< (nclusters x dimension)-array of centroid variables */
   SCIP_VAR*             objvar;             /**< variable modeling objective value */
   SCIP_VAR***           evars;              /**< (ndatapoints x nclusters)-array of objective variables */
   int**                 pointorder;         /**< (dimension x ndatapoints)-array containing ordering
                                              *   or datapoints per coordinate in non-decreasing order */

   SCIP_Bool             usestrongcprop;      /**< whether strong centroid propagation shall be used */
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
   int i;

   assert( scip != NULL );
   assert( propdata != NULL );

   for (i = (*propdata)->dimension - 1; i >= 0; --i)
   {
      SCIPfreeBlockMemoryArray(scip, &(*propdata)->pointorder[i], (*propdata)->ndatapoints);
   }
   SCIPfreeBlockMemoryArray(scip, &(*propdata)->pointorder, (*propdata)->dimension);
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
   SCIP_VAR***           cvars,              /**< array of cluster-centroid-variables */
   SCIP_VAR*             objvar,             /**< variable modeling objective value */
   SCIP_VAR***           evars               /**< (ndatapoints x nclusters)-array of objective variables */
   )
{
   SCIP_Real* coordinates;
   int i;
   int j;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( datapoints != NULL );
   assert( nclusters > 0 );
   assert( ndatapoints > 0 );
   assert( dimension > 0 );
   assert( xvars != NULL );
   assert( cvars != NULL );
   assert( objvar != NULL || evars != NULL );

   propdata->datapoints = datapoints;
   propdata->nclusters = nclusters;
   propdata->ndatapoints = ndatapoints;
   propdata->dimension = dimension;
   propdata->xvars = xvars;
   propdata->cvars = cvars;
   propdata->objvar = objvar;
   propdata->evars = evars;

   /* sort datapoints per coordinate and store ordering */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->pointorder, dimension) );
   for (i = 0; i < dimension; ++i)
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->pointorder[i], ndatapoints) );
      for (j = 0; j < ndatapoints; ++j)
         propdata->pointorder[i][j] = j;
   }
   SCIP_CALL( SCIPallocBufferArray(scip, &coordinates, ndatapoints) );

   /* sort datapoints */
   for (i = 0; i < dimension; ++i)
   {
      for (j = 0; j < ndatapoints; ++j)
         coordinates[j] = datapoints->points[j][i];

      SCIPsortRealInt(coordinates, propdata->pointorder[i], ndatapoints);
   }

   SCIPfreeBufferArray(scip, &coordinates);

   return SCIP_OKAY;
}


/** finds lower bounds on the objective variable(s) */
static
SCIP_RETCODE propagateObjLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of propagator */
   int*                  nfixings,           /**< pointer to store number of fixings found by propagator */
   SCIP_Bool*            infeasible          /**< pointer to store if infeasibility has been detected */
   )
{
   Datapoints* datapoints;
   SCIP_VAR*** xvars;
   SCIP_VAR*** evars;
   SCIP_VAR* objvar;
   SCIP_Real* barycenter;
   int npointsincluster;
   int nclusters;
   int ndatapoints;
   int dimension;
   int k;
   int p;
   int i;
   SCIP_Real lb = 0.0;
   int model;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( nfixings != NULL );
   assert( infeasible != NULL );

   SCIP_CALL( SCIPgetIntParam(scip, "clustering/model", &model) );

   *nfixings = 0;
   *infeasible = FALSE;

   datapoints = propdata->datapoints;
   xvars = propdata->xvars;
   evars = propdata->evars;
   objvar = propdata->objvar;
   nclusters = propdata->nclusters;
   ndatapoints = propdata->ndatapoints;
   dimension = propdata->dimension;

   assert( datapoints != NULL );
   assert( xvars != NULL );
   assert( objvar != NULL || evars != NULL );
   assert( nclusters > 0 );
   assert( ndatapoints > 0 );
   assert( dimension > 0 );

   SCIP_CALL( SCIPallocCleanBufferArray(scip, &barycenter, dimension) );

   /* iterate over all clusters and compute distance to current barycenter */
   for (k = 0; k < nclusters; ++k)
   {
      npointsincluster = 0;

      /* compute barycenter */
      for (p = 0; p < ndatapoints; ++p)
      {
         /* skip data points not (yet) assigned to cluster k */
         if ( SCIPvarGetLbLocal(xvars[p][k]) < 0.5 )
            continue;

         /* update barycenter for each coordinate */
         for (i = 0; i < dimension; ++i)
            barycenter[i] = (datapoints->points[p][i] + npointsincluster * barycenter[i]) / (npointsincluster + 1);
         ++npointsincluster;
      }

      /* compute distance of points to barycenter */
      for (p = 0; p < ndatapoints; ++p)
      {
         /* if the point is not assigned to this cluster, set its e-variable to 0 */
         if ( model == MODEL_QUADSOC && SCIPvarGetUbLocal(xvars[p][k]) < 0.5 )
         {
            if ( SCIPisGT(scip, SCIPvarGetLbLocal(evars[p][k]), 0.0) )
            {
               *infeasible = TRUE;

               for (i = 0; i < dimension; ++i)
                  barycenter[i] = 0.0;

               SCIPfreeCleanBufferArray(scip, &barycenter);
               return SCIP_OKAY;
            }
            else if ( SCIPisGT(scip, SCIPvarGetUbLocal(evars[p][k]), 0.0) )
            {
               SCIP_CALL( SCIPchgVarUb(scip, evars[p][k], 0.0) );
               *nfixings += 1;
            }
         }

         /* skip data points not (yet) assigned to cluster k */
         if ( SCIPvarGetLbLocal(xvars[p][k]) < 0.5 )
            continue;

         for (i = 0; i < dimension; ++i)
            lb += datapoints->points[p][i] * datapoints->points[p][i] - 2 * datapoints->points[p][i] * barycenter[i] + barycenter[i] * barycenter[i];
      }
      for (i = 0; i < dimension; ++i)
         barycenter[i] = 0.0;
   }

   SCIPfreeCleanBufferArray(scip, &barycenter);

   if ( model == MODEL_QUADRATIC )
   {
      /* possibly update lower bound on objective */
      if ( SCIPisGT(scip, lb, SCIPvarGetUbLocal(objvar)) )
         *infeasible = TRUE;
      else if ( SCIPisLT(scip, SCIPvarGetLbLocal(objvar), lb) )
      {
         SCIP_CALL( SCIPchgVarLb(scip, objvar, lb) );
         *nfixings += 1;
      }
   }

   return SCIP_OKAY;
}


/** strengthen upper/lower bounds of centroids based on maximum/minimum coordinate values
 *  of data points assigned to the corresponding cluster */
static
SCIP_RETCODE propagateCentroids(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of propagator */
   int*                  nprop,              /**< pointer to store number of bound tightenings */
   SCIP_Bool*            infeasible          /**< pointer to store if infeasibility has been detected */
   )
{
   Datapoints* datapoints;
   SCIP_VAR*** xvars;
   SCIP_VAR*** cvars;
   SCIP_Real ub;
   SCIP_Real lb;
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

   assert( datapoints != NULL );
   assert( xvars != NULL );
   assert( cvars != NULL );
   assert( nclusters > 0 );
   assert( ndatapoints > 0 );
   assert( dimension > 0 );

   /* for each dimension, find the minimum and maximum coordinate for a cluster k */
   for (k = 0; k < nclusters; ++k)
   {
      for (i = 0; i < dimension; ++i)
      {
         ub = - SCIPinfinity(scip);
         lb = SCIPinfinity(scip);

         for (p = 0; p < ndatapoints; ++p)
         {
            /* skip variables that cannot be part of cluster k */
            if ( SCIPvarGetUbLocal(xvars[p][k]) < 0.5 )
               continue;

            if ( datapoints->points[p][i] > ub )
               ub = datapoints->points[p][i];
            if ( datapoints->points[p][i] < lb )
               lb = datapoints->points[p][i];

            /* stop if ub and lb are at extreme points */
            if ( ub == datapoints->maxvals[i] && lb == datapoints->minvals[i] )
               break;
         }

         /* possibly tighten variable bounds */
         if ( (SCIPisGT(scip, ub, - SCIPinfinity(scip)) && SCIPisLT(scip, ub, SCIPvarGetLbLocal(cvars[k][i])))
            || (SCIPisLT(scip, lb, SCIPinfinity(scip)) &&  SCIPisGT(scip, lb, SCIPvarGetUbLocal(cvars[k][i]))) )
         {
            *infeasible = TRUE;

            return SCIP_OKAY;
         }

         if ( !SCIPisEQ(scip, ub, -SCIPinfinity(scip)) && SCIPisLT(scip, ub, SCIPvarGetUbLocal(cvars[k][i])) )
         {
            SCIP_CALL( SCIPchgVarUb(scip, cvars[k][i], ub) );
            *nprop += 1;
         }
         if ( !SCIPisEQ(scip, lb, SCIPinfinity(scip)) && SCIPisGT(scip, lb, SCIPvarGetLbLocal(cvars[k][i])) )
         {
            SCIP_CALL( SCIPchgVarLb(scip, cvars[k][i], lb) );
            *nprop += 1;
         }
      }
   }

   return SCIP_OKAY;
}


/** strengthen upper/lower bounds of centroids based on maximum/minimum possible value
 *  of centroids based on points that can be assigned to a cluster */
static
SCIP_RETCODE propagateCentroidsStrong(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of propagator */
   int*                  nprop,              /**< pointer to store number of bound tightenings */
   SCIP_Bool*            infeasible          /**< pointer to store if infeasibility has been detected */
   )
{
   Datapoints* datapoints;
   int** pointorder;
   SCIP_VAR*** xvars;
   SCIP_VAR*** cvars;
   SCIP_Real avg;
   SCIP_Real avginit;
   SCIP_Real avgold;
   int nclusters;
   int ndatapoints;
   int dimension;
   int k;
   int p;
   int q;
   int i;
   int npoints;
   int npointsinit;

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
   pointorder = propdata->pointorder;

   assert( datapoints != NULL );
   assert( xvars != NULL );
   assert( cvars != NULL );
   assert( nclusters > 0 );
   assert( ndatapoints > 0 );
   assert( dimension > 0 );
   assert( pointorder != NULL );

   /* for each dimension, find the minimum and maximum coordinate for a cluster k */
   for (k = 0; k < nclusters; ++k)
   {
      for (i = 0; i < dimension; ++i)
      {
         /* compute average of already assigned points */
         avginit = SCIPinfinity(scip);
         npointsinit = 0;

         for (p = 0; p < ndatapoints; ++p)
         {
            if ( SCIPvarGetLbLocal(xvars[p][k]) > 0.5 )
            {
               avginit = avginit * npointsinit + datapoints->points[p][i];
               avginit = avginit / (++npointsinit);
            }
         }

         /* compute lower bound on centroid */
         avg = avginit;
         npoints = npointsinit;
         avgold = avg;

         /* compute smallest possible value if unassigned points are assigned to cluster k */
         for (p = 0; p < ndatapoints; ++p)
         {
            q = pointorder[i][p];

            /* skip already treated or excluded points */
            if ( SCIPvarGetLbLocal(xvars[q][k]) > 0.5 || SCIPvarGetUbLocal(xvars[q][k]) < 0.5 )
               continue;

            avg = avg * npoints + datapoints->points[q][i];
            avg = avg / (++npoints);

            /* stop: avgold is minimum average */
            if ( SCIPisGT(scip, avg, avgold) )
               break;

            avgold = avg;
         }

         /* possibly tighten variable bounds */
         if ( (SCIPisLT(scip, avgold, SCIPinfinity(scip)) && SCIPisGT(scip, avgold, SCIPvarGetUbLocal(cvars[k][i]))) )
         {
            *infeasible = TRUE;

            return SCIP_OKAY;
         }
         else if ( SCIPisLT(scip, avgold, SCIPinfinity(scip)) && SCIPisGT(scip, avgold, SCIPvarGetLbLocal(cvars[k][i])) )
         {
            SCIP_CALL( SCIPchgVarLb(scip, cvars[k][i], avgold) );
            *nprop += 1;
         }

         /* compute largest possible value if unassigned points are assigned to cluster k */
         avg = avginit;
         if ( npointsinit == 0 )
            avg = -SCIPinfinity(scip);
         npoints = npointsinit;
         avgold = avg;

         for (p = ndatapoints - 1; p >= 0; --p)
         {
            q = pointorder[i][p];

            /* skip already treated or excluded points */
            if ( SCIPvarGetLbLocal(xvars[q][k]) > 0.5 || SCIPvarGetUbLocal(xvars[q][k]) < 0.5 )
               continue;

            avg = avg * npoints + datapoints->points[q][i];
            avg = avg / (++npoints);

            if ( SCIPisLT(scip, avg, avgold) )
               break;

            avgold = avg;
         }

         /* possibly tighten variable bounds */
         if ( (SCIPisGT(scip, avgold, -SCIPinfinity(scip)) && SCIPisLT(scip, avgold, SCIPvarGetLbLocal(cvars[k][i]))) )
         {
            *infeasible = TRUE;

            return SCIP_OKAY;
         }
         else if ( (SCIPisGT(scip, avgold, -SCIPinfinity(scip)) && SCIPisLT(scip, avgold, SCIPvarGetUbLocal(cvars[k][i]))) )
         {
            SCIP_CALL( SCIPchgVarUb(scip, cvars[k][i], avgold) );
            *nprop += 1;
         }
      }
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of propagator
 */


/** initialization method of propagator (called after problem was transformed) */
static
SCIP_DECL_PROPINIT(propInitBarycenter)
{  /*lint --e{715}*/
   SCIP_PROBDATA* probdata;
   SCIP_PROPDATA* propdata;
   Datapoints* datapoints;
   int nclusters;
   int ndatapoints;
   int dimension;
   SCIP_VAR*** xvars;
   SCIP_VAR*** cvars;
   SCIP_VAR* objvar;
   SCIP_VAR*** evars;
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
   objvar = SCIPprobdataGetObjvar(probdata, model);
   evars = SCIPprobdataGetEvars(probdata, model);

   /* create barycenter propagator data */
   SCIP_CALL( propdataSet(scip, propdata, datapoints, nclusters, ndatapoints, dimension, xvars, cvars, objvar, evars) );

   return SCIP_OKAY;
}


/** destructor of propagator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PROPFREE(propFreeBarycenter)
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
SCIP_DECL_PROPEXEC(propExecBarycenter)
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

   SCIP_CALL( propagateObjLb(scip, propdata, &nfixings, &infeasible) );

   if ( infeasible )
   {
      *result = SCIP_CUTOFF;

      return SCIP_OKAY;
   }
   else if ( nfixings > 0 )
      *result = SCIP_REDUCEDDOM;

   if ( propdata->usestrongcprop )
   {
      SCIP_CALL( propagateCentroidsStrong(scip, propdata, &nfixings, &infeasible) );
   }
   else
   {
      SCIP_CALL( propagateCentroids(scip, propdata, &nfixings, &infeasible) );
   }

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
SCIP_DECL_PROPRESPROP(propRespropBarycenter)
{  /*lint --e{715}*/
   assert( result != NULL );

   *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/*
 * propagator specific interface methods
 */

/** creates the barycenter propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropBarycenter(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROP* prop;
   SCIP_PROPDATA* propdata;

   SCIP_CALL( SCIPallocBlockMemory(scip, &propdata) );

   /* include propagator */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecBarycenter, propdata) );

   assert( prop != NULL );

   /* set optional callbacks via setter functions */
   SCIP_CALL( SCIPsetPropInit(scip, prop, propInitBarycenter) );
   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeBarycenter) );
   SCIP_CALL( SCIPsetPropResprop(scip, prop, propRespropBarycenter) );

   SCIP_CALL ( SCIPaddBoolParam(scip, "propagators/barycenter/usestrongcprop",
         "Shall strong centroid propagation be used?",
         &propdata->usestrongcprop, TRUE, DEFAULT_USESTRONGCPROP, NULL, NULL) );

   return SCIP_OKAY;
}
