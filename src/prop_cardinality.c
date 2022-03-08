/**@file   prop_cardinality.c
 * @brief  cardinality propagator
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "datapoints.h"
#include "getProbdata.h"
#include "prop_cardinality.h"

/* fundamental propagator properties */
#define PROP_NAME              "cardinality"
#define PROP_DESC              "updates bound on kappa variables"
#define PROP_PRIORITY                 0 /**< propagator priority */
#define PROP_FREQ                     -1 /**< propagator frequency */
#define PROP_DELAY                FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define PROP_TIMING             SCIP_PROPTIMING_BEFORELP/**< propagation timing mask */

/* optional propagator properties */
#define PROP_PRESOL_PRIORITY         -1 /**< priority of the presolving method (>= 0: before, < 0: after constraint handlers); combined with presolvers */
#define PROP_PRESOLTIMING       SCIP_PRESOLTIMING_MEDIUM /* timing of the presolving method (fast, medium, or exhaustive) */
#define PROP_PRESOL_MAXROUNDS        -1 /**< maximal number of presolving rounds the presolver participates in (-1: no
                                         *   limit) */


/*
 * Data structures
 */

/** propagator data */
struct SCIP_PropData
{
   int                   nclusters;          /**< number of clusters */
   int                   ndatapoints;        /**< number of data points */
   SCIP_VAR***           xvars;              /**< (ndatapoints x nclusters)-array of assignment var */
   SCIP_VAR**            kappavars;          /**< nclusters-array of kappa variables */
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


/** creates data structure of propagator */
static
SCIP_RETCODE propdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA**       propdata            /**< pointer to store propagator data */
   )
{
   assert( scip != NULL );
   assert( propdata != NULL );

   SCIP_CALL( SCIPallocBlockMemory(scip, propdata) );

   return SCIP_OKAY;
}


/** set data structure of propagator */
static
SCIP_RETCODE propdataSet(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< pointer to store propagator data */
   int                   nclusters,          /**< number of clusters */
   int                   ndatapoints,        /**< number of data points in datapoints */
   SCIP_VAR***           xvars,              /**< array of point-cluster-assignment-variables */
   SCIP_VAR**            kappavars           /**< nclusters-array of kappa variables */
   )
{
   assert( scip != NULL );
   assert( propdata != NULL );
   assert( nclusters > 0 );
   assert( ndatapoints > 0 );
   assert( xvars != NULL );
   assert( kappavars != NULL );

   propdata->nclusters = nclusters;
   propdata->ndatapoints = ndatapoints;
   propdata->xvars = xvars;
   propdata->kappavars = kappavars;

   return SCIP_OKAY;
}


/** performs the convex hull and special cone fixings */
static
SCIP_RETCODE propagate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of propagator */
   int*                  nreductions,        /**< pointer to store number of rerductions found by propagator */
   SCIP_Bool*            infeasible          /**< pointer to store if infeasibility has been detected */
   )
{
   SCIP_VAR*** xvars;
   SCIP_VAR** kappavars;
   int ndatapoints;
   int nclusters;
   int c;
   int i;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( nreductions != NULL );
   assert( infeasible != NULL );

   *nreductions = 0;
   *infeasible = FALSE;

   xvars = propdata->xvars;
   kappavars = propdata->kappavars;
   ndatapoints = propdata->ndatapoints;
   nclusters = propdata->nclusters;

   assert( xvars != NULL );
   assert( kappavars != NULL );
   assert( ndatapoints > 0 );
   assert( nclusters > 0 );

   /* iterate over clusters and derive lower and upper bounds on the corresponding kappa-variable */
   for (c = 0; c < nclusters; ++c)
   {
      int ub;
      int lb;
      int nassigned = 0;
      int nnotassigned = 0;

      for (i = 0; i < ndatapoints; ++i)
      {
         if ( SCIPvarGetLbLocal(xvars[i][c]) > 0.5 )
            ++nassigned;
         else if ( SCIPvarGetUbLocal(xvars[i][c]) < 0.5 )
            ++nnotassigned;
      }

      /* determine upper and lower bounds on kappa[c] */
      ub = ndatapoints - nassigned;
      lb = nnotassigned;

      if ( lb > ub )
      {
         *infeasible = TRUE;
         break;
      }

      if ( SCIPvarGetUbLocal(kappavars[c]) > ub )
      {
         SCIP_CALL( SCIPchgVarUb(scip, kappavars[c], ub) );
         (*nreductions)++;
      }
      else if ( SCIPvarGetLbLocal(kappavars[c]) < lb )
      {
         SCIP_CALL( SCIPchgVarLb(scip, kappavars[c], lb) );
         (*nreductions)++;
      }
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of propagator
 */


/** initialization method of propagator (called after problem was transformed) */
static
SCIP_DECL_PROPINIT(propInitCardinality)
{  /*lint --e{715}*/
   SCIP_PROBDATA* probdata;
   SCIP_PROPDATA* propdata;
   int nclusters;
   int ndatapoints;
   SCIP_VAR*** xvars;
   SCIP_VAR** kappavars;
   int model;

   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   SCIP_CALL( SCIPgetIntParam(scip, "clustering/model", &model) );

   nclusters = SCIPprobdataGetNClusters(probdata, model);
   ndatapoints = SCIPprobdataGetNDatapoints(probdata, model);
   xvars = SCIPprobdataGetXvars(probdata, model);
   kappavars = SCIPprobdataGetKappaVars(probdata, model);

   /* create convexity propagator data */
   SCIP_CALL( propdataSet(scip, propdata, nclusters, ndatapoints, xvars, kappavars) );

   return SCIP_OKAY;
}


/** destructor of propagator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PROPFREE(propFreeCardinality)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   assert( scip != NULL );
   assert( prop != NULL );

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   SCIP_CALL( propdataFree(scip, &propdata) );

   return SCIP_OKAY;
}


/** presolving method of propagator */
static
SCIP_DECL_PROPPRESOL(propPresolCardinality)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   int nreductions = 0;
   SCIP_Bool infeasible = FALSE;

   assert( scip != NULL );
   assert( prop != NULL );
   assert( result != NULL );
   assert( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING );

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( propagate(scip, propdata, &nreductions, &infeasible) );

   if ( infeasible )
      *result = SCIP_CUTOFF;
   else if ( nreductions > 0 )
   {
      *result = SCIP_SUCCESS;
      *nchgbds += nreductions;
   }

   return SCIP_OKAY;
}


/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecCardinality)
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
      *result = SCIP_CUTOFF;
   else if ( nfixings > 0 )
      *result = SCIP_REDUCEDDOM;

   return SCIP_OKAY;
}


/** propagation conflict resolving method of propagator */
static
SCIP_DECL_PROPRESPROP(propRespropCardinality)
{  /*lint --e{715}*/
   assert( result != NULL );

   *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/*
 * propagator specific interface methods
 */

/** creates the cardinality propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropCardinality(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROP* prop;
   SCIP_PROPDATA* propdata;

   SCIP_CALL( propdataCreate(scip, &propdata) );

   /* include propagator */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecCardinality, propdata) );

   assert( prop != NULL );

   /* set optional callbacks via setter functions */
   SCIP_CALL( SCIPsetPropInit(scip, prop, propInitCardinality) );
   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeCardinality) );
   SCIP_CALL( SCIPsetPropPresol(scip, prop, propPresolCardinality, PROP_PRESOL_PRIORITY,
         PROP_PRESOL_MAXROUNDS, PROP_PRESOLTIMING) );
   SCIP_CALL( SCIPsetPropResprop(scip, prop, propRespropCardinality) );

   return SCIP_OKAY;
}
