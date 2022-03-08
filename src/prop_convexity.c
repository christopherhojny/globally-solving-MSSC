/**@file   prop_convexity.c
 * @brief  convexity propagator
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "cddlib/setoper.h"
#include "cddlib/cddmp.h"
#include "cddlib/cdd.h"
#include "auxiliary_cdd.h"
#include "math.h"
#include "convex_hull.h"

#include "datapoints.h"
#include "getProbdata.h"
#include "prop_convexity.h"

/* fundamental propagator properties */
#define PROP_NAME              "convexity"
#define PROP_DESC              "ensures that clusters are convex"
#define PROP_PRIORITY                 0 /**< propagator priority */
#define PROP_FREQ                     -1 /**< propagator frequency */
#define PROP_DELAY                FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define PROP_TIMING             SCIP_PROPTIMING_BEFORELP/**< propagation timing mask */

/* optional propagator properties */
#define PROP_PRESOL_PRIORITY         -1 /**< priority of the presolving method (>= 0: before, < 0: after constraint handlers); combined with presolvers */
#define PROP_PRESOLTIMING       SCIP_PRESOLTIMING_MEDIUM /* timing of the presolving method (fast, medium, or exhaustive) */
#define PROP_PRESOL_MAXROUNDS        -1 /**< maximal number of presolving rounds the presolver participates in (-1: no
                                         *   limit) */

#define DEFAULT_USECONEPROP        TRUE /**< whether cone propagation shall be applied */
#define DEFAULT_PERFORMPRESOLVING FALSE /**< whether propagation shall be used in presolving */

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
   SCIP_Bool             useconeprop;        /**< whether cone propagation shall be applied */
   SCIP_Bool             performpresolving;  /**< whether propagation shall be used in presolving */
   SCIP_NODE*            lastnode;           /**< last node considered in propagation */
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
   SCIP_NODE*            lastnode            /**< last node considered in propagation */
   )
{
   assert( scip != NULL );
   assert( propdata != NULL );
   assert( datapoints != NULL );
   assert( nclusters > 0 );
   assert( ndatapoints > 0 );
   assert( dimension > 0 );
   assert( xvars != NULL );

   propdata->datapoints = datapoints;
   propdata->nclusters = nclusters;
   propdata->ndatapoints = ndatapoints;
   propdata->dimension = dimension;
   propdata->xvars = xvars;
   propdata->lastnode = lastnode;

   return SCIP_OKAY;
}


/** computes generators of the special cone of a point and a polytope */
static
dd_MatrixPtr computeSpecialConeGenerators(
   SCIP*                 scip,
   Datapoints*           datapoints,
   int*                  pointsincluster,
   int                   npointsincluster,
   int                   apexidx
   )
{
   dd_MatrixPtr generators;
   dd_rowrange ngenerators;
   dd_colrange dim;
   dd_rowrange i;
   dd_colrange j;

   assert( scip != NULL );
   assert( datapoints != NULL );
   assert( pointsincluster != NULL );
   assert( npointsincluster > 0 );
   assert( apexidx >= 0 );

   /* initialize generators matrix */
   ngenerators = npointsincluster + 1;
   dim = datapoints->dimension + 1;
   generators = dd_CreateMatrix(ngenerators, dim);

   /* add apex (we use homogeneous coordinates, i.e., first coordinate is 1) */
   dd_set_si(generators->matrix[0][0], 1);
   for (j = 1; j < dim; ++j)
      dd_set_si(generators->matrix[0][j], datapoints->points[apexidx][j-1]);

   /* add rays (we use homogeneous coordinates, i.e., first coordinate is 0) */
   for (i = 1; i < ngenerators; ++i)
   {
      dd_set_si(generators->matrix[i][0], 0);
      for (j = 1; j < dim; ++j)
         dd_set_si(generators->matrix[i][j], datapoints->points[apexidx][j-1] - datapoints->points[pointsincluster[i-1]][j-1]);
   }

   return generators;
}


/** performs the convex hull and special cone fixings */
static
SCIP_RETCODE propagate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of propagator */
   int*                  nfixings,           /**< pointer to store number of fixings found by propagator */
   SCIP_Bool*            infeasible          /**< pointer to store if infeasibility has been detected */
   )
{
   Datapoints* datapoints;
   SCIP_NODE* node;
   SCIP_BOUNDCHG* boundchg;
   SCIP_BDCHGIDX* bdchgidx;
   SCIP_DOMCHG* domchg;
   SCIP_VARDATA* vardata;
   SCIP_VAR* branchvar;
   SCIP_VAR*** xvars;
   int* pointsincluster;
   int npointsincluster;

   int nboundchgs;
   int ndatapoints;
   int nclusters;
   int dimension;
   int branchpoint;
   int p;
   int k;
   int model;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( nfixings != NULL );
   assert( infeasible != NULL );

   *nfixings = 0;
   *infeasible = FALSE;

   SCIP_CALL( SCIPgetIntParam(scip, "clustering/model", &model) );

   datapoints = propdata->datapoints;
   xvars = propdata->xvars;
   ndatapoints = propdata->ndatapoints;
   dimension = propdata->dimension;
   nclusters = propdata->nclusters;

   assert( datapoints != NULL );
   assert( xvars != NULL );
   assert( ndatapoints > 0 );
   assert( dimension > 0 );
   assert( nclusters > 0 );

   SCIP_CALL( SCIPallocBufferArray(scip, &pointsincluster, ndatapoints) );

   /* get last branching decision */
   node = SCIPgetCurrentNode(scip);
   assert( node != NULL );

   domchg = SCIPnodeGetDomchg(node);
   if ( domchg == NULL )
      goto FREEPOINTSINCLUSTER;

   nboundchgs = SCIPdomchgGetNBoundchgs(domchg);
   if ( nboundchgs == 0 )
      goto FREEPOINTSINCLUSTER;

   boundchg = SCIPdomchgGetBoundchg(domchg, 0);

   /* branching decisions have to be in the beginning of the bound change array */
   if ( SCIPboundchgGetBoundchgtype(boundchg) != SCIP_BOUNDCHGTYPE_BRANCHING )
      goto FREEPOINTSINCLUSTER;

   branchvar = SCIPboundchgGetVar(boundchg);
   bdchgidx = SCIPvarGetLastBdchgIndex(branchvar);

   /* skip non-binary branching variables */
   if ( SCIPvarGetType(branchvar) != SCIP_VARTYPE_BINARY )
      goto FREEPOINTSINCLUSTER;

   vardata = SCIPvarGetData(branchvar);
   if ( vardata == NULL )
      goto FREEPOINTSINCLUSTER;

   branchpoint = SCIPvardataGetPoint(vardata, model);

   /* iterate over clusters and check whether we have to perform propagation for branched data point */
   for (k = 0; k < nclusters; ++k)
   {
      SCIP_VAR* var;
      SCIP_Real lbbeforebdchg;
      SCIP_Real ubbeforebdchg;
      SCIP_Real lbafterbdchg;
      SCIP_Real ubafterbdchg;
      dd_rowrange rowidx;
      dd_colrange colidx;

      var = xvars[branchpoint][k];
      lbbeforebdchg = SCIPvarGetLbAtIndex(var, bdchgidx, FALSE);
      ubbeforebdchg = SCIPvarGetUbAtIndex(var, bdchgidx, FALSE);
      lbafterbdchg = SCIPvarGetLbLocal(var);
      ubafterbdchg = SCIPvarGetUbLocal(var);

      /* skip variables with the same bounds before and after the bound change */
      if ( SCIPisEQ(scip, lbbeforebdchg, lbafterbdchg) && SCIPisEQ(scip, ubbeforebdchg, ubafterbdchg) )
         continue;

      /* skip the branching variable if it has already been propagated */
      if ( var == branchvar && propdata->lastnode == node )
         continue;

      /* whether we perform convex hull fixings: fix points in the convex hull of cluster branchcluster */
      if ( lbafterbdchg > 0.5 && lbbeforebdchg < 0.5 )
      {
         /* compute convex hull of this cluster (cyclic order of its vertices) */
         dd_MatrixPtr generators;
         dd_MatrixPtr facetsconvexhull;
         SCIP_Bool success;

         npointsincluster = 0;
         for (p = 0; p < ndatapoints; ++p)
         {
            if ( SCIPvarGetLbLocal(xvars[p][k]) > 0.5 )
               pointsincluster[npointsincluster++] = p;
         }

         /* skip clusters being too small */
         if ( npointsincluster <= dimension )
            continue;

         /* set up matrix of generators (points in cluster) */
         dd_set_global_constants();
         generators = constructGeneratorMatrixPoints(datapoints, pointsincluster, npointsincluster);

         facetsconvexhull = computeConvexHullFacets(scip, generators, &success);

         dd_FreeMatrix(generators);

         if ( ! success )
         {
            dd_free_global_constants();

            goto FREEPOINTSINCLUSTER;
         }

         /* ignore convex sets that are not full-dimensional (lineality space is non-empty) */
         if ( set_card(facetsconvexhull->linset) > 0 )
         {
            dd_FreeMatrix(facetsconvexhull);
            dd_free_global_constants();

            goto FREEPOINTSINCLUSTER;
         }

         assert( SCIPisGE(scip, facetsconvexhull->rowsize, 3) );
         /* for each data point not in k, check whether it is contained in the convex hull */
         for (p = 0; p < ndatapoints && ! *infeasible; ++p)
         {
            /* skip points already contained in cluster k */
            if ( SCIPvarGetLbLocal(xvars[p][k]) > 0.5 )
               continue;

            /* iterate over inequalities and check whether one is violated*/
            for (rowidx = 0; rowidx < facetsconvexhull->rowsize; ++rowidx)
            {
               SCIP_Real value;
               value = getReal(facetsconvexhull->matrix[rowidx][0]);

               for (colidx = 1; colidx < facetsconvexhull->colsize; ++colidx)
                  value += getReal(facetsconvexhull->matrix[rowidx][colidx]) * datapoints->points[p][colidx-1];

               /* we have found a violated inequality */
               if ( SCIPisLT(scip, value, 0.0) )
                  break;
            }

            /* if point p contained in the convex hull, assign p to cluster k too */
            if ( rowidx == facetsconvexhull->rowsize )
            {
               if ( SCIPvarGetUbLocal(xvars[p][k]) < 0.5 )
               {
                  *infeasible = TRUE;
                  break;
               }
               else
               {
                  SCIP_CALL( SCIPchgVarLb(scip, xvars[p][k], 1.0) );
                  (*nfixings)++;
               }
            }
         }

         dd_FreeMatrix(facetsconvexhull);
         dd_free_global_constants();
      }
      else if ( propdata->useconeprop && ubafterbdchg < 0.5 && ubbeforebdchg > 0.5 )
      {
         dd_MatrixPtr generators;
         dd_MatrixPtr facetsconvexhull;
         SCIP_Bool success;

         /* perform cone propagation cluster branchcluster
          *
          * Compute convex hull of branch cluster and the special cone with branchpoint.
          * Assign all points in the special cone to be not contained in cluster branchcluster.
          */
         npointsincluster = 0;
         for (p = 0; p < ndatapoints; ++p)
         {
            if ( SCIPvarGetLbLocal(xvars[p][k]) > 0.5 )
               pointsincluster[npointsincluster++] = p;
         }

         /* skip points already contained in cluster k */
         if ( npointsincluster <= dimension )
            continue;

         /* compute special cone */
         dd_set_global_constants();  /* First, this must be called to use cddlib. */

         generators = computeSpecialConeGenerators(scip, datapoints, pointsincluster, npointsincluster, branchpoint);

         facetsconvexhull = computeConvexHullFacets(scip, generators, &success);
         dd_FreeMatrix(generators);

         if ( ! success )
         {
            dd_free_global_constants();

            goto FREEPOINTSINCLUSTER;
         }

         /* ignore convex sets that are not full-dimensional (lineality space is non-empty) */
         if ( set_card(facetsconvexhull->linset) > 0 )
         {
            dd_FreeMatrix(facetsconvexhull);
            dd_free_global_constants();

            goto FREEPOINTSINCLUSTER;
         }

         /* forbid points q in the special cone to be contained in cluster branchcluster */
         for (p = 0; p < ndatapoints && ! *infeasible; ++p)
         {
            /* skip points already contained in different cluster k */
            if ( SCIPvarGetUbLocal(xvars[p][k]) < 0.5 )
               continue;

            /* iterate over inequalities and check whether one is violated */
            for (rowidx = 0; rowidx < facetsconvexhull->rowsize; ++rowidx)
            {
               SCIP_Real value;
               value = getReal(facetsconvexhull->matrix[rowidx][0]);

               for (colidx = 1; colidx < facetsconvexhull->colsize; ++colidx)
                  value += getReal(facetsconvexhull->matrix[rowidx][colidx]) * datapoints->points[p][colidx-1];

               /* we have found a violated inequality (because of numerical inaccuracies in the coefficients,
                * ensure that contained points are sufficiently far away from the boundary)
                */
               if ( SCIPisLT(scip, value, 1e-4) )
                  break;
            }

            /* if point p contained in the convex hull, forbid p to be contained in cluster k */
            if ( rowidx == facetsconvexhull->rowsize )
            {
               if ( SCIPvarGetLbLocal(xvars[p][k]) > 0.5 )
               {
                  *infeasible = TRUE;
                  break;
               }
               else
               {
                  SCIP_CALL( SCIPchgVarUb(scip, xvars[p][k], 0.0) );
                  (*nfixings)++;
               }
            }
         }

         dd_FreeMatrix(facetsconvexhull);
         dd_free_global_constants();
      }
   }

   propdata->lastnode = node;

 FREEPOINTSINCLUSTER:
   SCIPfreeBufferArray(scip, &pointsincluster);

   return SCIP_OKAY;
}


/*
 * Callback methods of propagator
 */


/** initialization method of propagator (called after problem was transformed) */
static
SCIP_DECL_PROPINIT(propInitConvexity)
{  /*lint --e{715}*/
   SCIP_PROBDATA* probdata;
   SCIP_PROPDATA* propdata;
   Datapoints* datapoints;
   int nclusters;
   int ndatapoints;
   int dimension;
   SCIP_VAR*** xvars;
   int model;

   SCIP_CALL( SCIPgetIntParam(scip, "clustering/model", &model) );

   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   datapoints = SCIPprobdataGetDatapoints(probdata, model);
   nclusters = SCIPprobdataGetNClusters(probdata, model);
   ndatapoints = SCIPprobdataGetNDatapoints(probdata, model);
   dimension = SCIPprobdataGetDimension(probdata, model);
   xvars = SCIPprobdataGetXvars(probdata, model);

   /* create convexity propagator data */
   SCIP_CALL( propdataSet(scip, propdata, datapoints, nclusters, ndatapoints, dimension, xvars, NULL) );

   return SCIP_OKAY;
}


/** destructor of propagator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PROPFREE(propFreeConvexity)
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
SCIP_DECL_PROPPRESOL(propPresolConvexity)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   int nfixings = 0;
   SCIP_Bool infeasible = FALSE;

   assert( scip != NULL );
   assert( prop != NULL );
   assert( result != NULL );
   assert( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING );

   *result = SCIP_DIDNOTRUN;

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   if ( ! propdata->performpresolving )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( propagate(scip, propdata, &nfixings, &infeasible) );

   if ( infeasible )
      *result = SCIP_CUTOFF;
   else if ( nfixings > 0 )
   {
      *result = SCIP_SUCCESS;
      *nfixedvars += nfixings;
   }

   return SCIP_OKAY;
}


/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecConvexity)
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

   /* /\* do nothing if we are in a probing node *\/ */
   /* if ( SCIPinProbing(scip) ) */
   /*    return SCIP_OKAY; */

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
SCIP_DECL_PROPRESPROP(propRespropConvexity)
{  /*lint --e{715}*/
   assert( result != NULL );

   *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/*
 * propagator specific interface methods
 */

/** creates the convexity propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropConvexity(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROP* prop;
   SCIP_PROPDATA* propdata;

   SCIP_CALL( SCIPallocBlockMemory(scip, &propdata) );

   /* include propagator */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecConvexity, propdata) );

   assert( prop != NULL );

   /* set optional callbacks via setter functions */
   SCIP_CALL( SCIPsetPropInit(scip, prop, propInitConvexity) );
   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeConvexity) );
   SCIP_CALL( SCIPsetPropPresol(scip, prop, propPresolConvexity, PROP_PRESOL_PRIORITY,
         PROP_PRESOL_MAXROUNDS, PROP_PRESOLTIMING) );
   SCIP_CALL( SCIPsetPropResprop(scip, prop, propRespropConvexity) );

   SCIP_CALL( SCIPaddBoolParam(scip, "propagators/convexity/useconeprop",
         "whether cone propagation shall be applied",
         &propdata->useconeprop, TRUE, DEFAULT_USECONEPROP, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "propagators/convexity/performpresolving",
         "whether propagation shall be used in presolving",
         &propdata->performpresolving, TRUE, DEFAULT_PERFORMPRESOLVING, NULL, NULL) );

   return SCIP_OKAY;
}
