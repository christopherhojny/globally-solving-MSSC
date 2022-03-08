/**@file   prop_convexity_qhull.cpp
 * @brief  convexity propagator using qhull lib
 * @author Carina Costa
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <vector>
#include <string>

#include "libqhullcpp/RboxPoints.h"
#include "libqhullcpp/Qhull.h"
#include "libqhullcpp/QhullQh.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullVertexSet.h"
#include "libqhullcpp/QhullPoint.h"

#include "datapoints.h"
#include "getProbdata.h"
#include "prop_convexity_qhull.h"

/* fundamental propagator properties */
#define PROP_NAME              "convexity"
#define PROP_DESC              "ensures that clusters are convex"
#define PROP_PRIORITY                 0 //!< propagator priority
#define PROP_FREQ                     -1 //!< propagator frequency
#define PROP_DELAY                FALSE //!< should propagation method be delayed, if other propagators found reductions?
#define PROP_TIMING             SCIP_PROPTIMING_BEFORELP//!< propagation timing mask

/* optional propagator properties */
#define PROP_PRESOL_PRIORITY         -1 //!< priority of the presolving method (>= 0: before, < 0: after constraint handlers); combined with presolvers
#define PROP_PRESOLTIMING       SCIP_PRESOLTIMING_MEDIUM //!< timing of the presolving method (fast, medium, or exhaustive)
#define PROP_PRESOL_MAXROUNDS        -1 //!< maximal number of presolving rounds the presolver participates in (-1: no limit)

#define DEFAULT_USECONEPROP        TRUE //!< whether cone propagation shall be applied
#define DEFAULT_PERFORMPRESOLVING FALSE //!< whether propagation shall be used in presolving

/*
 * Data structures
 */

/** propagator data */
struct SCIP_PropData
{
   Datapoints*           datapoints;         //!< pointer to underlying data points
   int                   nclusters;          //!< number of clusters
   int                   ndatapoints;        //!< number of data points
   int                   dimension;          //!< dimension of data points
   SCIP_VAR***           xvars;              //!< (ndatapoints x dimension)-array of assignment var
   SCIP_Bool             useconeprop;        //!< whether cone propagation shall be applied
   SCIP_Bool             performpresolving;  //!< whether propagation shall be used in presolving
   SCIP_NODE*            lastnode;           //!< last node considered in propagation
};


/*
 * Local methods
 */

/** frees data of propagator */
static
SCIP_RETCODE propdataFree(
   SCIP*                 scip,               //!< SCIP data structure
   SCIP_PROPDATA**       propdata            //!< pointer to propdata
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
   SCIP*                 scip,               //!< SCIP data structure
   SCIP_PROPDATA*        propdata,           //!< pointer to store propagator data
   Datapoints*           datapoints,         //!< array of data points
   int                   nclusters,          //!< number of clusters
   int                   ndatapoints,        //!< number of data points in datapoints
   int                   dimension,          //!< dimension of points in data points
   SCIP_VAR***           xvars,              //!< array of point-cluster-assignment-variables
   SCIP_NODE*            lastnode            //!< last node considered in propagation
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


/** performs the convex hull and special cone fixings */
static
SCIP_RETCODE propagate(
   SCIP*                 scip,               //!< SCIP data structure
   SCIP_PROPDATA*        propdata,           //!< data of propagator
   int*                  nfixings,           //!< pointer to store number of fixings found by propagator
   SCIP_Bool*            infeasible          //!< pointer to store if infeasibility has been detected
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
   int npointsincluster;
   int countnumfacets;

   int nboundchgs;
   int ndatapoints;
   int nclusters;
   int dimension;
   int branchpoint;
   int p;
   int k;
   int i;
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

   // get last branching decision
   node = SCIPgetCurrentNode(scip);
   assert( node != NULL );

   domchg = SCIPnodeGetDomchg(node);
   if ( domchg == NULL )
      return SCIP_OKAY;

   nboundchgs = SCIPdomchgGetNBoundchgs(domchg);
   if ( nboundchgs == 0 )
      return SCIP_OKAY;

   boundchg = SCIPdomchgGetBoundchg(domchg, 0);

   // branching decisions have to be in the beginning of the bound change array
   if ( SCIPboundchgGetBoundchgtype(boundchg) != SCIP_BOUNDCHGTYPE_BRANCHING )
      return SCIP_OKAY;

   branchvar = SCIPboundchgGetVar(boundchg);
   bdchgidx = SCIPvarGetLastBdchgIndex(branchvar);

   // skip non-binary branching variables
   if ( SCIPvarGetType(branchvar) != SCIP_VARTYPE_BINARY )
      return SCIP_OKAY;

   vardata = SCIPvarGetData(branchvar);
   if ( vardata == NULL )
      return SCIP_OKAY;

   branchpoint = SCIPvardataGetPoint(vardata, model);

   // iterate over clusters and check whether we have to perform propagation for branched data point
   for (k = 0; k < nclusters; ++k)
   {
      SCIP_VAR* var;
      SCIP_Real lbbeforebdchg;
      SCIP_Real ubbeforebdchg;
      SCIP_Real lbafterbdchg;
      SCIP_Real ubafterbdchg;

      var = xvars[branchpoint][k];
      lbbeforebdchg = SCIPvarGetLbAtIndex(var, bdchgidx, FALSE);
      ubbeforebdchg = SCIPvarGetUbAtIndex(var, bdchgidx, FALSE);
      lbafterbdchg = SCIPvarGetLbLocal(var);
      ubafterbdchg = SCIPvarGetUbLocal(var);

      // skip variables with the same bounds before and after the bound change
      if ( SCIPisEQ(scip, lbbeforebdchg, lbafterbdchg) && SCIPisEQ(scip, ubbeforebdchg, ubafterbdchg) )
         continue;

      // skip the branching variable if it has already been propagated
      if ( var == branchvar && propdata->lastnode == node )
         continue;

      // whether we perform convex hull fixings: fix points in the convex hull of cluster branchcluster
      if ( lbafterbdchg > 0.5 && lbbeforebdchg < 0.5 )
      {
         // compute convex hull of this cluster
         orgQhull::RboxPoints rbox;
         rbox.setDimension(dimension);
         std::vector<SCIP_Real> generatorpoints;

         // set up vector of generators (points in cluster)
         npointsincluster = 0;
         for (p = 0; p < ndatapoints; ++p)
         {
            if ( SCIPvarGetLbLocal(xvars[p][k]) > 0.5 )
            {
               for (i = 0; i < dimension; ++i)
                  generatorpoints.push_back(datapoints->points[p][i]);
               ++npointsincluster;
            }
         }

         // skip clusters being too small
         if ( npointsincluster <= dimension )
            continue;

         rbox.append(generatorpoints);
         assert( SCIPisEQ(scip, npointsincluster, rbox.count()) );

         // compute facets of convex hull
         // option 'QJ': joggled input to avoid precision problems
         orgQhull::Qhull qhull(rbox, "QJ");

         orgQhull::QhullFacetList facets = qhull.facetList();
         int numfacets = qhull.facetList().count();

         // ignore convex sets that are not full-dimensional
         if ( qhull.hullDimension() < dimension )
            return SCIP_OKAY;

         // number of facets should be at least (dim + 1)
         assert( SCIPisGE(scip, numfacets, dimension + 1) );

         // for each data point not in k, check whether it is contained in the convex hull
         for (p = 0; p < ndatapoints && ! *infeasible; ++p)
         {
            // skip points already contained in cluster k
            if ( SCIPvarGetLbLocal(xvars[p][k]) > 0.5 )
               continue;

            // iterate over inequalities and check whether one is violated
            countnumfacets = 0;
            for (orgQhull::QhullFacetList::iterator f = facets.begin(); f != facets.end(); ++f)
            {
               SCIP_Real value;
               orgQhull::QhullHyperplane hyperplane = (*f).hyperplane();

               value = hyperplane.offset();
               for (i = 0; i < dimension; ++i)
               {
                  value += hyperplane[i] * datapoints->points[p][i];
               }

               // we have found a violated inequality
               if ( SCIPisGT(scip, value, 0.0) )
                  break;

               ++countnumfacets;
            }

            // if point p contained in the convex hull, assign p to cluster k too
            if ( countnumfacets == numfacets )
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
      }
      else if ( propdata->useconeprop && ubafterbdchg < 0.5 && ubbeforebdchg > 0.5 )
      {
         std::vector<SCIP_Real> generatorpoints;

         // perform cone propagation cluster branchcluster
         // Compute convex hull of branch cluster and the special cone with branchpoint.
         // Assign all points in the special cone to be not contained in cluster branchcluster.

         // add apex
         for (i = 0; i < dimension; ++i)
            generatorpoints.push_back(datapoints->points[branchpoint][i]);

         // add rays
         npointsincluster = 0;
         for (p = 0; p < ndatapoints; ++p)
         {
            if ( SCIPvarGetLbLocal(xvars[p][k]) > 0.5 )
            {
               for (i = 0; i < dimension; ++i)
                  generatorpoints.push_back(2 * datapoints->points[branchpoint][i] - datapoints->points[p][i]);
               ++npointsincluster;
            }
         }

         // skip clusters being too small
         if ( npointsincluster <= dimension )
            continue;

         // option 's': cospherical points (to consider points in a cone)
         // Option 'P': defines the apex
         orgQhull::RboxPoints rbox;
         rbox.appendComment(" " + std::to_string(npointsincluster));
         rbox.appendComment(" s P");
         rbox.appendComment(std::to_string(datapoints->points[branchpoint][0]));
         for (i = 1; i < dimension; ++i)
         {
            rbox.appendComment(",");
            rbox.appendComment(std::to_string(datapoints->points[branchpoint][i]));
         }

         rbox.setDimension(dimension);
         rbox.append(generatorpoints);

         assert( SCIPisEQ(scip, npointsincluster + 1, rbox.count()) );

         // compute the polyhedral cone from the point to its horizon facets.
         // option 'QJ': joggled input to avoid precision problems
         // option 'QVn': a facet is good if one of its vertices is point n. Here n = 0.
         orgQhull::Qhull qhull;
         qhull.runQhull(rbox, "QJ QV0");

         orgQhull::QhullFacetList facets = qhull.facetList();
         int numfacets= qhull.facetList().count();

         if ( SCIPisLT(scip, numfacets, dimension) )
            return SCIP_OKAY;

         // ignore convex sets that are not full-dimensional
         if ( qhull.hullDimension() < dimension )
            return SCIP_OKAY;

         // forbid points q in the special cone to be contained in cluster branchcluster
         for (p = 0; p < ndatapoints && ! *infeasible; ++p)
         {
            // skip points already contained in different cluster k
            if ( SCIPvarGetUbLocal(xvars[p][k]) < 0.5 )
               continue;

            countnumfacets = 0;
            // iterate over inequalities and check whether one is violated
            for (orgQhull::QhullFacetList::iterator f = facets.begin(); f != facets.end(); ++f)
            {
               if ( !(*f).isGood() )
                  continue;

               SCIP_Real value;
               orgQhull::QhullHyperplane hyperplane = (*f).hyperplane();

               value = hyperplane.offset();
               for (i = 0; i < dimension; ++i)
               {
                  value += hyperplane[i] * datapoints->points[p][i];
               }

               // we have found a violated inequality
               if ( SCIPisGT(scip, value, 0.0) )
                  break;

               ++countnumfacets;
            }

            // if point p contained in the cone, forbid p to be contained in cluster k
            if ( countnumfacets == numfacets )
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
      }
   }

   propdata->lastnode = node;

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
SCIP_RETCODE SCIPincludePropConvexityQhull(
   SCIP*                 scip                //!< SCIP data structure
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
