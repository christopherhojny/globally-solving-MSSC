/**@file   probdata_clustering_quadsoc.c
 * @brief  Problem data for clustering problem
 * @author Carina Costa
 *
 * This file handles the main problem data used in that project.
 *
 * @page CLUSTERING_PROBLEMDATA Main problem data
 *
 * The problem data is accessible in all plugins. The function SCIPgetProbData() returns the pointer to that
 * structure. We use this data structure to store all the information of the kidney exchange problem. Since this structure is
 * not visible in the other plugins, we implemented setter and getter functions to access this data. The problem data
 * structure SCIP_ProbData is shown below.
 *
 * A list of all interface methods can be found in probdata_clustering_quadsoc.h.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "datapoints.h"
#include "probdata_clustering_quadsoc.h"

#include "scip/scip.h"
#include "scip/cons_linear.h"
#include "scip/cons_orbitope.h"
#include "scip/cons_setppc.h"
#include "scip/cons_quadratic.h"

/** @brief Problem data which is accessible in all places
 *
 * This problem data is used to store the input of the clustering problem, all variables which are created, and all
 * constraints.
 */
struct SCIP_ProbData
{
   /* data of the problem */
   Datapoints*           datapoints;         /**< pointer to underlying data points */
   int                   nclusters;          /**< number of clusters */
   int                   ndatapoints;        /**< number of data points */
   int                   dimension;          /**< dimension of data points */

   /* data of the model */
   SCIP_VAR***           xvars;              /**< (ndatapoints x nclusters)-array of assignment variables */
   SCIP_VAR***           cvars;              /**< (nclusters x dimension)-array of centroid variables */
   SCIP_VAR***           evars;              /**< (ndatapoints x nclusters)-array of auxiliary variables */
   SCIP_VAR**            kappavars;          /**< ncluster-array to encode cardinality of clusters */
   SCIP_CONS***          objconss;           /**< (ndatapoints x nclusters)-array of partial objective constraints */
   SCIP_CONS**           partconss;          /**< partitioning constraints for xvars */
   SCIP_CONS*            orbitopecons;       /**< orbitope constraint to handle symmetries */
   SCIP_CONS**           cardinalityconss;   /**< constraints bounding the cardinality of clusters */
   SCIP_CONS**           linkkappaconss;     /**< constraints linking kappavars */
};


/** Variable data for xvars */
struct SCIP_VarData
{
   int                   datapointidx;       /**< index of corresponding data point */
   int                   clusteridx;         /**< index of corresponding cluster */
};


/**@name Local methods
 *
 * @{
 */


/** creates variable data */
static
SCIP_RETCODE vardataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA**        vardata,            /**< pointer to variable data */
   int                   datapointidx,       /**< index of variable's data point */
   int                   clusteridx          /**< index of variable's cluster */
   )
{
   assert( scip != NULL );
   assert( vardata != NULL );
   assert( datapointidx >= 0 );
   assert( clusteridx >= 0 );

   /* allocate memory */
   SCIP_CALL( SCIPallocBlockMemory(scip, vardata) );

   (*vardata)->datapointidx = datapointidx;
   (*vardata)->clusteridx = clusteridx;

   return SCIP_OKAY;
}

/** frees variable data */
static
SCIP_RETCODE vardataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA*         vardata             /**< variable data */
   )
{
   assert( scip != NULL );
   assert( vardata != NULL );

   SCIPfreeBlockMemory(scip, &vardata);

   return SCIP_OKAY;
}

/** creates problem data */
static
SCIP_RETCODE probdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata,           /**< pointer to problem data */
   Datapoints*           datapoints,         /**< pointer to underlying data points */
   int                   nclusters,          /**< number of clusters */
   SCIP_VAR***           xvars,              /**< (ndatapoints x nclusters)-array of assignment variables */
   SCIP_VAR***           cvars,              /**< (nclusters x dimension)-array of centroid variables */
   SCIP_VAR***           evars,              /**< (ndatapoints x nclusters)-array of auxiliary variables */
   SCIP_VAR**            kappavars,          /**< ncluster-array to encode cardinality of clusters */
   SCIP_CONS***          objconss,           /**< (ndatapoints x nclusters)-array of partial objective constraints */
   SCIP_CONS**           partconss,          /**< partitioning constraints for xvars */
   SCIP_CONS*            orbitopecons,       /**< orbitope constraint to handle symmetries */
   SCIP_CONS**           cardinalityconss,   /**< constraints bounding the cardinality of clusters */
   SCIP_CONS**           linkkappaconss      /**< constraints linking kappavars */
   )
{
   int ndatapoints;
   int dimension;
   int i;

   assert( scip != NULL );
   assert( probdata != NULL );
   assert( datapoints != NULL );
   assert( nclusters > 0 );

   /* allocate memory */

   (*probdata)->datapoints = datapoints;
   (*probdata)->nclusters = nclusters;
   (*probdata)->ndatapoints = datapoints->ndatapoints;
   (*probdata)->dimension = datapoints->dimension;

   ndatapoints = datapoints->ndatapoints;
   dimension = datapoints->dimension;

   /* possible copy variable arrays */

   if ( xvars != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->xvars, xvars, ndatapoints) );
      for (i = 0; i < ndatapoints; ++i)
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->xvars[i], xvars[i], nclusters) );
      }
   }
   else
      (*probdata)->xvars = NULL;
   if ( cvars != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->cvars, cvars, nclusters) );
      for (i = 0; i < nclusters; ++i)
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->cvars[i], cvars[i], dimension) );
      }
   }
   else
      (*probdata)->cvars = NULL;

   if ( evars != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->evars, evars, ndatapoints) );
      for (i = 0; i < ndatapoints; ++i)
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->evars[i], evars[i], nclusters) );
      }
   }
   else
      (*probdata)->evars = NULL;
   if ( kappavars != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->kappavars, kappavars, nclusters) );
   }
   else
      (*probdata)->kappavars = NULL;

   /* duplicate arrays */
   if ( partconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->partconss, partconss, ndatapoints) );
   }
   else
      (*probdata)->partconss = NULL;

   if ( objconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->objconss, objconss, ndatapoints) );
      for (i = 0; i < ndatapoints; ++i)
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->objconss[i], objconss[i], nclusters) );
      }
   }
   else
      (*probdata)->objconss = NULL;

   if ( orbitopecons != NULL )
      (*probdata)->orbitopecons = orbitopecons;
   else
      (*probdata)->orbitopecons = NULL;

   if ( cardinalityconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->cardinalityconss, cardinalityconss, nclusters) );
   }
   else
      (*probdata)->cardinalityconss = NULL;
   if ( linkkappaconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->linkkappaconss, linkkappaconss, nclusters) );
   }
   else
      (*probdata)->linkkappaconss = NULL;

   return SCIP_OKAY;
}

/** frees the memory of the given problem data */
static
SCIP_RETCODE probdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata,           /**< pointer to problem data */
   SCIP_Bool             isoriginal          /**< whether we want to free the original problem */
   )
{
   int i;
   int j;
   int nclusters;
   int ndatapoints;
   int dimension;
   SCIP_Bool uselocalizedcut;

   assert( scip != NULL );
   assert( probdata != NULL );
   assert( (*probdata)->datapoints != NULL );
   assert( (*probdata)->nclusters > 0 );

   ndatapoints = (*probdata)->ndatapoints;
   dimension = (*probdata)->dimension;
   nclusters = (*probdata)->nclusters;

   assert( ndatapoints > 0 );
   assert( nclusters > 0 );

   SCIP_CALL( SCIPgetBoolParam(scip, "clustering/uselocalizedcardcut", &uselocalizedcut) );

   /* release all variables */
   for (i = 0; i < ndatapoints; ++i)
   {
      for (j = 0; j < nclusters; ++j)
      {
         if ( isoriginal )
         {
            SCIP_VARDATA* vardata;

            vardata = SCIPvarGetData((*probdata)->xvars[i][j]);
            SCIP_CALL( vardataFree(scip, vardata) );
         }
         SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->xvars[i][j]) );
      }
   }
   for (i = 0; i < nclusters; ++i)
   {
      for (j = 0; j < dimension; ++j)
      {
         SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->cvars[i][j]) );
      }
   }

   for (i = 0; i < ndatapoints; ++i)
   {
      for (j = 0; j < nclusters; ++j)
      {
         SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->evars[i][j]) );
      }
   }

   for (i = 0; i < nclusters; ++i)
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->kappavars[i]) );
   }

   /* release all constraints */
   if ( uselocalizedcut )
   {
      for (i = nclusters - 1; i >=0; --i)
      {
         SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->linkkappaconss[i]) );
      }
      SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->linkkappaconss, nclusters);
   }

   for (i = nclusters - 1; i >=0; --i)
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->cardinalityconss[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->cardinalityconss, nclusters);

   for (i = ndatapoints - 1; i >= 0; --i)
   {
      for (j = 0; j < nclusters; ++j)
      {
         SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->objconss[i][j]) );
      }
         SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->objconss[i], nclusters);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->objconss, ndatapoints);

   for (i = ndatapoints - 1; i >= 0; --i)
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->partconss[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->partconss, ndatapoints);

   SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->orbitopecons) );

   /* free memory of arrays */
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->kappavars, nclusters);

   for (i = ndatapoints - 1; i >= 0; --i)
   {
      SCIPfreeBlockMemoryArray(scip, &(*probdata)->evars[i], nclusters);
   }
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->evars, ndatapoints);

   for (i = nclusters - 1; i >= 0; --i)
   {
      SCIPfreeBlockMemoryArray(scip, &(*probdata)->cvars[i], dimension);
   }
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->cvars, nclusters);

   for (i = ndatapoints - 1; i >= 0; --i)
   {
      SCIPfreeBlockMemoryArray(scip, &(*probdata)->xvars[i], nclusters);
   }
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->xvars, ndatapoints);

   /* free probdata */
   SCIPfreeBlockMemory(scip, probdata);

   return SCIP_OKAY;
}

/** creates the initial variables of the problem */
static
SCIP_RETCODE SCIPcreateVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_Real* minvals;
   SCIP_Real* maxvals;
   int ndatapoints;
   int dimension;
   int nclusters;
   int i;
   int j;
   SCIP_Bool uselocalizedcut;

   assert( scip != NULL );
   assert( probdata != NULL );

   assert( probdata->datapoints != NULL );
   assert( probdata->datapoints->ndatapoints > 0 );
   assert( probdata->datapoints->dimension > 0 );
   assert( probdata->datapoints->minvals != NULL );
   assert( probdata->datapoints->maxvals != NULL );
   assert( probdata->nclusters > 0 );

   ndatapoints = probdata->datapoints->ndatapoints;
   dimension = probdata->datapoints->dimension;
   minvals = probdata->datapoints->minvals;
   maxvals = probdata->datapoints->maxvals;
   nclusters = probdata->nclusters;

   SCIP_CALL( SCIPgetBoolParam(scip, "clustering/uselocalizedcardcut", &uselocalizedcut) );

   /* create objective variables evars */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->evars, ndatapoints) );
   for (i = 0; i < ndatapoints; ++i)
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->evars[i], nclusters) );
      for (j = 0; j < nclusters; ++j)
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "e_%d_%d", i, j);

         SCIP_CALL( SCIPcreateVar(scip, &probdata->evars[i][j], name, 0.0, SCIPinfinity(scip), 1.0,
               SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(scip, probdata->evars[i][j]) );
      }
   }

   /* create assignment variables */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->xvars, ndatapoints) );
   for (i = 0; i < ndatapoints; ++i)
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->xvars[i], nclusters) );
      for (j = 0; j < nclusters; ++j)
      {
         SCIP_VARDATA* vardata;

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x_%d_%d", i, j);

         SCIP_CALL( vardataCreate(scip, &vardata, i, j) );

         SCIP_CALL( SCIPcreateVar(scip, &probdata->xvars[i][j], name, 0.0, 1.0, 0.0,
               SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, vardata) );
         SCIP_CALL( SCIPaddVar(scip, probdata->xvars[i][j]) );
      }
   }

   /* create centroid variables */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->cvars, nclusters) );
   for (i = 0; i < nclusters; ++i)
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->cvars[i], dimension) );
      for (j = 0; j < dimension; ++j)
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "c_%d_%d", i, j);

         SCIP_CALL( SCIPcreateVar(scip, &probdata->cvars[i][j], name, minvals[j], maxvals[j], 0.0,
               SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(scip, probdata->cvars[i][j]) );
      }
   }

   /* create cardinality variables */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->kappavars, nclusters) );
   for (i = 0; i < nclusters; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "kappa_%d", i);

      if ( uselocalizedcut )
      {
         SCIP_VARTYPE vtype;
         SCIP_Bool decision;

         SCIP_CALL( SCIPgetBoolParam(scip, "clustering/kappaisintegral", &decision) );
         vtype = decision ? SCIP_VARTYPE_INTEGER : SCIP_VARTYPE_IMPLINT;

         SCIP_CALL( SCIPcreateVar(scip, &probdata->kappavars[i], name, nclusters - 1, ndatapoints - 1, 0.0,
               vtype, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );

         SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, probdata->kappavars[i]) );
      }
      else
      {
         SCIP_CALL( SCIPcreateVar(scip, &probdata->kappavars[i], name, 0.0, 0.0, 0.0,
               SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      }
      SCIP_CALL( SCIPaddVar(scip, probdata->kappavars[i]) );
   }

   return SCIP_OKAY;
}

/** creates the constraints of the problem */
static
SCIP_RETCODE SCIPcreateConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_Real** points;
   SCIP_Real* minvals;
   SCIP_Real* maxvals;
   SCIP_Real val;
   SCIP_Real M;
   int ndatapoints;
   int dimension;
   int nclusters;
   int i;
   int j;
   int l;
   int c;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_Bool uselocalizedcut;

   assert( scip != NULL );
   assert( probdata != NULL );

   assert( probdata->datapoints != NULL );
   assert( probdata->datapoints->points != NULL );
   assert( probdata->datapoints->ndatapoints > 0 );
   assert( probdata->datapoints->dimension > 0 );
   assert( probdata->nclusters > 0 );
   assert( probdata->datapoints->minvals != NULL );
   assert( probdata->datapoints->maxvals != NULL );

   points = probdata->datapoints->points;
   ndatapoints = probdata->datapoints->ndatapoints;
   dimension = probdata->datapoints->dimension;
   nclusters = probdata->nclusters;
   minvals = probdata->datapoints->minvals;
   maxvals = probdata->datapoints->maxvals;

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, MAX(ndatapoints + 1, nclusters)) );

   /* create partitioning constraints */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->partconss), ndatapoints) );
   for (i = 0; i < ndatapoints; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "partcons_%d", i);

      for (j = 0; j < nclusters; ++j)
         vars[j] = probdata->xvars[i][j];

      SCIP_CALL( SCIPcreateConsBasicSetpart(scip, &probdata->partconss[i], name, nclusters, vars) );
      SCIP_CALL( SCIPaddCons(scip, probdata->partconss[i]) );
   }

   /* compute big M value */
   M = 0.0;
   for (l = 0; l < dimension; ++l)
      M += (maxvals[l] - minvals[l]) * (maxvals[l] - minvals[l]);

   /* create objective epigraph constraint */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->objconss), ndatapoints) );
   for (i = 0; i < ndatapoints; ++i)
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->objconss[i]), nclusters) );
      for (j = 0; j < nclusters; ++j)
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "objconss_%d_%d", i, j);

         SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &probdata->objconss[i][j], name, 0, NULL, NULL, 0, NULL, NULL, NULL, 0.0, SCIPinfinity(scip)) );

         /* add constants */
         val = 0.0;
         for (l = 0; l < dimension; ++l)
            val += points[i][l] * points[i][l];
         SCIPaddConstantQuadratic(scip, probdata->objconss[i][j], M - val);

         /* add linear part in evars */
         SCIP_CALL( SCIPaddLinearVarQuadratic(scip, probdata->objconss[i][j], probdata->evars[i][j], 1.0) );

         /* add linear part in xvars */
         SCIP_CALL( SCIPaddLinearVarQuadratic(scip, probdata->objconss[i][j], probdata->xvars[i][j], -M) );

         for (l = 0; l < dimension; ++l)
         {
            SCIP_CALL( SCIPaddQuadVarQuadratic(scip, probdata->objconss[i][j], probdata->cvars[j][l], 2.0 * points[i][l], -1.0) );
         }

         SCIP_CALL( SCIPaddCons(scip, probdata->objconss[i][j]) );
      }
   }

   /* add orbitope constraint to handle symmetries */
#if SCIP_VERSION <= 703
   SCIP_CALL( SCIPcreateConsBasicOrbitope(scip, &probdata->orbitopecons, "orbitope", probdata->xvars,
         SCIP_ORBITOPETYPE_PARTITIONING, probdata->ndatapoints, probdata->nclusters, TRUE, FALSE) );
#else
   SCIP_CALL( SCIPcreateConsBasicOrbitope(scip, &probdata->orbitopecons, "orbitope", probdata->xvars,
         SCIP_ORBITOPETYPE_PARTITIONING, probdata->ndatapoints, probdata->nclusters, FALSE, TRUE, FALSE, FALSE) );
#endif
   SCIP_CALL( SCIPaddCons(scip, probdata->orbitopecons) );

   /* add cardinality constraints */
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, ndatapoints + 1) );
   for (i = 0; i <= ndatapoints; ++i)
      vals[i] = 1.0;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->cardinalityconss), nclusters) );
   for (c = 0; c < nclusters; ++c)
   {
      for (i = 0; i < ndatapoints; ++i)
         vars[i] = probdata->xvars[i][c];

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cardcons_%d", c);

      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &probdata->cardinalityconss[c], name,
            ndatapoints, vars, vals, 1.0, ndatapoints - nclusters + 1) );
      SCIP_CALL( SCIPaddCons(scip, probdata->cardinalityconss[c]) );
   }

   /* add link-kappa constraints */
   SCIP_CALL( SCIPgetBoolParam(scip, "clustering/uselocalizedcardcut", &uselocalizedcut) );

   if ( uselocalizedcut )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->linkkappaconss), nclusters) );
      for (c = 0; c < nclusters; ++c)
      {
         for (i = 0; i < ndatapoints; ++i)
            vars[i] = probdata->xvars[i][c];
         vars[ndatapoints] = probdata->kappavars[c];

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cardcons_%d", c);

         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &probdata->linkkappaconss[c], name,
               ndatapoints  + 1, vars, vals, ndatapoints, ndatapoints) );
         SCIP_CALL( SCIPaddCons(scip, probdata->linkkappaconss[c]) );
      }
   }

   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}


/**@} */

/**@name Callback methods of problem data
 *
 * @{
 */

/** frees user data of original problem (called when the original problem is freed) */
static
SCIP_DECL_PROBDELORIG(probdelorigClustering)
{
   SCIPdebugMsg(scip, "free original problem data\n");

   SCIP_CALL( probdataFree(scip, probdata, TRUE) );

   return SCIP_OKAY;
}

/** creates user data of transformed problem by transforming the original user problem data
 *  (called after problem was transformed) */
static
SCIP_DECL_PROBTRANS(probtransClustering)
{
   int ndatapoints;
   int dimension;
   int nclusters;
   int i;
   SCIP_Bool uselocalizedcut;

   nclusters = sourcedata->nclusters;
   dimension = sourcedata->datapoints->dimension;
   ndatapoints = sourcedata->datapoints->ndatapoints;

   /* create transform probdata */
   SCIP_CALL( SCIPallocBlockMemory(scip, targetdata) );
   SCIP_CALL( probdataCreate(scip, targetdata, sourcedata->datapoints, sourcedata->nclusters,
         sourcedata->xvars, sourcedata->cvars, sourcedata->evars,
         sourcedata->kappavars, sourcedata->objconss, sourcedata->partconss,
         sourcedata->orbitopecons, sourcedata->cardinalityconss,
         sourcedata->linkkappaconss) );

   SCIP_CALL( SCIPgetBoolParam(scip, "clustering/uselocalizedcardcut", &uselocalizedcut) );

   /* transform all constraints */
   SCIP_CALL( SCIPtransformConss(scip, ndatapoints, (*targetdata)->partconss, (*targetdata)->partconss) );
    for (i = 0; i < ndatapoints; ++i)
   {
      SCIP_CALL( SCIPtransformConss(scip, nclusters, (*targetdata)->objconss[i], (*targetdata)->objconss[i]) );
   }
   SCIP_CALL( SCIPtransformCons(scip, (*targetdata)->orbitopecons, &(*targetdata)->orbitopecons) );
   SCIP_CALL( SCIPtransformConss(scip, nclusters, (*targetdata)->cardinalityconss, (*targetdata)->cardinalityconss) );
   if ( uselocalizedcut )
   {
      SCIP_CALL( SCIPtransformConss(scip, nclusters, (*targetdata)->linkkappaconss, (*targetdata)->linkkappaconss) );
   }

   /* transform all variables */
   for (i = 0; i < ndatapoints; ++i)
   {
      SCIP_CALL( SCIPtransformVars(scip, nclusters, (*targetdata)->xvars[i], (*targetdata)->xvars[i]) );
   }
   for (i = 0; i < nclusters; ++i)
   {
      SCIP_CALL( SCIPtransformVars(scip, dimension, (*targetdata)->cvars[i], (*targetdata)->cvars[i]) );
   }

   for (i = 0; i < ndatapoints; ++i)
   {
      SCIP_CALL( SCIPtransformVars(scip, nclusters, (*targetdata)->evars[i], (*targetdata)->evars[i]) );
   }
   SCIP_CALL( SCIPtransformVars(scip, nclusters, (*targetdata)->kappavars, (*targetdata)->kappavars) );

   return SCIP_OKAY;
}

/** frees user data of transformed problem (called when the transformed problem is freed) */
static
SCIP_DECL_PROBDELTRANS(probdeltransClustering)
{
   SCIPdebugMsg(scip, "free transformed problem data\n");

   SCIP_CALL( probdataFree(scip, probdata, FALSE) );

   return SCIP_OKAY;
}

/** solving process initialization method of transformed data (called before the branch and bound process begins) */
static
SCIP_DECL_PROBINITSOL(probinitsolClustering)
{
   assert(probdata != NULL);

   return SCIP_OKAY;
}

/** solving process deinitialization method of transformed data (called before the branch and bound data is freed) */
static
SCIP_DECL_PROBEXITSOL(probexitsolClustering)
{  /*lint --e{715}*/
   assert(probdata != NULL);

   return SCIP_OKAY;
}

/**@} */


/**@name Interface methods
 *
 * @{
 */

/** sets up the problem data */
SCIP_RETCODE SCIPprobdataCreateQuadSOC(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           probname,           /**< problem name */
   Datapoints*           datapoints,         /**< pointer to underlying data points */
   int                   nclusters           /**< number of clusters */
   )
{
   SCIP_PROBDATA* probdata;

   assert( scip != NULL );
   assert( datapoints != NULL );
   assert( nclusters > 0 );

   /* create problem in SCIP and add non-NULL callbacks via setter functions */
   SCIP_CALL( SCIPcreateProbBasic(scip, probname) );

   SCIP_CALL( SCIPsetProbDelorig(scip, probdelorigClustering) );
   SCIP_CALL( SCIPsetProbTrans(scip, probtransClustering) );
   SCIP_CALL( SCIPsetProbDeltrans(scip, probdeltransClustering) );
   SCIP_CALL( SCIPsetProbInitsol(scip, probinitsolClustering) );
   SCIP_CALL( SCIPsetProbExitsol(scip, probexitsolClustering) );

   /* set objective sense */
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );

   /* create problem data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &probdata) );
   SCIP_CALL( probdataCreate(scip, &probdata, datapoints, nclusters,
         NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   SCIP_CALL( SCIPcreateVariables(scip, probdata) );
   SCIP_CALL( SCIPcreateConstraints(scip, probdata) );

   /* set user problem data */
   SCIP_CALL( SCIPsetProbData(scip, probdata) );

   return SCIP_OKAY;
}



/** returns array pointer to data points */
Datapoints* SCIPprobdataGetDatapointsQuadSOC(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert( probdata != NULL );

   return probdata->datapoints;
}

/** returns number of clusters */
int SCIPprobdataGetNClustersQuadSOC(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert( probdata != NULL );

   return probdata->nclusters;
}

/** returns dimension */
int SCIPprobdataGetDimensionQuadSOC(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert( probdata != NULL );

   return probdata->dimension;
}

/** returns number of data points */
int SCIPprobdataGetNDatapointsQuadSOC(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert( probdata != NULL );

   return probdata->ndatapoints;
}

/** returns evars variables */
SCIP_VAR*** SCIPprobdataGetEvarsQuadSOC(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert( probdata != NULL );

   return probdata->evars;
}

/** returns kappavars variables */
SCIP_VAR** SCIPprobdataGetKappaVarsQuadSOC(
   SCIP_PROBDATA*        probdata            /**< problem data */
)
{
   assert( probdata != NULL );

   return probdata->kappavars;
}

/** returns xvars variables */
SCIP_VAR*** SCIPprobdataGetXvarsQuadSOC(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert( probdata != NULL );

   return probdata->xvars;
}

/** returns cvars variables */
SCIP_VAR*** SCIPprobdataGetCvarsQuadSOC(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert( probdata != NULL );

   return probdata->cvars;
}

/** returns index of data point associated with a x-variable */
int SCIPvardataGetPointQuadSOC(
   SCIP_VARDATA*         vardata             /**< variable data */
   )
{
   assert( vardata != NULL );

   return vardata->datapointidx;
}

/** returns index of cluster associated with a x-variable */
int SCIPvardataGetClusterQuadSOC(
   SCIP_VARDATA*         vardata             /**< variable data */
   )
{
   assert( vardata != NULL );

   return vardata->clusteridx;
}

/**@} */
