/**@file   clusteringParams.c
 * @brief  set up parameters for the clustering problem
 * @author Christopher Hojny
 */

#include "clusteringParams.h"
#include "typedefs.h"


/* default settings for parameters - for documentation see below */

#define DEFAULT_PRINTDATAINFO                FALSE
#define DEFAULT_USELOCALIZEDCARDCUT          FALSE
#define DEFAULT_MODEL                            2
#define DEFAULT_DISABLEDIMENSION                 9
#define DEFAULT_SHIFTDATAPOINTS              FALSE
#define DEFAULT_CHECKIFRESCALEDATA           FALSE
#define DEFAULT_SCALERANGE                     1e3
#define DEFAULT_USEWARMSTARTSOL               TRUE
#define DEFAULT_KAPPAISINTEGRAL               TRUE
#define DEFAULT_BRANCHPAIRSMOSTFRACTIONAL     TRUE
#define DEFAULT_BRANCHPAIRSNMIDPOINTS           10
#define DEFAULT_USEMINENTROPY                FALSE
#define DEFAULT_LIBCONVEXHUll            QHULL_LIB

/** Set basic SCIP parameters that are relevant for the clustering problem */
SCIP_RETCODE setSCIPParameters(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert( scip != NULL );

   /* disable some heuristics as they may require too much time */
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/alns/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/mpec/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/subnlp/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/rens/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/locks/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/objpscostdiving/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/distributiondiving/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/clique/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/nlpdiving/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/rins/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/linesearchdiving/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/conflictdiving/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/crossover/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/fracdiving/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/guideddiving/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/pscostdiving/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/randrounding/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/veclendiving/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/adaptivediving/freq", -1) );

   return SCIP_OKAY;
}


/** Set parameters based on problem data */
SCIP_RETCODE setDatabasedParameters(
   SCIP*                 scip,               /**< SCIP data structure */
   Datapoints*           datapoints          /**< data points */
   )
{
   int disabledimension;

   assert( scip != NULL );

   SCIP_CALL( SCIPgetIntParam(scip, "clustering/disabledimconvexityprop", &disabledimension) );

   if ( datapoints->dimension >= disabledimension )
   {
      SCIP_CALL( SCIPsetIntParam(scip, "propagating/convexity/freq", -1) );
   }

   return SCIP_OKAY;
}


/** Introduce parameters that are relevant for the clustering problem */
SCIP_RETCODE addClusteringParameters(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert( scip != NULL );

   /* information regarding problem */
   SCIP_CALL( SCIPaddBoolParam(scip, "clustering/printdatainfo", "Shall information on the data points be printed?",
         NULL, TRUE, DEFAULT_PRINTDATAINFO, NULL, NULL) );

   /* whether we use a localized version of the cardinality cut  */
   SCIP_CALL( SCIPaddBoolParam(scip, "clustering/uselocalizedcardcut",
         "Shall we use a localized version of the cardinality cut?",
         NULL, TRUE, DEFAULT_USELOCALIZEDCARDCUT, NULL, NULL) );

   /* whether cardinality variables are integral  */
   SCIP_CALL( SCIPaddBoolParam(scip, "clustering/kappaisintegral",
         "Whether cardinality variables are integral?",
         NULL, TRUE, DEFAULT_KAPPAISINTEGRAL, NULL, NULL) );

   /* information regarding problem */
   SCIP_CALL( SCIPaddIntParam(scip, "clustering/model", "Model that is used (0: quadratic; 1: SOC; 2: quadSOC)",
         NULL, TRUE, DEFAULT_MODEL, 0, 2, NULL, NULL) );

   /* when we disable prop_convexity */
   SCIP_CALL( SCIPaddIntParam(scip, "clustering/disabledimconvexityprop",
         "If the dimension of datapoints is at least this large, convexity propagator gets disabled.",
         NULL, TRUE, DEFAULT_DISABLEDIMENSION, 1, INT_MAX, NULL, NULL) );

   /* whether data points should be shifted to be centered at the origin */
   SCIP_CALL( SCIPaddBoolParam(scip, "clustering/shiftdatapoints", "Shall datapoints to be centered at the origin?",
         NULL, TRUE, DEFAULT_SHIFTDATAPOINTS, NULL, NULL) );

   /* whether data points should be shifted and if needed rescaled */
   SCIP_CALL( SCIPaddBoolParam(scip, "clustering/checkifrescaledata", "Shall data be rescaled if shifting is not enough?",
         NULL, TRUE, DEFAULT_CHECKIFRESCALEDATA, NULL, NULL) );

   /* in which range the data should be rescaled */
   SCIP_CALL( SCIPaddRealParam(scip, "clustering/scalerange", "The range in which the data should be.",
         NULL, TRUE, DEFAULT_SCALERANGE, 1, 10e6, NULL, NULL) );

   /* whether we use a warm start solution */
   SCIP_CALL( SCIPaddBoolParam(scip, "clustering/usewarmstartsol", "Shall we use a warm start solution?",
         NULL, TRUE, DEFAULT_USEWARMSTARTSOL, NULL, NULL) );

   /* perform branch on most fractional part in branch-pairs-midmost */
   SCIP_CALL( SCIPaddBoolParam(scip, "clustering/branchpairsmostfractional", "Branch on most fractional part in branch pairs?",
         NULL, TRUE, DEFAULT_BRANCHPAIRSMOSTFRACTIONAL, NULL, NULL) );

   /* number of middle points to compute in branch-pairs-midmost */
   SCIP_CALL( SCIPaddIntParam(scip, "clustering/branchpairsnmidpoints", "Number of additional middle points to compute (1 or 10).",
         NULL, TRUE, DEFAULT_BRANCHPAIRSNMIDPOINTS, 1, 10, NULL, NULL) );

   /* whether we use the min (least uncertainty) version of entropy branching rule */
   SCIP_CALL( SCIPaddBoolParam(scip, "clustering/useminentropy", "Shall we use the min version of entropy branching rule?",
         NULL, TRUE, DEFAULT_USEMINENTROPY, NULL, NULL) );

   /* which lib to use to compute convex hull */
   SCIP_CALL( SCIPaddIntParam(scip, "clustering/libconvexhull", "Lib that is used (0: CDD_LIB; 1: QHULL_LIB)",
         NULL, TRUE, DEFAULT_LIBCONVEXHUll, 0, 1, NULL, NULL) );

   return SCIP_OKAY;
}
