/**@file   main.cpp
 * @brief  main file for solving the k-means problem
 * @author Christopher Hojny
 */

#include "clusteringParams.h"
#include "clusteringPlugins.h"
#include "datapoints.h"
#include "getProblemName.h"
#include "parseOptions.h"
#include "problem_clustering.h"
#include "readInstance.h"
#include "warmStart.h"

// #include <scip/scipdefplugins.h>

#include <scip/scip.h>

#include <iostream>
#include <limits>

/**@brief Solve partitioning problem given in file @a filename */
static
SCIP_RETCODE solveKMeansClusteringProblem(
   std::string           filename,           //!< name of instance file
   int                   nclusters,          //!< number of clusters
   const char*           settings,           //!< Possible name of setting file
   double                timeLimit,          //!< time limit
   double                memLimit,           //!< memory limit
   SCIP_Longint          nodeLimit,          //!< node limit
   int                   displayFreq         //!< display frequency
   )
{  /*lint --e{429}*/

   // initialize SCIP
   SCIP* scip;
   SCIP_CALL( SCIPcreate(&scip) );

   // output SCIP banner
#ifdef SCIP_THREADSAFE_MESSAGEHDLRS
   SCIPprintVersion(scip, NULL);
#else
   SCIPprintVersion(NULL);
#endif
   std::cout << "\n" << std::endl;   // (to force flush)

   // add our own parameters
   SCIP_CALL( addClusteringParameters(scip) );

   // load basic plugins
   SCIP_CALL( includeClusteringPlugins(scip, nclusters) );

   SCIP_CALL( setSCIPParameters(scip) );

   if ( settings != 0 )
   {
      if ( ! SCIPfileExists(settings) )
      {
         SCIPerrorMessage("Setting file <%s> does not exist.\n\n", settings);
      }
      else
      {
         SCIPinfoMessage(scip, 0, "Reading parameters from <%s>.\n\n", settings);
         SCIP_CALL( SCIPreadParams(scip, settings) );
      }
   }

   // initialize problem class
   std::string problemName = getProblemName(filename.c_str());

   /* output changed parameters */
   SCIPinfoMessage(scip, 0, "Changed settings:\n");
   SCIP_CALL( SCIPwriteParams(scip, 0, FALSE, TRUE) );
   SCIPinfoMessage(scip, 0, "\n");

   // set limits
   if ( timeLimit < 1e20 )
   {
      SCIPinfoMessage(scip, 0, "Setting time limit to %g.\n", timeLimit);
      SCIP_CALL( SCIPsetRealParam(scip, "limits/time", timeLimit) );
   }
   if ( memLimit < 1e20 )
   {
      SCIPinfoMessage(scip, 0, "Setting memory limit to %g.\n", memLimit);
      SCIP_CALL( SCIPsetRealParam(scip, "limits/memory", memLimit) );
   }
   if ( nodeLimit < SCIP_LONGINT_MAX )
   {
      SCIPinfoMessage(scip, 0, "Setting node limit to %ld.\n", nodeLimit);
      SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", nodeLimit) );
   }
   if ( displayFreq < INT_MAX )
   {
      SCIPinfoMessage(scip, 0, "Setting display frequency to %d.\n", displayFreq);
      SCIP_CALL( SCIPsetIntParam(scip, "display/freq", displayFreq) );
   }
   SCIPinfoMessage(scip, 0, "\n");

   // create problem
   SCIP_Bool printinfo;
   Datapoints* datapoints;
   SCIP_CALL( SCIPgetBoolParam(scip, "clustering/printdatainfo", &printinfo) );
   SCIP_CALL( SCIPallocBlockMemory(scip, &datapoints) );
   SCIP_CALL( readInstance(scip, filename, &datapoints, SCIPinfinity(scip), printinfo) );
   SCIP_CALL( SCIPcreateModel(scip, datapoints, nclusters) );

   // disable some plug-ins based on data sepcifications
   SCIP_CALL( setDatabasedParameters(scip, datapoints) );

   // create and add k-means warmstart solution
   SCIP_Bool usewarmstartsol;
   SCIP_CALL( SCIPgetBoolParam(scip, "clustering/usewarmstartsol", &usewarmstartsol) );
   if ( usewarmstartsol )
   {
      SCIP_CALL( addClusteringWarmStart(scip, datapoints, nclusters) );
   }

   // solve the problem ...
   // SCIP_CALL( SCIPwriteOrigProblem(scip, NULL, "lp", FALSE) );
   SCIP_CALL( SCIPsolve(scip) );

   // output statistics
   SCIPinfoMessage(scip, NULL, "\n");
   SCIP_CALL( SCIPprintStatistics(scip, NULL) );

   SCIP_CALL( SCIPprintSol(scip, SCIPgetBestSol(scip), NULL, FALSE) );

   // free data
   SCIP_CALL( SCIPfreeModel(scip, &datapoints, nclusters) );

   SCIP_CALL( SCIPfreeTransform(scip) );
   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   //SCIPprintMemoryDiagnostic(scip);
   return SCIP_OKAY;
}


/** main function for solving the partitioning problem */
int main(int argc, const char** argv)
{
   // check parameters
   parseOption options[] = { { "-s", true, 0 }, { "-t", true, 0 }, { "-m", true, 0 }, { "-n", true, 0 }, { "-d", true, 0 } };
   std::vector<std::string> otherArgs;
   if ( ! extractOptions(argc, argv, options, 2, 2, otherArgs) ) /*lint !e1025*/
   {
      std::cerr << "usage: " << argv[0] << " <file> <k> [-s <settings>] [-t <time limit>] [-m <memory limit>] [-n <node limit>] ";
      std::cerr << "[-d <disp. freq>]" << std::endl;
      exit(1);
   }

   // check for setting file
   const char* settings = 0;
   if ( options[0].value != 0 )
      settings = options[0].value;

   // check for time, memory, and node limit
   double timeLimit = 1e20;
   double memLimit = 1e20;
   SCIP_Longint nodeLimit = SCIP_LONGINT_MAX;
   int dispFreq = INT_MAX;
   if ( options[1].value != 0 )
      timeLimit = atof(options[1].value);
   if ( options[2].value != 0 )
      memLimit = atof(options[2].value);
   if ( options[3].value != 0 )
      nodeLimit = atol(options[3].value);
   if ( options[4].value != 0 )
      dispFreq = atoi(options[4].value);
   else
      dispFreq = 1000;

   // check for number of clusters
   int nclusters = -1;
   if ( otherArgs[1] == "" )
   {
      std::cerr << "Need to specify numbers of clusters <k>." << std::endl;
      exit(1);
   }
   else
      nclusters = atoi(otherArgs[1].c_str());


   // run clustering code
   SCIP_RETCODE retcode;
   retcode = solveKMeansClusteringProblem(otherArgs[0], nclusters, settings, timeLimit, memLimit, nodeLimit, dispFreq);

   if ( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode);
      return -1;
   }
   return 0;
}
