
#include <iostream>
#include <fstream>
#include <sstream>

#include "datapoints.h"
#include "readInstance.h"

using namespace std;

/** reads a directed graph from a file */
SCIP_RETCODE readInstance(
   SCIP*                 scip,               /**< SCIP data structure */
   string                filename,           /**< name of file encoding data points */
   Datapoints**          datapoints,         /**< pointer to store data points */
   SCIP_Real             infinity,           /**< real value interpreted as infinity */
   SCIP_Bool             printinfo           /**< Should information on instance be printed on screen? */
   )
{
   int ndatapoints;
   int dimension;
   SCIP_Real val;
   SCIP_Bool shiftdatapoints;
   SCIP_Bool checkifrescaledata;
   SCIP_Real scalerange;

   assert( scip != NULL );
   assert( datapoints != NULL );

   cout << "Reading data points from file" << endl;
   ifstream input(filename);
   if (!input)
   {
      cout << "No input file found." << endl;
   }

   // read number of data points and dimension of data points
   input >> ndatapoints;
   input >> dimension;

   (*datapoints)->ndatapoints = ndatapoints;
   (*datapoints)->dimension = dimension;

   cout << "Read file containing " << ndatapoints << " data points of dimension " << dimension << endl;

   // allocate memory for data points
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*datapoints)->points, ndatapoints) );
   for (int i = 0; i < ndatapoints; ++i)
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*datapoints)->points[i], dimension) );
   }

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*datapoints)->minvals, dimension) );
   for (int i = 0; i < dimension; ++i)
      (*datapoints)->minvals[i] = infinity;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*datapoints)->maxvals, dimension) );
   for (int i = 0; i < dimension; ++i)
      (*datapoints)->maxvals[i] = -infinity;

   // read information from file
   for (int i = 0; i < ndatapoints; ++i)
   {
      for (int j = 0; j < dimension; ++j)
      {
         input >> val;
         (*datapoints)->points[i][j] = val;

         if ( val > (*datapoints)->maxvals[j] )
            (*datapoints)->maxvals[j] = val;
         if ( val < (*datapoints)->minvals[j] )
            (*datapoints)->minvals[j] = val;
      }
   }

   if ( printinfo )
   {
      printf("information on coordinates of data points\n");
      printf("coordinate |  minvalue  | maxvalue\n");
      printf("------------------------------------\n");
      for (int i = 0; i < dimension; ++i)
         printf("%10d | %10.2f | %10.2f\n", i, (*datapoints)->minvals[i], (*datapoints)->maxvals[i]);
   }

   SCIP_CALL( SCIPgetBoolParam(scip, "clustering/shiftdatapoints", &shiftdatapoints) );
   SCIP_CALL( SCIPgetBoolParam(scip, "clustering/checkifrescaledata", &checkifrescaledata) );
   SCIP_CALL( SCIPgetRealParam(scip, "clustering/scalerange", &scalerange) );

   // shift datapoints to be centered at the origin, if shifted data is not in [-scalerange, scalerange] we rescale
   if ( shiftdatapoints )
   {
      SCIP_Bool rescaledata = FALSE;

      // for each coordinate, shift datapoints
      for (int i = 0; i < dimension; ++i)
      {
         SCIP_Real midpoint = 0.5 * ((*datapoints)->maxvals[i] + (*datapoints)->minvals[i]);
         SCIP_Real scaleaftershift = (*datapoints)->maxvals[i] - midpoint;

         (*datapoints)->maxvals[i] = scaleaftershift;
         (*datapoints)->minvals[i] = -scaleaftershift;

         if ( scalerange < scaleaftershift )
         {
            rescaledata = TRUE;
         }

         // shift datapoints
         for (int p = 0; p < ndatapoints; ++p)
         {
            SCIP_Real shiftedval = (*datapoints)->points[p][i] - midpoint;
            (*datapoints)->points[p][i] = shiftedval;
         }
      }
      cout << "Data was shifted." << endl;

      // rescale data to be within the [-scalerange, scalerange]-cube
      if ( rescaledata && checkifrescaledata )
      {
         // find the minimum and maximum values among all coordinates
         SCIP_Real minshiftedrange = infinity;
         SCIP_Real maxshiftedrange = -infinity;
         for (int i = 0; i < dimension; ++i)
         {
            if ( (*datapoints)->maxvals[i] > maxshiftedrange )
               maxshiftedrange = (*datapoints)->maxvals[i];
            if ( (*datapoints)->minvals[i] < minshiftedrange )
               minshiftedrange = (*datapoints)->minvals[i];
         }

         // reset minvals and maxvals
         for (int i = 0; i < dimension; ++i)
            (*datapoints)->minvals[i] = infinity;

         for (int i = 0; i < dimension; ++i)
            (*datapoints)->maxvals[i] = -infinity;

         // compute new scaled data
         for (int i = 0; i < ndatapoints; ++i)
         {
            for (int j = 0; j < dimension; ++j)
            {
               SCIP_Real scaledval = -scalerange + (((*datapoints)->points[i][j] - minshiftedrange) * (2.0 * scalerange)) / (maxshiftedrange - minshiftedrange);
               (*datapoints)->points[i][j] = scaledval;

               if ( scaledval > (*datapoints)->maxvals[j] )
                  (*datapoints)->maxvals[j] = scaledval;
               if ( scaledval < (*datapoints)->minvals[j] )
                  (*datapoints)->minvals[j] = scaledval;
            }
         }
         cout << "Data was rescaled." << endl;
      }
   }

   return SCIP_OKAY;
}
