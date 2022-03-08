// -*- C++ -*-

#ifndef READ_INSTANCE_H
#define READ_INSTANCE_H

#include <string>
#include "datapoints.h"

using namespace std;

/** reads a directed graph from a file */
extern
SCIP_RETCODE readInstance(
   SCIP*                 scip,               /**< SCIP data structure */
   string                filename,           /**< name of file encoding data points */
   Datapoints**          datapoints,         /**< pointer to store data points */
   SCIP_Real             infinity,           /**< real value interpreted as infinity */
   SCIP_Bool             printinfo           /**< Should information on instance be printed on screen? */
   );

#endif
