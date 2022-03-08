/**@file   branch_maxdist.h
 * @brief  braching rule based on maximum distance points
 * @author
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_MAXDIST_H__
#define __SCIP_BRANCH_MAXDIST_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the maxdist branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleMaxdist(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
