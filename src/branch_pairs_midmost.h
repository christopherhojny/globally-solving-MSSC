/**@file   branch_pairs_midmost.h
 * @brief  braching rule based on pairs of points (midmost) and clusters
 * @author Carina Costa
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_PAIRS_MIDMOST_H__
#define __SCIP_BRANCH_PAIRS_MIDMOST_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the pairs branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchrulePairsMidmost(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif