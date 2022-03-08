/**@file   branch_entropy.h
 * @brief  braching rule based on entropy
 * @author
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_ENTROPY_H__
#define __SCIP_BRANCH_ENTROPY_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the entropy branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleEntropy(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
