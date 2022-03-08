
/**@file   heur_improve.h
 * @brief  improvement primal heuristic
 * @author Carina Moreira Costa
 *
 * An improvement heuristic for the k-means clustering problem.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_IMPROVE_H__
#define __SCIP_HEUR_IMPROVE_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the improvement primal heuristic and includes it in SCIP */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeHeurImprove(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
