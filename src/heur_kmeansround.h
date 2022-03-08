/**@file    heur_kmeansround.h
 * @brief   kmeans rounding primal heuristic
 * @author  Carina Moreira Costa
 *
 * A rounding heuristic for the kmeans clustering problem.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_KMEANSROUND_H__
#define __SCIP_HEUR_KMEANSROUND_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the rounding primal heuristic and includes it in SCIP */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeHeurKmeansround(
   SCIP*                 scip                //!< SCIP data structure
   );

#ifdef __cplusplus
}
#endif

#endif
