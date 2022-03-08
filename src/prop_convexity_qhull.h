/**@file   prop_convexity_qhull.h
 * @brief  convexity propagator using qhull lib
 * @author Carina Costa
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROP_CONVEXITY_QHULL_H__
#define __SCIP_PROP_CONVEXITY_QHULL_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the convexity propagator and includes it in SCIP */
SCIP_EXPORT
SCIP_RETCODE SCIPincludePropConvexityQhull(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
