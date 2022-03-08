/**@file   prop_barycenter.h
 * @brief  barycenter propagator
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROP_BARYCENTER_H__
#define __SCIP_PROP_BARYCENTER_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the braycenter propagator and includes it in SCIP */
SCIP_EXPORT
SCIP_RETCODE SCIPincludePropBarycenter(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
