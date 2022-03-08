/**@file   sepa_gradient.h
 * @brief  gradient cut separator
 * @author
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SEPA_GRADIENT_H__
#define __SCIP_SEPA_GRADIENT_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the gradient cut separator and includes it in SCIP */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeSepaGradient(
   SCIP*                 scip                /**< SCIP data structure */
   );


#ifdef __cplusplus
}
#endif

#endif
