#pragma once

/*----------------------------------------------------------------------------*/
/**********   INCLUDES                                               **********/
/*----------------------------------------------------------------------------*/
#include "Naomini/Forward.hpp"
#include <set>

/*----------------------------------------------------------------------------*/
/**********   NAMESPACE                                              **********/
/*----------------------------------------------------------------------------*/
namespace Naomini {

/** This function calculates an Extended Connectivity Fingerprint of a
 * molecule m. The number of iterations (nofIterations) defines how often
 * the iterative update procedure of the ecfp algorithm is executed.
 */
  ECFP getECFP(MoleculePtr m, unsigned nofIterations);

/*----------------------------------------------------------------------------*/
/**********   NAMESPACE END                                          **********/
/*----------------------------------------------------------------------------*/
}
