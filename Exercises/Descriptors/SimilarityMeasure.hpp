#pragma once

/*----------------------------------------------------------------------------*/
/**********   INCLUDES                                               **********/
/*----------------------------------------------------------------------------*/
#include "Naomini/Forward.hpp"
#include "ExtendedConnectivityFingerPrint.hpp"

/*----------------------------------------------------------------------------*/
/**********   NAMESPACE                                              **********/
/*----------------------------------------------------------------------------*/
namespace Naomini {

/** This function returns the tanimoto coefficient of two ECFP */
double getTanimotoCoefficient(const ECFP &a, const ECFP &b);

/** This function returns the cosine coefficient of two ECFP */
double getCosineCoefficient(const ECFP &a, const ECFP &b);

/** This function returns the hamming distance of two ECFP */
double getHammingDistance(const ECFP &a, const ECFP &b);

/*----------------------------------------------------------------------------*/
/**********   NAMESPACE END                                          **********/
/*----------------------------------------------------------------------------*/
}
