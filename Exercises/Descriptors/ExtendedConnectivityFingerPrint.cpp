/**
 @file
 @brief     
 @author    Bietz
 @date      Mar 29, 2011
 @if DoNotInclude
 Copyright ZBH  Center for Bioinformatics
 University of Hamburg
 Bundesstrasse 43, 20146 Hamburg, Germany
 ================================================================================
 This module is part of the molecule software library Naomi,
 developed at the Center for Bioinformatics Hamburg (ZBH).
 It is not allowed to use, modify, copy or redistribute this software without
 written approval of the ZBH. For further information contact

 ZBH  Center for Bioinformatics, University of Hamburg
 Research Group for Computational Molecular Design

 Voice:  +49 (40) 42838 7350
 Fax:    +49 (40) 42838 7352
 E-Mail: info@zbh.uni-hamburg.de
 ==============================================================================
 @endif
 */

#include "ExtendedConnectivityFingerPrint.hpp"

#include "Naomini/Atom.hpp"
#include "Naomini/Bond.hpp"
#include "Naomini/Molecule.hpp"
#include "Naomini/Helpers.hpp"

using namespace Naomini;

ECFP Naomini::getECFP(MoleculePtr m, unsigned nofIterations) {
  ECFP fingerprint;
  throw "Insert your code in ExtendedConnectivityFingerPrint.cpp";
  return fingerprint;
}
