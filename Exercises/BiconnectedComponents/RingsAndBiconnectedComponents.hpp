#pragma once

/*----------------------------------------------------------------------------*/
/**********   INCLUDES                                               **********/
/*----------------------------------------------------------------------------*/

#include <map>
#include <vector>
#include "Naomini/Forward.hpp"
#include "Naomini/MoleculeDrawer.hpp"
/*----------------------------------------------------------------------------*/
/**********   NAMESPACE                                              **********/
/*----------------------------------------------------------------------------*/

namespace Naomini {

/** Calculates a non-unique, but all-encompassing set of rings
 * of the molecule. Cycles (rings) are stored as atom sets (Ring). Returns a
 * vector of atom sets (RingVector).
 *
 * @brief calculates rings of a molecule
 */
RingVector moleculeGetRings(MoleculePtr mol);

/** Calculates a set of rings of the molecule (@see moleculeGetRings) and
 * extends every ring by the neighbor atoms that are exclusively connected to
 * the ring.
 *
 * @brief calculates rings of a molecule including one-atom ring substituents
 */
RingVector moleculeGetExtendedRings(RingVector rings);

/** Returns all biconnected components of a molecule.
 *
 * @brief calculates biconnected components
 */
BCCVector moleculeGetBiconnectedComponents(MoleculePtr mol);

void moleculeDyeComponents(MoleculePtr mol, RingVector rings, BCCVector bccs,
						AtomVector &dyed_rings, AtomVector &dyed_linkers, AtomVector &dyed_misc);

void moleculeDyeBonds(MoleculePtr mol, MoleculeDrawer drawer, RingVector rings, BCCVector bccs);

/** Wanders through an atom's neighbours recursively
 *  and marks all bonds that are not cyclic.
 *
 * @brief recursive depth first search
 */
void DFS_Visit(AtomPtr atom, AtomPtr parent,
		std::vector<unsigned> &discovery, std::map<AtomPtr, unsigned> &low,
		unsigned &time, std::map<BondPtr,bool> &cyclic);

void DFS_Linker(AtomPtr predecessor, AtomPtr atom, RingVector rings, AtomSet &linker);

bool isCyclic(AtomPtr candidate, RingVector rings);

bool isRingAtom(AtomPtr atom, RingVector rings);

bool isLinkerAtom(AtomPtr atom, BCCVector bccs);

bool hasHydrogen(BondPtr bond);

/*----------------------------------------------------------------------------*/
/**********   NAMESPACE END                                          **********/
/*----------------------------------------------------------------------------*/
}
