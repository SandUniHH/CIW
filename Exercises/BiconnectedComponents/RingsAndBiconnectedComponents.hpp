/* Edna Teich, Claudius Sandmeier
 * CIW Blatt 3
 * Abgabe 18.05.2017
 */

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
 * vector of atom sets (RingVector). Since the atoms are evaluated using their
 * lowpoint values, several rings but not all of a poly-ring structure may be
 * classified as one ring. Since all rings will be returned, the end result
 * are still all rings of the input molecule.
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

/** Returns the linkers of a molecule. Removing any bond of a linker results in
 *  splitting the molecule in two biconnected components.
 *
 *	Note: Due to the way this method was implemented, the method name and return
 *	type is actually a bit of a misnomer. To reduce redundancy not the whole
 *	component will be returned, only the AtomSets of linkers. Since we colour
 *	the rest of the molecule later, we would otherwise colour some atoms doubly
 *	or would need to split off the linkers from the rest of
 *	the biconnected components first. The result remains the same.
 *
 * @brief calculates linkers of a molecule
 */
BCCVector moleculeGetBiconnectedComponents(MoleculePtr mol);

/* Fills the respective AtomVectors depending of whether the atom is
 * part of a ring, a linker or single chain.
 *
 * @brief sorts atoms into their respective group
 */
void moleculeDyeComponents(MoleculePtr mol, RingVector rings, BCCVector bccs,
							AtomVector &dyed_rings, AtomVector &dyed_linkers,
							AtomVector &dyed_misc);

/* Similar to moleculeDyeComponents, but colours the bonds instead.
 * Also colours directly and doesn't simply return the respective AtomVectors
 *
 * @brief colours bond according to their respective group
 */
void moleculeDyeBonds(MoleculePtr mol, MoleculeDrawer drawer, RingVector rings,
						BCCVector bccs);


/*----------------------------------------------------------------------------*/
/**********   NAMESPACE END                                          **********/
/*----------------------------------------------------------------------------*/
}
