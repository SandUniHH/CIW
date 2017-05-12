/* Edna Teich, Claudius Sandmeier
 * CIW Blatt 2
 * Abgabe 04.05.2017
 */

#include "CompatibilityGraph.hpp"

#include <stdexcept>
#include <algorithm>			// to find Nodes in a vector

#include "Naomini/Forward.hpp"
#include "Naomini/Molecule.hpp"
#include "Naomini/Bond.hpp"
#include "Naomini/Helpers.hpp" 	// for atomIsHydrogen method


namespace Naomini {

/* Helper Function ************************************************************/
namespace {

bool moleculeHasBond(MoleculePtr mol, AtomPtr atom1, AtomPtr atom2)
{
  for (BondPtr bond : mol->getBonds()) {
    AtomPair atoms = bond->getAtoms();
    if ((atoms.first == atom1 && atoms.second == atom2) ||
        (atoms.first == atom2 && atoms.second == atom1)) {
      return true;
    }
  }
  return false;
}

} // end anonymouse namespace
/* Helper Function END ********************************************************/

/* Compatibility Class Functions **********************************************/

CompatibilityGraph::CompatibilityGraph(MoleculePtr mol1, MoleculePtr mol2)
{
	/* Befuellen des Nodevektors mit jeder moeglichen Atomkombination der beiden Molekuele.
	 * Wasserstoffe und ungleiche Element werden nicht zusammen gespeichert.
	 */
	for (AtomPtr atom_i : mol1->getAtoms())
	{
		if(atomIsHydrogen(atom_i)) continue;

		for (AtomPtr atom_j : mol2->getAtoms())
		{
			if(atom_i->getAtomicNumber() == atom_j->getAtomicNumber())
			{
				nodevector.push_back(Node(atom_i->getAtomicNumber(), atom_i, atom_j));
			}
		}
	}

	/* Befuellen des Edgevektors, wenn bei beiden Nodes moleculeHasBond zutrifft (c-Kanten)
	 * oder bei beiden nicht (d-Kanten).
	 */
	for (Node i : nodevector)
	{
		for (Node j : nodevector)
		{
			/* Ueberspringen, wenn das erste bzw. das zweite Atom der beides Nodes gleich sind,
			 * um horizontale bzw. vertikale Edges zu verhindern. */
			if (i.atom1 == j.atom1 || i.atom2 == j.atom2)
				continue;

			/* c-Kanten, bei beiden Molekuelen eine Kante vorhanden und
			 * d-Kanten, bei beiden Molekuelen KEINE Kante vorhanden */
			if( moleculeHasBond(mol1, i.atom1, j.atom1) ==
				moleculeHasBond(mol2, i.atom2, j.atom2))
			{
				edgevector.push_back(Edge(i, j));
			}
		}
	}
}

/* Destructor */
CompatibilityGraph::~CompatibilityGraph()
{}

size_t CompatibilityGraph::getNofNodes()
{
	return nodevector.size();
}

std::vector<Node> CompatibilityGraph::getNodes() const
{
	return nodevector;
}

Node CompatibilityGraph::getNode(unsigned i) const
{
	return nodevector[i];
}

bool CompatibilityGraph::hasEdge(const Node &a, const Node &b) const
{
	/* Iteriert durch den Edgevektor, bis eine Edge zwischen den beiden Nodes gefunden wird
	 *  und bricht dann ab. */
	return (std::find(edgevector.begin(), edgevector.end(), Edge(a, b)) != edgevector.end());
}
}

