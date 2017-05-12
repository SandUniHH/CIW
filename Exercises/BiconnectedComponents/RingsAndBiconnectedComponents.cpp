#include "RingsAndBiconnectedComponents.hpp"

#include "Naomini/Forward.hpp"
#include "Naomini/Molecule.hpp"
#include "Naomini/Atom.hpp"
#include "Naomini/Bond.hpp"
#include "Naomini/Helpers.hpp"

using namespace Naomini;

RingVector Naomini::moleculeGetRings(MoleculePtr mol){

	// Insert your code here but do not use this function for your solution!
	//RingVector rings = getRingsOfMolecule(mol);

	RingVector rings;
	AtomSet atoms;

	std::vector<unsigned> discovery; // storing the discovery time
	std::map<AtomPtr, unsigned> low; // storing the lowpoint value
	std::map<BondPtr, bool> cyclic;
	unsigned time = 0;

	for (AtomPtr atom : mol->getAtoms())
	{
		if(atomIsHydrogen(atom)) continue;

		/* Since the molecules are sorted by decreasing atomic number,
		 * we can simply push 0s in the vectors and know that the index will be correct later on.
		 * This is not relevant yet, but later on hydrogen getID() returns
		 * would be problematic if their IDs weren't higher than all the other atoms'.
		 */
		discovery.push_back(0);
		low[atom] = 0;

		atoms.insert(atom);
	}

	/* set all bonds as cyclic for now */
	for (BondPtr bond : mol->getBonds())
	{
		AtomPair atom_pair = bond->getAtoms();

		if (atomIsHydrogen(atom_pair.first) ||
			atomIsHydrogen(atom_pair.second)) continue;

		cyclic[bond] = true;
	}

	/* recursively call DFS_Visit and set bonds as acyclic if applicable */
	for (AtomPtr atom : atoms)
	{
		if (discovery[atom->getID()] == 0)
			Naomini::DFS_Visit(atom, NULL, discovery, low, time, cyclic);
	}

	/* Sort and collect the cyclic atoms with equal lowpoints */
	for (unsigned lowpoint = 0; lowpoint <= time; lowpoint++)
	{
		Ring ring;
		ring.clear();

		/* go through each entry in the lowpoint map */
		for (const auto &entry : low)
		{
			/* check all the bonds for cyclic property */
			for (BondPtr bond : entry.first->getBonds())
			{
				if (cyclic[bond])
				{
					AtomPair atom_pair = bond->getAtoms();

					if (entry.second == lowpoint)
					{
						ring.insert(atom_pair.first);
						ring.insert(atom_pair.second);
					}
					//current_set.insert(atom_pair.first);
					//current_set.insert(atom_pair.second);

					//cyclic_atoms.insert(atom_pair.first);
					//cyclic_atoms[lowpoint] = current_set;
				}
			}
		}

		if(!ring.empty())
			rings.push_back(ring);
	}

	return rings;
}

/*----------------------------------------------------------------------------*/

RingVector Naomini::moleculeGetExtendedRings(MoleculePtr mol){

	RingVector rings;
	RingVector rings_to_extend = Naomini::moleculeGetRings(mol);
	unsigned neighbour_counter;

	for (Ring ring : rings_to_extend)
	{
		for (AtomPtr atom : ring)
		{
			Ring extRing = ring;

			for (AtomPtr neighbour : atom->getNeighborAtoms())
			{
				neighbour_counter = 0;

				if (atomIsHydrogen(neighbour)) continue;

				for (AtomPtr next_neighbour : neighbour->getNeighborAtoms())
				{
					if (!atomIsHydrogen(next_neighbour)) neighbour_counter++;
				}
				if (neighbour_counter == 1) extRing.insert(neighbour);
			}
			rings.push_back(extRing);
		}
	}


	return rings;
}

/*----------------------------------------------------------------------------*/

BCCVector Naomini::moleculeGetBiconnectedComponents(MoleculePtr mol){
	//throw "Insert your code in RingsAndBiconnectedComponents.cpp";

	BCCVector allBCCs;
	RingVector rings = Naomini::moleculeGetRings(mol);
/*
	for (Ring ring : rings)
	{
		for (AtomPtr atom : ring)
		{
			AtomSet linkers;
			linkers.clear();

			for (AtomPtr neighbour : atom->getNeighborAtoms())
			{
				linkers.insert(findOtherRing(atom, neighbour, rings));
			}

			if(!linkers.empty())
				allBCCs.push_back(linkers);
		}
	}
*/

	return allBCCs;
}

/*----------------------------------------------------------------------------*/
/* Helper functions															  */
/*----------------------------------------------------------------------------*/

void Naomini::DFS_Visit(AtomPtr atom, AtomPtr parent,
		std::vector<unsigned> &discovery, std::map<AtomPtr, unsigned> &low,
		unsigned &time, std::map<BondPtr,bool> &cyclic)
{
	unsigned atom_id, adj_id;

	atom_id = atom->getID();
	discovery[atom_id] = ++time;
	low[atom] = discovery[atom_id];

	for (AtomPtr adj_atom : atom->getNeighborAtoms())
	{
		if(atomIsHydrogen(adj_atom)) continue;

		adj_id = adj_atom->getID();

		if (discovery[adj_id] == 0)
		{
			DFS_Visit(adj_atom, atom, discovery, low, time, cyclic);
			low[atom] = std::min(low[atom], low[adj_atom]);

			if(low[adj_atom] > discovery[atom_id])
			{
				for (BondPtr bond : adj_atom->getBonds())
				{
					AtomPair atom_pair = bond->getAtoms();
					if ((atom_pair.first == adj_atom && atom_pair.second == atom) ||
						(atom_pair.first == atom && atom_pair.second == adj_atom))
					{
						cyclic[bond] = false;
						break;
					}
				}
			}
		}
		else if (parent != NULL && adj_id != parent->getID())
			low[atom] = std::min(low[atom], discovery[adj_id]);
	}
}


