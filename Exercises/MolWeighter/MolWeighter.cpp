/*
 * Edna Teich, Claudius Sandmeier
 * Praesenzuebung
 * 04.04.2017
 */

#include "MolWeighter.hpp"

#include "Naomini/Helpers.hpp"
#include "Naomini/Molecule.hpp"
#include "Naomini/Atom.hpp"
#include "Naomini/Bond.hpp"

using namespace Naomini;

MolWeighter::MolWeighter(MoleculeVector mols)
{
	molecules = mols;
}

MoleculeVector MolWeighter::getMoleculesWithWeight(Operator op, double weight_value)
{
	MoleculeVector filtered_molecules;

	for (MoleculePtr molecule : molecules) // iterate over each molecule
	{
		double mol_weight = 0;

		for (auto iter = molecule->getAtoms().begin(); iter != molecule->getAtoms().end(); iter++)
		{
			const AtomPtr &aptr = *iter;
			mol_weight += aptr->getAtomicWeight(); // add each atomic weight to total molecule weight
		}

		/* push molecule to vector if it fits criteria */
		if (moleculeFitsCriteria(mol_weight, op, weight_value))
			filtered_molecules.push_back(molecule);
	}

	return filtered_molecules;
}


// private

bool MolWeighter::moleculeFitsCriteria(double mol_weight, Operator op, double weight_value)
{
	return ((op == Equal && mol_weight == weight_value)   ||
			(op == MoreThan && mol_weight > weight_value) ||
			(op == LessThan && mol_weight < weight_value));
}
