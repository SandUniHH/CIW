/* Edna Teich, Claudius Sandmeier
 * CIW Blatt 3
 * Abgabe 18.05.2017
 */

#include <iostream>

#include "Naomini/Atom.hpp"
#include "Naomini/Bond.hpp"
#include "Naomini/Molecule.hpp"
#include "Naomini/MoleculeFactory.hpp"
#include "Naomini/MoleculeDrawer.hpp"
#include "Naomini/Helpers.hpp"

#include "RingsAndBiconnectedComponents.hpp"
namespace Naomini{

int main(int argc, char *argv[]) {
  // check the number of arguments
  if (argc != 2){
    std::cout << "Usage: " << argv[0] << " <mol2-file>" << std::endl;
    return 1;
  }
  // get a molecule for each molecule entry in the provided file
  MoleculeVector mols = MoleculeFactory::getAllMolecules(std::string(argv[1]));

  try{
    // calculate and draw all rings of all molecules
    for (MoleculePtr mol : mols) {

		RingVector rings = moleculeGetRings(mol);
		BCCVector bccs = moleculeGetBiconnectedComponents(mol);

		/* Mark all rings */
		Naomini::MoleculeDrawer drawer1(mol);
		drawer1.markSubstructures(rings);

		/* Mark all rings and single terminal atoms connected to rings */
		Naomini::MoleculeDrawer drawer2(mol);
		//RingVector extendedRings = moleculeGetExtendedRings(mol);
		RingVector extendedRings = moleculeGetExtendedRings(rings);
		drawer2.markSubstructures(extendedRings);

		Naomini::MoleculeDrawer drawer3(mol);

		/* dye all biconnected components (linkers) */
		//drawer3.markSubstructures(bccs);

		/* dye all bonds according to whether they are rings, linkers or other. */
		moleculeDyeBonds(mol, drawer3, rings, bccs);

		/* dye all atoms instead */
		/*
		AtomVector dyed_rings;
		AtomVector dyed_linkers;
		AtomVector dyed_misc;
		moleculeDyeComponents(mol, rings, bccs, dyed_rings, dyed_linkers, dyed_misc);
		drawer3.markSubstructure(dyed_rings, MoleculeDrawer::RED);
		drawer3.markSubstructure(dyed_linkers, MoleculeDrawer::GREEN);
		drawer3.markSubstructure(dyed_misc, MoleculeDrawer::YELLOW);
		*/
    }
  }
  catch (const char *err){
    std::cerr << err << std::endl;
  }

  // delete all molecules
  for (MoleculePtr molecule : mols) {
    delete molecule;
  }
  return 0;
}

} // end namespace Naomini
