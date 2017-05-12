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

		Naomini::MoleculeDrawer drawer1(mol);
		//RingVector rings = moleculeGetRings(mol);
		drawer1.markSubstructures(rings);

		Naomini::MoleculeDrawer drawer2(mol);
		//RingVector extendedRings = moleculeGetExtendedRings(mol);
		RingVector extendedRings = moleculeGetExtendedRings(rings);
		drawer2.markSubstructures(extendedRings);

		Naomini::MoleculeDrawer drawer3(mol);

		moleculeDyeBonds(mol, drawer3, rings, bccs);

		/*
		AtomVector dyed_rings;
		AtomVector dyed_linkers;
		AtomVector dyed_misc;
		moleculeDyeComponents(mol, rings, bccs, dyed_rings, dyed_linkers, dyed_misc);
		drawer3.markSubstructure(dyed_rings, MoleculeDrawer::RED);
		drawer3.markSubstructure(dyed_linkers, MoleculeDrawer::GREEN);
		drawer3.markSubstructure(dyed_misc, MoleculeDrawer::YELLOW);
		*/

		//BCCVector bccs = moleculeGetBiconnectedComponents(mol);
		//drawer3.markSubstructures(bccs);

		//Naomini::MoleculeDrawer drawer4(mol);
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
