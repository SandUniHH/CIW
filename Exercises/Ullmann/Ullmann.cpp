/**
 @file
 @brief     
 @author    Bietz
 @date      Apr 10, 2012
 @if DoNotInclude
 Copyright ZBH  Center for Bioinformatics
 University of Hamburg
 Bundesstrasse 43, 20146 Hamburg, Germany
 ===============================================================================
 This module is part of the molecule software library Naomi,
 developed at the Center for Bioinformatics Hamburg (ZBH).
 It is not allowed to use, modify, copy or redistribute this software without
 written approval of the ZBH. For further information contact

 ZBH  Center for Bioinformatics, University of Hamburg
 Research Group for Computational Molecular Design

 Voice:  +49 (40) 42838 7350
 Fax:    +49 (40) 42838 7352
 E-Mail: info@zbh.uni-hamburg.de
 ===============================================================================
 @endif
 */
#include <iomanip>
#include <stdio.h>

#include "Ullmann.hpp"
#include "Naomini/Forward.hpp"
#include "Naomini/Atom.hpp"
#include "Naomini/Bond.hpp"
#include "Naomini/Molecule.hpp"
#include "Naomini/Helpers.hpp"
#include "Naomini/MoleculeDrawer.hpp"

/*----------------------------------------------------------------------------*/
/**********   NAMESPACE NAOMINI                                      **********/
/*----------------------------------------------------------------------------*/
namespace Naomini {

/*----------------------------------------------------------------------------*/
/**********   NAMELESS NAMESPACE                                     **********/
/*----------------------------------------------------------------------------*/
namespace{

/*------------*/
/** HELPERS  **/
/*------------*/

bool atoms_are_adjacent(AtomPtr a, AtomPtr b){
  for (AtomPtr n : a->getNeighborAtoms()) {
    if (b == n){
      return true;
    }
  }
  return false;
}

/*----------------------------------------------------------------------------*/

BoolMatrix get_adjacency_matrix(const AtomVector& atoms){
  unsigned nofAtoms = atoms.size();
  BoolMatrix adj(nofAtoms, nofAtoms);
  for (unsigned i = 0; i < nofAtoms - 1; ++i) {
    adj(i,i) = false;
    for (unsigned j = i+1 ; j < nofAtoms; ++j) {
      adj(i,j) = adj(j,i) = atoms_are_adjacent(atoms.at(i), atoms.at(j));
    }
  }
  return adj;
}

/*----------------------------------------------------------------------------*/

unsigned nofHeavyNeighbors(AtomPtr atom){
  unsigned count = 0;
  for (AtomPtr n : atom->getNeighborAtoms()) {
    if (!atomIsHydrogen(n)) count++;
  }
  return count;
}

/*----------------------------------------------------------------------------*/

AtomVector getNeighborHeavyAtoms(AtomPtr atom)
{
  AtomVector atom_vec = atom->getNeighborAtoms();
  AtomVector output_atom_vec;
  unsigned nofAtoms = atom_vec.size();
  for(unsigned i = 0; i < nofAtoms; ++i)
  {
    if(!atomIsHydrogen(atom_vec.at(i)) )
    {
      output_atom_vec.push_back(atom_vec.at(i));
    }
  }
  return output_atom_vec;
}

/*----------------------------------------------------------------------------*/
/**********   END NAMELESS NAMESPACE                                 **********/
/*----------------------------------------------------------------------------*/
}


/*---------------*/
/** CONSTRUCTOR **/
/*---------------*/

Ullmann::Ullmann(MoleculePtr substructure, MoleculePtr molecule)
:m_Molecule(molecule){
  for (AtomPtr sAtom : substructure->getAtoms()) {
    if (!atomIsHydrogen(sAtom)){
      m_SubstructureAtoms.push_back(sAtom);
    }
  }
  for (AtomPtr mAtom : molecule->getAtoms()) {
    if (!atomIsHydrogen(mAtom)){
      m_MoleculeAtoms.push_back(mAtom);
    }
  }

  m_MoleculeAdjacencyMatrix = get_adjacency_matrix(m_MoleculeAtoms);
  m_SubstructureAdjacencyMatrix = get_adjacency_matrix(m_SubstructureAtoms);
}

/*------------*/
/** METHODS  **/
/*------------*/

void Ullmann::run(){
  BoolMatrix initMatrix = getInitialMatrix();
  printMatrix(initMatrix, 0, "Initialmatrix");

  BoolVector F(m_MoleculeAtoms.size(), false);

  if (refine(initMatrix, 0)) // optional implementation according to the Ullmann paper, but time difference is negligible
	  enumerateSubgraph(initMatrix, 0, F);
}

/*----------------------------------------------------------------------------*/

bool Ullmann::checkSubgraph(const BoolMatrix &M){
  using namespace boost::numeric::ublas;

  BoolMatrix matchingAdj = prod(M, trans(prod(M, m_MoleculeAdjacencyMatrix)));

  for (unsigned i = 0; i < m_SubstructureAdjacencyMatrix.size1()-1; ++i ) {
    for (unsigned j = i+1; j < m_SubstructureAdjacencyMatrix.size2(); ++j ) {
      if (m_SubstructureAdjacencyMatrix(i,j) && !matchingAdj(i,j)){
        return false;
      }
    }
  }
  return true;
}

/*----------------------------------------------------------------------------*/

void Ullmann::drawerSubstructureMatching(const BoolMatrix &M){
  AtomVector substructure;

  for (unsigned i = 0; i < M.size1(); ++i ) {
    for (unsigned j = 0; j < M.size2(); ++j ) {
      if (M(i,j)) {
        substructure.push_back(m_MoleculeAtoms.at(j));
        break;
      }
    }
  }

  MoleculeDrawer drawer(m_Molecule);
  drawer.markSubstructure(substructure, MoleculeDrawer::GREEN);

}

/*----------------------------------------------------------------------------*/

void Ullmann::enumerateSubgraph(BoolMatrix &M, unsigned k, BoolVector F){

  if (k == m_SubstructureAtoms.size()){
    // if (checkSubgraph(M)){ // not needed due to refinement
      drawerSubstructureMatching(M);
    // }
  }
  else{
    for (unsigned l = 0; l < M.size2(); ++l) {
      if (M(k,l) && !F.at(l)){
        BoolMatrix copyOfM = M;
        for (unsigned j = 0; j < M.size2(); ++j) {
          M(k,j) = false;
        }
        M(k,l) = true;
        F.at(l) =  true;
//        printMatrix(M, k+1, "");
        if(refine(M, k)) // refinement procedure
        	enumerateSubgraph(M, k+1, F);

        F.at(l) =  false;
        M = copyOfM;
      }
    }
  }
}

/*----------------------------------------------------------------------------*/

BoolMatrix Ullmann::getInitialMatrix(){
  unsigned nofSubstrAtoms = m_SubstructureAtoms.size();
  unsigned nofMolAtoms = m_MoleculeAtoms.size();

  BoolMatrix initMatrix(nofSubstrAtoms,nofMolAtoms);
  for (unsigned i = 0; i < nofSubstrAtoms; ++i) {
    AtomPtr subAtom =  m_SubstructureAtoms.at(i);
    for (unsigned j = 0; j < nofMolAtoms; ++j) {
      AtomPtr molAtom = m_MoleculeAtoms.at(j);
      if ((subAtom->getAtomicNumber() == molAtom->getAtomicNumber()) &&
          nofHeavyNeighbors(subAtom) <= nofHeavyNeighbors(molAtom)){
        initMatrix(i, j) = true;
      }
      else{
        initMatrix(i, j) = false;
      }
    }
  }
  return initMatrix;
}

/*----------------------------------------------------------------------------*/

void Ullmann::printMatrix(
    BoolMatrix &M, unsigned k, const std::string& name){
  std::cout << name << std::endl;
  std::cout << "  "<< k <<" ║";
  for (unsigned j = 0; j < M.size2(); j++) {
    std::cout << std::setw( 4 ) << std::setprecision(2)
    << m_MoleculeAtoms.at(j)->getName() << " │";
  }
  std::cout << "\n════╬═════";
  for (unsigned j = 1; j < M.size2(); j++) {
    std::cout << "╪═════";
  }
  std::cout << "╡\n";
  for (unsigned i = 0; i < M.size1(); i++) {
    std::cout << " "  << std::setw( 2 )
    << m_SubstructureAtoms.at(i)->getName() <<" ║";
    for (unsigned j = 0; j < M.size2(); j++) {
      if (M(i,j)){printf("  1  ");}
      else       {printf("     ");}
      std::cout << "│";
    }
    std::cout << "\n────╫─────";
    for (unsigned j = 1; j < M.size2(); j++) {
      std::cout << "┼─────";
    }
    std::cout << "┤\n";
  }
  std::cout << std::endl;
}

/*----------------------------------------------------------------------------*/


bool Ullmann::refine(BoolMatrix &M, unsigned k)
{
	bool valid_index, found_match, possible_substructure, matrix_changed;

	do // as long as there are changes done to the matrix
	{
		matrix_changed = false; // becomes true if 1s have been nullified
		for (unsigned i = k + 1; i < M.size1(); i++) // size1: number of rows
		{
			for (unsigned j = 0; j < M.size2(); j++) // size2: number of columns
			{
				if(M(i,j)) // only the indices with value 1 are relevant
				{
					valid_index = true;
					/* AtomVector a_molecule_neighbours = m_MoleculeAtoms.at(j)->getNeighborAtoms();
					AtomVector a_subgraph_neighbours = m_SubstructureAtoms.at(i)->getNeighborAtoms();
					*/
					AtomVector a_molecule_neighbours = getNeighborHeavyAtoms(m_MoleculeAtoms.at(j));
					AtomVector a_subgraph_neighbours = getNeighborHeavyAtoms(m_SubstructureAtoms.at(i));

					for (unsigned x = 0; x < a_subgraph_neighbours.size(); x++)
					{
						found_match = false;

						for (unsigned y = 0; y < a_molecule_neighbours.size(); y++)
						{
							if(M(a_subgraph_neighbours.at(x)->getID(),
									a_molecule_neighbours.at(y)->getID()))
							{
								found_match = true;
								break;
							}
						}
						if (!found_match)
						{
							valid_index = false;
							break;
						}
					}
					if (!valid_index)
					{
						M(i,j) = 0;
						matrix_changed = true;
						possible_substructure = false;
						for (unsigned h = 0; h < M.size2(); h++)
						{
							if(M(i,h))
							{
								possible_substructure = true;
								break;
							}
						}
						if(!possible_substructure) return false; // if there is not a single 1 left in the line
					}
				}
			}
		}
	}
	while (matrix_changed);

	return true;
}

/*----------------------------------------------------------------------------*/
// REFINEMENT

// parameter:
// M ist the bool matrix of the current permutation
// k is the last row permutation is set
//
//bool Ullmann::refine(BoolMatrix &M, unsigned k)
//{
//  bool change, valid, found, allocation;
//  do
//  {
//    change = false;
//    // iterate through all rows with i>k
//    for(unsigned i = k+1; i < M.size1(); i++)
//    {
//      // iterate through all columns
//      for(unsigned j = 0; j < M.size2(); j++)
//      {
//	// search for all matches in that row
//	if( M(i, j) )
//	{
//	  valid = true;
//	  // get all neigbor atoms of molecule and substructure of the current
//	  // allocaed pair
//	  AtomVector sub_neighbor_atoms =
//	  getNeighborHeavyAtoms(m_SubstructureAtoms.at(i));
//	  AtomVector mol_neighbor_atoms =
//	  getNeighborHeavyAtoms(m_MoleculeAtoms.at(j));
//	  //iterate through all neighbor atoms of the substructure
//	  for(unsigned l = 0; l < sub_neighbor_atoms.size(); l++)
//	  {
//	    found = false;
//	    // iterate through all neighbor atoms of the molecule
//	    for(unsigned m = 0; m < mol_neighbor_atoms.size(); m++)
//	    {
//	      // check if each neighbor of the current atom pair can be
//	      // allocated to a neighbor of the partner atom
//	      if(M(sub_neighbor_atoms.at(l)->getID(),
//		mol_neighbor_atoms.at(m)->getID()))
//	      {
//		found = true;
//		break;
//	      }
//	    }
//	    // stop iteration if a neighbor atom in the substructure could not
//	    // be allocated to any neighbor of the partner atom in the molecule
//	    if(!found)
//	    {
//	      valid = false;
//	      break;
//	    }
//	  }
//	  if(!valid)
//	  {
//	    // erase allcoation in the bool matrix
//	    // and mark permutation as changed
//	    M(i,j) = 0;
//	    change = true;
//	    allocation = false;
//	    // check if there is any valid allocation of the current atom of the
//	    // substructure
//	    for(unsigned h = 0; h < M.size2(); h++)
//	    {
//	      if(M(i,h))
//	      {
//		allocation = true;
//	      }
//	    }
//	    // if no allocation was found, return false
//	    if(!allocation)
//	    {
//	      return false;
//	    }
//	  }
//	}
//      }
//    }
//  // go on while Matrix was changed
//  }while(change);
//  return true;
//}


/*----------------------------------------------------------------------------*/
/**********   NAMESPACE END                                          **********/
/*----------------------------------------------------------------------------*/
}
