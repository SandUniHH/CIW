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

#ifndef NAOMINI_ULLMANN_HPP_
#define NAOMINI_ULLMANN_HPP_


/*----------------------------------------------------------------------------*/
/**********   INCLUDES                                               **********/
/*----------------------------------------------------------------------------*/

#include <boost/numeric/ublas/matrix.hpp>
#include "Naomini/Forward.hpp"

/*----------------------------------------------------------------------------*/
/**********   NAMESPACE                                              **********/
/*----------------------------------------------------------------------------*/

namespace Naomini {

typedef boost::numeric::ublas::matrix<bool> BoolMatrix;
typedef std::vector<bool> BoolVector;

/*----------------------------------------------------------------------------*/
/**********   CLASS                                                  **********/
/*----------------------------------------------------------------------------*/

/** Claculates all substructure matches of a query substructure in a molecule.
 *
 * @brief class for substructure matching
 */

class Ullmann {

public:
  /*---------------*/
   /** CONSTRUCTOR **/
   /*---------------*/

   /** Constructs an object of class Ullmann.
   *
   * @brief constructor  */
  Ullmann(MoleculePtr substructure, MoleculePtr molecule);


  /*------------*/
  /** METHODS  **/
  /*------------*/

  /** Triggers the substructure calculation.
   * @brief run ullman algorithm */
  void run();

private:

  /** Checks if a permutation matrix meets the subgraph isomorphism criterion.
    * @brief checks subgraph isomorphism criterion*/
  bool checkSubgraph(const BoolMatrix &M);

  /** Draws the molecule and marks the substructure matching which is encoded
   * in the matrix.
   * @brief draws substructure matching*/
  void drawerSubstructureMatching(const BoolMatrix &M);

  /** Calculates the initial matrix on basis of element type and number of heavy
   * atom neighbors.
   * @brief get initial matrix*/
  void enumerateSubgraph(BoolMatrix &M, unsigned k, BoolVector F);

  /** Calculates the initial matrix on basis of element type and number of heavy
   * atom neighbors.
   * @brief get initial matrix*/
  BoolMatrix getInitialMatrix();

  /** Prints a matrix to standard output.
   * @brief clear atom selection*/
  void printMatrix(BoolMatrix &M, unsigned k, const std::string& debugInfo);

  /* Speed up the algorithm by adjusting the matrix from line k on downward.
   *
   * time without refinement, Intel(R) Core(TM) i5-4590 CPU @ 3.30GHz, 1 core:
   *
   * user    21m22.971s
   * sys     0m0.033s
   *
   * user    21m28.062s
   * sys     0m0.033s
   *
   * with refinement:
   *
   * user    0m0.999s
   * sys     0m0.034s

   *
   * user    0m1.016s
   * sys     0m0.027s
   *
   * @brief refine the current matrix */
  bool refine(BoolMatrix &M, unsigned k);

  /*--------------------*/
  /** MEMBER VARIABLES **/
  /*--------------------*/

  MoleculePtr m_Molecule; ///< target molecule
  AtomVector m_MoleculeAtoms; ///< target atoms
  BoolMatrix m_MoleculeAdjacencyMatrix; ///< target adjacency matrix

  AtomVector m_SubstructureAtoms; ///< substructure query atom
  BoolMatrix m_SubstructureAdjacencyMatrix; ///< substructure query adjacency matrix
};

/*----------------------------------------------------------------------------*/
/**********   NAMESPACE END                                          **********/
/*----------------------------------------------------------------------------*/
}
#endif /* NAOMINI_ULLMANN_HPP_ */
