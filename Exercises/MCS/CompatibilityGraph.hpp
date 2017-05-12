/* Edna Teich, Claudius Sandmeier
 * CIW Blatt 2
 * Abgabe 04.05.2017
 */

#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <map>

#include "Naomini/Forward.hpp"
#include "Naomini/Atom.hpp"

namespace Naomini {

class Node {

public:
  Node(unsigned element, AtomPtr atom1, AtomPtr atom2)
    :element(element), atom1(atom1), atom2(atom2)
  {}

  unsigned element;
  AtomPtr atom1;
  AtomPtr atom2;

  bool operator== (const Node& n) const
  {
    return (element == n.element)
        && (atom1->getID() == n.atom1->getID())
        && (atom2->getID() == n.atom2->getID());
  }

  bool operator< (const Node& n) const
  {
    if (element != n.element) {
      return element < n.element;
    }
    else if (atom1->getID() != n.atom1->getID())
    {
      return atom1->getID() < n.atom1->getID();
    }
    else {
      return atom2->getID() < n.atom2->getID();
    }
  }

  friend std::ostream& operator<< (std::ostream& os, const Node& n)
  {
    os << n.element << "(" << n.atom1->getID() << "," << n.atom2->getID() << ")";
    return os;
  }
};

/* Edgeklasse analog zur Nodeklasse.
 * Ermoeglicht die gemeinsame Speicherung zweier Nodes, die durch eine Edge verbunden sind.
 */
class Edge
{
public:
	/* Konstruktor mit Instanzierung, analog zur Nodeklasse */
	Edge(Node node1, Node node2)
	:node1(node1), node2(node2)
	{}

	Node node1;
	Node node2;

	bool operator== (const Edge &edge) const
	{
		return (((node1 == edge.node1) && (node2 == edge.node2)) ||
				((node1 == edge.node2) && (node2 == edge.node1)));
	}
};

/** Graph Klasse zur Nutzung in BronKerbosch.hpp */

class CompatibilityGraph {

public:

	/* Konstruktor des Compatibility-Graphen zweier Molekuele.
	 * Erzeugt den Node-Vektor Vp und den Edgevektor Ep und befuellt sie.
	 */
	CompatibilityGraph(MoleculePtr mol1, MoleculePtr mol2);

	/* Destruktor */
	~CompatibilityGraph();

	/* Liefert die Anzahl der Nodes im Nodevektor */
	size_t getNofNodes();

	/* Liefert den Nodevektor zurueck. */
	std::vector<Node> getNodes() const;

	/* Liefert einen bestimmten Node im Nodevektor. */
	Node getNode(unsigned i) const;

	/* Ueberprueft, ob eine Edge zwischen zwei Nodes besteht. */
	bool hasEdge(const Node &a, const Node &b) const;

private:
/* Hier bitte Membervariablen zum Speichern von Knoten und Kanten einfuegen. */

  std::vector<Node> nodevector; // Vp
  std::vector<Edge> edgevector; // Ep
};

}
