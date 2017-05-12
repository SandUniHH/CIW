#include <set>

#include "BronKerbosch.hpp"
#include "Naomini/Atom.hpp"
#include "Naomini/Bond.hpp"
#include "Naomini/Molecule.hpp"
#include "Naomini/Helpers.hpp"

#include "CompatibilityGraph.hpp"

namespace Naomini {

namespace {

std::set<Node> get_adj(const CompatibilityGraph &comp, Node n)
{
  std::vector<Node> nodes = comp.getNodes();
  std::set<Node> adj;

  for (unsigned i = 0; i < nodes.size(); i++)
  {
    Node &node = nodes.at(i);
    if (comp.hasEdge(node, n))
    {
      adj.insert(node);
    }
  }
  return adj;
}

std::set<Node> add_to_set(std::set<Node> set, Node x)
{
  std::set<Node> new_set (set);
  new_set.insert(x);
  return new_set;
}

std::set<Node> cut_set(const CompatibilityGraph &comp, std::set<Node> set, Node x)
{
  std::set<Node> new_set;
  std::set<Node> X_adj = get_adj(comp,x);

  std::set<Node>::iterator iter;
  for (iter = set.begin(); iter != set.end(); iter++)
  {
    const Node &node = *iter;
    if (X_adj.find(node) != X_adj.end())
    {
      new_set.insert(node);
    }
  }
  return new_set;
}

void output_clique(std::set<Node> &set)
{
  std::set<Node>::iterator it;
  std::cout << "clique:";

  for (it = set.begin(); it != set.end(); it++)
  {
    std::cout << " " << *it;
  }
  std::cout << std::endl;
}

}

/*----------------------------------------------------------------------------*/

void clique_bk(const CompatibilityGraph &comp, std::set<Node> &set_i, std::set<Node> &set_x, std::set<Node> &set_n)
{
  if (set_x.empty()) {
    if (set_n.empty()) {
      output_clique(set_i);
    } else {
//      std::cerr << "  Break: X is empty but N is not\n";
      return;
    }
  }
  else
  {
    std::set<Node>::iterator it;
    for (it = set_n.begin(); it != set_n.end(); it++)
    {
      Node n = *it;
      std::set<Node> N_adj = get_adj(comp,n);
      bool breakup = true;

      std::set<Node>::iterator jt;
      for (jt = set_x.begin(); jt != set_x.end(); jt++)
      {
        Node x = *jt;
        std::set<Node>::iterator find = N_adj.find(x);

        if (find == N_adj.end())
        {
          breakup = false;
        }
      }
      if (breakup) {
//        std::cout << "  Break: " << n << " is adjacent to all elements of X\n";
        return;
      }
    }

    std::set<Node>::iterator nodeIter;
    for (nodeIter = set_x.begin(); nodeIter != set_x.end(); nodeIter++)
    {
      Node x (*nodeIter);

      std::set<Node> new_i (add_to_set(set_i, x));
      std::set<Node> new_x (cut_set(comp, set_x, x));
      std::set<Node> new_n (cut_set(comp, set_n, x));

      clique_bk(comp, new_i, new_x, new_n);
      set_n.insert(x);
      set_x.erase(set_x.find(x));
    }
  }
}


void calculateBronKerbosch(Naomini::MoleculePtr mol1, MoleculePtr mol2)
{
  CompatibilityGraph comp(mol1,mol2);

  std::set<Node> set_i;
  std::set<Node> set_x;
  for (unsigned i = 0 ; i < comp.getNofNodes(); i++)
  {
    set_x.insert(comp.getNode(i));
  }
  std::set<Node> set_n;

  clique_bk(comp, set_i, set_x, set_n);
}

}








