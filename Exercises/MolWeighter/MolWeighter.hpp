/*
 * Edna Teich, Claudius Sandmeier
 * Praesenzuebung
 * 04.04.2017
 */

/*
 * - Warum macht es Sinn, manche Funktionen als "private" und nicht als "public" in der MolWeighter-Klasse zu definieren? -
 *
 * Private Funktionen sind nicht fuer den Zugriff von ausserhalb der Klasse bestimmt und liefern
 * im schlimmsten Fall falsche Ergebnisse, wenn sie nicht an der korrekten Stelle innerhalb einer anderen Methode aufgerufen werden.
 *
 * Sie sind zudem nur Hilfsmethoden fuer die Methoden der Klasse, auf die von aussen zugegriffen werden kann ( = public).
 * Daher ist es auch leichter zu erkennen, welche Operationen eine Klasse ueber ein Interface anbietet,
 * wenn nur die von ausserhalb nutzbaren Methoden gelistet sind.
 *
 */

#pragma once

#include "Naomini/Forward.hpp"

namespace Naomini{

enum Operator {
  Equal = 0,
  MoreThan,
  LessThan
};

class MolWeighter{

public:

	MolWeighter(MoleculeVector mols);
	MoleculeVector getMoleculesWithWeight (Operator op, double weight_value);

private:

	MoleculeVector molecules; // Klassenvariable
	bool moleculeFitsCriteria(double mol_value, Operator op, double value);

};

}
