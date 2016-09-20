#include <list>
#include <vector>
#include <memory>
#include "AtomSet.h"
#include "Molecule.hpp"

namespace MoleculeManip {

    class MoleculeCollection {
        /* Typedefs */
        using unsigned_type = unsigned;

    private:
        std::list<
            std::shared_ptr<
                Molecule
            >
        > _molecule_instances;

    public:
        /* Additional constructors:
         *
         * Detect a molecular graph from an AtomSet's spatial coordinates 
         * \param[in] atomset The Delib::AtomSet to initialize from
         * MoleculeCollection(const Delib::AtomSet& atomset);
         */

        /* Destructor
         */
        virtual ~I_MoleculeCollection() {};

        /* Public member functions */
        /*!
         * Get the number of atoms in the MoleculeCollection
         * \returns The number of atoms in the MoleculeCollection
         */
        virtual unsigned_type get_num_atoms() const = 0;

        /*!
         * Get the number of molecules in the MoleculeCollection
         * \returns The number of molecules in the MoleculeCollection
         */
        virtual unsigned_type get_num_molecules() const = 0;

        /*!
         * Get the Molecules in the MoleculeCollection
         * \returns A list of groups in the MoleculeCollection
         */
        virtual std::list<
            std::shared_ptr<
                Molecule
            >
        > get_molecules() const = 0;

        /*!
         * Add a Molecule to the MoleculeCollection
         * \param molecule The molecule to add
         */
        virtual void add_molecule(const Molecule& molecule) = 0;

        /*!
         * Get a Delib::AtomSet 
         * \returns an AtomSet constructed from all Molecules
         */
        virtual Delib::AtomSet get_atom_set() const = 0;

    };

}
