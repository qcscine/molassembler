#include <vector>
#include <algorithm>
#include <tuple>

using namespace std;

namespace MoleculeManip {

    enum BondType : uint8_t {
        Single,
        Double,
        Triple,
        Quadruple,
        Quintuple,
        Sextuple
    };

    class AdjacencyList {
    public:
    /* Typedefs */
        using unsigned_type = unsigned;
        using AtomIndexType = uint16_t; // 65k max is sufficient

    private:
    /* Private members */
        

    /* Private member functions */

    public:
    /* Public member functions */
        /*!
         * Checks whether two atoms are bonded.
         * \param[in] a The first atom index
         * \param[in] b The second atom index
         * \returns Whether two atoms are bonded or not.
         */
        virtual bool bond_exists(
            const AtomIndexType& a,
            const AtomIndexType& b
        ) const noexcept = 0;

        /*!
         * Checks whether two atoms are bonded by the specified BondType.
         * \param[in] a The first atom index
         * \param[in] b The second atom index
         * \param[in] bond_type The bond type to check
         * \returns Whether the two atoms are bonded as specified.
         */
        virtual bool bond_exists(
            const AtomIndexType& a,
            const AtomIndexType& b,
            const BondType& bond_type
        ) const noexcept = 0;

        /*!
         * Adds a bond to the AdjacencyList. No checks are made, and no
         * exceptions are thrown.
         * \param[in] a The first atom index
         * \param[in] b The second atom index
         * \param[in] bond_type The bond type to register
         */
        virtual void add_bond(
            const AtomIndexType& a,
            const AtomIndexType& b,
            const BondType& bond_type
        ) noexcept = 0;

        /*!
         * Adds an entry for an atom to the AdjacencyList. No checks are made
         * and no exceptions are thrown.
         */
        virtual AtomIndexType add_atom() noexcept = 0;

        /*!
         * Removes a bond from the AdjacencyList. If the bond doesn't exist,
         * the function call is ignored.
         * \param[in] a The first atom index
         * \param[in] b The second atom index
         */
        virtual void remove_bond(
            const AtomIndexType& a,
            const AtomIndexType& b
        ) noexcept = 0;

        /*!
         * Gets the BondType between two atoms.
         * \param[in] a The first atom index
         * \param[in] b The second atom index
         * \returns The bondtype between the atoms.
         */
        virtual BondType get_bond_type(
            const AtomIndexType& a,
            const AtomIndexType& b
        ) noexcept = 0;

        /*!
         * Get a vector of all bonds outwards from an atom.
         * \param[in] a The atom index
         * \returns A vector of all bonds outwards from an atom as std::pairs,
         * first is the BondType, second the atom index.
         */
        virtual vector<
            pair<
                BondType,
                AtomIndexType
            >
        > get_bonds(const AtomIndexType& a) const = 0;
    };

    template<typename T>
    class MinimalAdjacencyList : public AdjacencyList {
    private:
    /* Private members */
        T _adjacencies;

    public:
        bool bond_exists(
            const AtomIndexType& a,
            const AtomIndexType& b
        ) const noexcept override final;
    };

    template<typename T>
    bool MinimalAdjacencyList<T>::bond_exists(
        const AtomIndexType& a,
        const AtomIndexType& b
    ) const noexcept {
    }
}
