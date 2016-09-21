#include <vector>
#include <algorithm>
#include <tuple>
#include <cmath>
#include <functional>
#include <iostream>

using namespace std;

template<typename T, typename U>
ostream& operator << (
    ostream& os,
    const std::pair<T, U>& pair
) {
    os << "(" << pair.first << ", " << pair.second << ")";
    return os;
}

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
        using TupleType = tuple<
            AtomIndexType,
            AtomIndexType,
            BondType
        >;

    private:
    /* Private members */
        

    /* Private member functions */

    public:
        virtual ~AdjacencyList() {};

    /* Public member functions */
        /*!
         * Checks whether the list is ordered.
         * \returns whether the list is ordered.
         */
        virtual bool is_ordered() const noexcept = 0;

        /*!
         * Checks whether two atoms are bonded.
         * \param[in] a The first atom index
         * \param[in] b The second atom index
         * \returns Whether two atoms are bonded or not.
         */
        /*virtual bool bond_exists(
            const AtomIndexType& a,
            const AtomIndexType& b
        ) const noexcept = 0;*/

        /*!
         * Checks whether two atoms are bonded by the specified BondType.
         * \param[in] a The first atom index
         * \param[in] b The second atom index
         * \param[in] bond_type The bond type to check
         * \returns Whether the two atoms are bonded as specified.
         */
        /*virtual bool bond_exists(
            const AtomIndexType& a,
            const AtomIndexType& b,
            const BondType& bond_type
        ) const noexcept = 0;*/

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
        /*virtual AtomIndexType add_atom() noexcept = 0;*/

        /*!
         * Removes a bond from the AdjacencyList. If the bond doesn't exist,
         * the function call is ignored.
         * \param[in] a The first atom index
         * \param[in] b The second atom index
         */
        /*virtual void remove_bond(
            const AtomIndexType& a,
            const AtomIndexType& b
        ) noexcept = 0;*/

        /*!
         * Gets the BondType between two atoms.
         * \param[in] a The first atom index
         * \param[in] b The second atom index
         * \returns The bondtype between the atoms.
         */
        /*virtual BondType get_bond_type(
            const AtomIndexType& a,
            const AtomIndexType& b
        ) noexcept = 0;*/

        /*!
         * Get a vector of all bonds outwards from an atom.
         * \param[in] a The atom index
         * \returns A vector of all bonds outwards from an atom as std::pairs,
         * first is the BondType, second the atom index.
         */
        /*virtual vector<
            pair<
                BondType,
                AtomIndexType
            >
        > get_bonds(const AtomIndexType& a) const = 0;*/
    };

    class MinimalAdjacencyList : public AdjacencyList {
    public:
        using EdgeIndexType = uint32_t;
    private:
    /* Private members */
        std::vector<TupleType> _adjacencies;

        /* This function assumes a < b! */
        EdgeIndexType _get_insert_pos(
            const AtomIndexType& a, 
            const AtomIndexType& b
        ) const noexcept {
            if(_adjacencies.size() == 0) {
                return 0;
            }

            auto tuple_lt = [&a, &b](const TupleType& tuple) -> bool {
                return (
                    std::get<0>(tuple) < a
                    || (
                        std::get<0>(tuple) == a
                        && std::get<1>(tuple) < b 
                    )
                );
            };
            auto tuple_gt = [&a, &b](const TupleType& tuple) -> bool {
                return (
                    std::get<0>(tuple) > a 
                    || (
                        std::get<0>(tuple) == a 
                        && std::get<1>(tuple) > b 
                    )
                );
            };

            auto print_tuple = [](const TupleType& tuple) -> void {
                std::cout << "(" << std::get<0>(tuple) << ", ";
                std::cout << std::get<1>(tuple) << ", ";
                std::cout << std::get<2>(tuple) << ")" << std::endl;
            };
            

            std::cout << "Table is:" << std::endl;
            for(const auto& tuple : _adjacencies) {
                print_tuple(tuple);
            }

            std::cout << "Inserting (" << a << ", " << b << ")" << std::endl;

            unsigned_type left = 0;
            unsigned_type right = _adjacencies.size() - 1;
            unsigned_type middle;

            unsigned_type debug_counter = 0;

            auto print_state = [&](const unsigned_type& counter) -> void {
                std::cout << counter << ": " << left << " " << middle << " " << right << std::endl;
            };

            while(left <= right) {
                debug_counter += 1;
                if(debug_counter == 100) {
                    std::cout << "Broken out due to max cycles!" << std::endl;
                    break;
                }

                middle = std::floor(
                    ((double) left + (double) right) / 2.0
                );
                print_state(debug_counter);

                if(left == right) {
                    std::cout << "-> L == R, middle is ";
                    print_tuple(
                        _adjacencies[middle]
                    );
                    std::cout << std::endl;

                    if(tuple_lt(_adjacencies[middle])) {
                        return middle + 1;
                    } else {
                        return middle;
                    }
                }

                if(tuple_lt(_adjacencies[middle])) {
                    left = middle + 1;
                    print_state(debug_counter);
                    continue;
                }

                if(tuple_gt(_adjacencies[middle])) {
                    if(middle == 0) right = middle;
                    else right = middle - 1; 
                    print_state(debug_counter);
                    continue;
                }
            }


            std::cout << "This is a case without a stage where L==R" << std::endl;
            return 0;
        }

    public:
        bool is_ordered() const noexcept override final {
            auto lt = [](
                const TupleType& a,
                const TupleType& b
            ) {
                return (
                    std::get<0>(a) < std::get<0>(b)
                    || (
                        std::get<0>(a) == std::get<0>(b)
                        && std::get<1>(a) < std::get<1>(b)
                    )
                );
            };

            /* TODO
             * - Can't think of a way to put this in functional paradigm
             * - Not allowed to use indices if you want to be able to use list 
             *   in template instantiation
             */
            for(EdgeIndexType i = 0; i < _adjacencies.size() - 1; i++) {
                if(!lt(
                    _adjacencies[i],
                    _adjacencies[i+1]
                )) {
                    return false;
                }
            }
            return true;
        }

        /* No performance sacrifices are made for stability.
         * a, b must fulfill: a < b.
         */
        void add_bond(
            const AtomIndexType& a,
            const AtomIndexType& b,
            const BondType& bond_type
        ) noexcept override final {
            auto pos = _get_insert_pos(a, b);

            std::cout << "Inserting before position " << pos << std::endl;

            _adjacencies.insert(
                _adjacencies.begin() + pos,
                make_tuple(
                    a,
                    b,
                    bond_type
                )
            );
        }
    };
}
