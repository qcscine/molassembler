#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE AdjacencyListTests
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>
#include <memory>
#include <algorithm>
#include <iostream>
#include <random>

#include "AdjacencyList.hpp"

using namespace MoleculeManip;

using PairOfAtomIndices = std::pair<
    AtomIndexType,
    AtomIndexType
>;

std::vector<PairOfAtomIndices> all_valid_combinations(
    const AtomIndexType& n_atoms
) {
    using unsigned_t = AtomIndexType;
    std::vector<PairOfAtomIndices> return_vector;

    for(unsigned_t i = 0; i < n_atoms; i++) {
        for(unsigned_t j = i + 1; j < n_atoms; j++) {
            return_vector.emplace_back(i, j);
        }
    }

    return return_vector;
}

template<typename T>
ostream& operator << (
    ostream& os,
    const std::vector<T>& vector
) {
    for(const auto& element: vector) {
        os << element << std::endl;
    }
    return os;
}

ostream& operator << (
    ostream& os,
    const PairOfAtomIndices& pair
) {
    std::cout << "(" << pair.first << ", " << pair.second << ")";
    return os;
}

typedef boost::mpl::list<
    MinimalAdjacencyList,
    FastAdjacencyList
> test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE(insert_stability, T, test_types) {
    std::shared_ptr<AdjacencyList> test_object = std::make_shared<T>();

    AtomIndexType n_atoms = 50;

    // generate all valid combinations if insertable atoms
    std::vector<PairOfAtomIndices> combinations = all_valid_combinations(
        n_atoms
    );

    // make a random sequence
    std::shuffle(
        combinations.begin(),
        combinations.end(),
        std::mt19937{
            std::random_device{}()
        }
    );

    // add the atoms
    for(AtomIndexType i = 0; i < n_atoms; i++) {
        test_object -> add_atom();
    }
    
    // add the bonds
    for(const auto& combination : combinations) {
        test_object -> add_bond(
            combination.first,
            combination.second,
            BondType::Single
        );

        auto valid_and_explanation_pair = test_object -> validate();
        if(!valid_and_explanation_pair.first) {
            std::cout << "Failed on adding bond " << combination << std::endl;
            BOOST_FAIL(valid_and_explanation_pair.second);
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(get_bond_correctness, T, test_types) {

    std::shared_ptr<AdjacencyList> test_object = std::make_shared<T>();
    AtomIndexType n_atoms = 50;

    // generate all valid combinations of insertable atoms
    std::vector<PairOfAtomIndices> combinations = all_valid_combinations(
        n_atoms
    );

    // make a random sequence
    std::shuffle(
        combinations.begin(),
        combinations.end(),
        std::mt19937{
            std::random_device{}()
        }
    );

    // add the atoms
    for(AtomIndexType i = 0; i < n_atoms; i++) {
        test_object -> add_atom();
    }
    
    // add the bonds
    for(const auto& combination : combinations) {
        test_object -> add_bond(
            combination.first,
            combination.second,
            BondType::Single
        );
    }

    // re-shuffle the combinations
    std::shuffle(
        combinations.begin(),
        combinations.end(),
        std::mt19937{
            std::random_device{}()
        }
    );

    // test getter functions
    for(const auto& combination : combinations) {
        // bond_exists (1)
        if(
            !test_object -> bond_exists(
                combination.first,
                combination.second
            )
        ) {
            BOOST_FAIL(
                std::string("The added edge (")
                    + std::to_string(combination.first)
                    + ", "
                    + std::to_string(combination.second)
                    + ") is reported non-existent by bond_exists (1)."
            );
        }

        // bond_exists (2)
        if(
            !test_object -> bond_exists(
                combination.first,
                combination.second,
                BondType::Single
            )
        ) {
            BOOST_FAIL(
                std::string("The added edge (")
                    + std::to_string(combination.first)
                    + ", "
                    + std::to_string(combination.second)
                    + ", Single) is reported non-existent by bond_exists (2)."
            );
        }

        // get_bond_type
        if(test_object -> get_bond_type(
                combination.first,
                combination.second
            ) != BondType::Single
        ) {
            BOOST_FAIL(
                std::string("The added edge (")
                    + std::to_string(combination.first)
                    + ", "
                    + std::to_string(combination.second)
                    + ") is not BondType::Single, says get_bond_type"
            );
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(remove_bond_stability, T, test_types) {
    std::shared_ptr<AdjacencyList> test_object = std::make_shared<T>();
    AtomIndexType n_atoms = 50;

    // generate all valid combinations of insertable atoms
    std::vector<PairOfAtomIndices> combinations = all_valid_combinations(
        n_atoms
    );

    // make a random sequence
    std::shuffle(
        combinations.begin(),
        combinations.end(),
        std::mt19937{
            std::random_device{}()
        }
    );

    // add the atoms
    for(AtomIndexType i = 0; i < n_atoms; i++) {
        test_object -> add_atom();
    }
    
    // add the bonds
    for(const auto& combination : combinations) {
        test_object -> add_bond(
            combination.first,
            combination.second,
            BondType::Single
        );
    }

    // re-shuffle the combinations
    std::shuffle(
        combinations.begin(),
        combinations.end(),
        std::mt19937{
            std::random_device{}()
        }
    );

    // remove bonds, check validity
    for(const auto& combination : combinations) {
        test_object -> remove_bond(
            combination.first,
            combination.second
        );

        auto valid_and_explanation_pair = test_object -> validate();
        if(!valid_and_explanation_pair.first) {
            std::cout << "Failed on removing bond " << combination << std::endl;
            BOOST_FAIL(valid_and_explanation_pair.second);
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(remove_atom_stability, T, test_types) {
    std::shared_ptr<AdjacencyList> test_object = std::make_shared<T>();
    AtomIndexType n_atoms = 50;

    // generate all valid combinations of insertable atoms
    std::vector<PairOfAtomIndices> combinations = all_valid_combinations(
        n_atoms
    );

    // make a random sequence
    std::shuffle(
        combinations.begin(),
        combinations.end(),
        std::mt19937{
            std::random_device{}()
        }
    );

    // add the atoms
    for(AtomIndexType i = 0; i < n_atoms; i++) {
        test_object -> add_atom();
    }
    
    // add the bonds
    for(const auto& combination : combinations) {
        test_object -> add_bond(
            combination.first,
            combination.second,
            BondType::Single
        );
    }

    std::vector<AtomIndexType> atom_indices(n_atoms);
    std::iota(
        atom_indices.begin(),
        atom_indices.end(),
        0
    );
    std::shuffle(
        atom_indices.begin(),
        atom_indices.end(),
        std::mt19937{
            std::random_device{}()
        }
    );

    // remove atoms, check validity
    for(const auto& atom_index : atom_indices) {
        test_object -> remove_atom(atom_index);

        auto valid_and_explanation_pair = test_object -> validate();
        if(!valid_and_explanation_pair.first) {
            std::cout << "Failed on removing atom " << atom_index << std::endl;
            BOOST_FAIL(valid_and_explanation_pair.second);
        }
    }
}
