#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE AdjacencyListTests
#include <boost/test/unit_test.hpp>
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
    const AtomIndexType& lower_bound,
    const AtomIndexType& upper_bound
) {
    using unsigned_t = AtomIndexType;
    std::vector<PairOfAtomIndices> return_vector;

    for(unsigned_t i = lower_bound; i <= upper_bound; i++) {
        for(unsigned_t j = i + 1; j <= upper_bound; j++) {
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

/* MinimalAdjacencyList tests */

BOOST_AUTO_TEST_CASE( minimal_operation_stability ) {
    // instantiate an object
    std::shared_ptr<AdjacencyList> test_object = std::make_shared<
        MinimalAdjacencyList
    >();

    AtomIndexType lower = 1;
    AtomIndexType upper = 50;

    // generate all valid combinations if insertable atoms
    std::vector<PairOfAtomIndices> combinations = all_valid_combinations(
        lower,
        upper
    );

    // make a random sequence
    std::shuffle(
        combinations.begin(),
        combinations.end(),
        std::mt19937{
            std::random_device{}()
        }
    );

    // add the bonds
    for(const auto& combination : combinations) {
        test_object->add_bond(
            combination.first,
            combination.second,
            BondType::Single
        );

        if(!test_object->is_ordered()) {
            std::cout << "Failed on " << combination << std::endl;
            BOOST_FAIL(
                "The adjacency list is not ordered anymore!"
            );
        }
    }
}
