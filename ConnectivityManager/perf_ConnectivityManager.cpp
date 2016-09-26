#include <memory>
#include <algorithm>
#include <iostream>
#include <random>
#include <chrono>
#include <iomanip>
#include <cmath>
#include <type_traits>

#include "ConnectivityManager.hpp"
#include "numeric/numeric.hpp"

/* TODO
 */

using namespace MoleculeManip;

using PairOfAtomIndices = std::pair<
    AtomIndexType,
    AtomIndexType
>;

std::vector<PairOfAtomIndices> all_valid_combinations(
    const AtomIndexType& lower,
    const AtomIndexType& upper
) {
    using unsigned_t = AtomIndexType;
    std::vector<PairOfAtomIndices> return_vector;

    for(unsigned_t i = lower; i < upper; i++) {
        for(unsigned_t j = i + 1; j < upper; j++) {
            return_vector.emplace_back(i, j);
        }
    }

    return return_vector;
}

template<typename T>
ostream& operator << (
    ostream& os,
    const std::vector<T>& container
) {
    for(unsigned i = 0; i < container.size(); i++) {
        os << container[i];
        if(i != container.size() - 1) {
            os << ",";
        }
    }
    return os;
}


template<typename T>
unsigned time_wrap_nullary_callable(
    const T& nullary_callable
) {
    using namespace std::chrono;

    time_point<system_clock> start, end;
    start = system_clock::now();

    nullary_callable();
    
    end = system_clock::now();
    duration<double> elapsed = end - start;

    return elapsed.count() * 1e9;
}

template<typename T>
std::vector<unsigned> performance_test(
    T& test_instance,
    const AtomIndexType& n_atoms
) {
    test_instance -> reset();
    std::vector<unsigned> exec_times;

    // add the required number of atoms
    std::vector<AtomIndexType> atom_indices;
    std::vector<unsigned> current_exec_times;
    
    for(AtomIndexType i = 0; i < n_atoms; i++) {
        AtomIndexType temp; 
        current_exec_times.push_back(time_wrap_nullary_callable(
            [&]() {
                temp = test_instance -> add_atom();
            }
        ));
        atom_indices.push_back(
            temp
        );
    }

    exec_times.push_back(
        numeric::average(current_exec_times)
    );

    // generate the combinations from the lower and upper bound on atom indices
    std::vector<PairOfAtomIndices> combinations = all_valid_combinations(
        atom_indices[0],
        atom_indices[atom_indices.size() - 1]
    );

    // shuffle the atom indices and combinations
    std::shuffle(
        atom_indices.begin(),
        atom_indices.end(),
        std::mt19937{
            std::random_device{}()
        }
    );
    std::shuffle(
        combinations.begin(),
        combinations.end(),
        std::mt19937{
            std::random_device{}()
        }
    );

    // TODO
    // - currently, we get V² edges, this is not a typical molecule situation
    //   just truncating the combinations to a rough estimation of the number
    //   of edges in a real molecule (like e.g. 4 * n_atoms) does not
    //   necessarily lead to a connected graph (which we want)
    //   -> speaking of which, don't we need to ensure that the 
    //      ConnectivityManager needs to stay connected? It doesn't, currently.
    //      Perhaps do this by storing chords

    // oh well, for now, let's just truncate the damn thing
    if( 4 * n_atoms < combinations.size()) combinations.resize(4 * n_atoms);
    
    // add bonds
    current_exec_times.clear();
    for(const auto& combination : combinations) {
        current_exec_times.push_back(time_wrap_nullary_callable(
            [&]() {
                test_instance -> add_bond(
                    combination.first,
                    combination.second,
                    BondType::Single
                );
            }
        ));
    }
    exec_times.push_back(numeric::average(current_exec_times));

    // bond_exists (1)
    current_exec_times.clear();
    for(const auto& combination: combinations) {
        current_exec_times.push_back(time_wrap_nullary_callable(
            [&]() {
                test_instance -> bond_exists(
                    combination.first,
                    combination.second
                );
            }
        ));
    }
    exec_times.push_back(numeric::average(current_exec_times));

    // bond_exists (2)
    current_exec_times.clear();
    for(const auto& combination: combinations) {
        current_exec_times.push_back(time_wrap_nullary_callable(
            [&]() {
                test_instance -> bond_exists(
                    combination.first,
                    combination.second,
                    BondType::Single
                );
            }
        ));
    }
    exec_times.push_back(numeric::average(current_exec_times));

    // get_bond_type
    current_exec_times.clear();
    for(const auto& combination: combinations) {
        current_exec_times.push_back(time_wrap_nullary_callable(
            [&]() {
                test_instance -> bond_exists(
                    combination.first,
                    combination.second
                );
            }
        ));
    }
    exec_times.push_back(numeric::average(current_exec_times));

    // get_bond_pairs
    current_exec_times.clear();
    for(const auto& atom_index: atom_indices) {
        current_exec_times.push_back(time_wrap_nullary_callable(
            [&]() {
                test_instance -> get_bond_pairs(atom_index);
            }
        ));
    }
    exec_times.push_back(numeric::average(current_exec_times));

    // get_bonded_atom_indices
    current_exec_times.clear();
    for(const auto& atom_index: atom_indices) {
        current_exec_times.push_back(time_wrap_nullary_callable(
            [&]() {
                test_instance -> get_bonded_atom_indices(atom_index);
            }
        ));
    }
    exec_times.push_back(numeric::average(current_exec_times));

    // validate
    exec_times.push_back(time_wrap_nullary_callable(
        [&]() {
            test_instance -> validate();
        }
    ));

    // remove_bond
    current_exec_times.clear();
    for(const auto& combination: combinations) {
        current_exec_times.push_back(time_wrap_nullary_callable(
            [&]() {
                test_instance -> remove_bond(
                    combination.first,
                    combination.second
                );
            }
        ));
    }
    exec_times.push_back(numeric::average(current_exec_times));

    // re-add all bonds
    for(const auto& combination : combinations) {
        test_instance -> add_bond(
            combination.first,
            combination.second,
            BondType::Single
        );
    }

    // remove_all atoms
    current_exec_times.clear();
    for(const auto& atom_index: atom_indices) {
        current_exec_times.push_back(time_wrap_nullary_callable(
            [&]() {
                test_instance -> remove_atom(atom_index);
            }
        ));
    }
    exec_times.push_back(numeric::average(current_exec_times));

    return exec_times;
}

int main() {
    std::vector<
        std::pair<
            std::string,
            std::shared_ptr<
                AbstractConnectivityManager
            >
        >
    > test_pairs = {
        std::make_pair(
            "MinimalConnectivityManager",
            make_shared<MinimalConnectivityManager>()
        ),
        std::make_pair(
            "FastConnectivityManager",
            make_shared<FastConnectivityManager>()
        )
    };

    const std::vector<std::string> function_aliases = {
        "add_atom",
        "add_bond",
        "bond_exists_1",
        "bond_exists_2",
        "get_bond_type",
        "get_bond_pairs",
        "get_bonded_atom_indices",
        "validate",
        "remove_bond",
        "remove_atom"
    };

    const unsigned repeat = 100;
    const unsigned max_n_atoms = 1000;

    for(auto& test_pair : test_pairs) {
        std::cout << "For class " << test_pair.first << std::endl;
        std::vector<
            std::vector<unsigned>
        > function_averages (function_aliases.size(), std::vector<unsigned>());
        std::vector<
            std::vector<unsigned>
        > function_stddevs (function_aliases.size(), std::vector<unsigned>());

        for(unsigned n_atoms = 10; n_atoms <= max_n_atoms; n_atoms += 55) {
            std::vector<
                std::vector<unsigned>
            > microsecond_exec_times (function_aliases.size(), std::vector<unsigned>());

            for(unsigned i = 0; i < repeat; i++) {
                auto iteration = performance_test(
                    test_pair.second,
                    n_atoms
                );
                for(unsigned j = 0; j < iteration.size(); j++) {
                    microsecond_exec_times[j].push_back(
                        iteration[j]
                    );
                }
            }

            /* due to the odd way microsecond_exec_times is constructed,
             * these must be constructed differently
             */
            auto averages = numeric::map_average(microsecond_exec_times);
            auto stddevs = numeric::map_population_stddev(
                microsecond_exec_times
            );

            /* now I have 
             * (add_atom for n_atoms) <-> averages[0] +- stddevs[0]
             * ...
             *
             * want
             * function_averages[function_index] = [average_2, average_3, ...]
             * function_stddevs[function_index] = [stddev_2, stddev_3, ...]
             */

            assert(averages.size() == function_averages.size());
            for(unsigned i = 0; i < averages.size(); i++) {
                function_averages[i].push_back(averages[i]);
                function_stddevs[i].push_back(stddevs[i]);
            }
        }

        // write it all to stdout
        for(unsigned i = 0; i < function_averages.size(); i++) {
            std::cout << "av_" << function_aliases[i] << ", ";
            std::cout << function_averages[i] << std::endl;
            
            std::cout << "sd_" << function_aliases[i] << ", ";
            std::cout << function_stddevs[i] << std::endl;
        }

        std::cout << std::endl;
    }
    
    return 0;
}
