#include <memory>
#include <algorithm>
#include <iostream>
#include <random>
#include <chrono>
#include <iomanip>
#include <cmath>
#include <type_traits>

#include "ConnectivityManager.hpp"

/* TODO
 * - fails! not enough testing done?
 */

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

unsigned average(
    const std::vector<unsigned>& values
) {
    if(values.size() == 0) return 0;
    else return std::floor(
            (double) std::accumulate(
            values.begin(),
            values.end(),
            0
        ) / (double) values.size()
    );
}

unsigned stddev(
    const std::vector<unsigned>& values,
    const unsigned& average
) {
    unsigned intermediate = std::accumulate(
        values.begin(),
        values.end(),
        0,
        [&average](const unsigned& carry, const unsigned& value) {
            return (
                carry 
                + std::pow(
                    value - average,
                    2
                )
            );
        }
    );
    return std::sqrt( (double) intermediate / (double) values.size() );
}


vector<unsigned> map_average(
    const std::vector<
        std::vector<unsigned>
    >& list_values
) {
    vector<unsigned> mapped;
    for(const auto& values: list_values) {
        mapped.push_back(
            average(values)
        );
    }
    return mapped;
}

vector<unsigned> map_stddev(
    const std::vector<
        std::vector<unsigned>
    >& list_values,
    const std::vector<unsigned>& averages
) {
    assert(list_values.size() == averages.size());
    vector<unsigned> mapped;
    for(unsigned i = 0; i < averages.size(); i++) {
        mapped.push_back(
            stddev(
                list_values[i],
                averages[i]
            )
        );
    }

    return mapped;
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

    return elapsed.count() * 1e6;
}

template<typename T>
std::vector<unsigned> performance_test(
    T& test_instance,
    const AtomIndexType& n_atoms
) {
    std::vector<unsigned> exec_times;

    // generate all valid combinations if insertable atoms
    std::vector<PairOfAtomIndices> combinations = all_valid_combinations(
        n_atoms
    );

    // make a list of AtomIndexTypes
    std::vector<AtomIndexType> atom_indices(n_atoms);
    std::iota(
        atom_indices.begin(),
        atom_indices.end(),
        0
    );

    // make a random sequence
    std::shuffle(
        combinations.begin(),
        combinations.end(),
        std::mt19937{
            std::random_device{}()
        }
    );
    std::shuffle(
        atom_indices.begin(),
        atom_indices.end(),
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
    combinations.resize(4 * n_atoms);

    // add the atoms
    exec_times.push_back(time_wrap_nullary_callable(
        [&]() {
            for(AtomIndexType i = 0; i < n_atoms; i++) {
                test_instance -> add_atom();
            }
        }
    ));
    
    // add bonds
    exec_times.push_back(time_wrap_nullary_callable(
        [&]() {
            for(const auto& combination : combinations) {
                test_instance -> add_bond(
                    combination.first,
                    combination.second,
                    BondType::Single
                );
            }
        }
    ));

    // bond_exists (1)
    exec_times.push_back(time_wrap_nullary_callable(
        [&]() {
            for(const auto& combination : combinations) {
                test_instance -> bond_exists(
                    combination.first,
                    combination.second
                );
            }
        }
    ));

    // bond_exists (2)
    exec_times.push_back(time_wrap_nullary_callable(
        [&]() {
            for(const auto& combination : combinations) {
                test_instance -> bond_exists(
                    combination.first,
                    combination.second,
                    BondType::Single
                );
            }
        }
    ));

    // get_bond_type
    exec_times.push_back(time_wrap_nullary_callable(
        [&]() {
            for(const auto& combination : combinations) {
                test_instance -> bond_exists(
                    combination.first,
                    combination.second
                );
            }
        }
    ));

    // get_bond_pairs
    exec_times.push_back(time_wrap_nullary_callable(
        [&]() {
            for(const auto& atom_index : atom_indices) {
                test_instance -> get_bond_pairs(atom_index);
            }
        }
    ));

    // get_bonded_atom_indices
    exec_times.push_back(time_wrap_nullary_callable(
        [&]() {
            for(const auto& atom_index : atom_indices) {
                test_instance -> get_bonded_atom_indices(atom_index);
            }
        }
    ));

    // validate
    exec_times.push_back(time_wrap_nullary_callable(
        [&]() {
            test_instance -> validate();
        }
    ));

    // remove_bond
    exec_times.push_back(time_wrap_nullary_callable(
        [&]() {
            for(const auto& combination : combinations) {
                test_instance -> remove_bond(
                    combination.first,
                    combination.second
                );
            }
        }
    ));

    // re-add all bonds
    for(const auto& combination : combinations) {
        test_instance -> add_bond(
            combination.first,
            combination.second,
            BondType::Single
        );
    }

    // remove_atom
    exec_times.push_back(time_wrap_nullary_callable(
        [&]() {
            for(const auto& atom_index : atom_indices) {
                test_instance -> remove_atom(atom_index);
            }
        }
    ));

    std::cout << std::endl;

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

    const unsigned repeat = 10;
    const unsigned max_n_atoms = 100;

    for(auto& test_pair : test_pairs) {
        std::cout << "For class " << test_pair.first << std::endl;
        std::vector<
            std::vector<unsigned>
        > function_averages (function_aliases.size(), std::vector<unsigned>());
        std::vector<
            std::vector<unsigned>
        > function_stddevs (function_aliases.size(), std::vector<unsigned>());

        for(unsigned n_atoms = 10; n_atoms <= max_n_atoms; n_atoms++) {
            std::vector<
                std::vector<unsigned>
            > microsecond_exec_times;

            for(unsigned i = 0; i < repeat; i++) {
                microsecond_exec_times.push_back(
                    performance_test(
                        test_pair.second,
                        n_atoms
                    )
                );
            }

            auto averages = map_average(microsecond_exec_times);
            auto stddevs = map_stddev(
                microsecond_exec_times,
                averages
            );

            /* now I have 
             * (add_atom for n_atoms) <-> averages[0] +- stddevs[0]
             * ...
             *
             * want
             * function_averages[function_index] = [average_2, average_3, ...]
             * function_stddevs[function_index] = [stddev_2, stddev_3, ...]
             */

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
