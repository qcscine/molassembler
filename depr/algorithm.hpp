#include <algorithm>
#include <vector>
#include <iostream>
#include <cassert>
#include <functional>

using Assignment = std::vector<char>;

struct LigandConnection {
    /* Data */
    std::vector<bool> connection;

    /* Constructor */
    LigandConnection(
        const std::vector<bool>& pass_connection
    ) : connection(pass_connection) {
        assert(pass_connection.size() == 6);
    };

    /* Operators */
    bool operator [] (const unsigned& pos) const {
        return connection[pos];
    }
    bool operator < (const LigandConnection& other) {
        for(unsigned i = 0; i < 6; i++) {
            if(
                !this -> connection[i] // false at i in this
                && other[i] // AND true at i in other
            ) return true; // then this is smaller than other
        }

        return false;
    }

    /* Iterators */
    std::vector<bool>::const_iterator begin() {
        return connection.begin();
    }
    std::vector<bool>::const_iterator end() {
        return connection.end();
    }
};

struct AssignmentColumn {
    char character;
    std::vector<bool> groups;

    AssignmentColumn(
        const char& pass_character,
        const std::vector<bool> pass_groups
    ) : 
        character(pass_character),
        groups(pass_groups) 
    {};

    bool operator < (const AssignmentColumn& other) const {
        assert(this -> groups.size() == other.groups.size());
        if(this -> character < other.character) return true;

        for(unsigned i = 0; i < this -> groups.size(); i++) {
            if((
                !(this -> groups[i]) // false at i in this
                && other.groups[i] // and true at i in other
            )) return true; // then this is smaller
        }

        return false;
    }
};

// TODO rename after refactor
struct AssignmentStruct {
    /* Members */
    std::vector<char> characters;
    std::vector<
        std::vector<bool>
    > ligand_connections;

    /* Constructor */
    AssignmentStruct(
        const std::vector<char>& pass_characters,
        const std::vector<
            std::vector<bool>
        >& pass_ligand_connections
    ) : 
        characters(pass_characters),
        ligand_connections(pass_ligand_connections)
    {};

    /* Information */
    std::vector<AssignmentColumn> get_columns() const {
        std::vector<AssignmentColumn> columns;
        for(unsigned i = 0; i < 6; i++) {
            // collect group information
            std::vector<bool> groups;
            for(const auto& connection : ligand_connections) {
                groups.push_back(connection[i]);
            }

            // add to columns
            columns.emplace_back(
                characters[i],
                groups
            );
        }

        return columns;
    }

    /* Operator */
    bool operator == (
        const AssignmentStruct& other
    ) {
        // character permutations equal
        if(!std::equal(
            this -> characters.begin(),
            this -> characters.end(),
            other.characters.begin()
        )) {
            return false;
        }

        // number of connected ligand vectors
        if(this -> ligand_connections.size() != other.ligand_connections.size()) {
            return false;
        }

        // ligand connection vectors
        for(const auto& ligand_connection_vector : this -> ligand_connections) {
            // every connection vector must be contained in the other's set
            if(!std::accumulate(
                other.ligand_connections.begin(),
                other.ligand_connections.end(),
                false,
                [&ligand_connection_vector](
                    const bool& carry,
                    const std::vector<bool>& other_ligand_connection_vector
                ) {
                    return (
                        carry
                        || std::equal(
                            ligand_connection_vector.begin(),
                            ligand_connection_vector.end(),
                            other_ligand_connection_vector.begin()
                        )
                    );
                }
            )) {
                // accumulate returns whether the current vector is in the other
                // if it's not, then the AssignmentStructs are unequal, so:
                return false;
            }
        }

        return true;
    }
};

struct EquivalentAssignmentGenerationInstructions {
    Assignment rotate(
        const Assignment& a,
        const std::function<Assignment(const Assignment&)> rotation,
        const unsigned& n_times
    ) {
        Assignment return_value = a;
        for(unsigned i = 0; i < n_times; i++) {
            return_value = rotation(a);
        }
        return return_value;
    }

    static AssignmentStruct rotate(
        const AssignmentStruct& a,
        const unsigned& rotation_index
    ) {
        assert(rotation_index < 3);
        // rotate the characters
        auto rotated_characters = 
            EquivalentAssignmentGenerationInstructions::Oh_C4_rotations<
                char
            >[rotation_index](a.characters);

        // rotate ligand connections
        std::vector<
            std::vector<bool>
        > rotated_ligand_connections;
        for(const auto& ligand_connection_vector : a.ligand_connections) {
            rotated_ligand_connections.emplace_back(
                EquivalentAssignmentGenerationInstructions::Oh_C4_rotations<
                    bool
                >[rotation_index](ligand_connection_vector)
            );
        }

        return AssignmentStruct(
            rotated_characters,
            rotated_ligand_connections
        );
    }

    static const std::vector<
        std::function<Assignment(const Assignment&)>
    > c4_rotations;

    template<typename T>
    static std::vector<
        std::function<
            std::vector<T>(
                const std::vector<T>&
            )
        >
    > Oh_C4_rotations;
};

const std::vector<
    std::function<Assignment(const Assignment&)>
> EquivalentAssignmentGenerationInstructions::c4_rotations = {
    [](const Assignment& a) -> Assignment { 
        // this is for {1, 2, 3, 4}
        return {
            a[3], // 1 <- 4 = [3]
            a[0], // 2 <- 1 = [0]
            a[1], // 3 <- 2 = [1]
            a[2], // 4 <- 3 = [2]
            a[4], // 5
            a[5]  // 6
        };
    },
    [](const Assignment& a) -> Assignment { 
        // this is for {2, 5, 4, 6}
        return {
            a[0], // 1
            a[5], // 2 <- 6 = [5]
            a[2], // 3
            a[4], // 4 <- 5 = [4]
            a[1], // 5 <- 2 = [1]
            a[3]  // 6 <- 4 = [3]
        };
    },
    [](const Assignment& a) -> Assignment { 
        // this is for {1, 6, 3, 5}
        return {
            a[4], // 1 <- 5 = [4]
            a[1], // 2
            a[5], // 3 <- 6 = [5]
            a[3], // 4 
            a[2], // 5 <- 3 = [2]
            a[0]  // 6 <- 1 = [0]
        };
    },
};


template<typename T>
std::vector<
    std::function<
        std::vector<T>(
            const std::vector<T>&
        )
    >
> EquivalentAssignmentGenerationInstructions::Oh_C4_rotations = {
    [](const std::vector<T>& a) -> std::vector<T> { 
        // this is for {1, 2, 3, 4}
        return {
            a[3], // 1 <- 4 = [3]
            a[0], // 2 <- 1 = [0]
            a[1], // 3 <- 2 = [1]
            a[2], // 4 <- 3 = [2]
            a[4], // 5
            a[5]  // 6
        };
    },
    [](const std::vector<T>& a) -> std::vector<T> { 
        // this is for {2, 5, 4, 6}
        return {
            a[0], // 1
            a[5], // 2 <- 6 = [5]
            a[2], // 3
            a[4], // 4 <- 5 = [4]
            a[1], // 5 <- 2 = [1]
            a[3]  // 6 <- 4 = [3]
        };
    },
    [](const std::vector<T>& a) -> std::vector<T> { 
        // this is for {1, 6, 3, 5}
        return {
            a[4], // 1 <- 5 = [4]
            a[1], // 2
            a[5], // 3 <- 6 = [5]
            a[3], // 4 
            a[2], // 5 <- 3 = [2]
            a[0]  // 6 <- 1 = [0]
        };
    },
};

bool is_equivalent_assignment(
    const Assignment& a,
    const Assignment& b
) {
    Assignment variant = a;
    unsigned rotation_function_index;
    // try all combinations of c4 rotations
    for(unsigned i = 0; i < 6; i++) {
        /* spin
         *
         * but which function?
         *
         * sequence in dependence of i is
         * index i 0 1 2 3 4 5 
         * group g 1 3 2 1 3 2 = 3 - ( (i + 2) % 3 )
         * index [] = 2 - ( (i + 2) % 3)
         */
        rotation_function_index = 2 - ( (i + 2) % 3 );

        for(unsigned j = 0; j < 4; j++) {
            variant = EquivalentAssignmentGenerationInstructions::c4_rotations[
                rotation_function_index
            ](variant);
            if(std::equal(
                variant.begin(),
                variant.end(),
                b.begin()
            )) {
                return true;
            }
        }

        /* rotate to next top position
         *
         * but which rotation to use?
         *
         * sequence in dependence of i
         * index i 0 1 2 3 4 5
         * group g 2 1 3 2 1 3 = 3 - ( (i + 1) % 3)
         * index [] = 2 - ( (i + 1) % 3
         */
        
        rotation_function_index = 2 - ( (i + 1) % 3 );
        variant = EquivalentAssignmentGenerationInstructions::c4_rotations[
            rotation_function_index
        ](variant);
    }

    return false;
}

bool is_equivalent_assignment(
    const AssignmentStruct& a,
    const AssignmentStruct& b
) {
    AssignmentStruct variant = a;
    unsigned rotation_function_index;
    // try all combinations of c4 rotations
    for(unsigned i = 0; i < 6; i++) {
        /* spin
         *
         * but which function?
         *
         * sequence in dependence of i is
         * index i 0 1 2 3 4 5 
         * group g 1 3 2 1 3 2 = 3 - ( (i + 2) % 3 )
         * index [] = 2 - ( (i + 2) % 3)
         */
        rotation_function_index = 2 - ( (i + 2) % 3 );

        for(unsigned j = 0; j < 4; j++) {
            variant = EquivalentAssignmentGenerationInstructions::rotate(
                variant,
                rotation_function_index
            );
            
            if(variant == b) return true;
        }

        /* rotate to next top position
         *
         * but which rotation to use?
         *
         * sequence in dependence of i
         * index i 0 1 2 3 4 5
         * group g 2 1 3 2 1 3 = 3 - ( (i + 1) % 3)
         * index [] = 2 - ( (i + 1) % 3
         */
        
        rotation_function_index = 2 - ( (i + 1) % 3 );
        variant = EquivalentAssignmentGenerationInstructions::rotate(
            variant,
            rotation_function_index
        );
    }

    return false;
}

std::vector<Assignment> unique_assignments(const Assignment& initial) {
    std::vector<
        std::vector<char>
    > unique_assignments;

    std::vector<char> assignment = initial;

    do {
        bool current_assignment_is_distinct = std::accumulate(
            unique_assignments.begin(),
            unique_assignments.end(),
            true,
            [&assignment](
                const bool& carry,
                const Assignment& unique_assignment
            ) {
                if(carry) { 
                    // only bother to compute if the accumulation is still true
                    return (
                        carry
                        && !is_equivalent_assignment(
                            unique_assignment,
                            assignment
                        )
                    );
                } else {
                    return carry;
                }
            }
        );
        if(current_assignment_is_distinct) {
            unique_assignments.push_back(assignment);
        }
    } while(std::next_permutation(
        assignment.begin(),
        assignment.end()
    ));

    return unique_assignments;
}

std::vector<AssignmentStruct> unique_assignments(const AssignmentStruct& initial) {
    std::vector<
        AssignmentStruct
    > unique_assignments;

    AssignmentStruct assignment = initial;

    do {
        bool current_assignment_is_distinct = std::accumulate(
            unique_assignments.begin(),
            unique_assignments.end(),
            true,
            [&assignment](
                const bool& carry,
                const AssignmentStruct& unique_assignment
            ) {
                if(carry) { 
                    // only bother to compute if the accumulation is still true
                    return (
                        carry
                        && !is_equivalent_assignment(
                            unique_assignment,
                            assignment
                        )
                    );
                } else {
                    return carry;
                }
            }
        );
        if(current_assignment_is_distinct) {
            unique_assignments.push_back(assignment);
        }
    } while(assignment.next_permutation());

    return unique_assignments;
}
