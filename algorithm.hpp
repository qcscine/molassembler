#include <vector>
#include <algorithm>
#include <functional>
#include <cassert>
#include <iostream>
#include <sstream>

template<typename T>
std::vector<
    std::function<
        std::vector<T>(
            const std::vector<T>&
        )
    >
> Oh_C4_rotations = {
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

bool Oh_indices_are_trans(
    const unsigned& a,
    const unsigned& b
) {
    if(
        (
            std::min(a, b) == 0
            && std::max(a, b) == 2
        ) || (
            std::min(a, b) == 1
            && std::max(a, b) == 3
        ) || (
            std::min(a, b) == 4
            && std::max(a, b) == 5
        )
    ) return true;
    else return false;
}

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
        else if(this -> character > other.character) return false;
        else {
            for(unsigned i = 0; i < this -> groups.size(); i++) {
                if(this -> groups[i] < other.groups[i]) return true; 
                else if(this -> groups[i] > other.groups[i]) return false; 
                else continue;
            }

            return false;
        }
    }

    bool operator == (const AssignmentColumn& other) const {
        if(this -> character != other.character) return false;

        for(unsigned i = 0; i < this -> groups.size(); i++) {
            if(this -> groups[i] != other.groups[i]) return false;
        }

        return true;
    }

    bool operator != (const AssignmentColumn& other) const {
        return !(
            (*this) == other
        );
    }
};

// TODO temporary declarations
struct Assignment;
std::ostream& operator << (
    std::ostream& os,
    const Assignment& a
);

struct Assignment {
private:
    /* TODO
     * - reduce_groups needs serious testing
     */
    std::vector<
        std::vector<unsigned>
    > _reduce_groups() const {
        std::vector<
            std::vector<unsigned>
        > group_reduction;

        for(unsigned i = 0; i < 6; i++) {
            if(group_reduction.size() < position_occupations[i].groups.size()) {
                group_reduction.resize(
                    position_occupations[i].groups.size()
                );
            }

            for(unsigned j = 0; j < position_occupations[i].groups.size(); j++) {
                if(position_occupations[i].groups[j]) {
                    group_reduction[j].push_back(i);
                }
            }
        }

        // sort! this is required from _reduced_groups_are_equal
        for(auto& reduced_group : group_reduction) {
            std::sort(
                reduced_group.begin(),
                reduced_group.end()
            );
        }

        return group_reduction;
    }

    bool _reduced_groups_are_equal(
        const std::vector<
            std::vector<unsigned>
        >& a,
        const std::vector<
            std::vector<unsigned>
        >& b
    ) const {
        // every element in a must be in b
        if(a.size() != b.size()) return false;
        
        for(const auto& a_vector_element : a) {
            if(!std::accumulate(
                b.begin(),
                b.end(),
                false,
                [&a_vector_element](
                    const bool& carry,
                    const std::vector<unsigned>& b_vector_element
                ) {
                    if(carry) return carry;
                    else return (
                        carry
                        || std::equal(
                            a_vector_element.begin(),
                            a_vector_element.end(),
                            b_vector_element.begin()
                        )
                    );
                }
            )) {
                return false;
            }
        }

        return true;
    }


public:
    std::vector<AssignmentColumn> position_occupations;

    /* Constructors */
    Assignment() = delete;
    Assignment(
        const std::vector<char>& characters
    ) {
        assert(characters.size() == 6);

        for(unsigned i = 0; i < 6; i++) {
            position_occupations.emplace_back(
                characters[i],
                std::vector<bool>()
            );
        }

        sort_occupations();
    }

    Assignment(
        const std::vector<char>& characters,
        const std::vector<
            std::vector<bool>
        >& ligand_connections
    ) {
        assert(characters.size() == 6);
        for(const auto& ligand_connection : ligand_connections) {
            assert(ligand_connection.size() == 6);
        }

        for(unsigned i = 0; i < 6; i++) {
            std::vector<bool> groups;
            for(const auto& ligand_connection: ligand_connections) {
                groups.push_back(ligand_connection[i]);
            }

            position_occupations.emplace_back(
                characters[i],
                groups
            );
        }

        sort_occupations();
    }

    /* Modification */
    void sort_occupations() {
        std::sort(
            position_occupations.begin(),
            position_occupations.end()
        );
    }

    bool next_permutation() {
        while(std::next_permutation(
            position_occupations.begin(),
            position_occupations.end()
        )) {
            if(ligand_connections_are_ordered()) {
                return true;
            }
        }
        return false;
    }

    /* doesn't work! an issue with operator < ?
    bool prev_permutation() {
        while(std::prev_permutation(
            position_occupations.begin(),
            position_occupations.end()
        )) {
            std::cout << "In std::prev_perm loop: " << std::endl;
            std::cout << *this << std::endl;
            if(ligand_connections_are_ordered()) {
                return true;
            }
        }
        return false;
    }*/

    void rotate(const unsigned& rotation_function_index) {
        position_occupations = Oh_C4_rotations<AssignmentColumn>[
            rotation_function_index
        ](position_occupations);
    }

    /* Information */
    bool occupations_are_ordered() {
        return std::is_sorted(
            position_occupations.begin(),
            position_occupations.end()
        );
    }

    std::vector<
        std::vector<bool>
    > get_groups_row_view() const {
        std::vector<
            std::vector<bool>
        > row_view;

        for(const auto& occupation : position_occupations) {
            if(row_view.size() < occupation.groups.size()) {
                row_view.resize(occupation.groups.size());
            }

            for(unsigned i = 0; i < occupation.groups.size(); i++) {
                row_view[i].push_back(occupation.groups[i]);
            }
        }

        return row_view;
    }

    bool ligand_connections_are_ordered() const {
        // shortcut if row_view will be empty
        if(position_occupations[0].groups.size() == 0) return true;

        std::vector<
            std::vector<bool>
        > row_view = get_groups_row_view();

        for(unsigned i = 0; i < row_view.size() - 1; i++) {
            // compare row i with row i+1
            for(unsigned j = 0; j < 6; j++) {
                if(row_view[i][j] < row_view[i + 1][j]) break;
                if(row_view[i][j] > row_view[i + 1][j]) return false;
            }
        }

        return true;
    }

    bool is_rotationally_superimposable(
        const Assignment& other
    ) const {
        Assignment twister = *this;
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
                twister.rotate(rotation_function_index);
                if(twister == other) {
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
            twister.rotate(rotation_function_index);
        }

        return false;
    }

    /* Operators */
    bool operator == (
        const Assignment& other
    ) const {
        // compare characters
        for(unsigned i = 0; i < 6; i++) {
            if(
                this -> position_occupations[i].character != 
                    other.position_occupations[i].character
            ) {
                return false;
            }
        }

        // compare reduced groups
        if(!_reduced_groups_are_equal(
            _reduce_groups(),
            other._reduce_groups()
        )) {
            return false;
        }
        
        return true;
    }
};

std::ostream& operator << (
    std::ostream& os,
    const Assignment& a
) {
    // make group.size + 1 stringstreams
    std::vector<
        std::stringstream
    > streams (
        1 + a.position_occupations[0].groups.size()
    );

    for(const auto& column : a.position_occupations ) {
        streams[0] << column.character << " ";
        for(unsigned i = 0; i < column.groups.size(); i++) {
            streams[i + 1] << (column.groups[i] ? "T" : "F") << " ";
        }
    }

    for(const auto& stream : streams) {
        os << stream.str() << std::endl;
    }

    return os;
}

std::vector<Assignment> map_remove_trans_spanning_groups(
    const std::vector<Assignment>& assignments
) {
    std::vector<Assignment> map;

    for(const auto& assignment : assignments) {
        bool passes_checks = true;

        auto row_view = assignment.get_groups_row_view();
        for(unsigned i = 0; i < row_view.size(); i++) {
            unsigned num_members = std::accumulate(
                row_view[i].begin(),
                row_view[i].end(),
                0,
                [](const unsigned& carry, const bool& is_member) {
                    if(is_member) return carry + 1;
                    else return carry;
                }
            );
            /* if a group has two members and they are trans, remove the
             * assignment
             */
            if(num_members == 2) {
                // are they trans to each other?
                // collect the indices
                std::vector<unsigned> pair;
                for(unsigned j = 0; j < 6; j++) {
                    if(row_view[i][j]) pair.push_back(j);
                }
                // check
                if(Oh_indices_are_trans(
                    pair[0],
                    pair[1]
                )) {
                    passes_checks = false;
                    break;
                }
            }
        }

        if(passes_checks) {
            map.push_back(assignment);
        }
    }

    return map;
}

/* Gives NO guarantees as to satisfiability.
 * E.g. M (A-A)3 generates a trans-trans-trans assignment, which is extremely 
 *  hard to find actual ligands for that work.
 * The satisfiability of assignments must be checked before trying to embed 
 *  structures with completely nonsensical constraints. Perhaps restrict A-A 
 *  ligands with bridge length 4 (chelating atoms included), maybe even up to 6
 *  to cis arrangements. Xantphos (with bridge length 7 is the smallest 
 *  trans-spanning ligand mentioned in Wikipedia).
 */
std::vector<Assignment> unique_assignments(
    const Assignment& initial,
    const bool& remove_trans_spanning_groups = true
) {
    std::vector<
        Assignment
    > unique_assignments;

    Assignment assignment = initial;

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
                        && !unique_assignment.is_rotationally_superimposable(
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

    if(remove_trans_spanning_groups) {
        return map_remove_trans_spanning_groups(
            unique_assignments
        );
    } else return unique_assignments;
}

