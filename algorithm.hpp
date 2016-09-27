#include <algorithm>
#include <vector>
#include <iostream>
#include <cassert>
#include <functional>

using Assignment = std::vector<char>;

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

    static const std::vector<
        std::function<Assignment(const Assignment&)>
    > c4_rotations;
    static const std::vector<
        std::function<Assignment(const Assignment&)>
    > c2_rotations;
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
