#include <algorithm>
#include <vector>
#include <iostream>
#include <cassert>
#include <functional>

template<typename T>
std::ostream& operator << (std::ostream& os, const std::vector<T>& vector) {
    for(unsigned i = 0; i < vector.size(); i++) {
        os << vector[i];
        if(i != vector.size() - 1) {
            os << ", ";
        }
    }

    return os;
}

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
const std::vector<
    std::function<Assignment(const Assignment&)>
> EquivalentAssignmentGenerationInstructions::c2_rotations  = {
    /* These are three pairs of switches, so e.g.
     * 13 25 46 means 1 <-> 3, 2 <-> 5, 4 <-> 6
     */
    [](const Assignment& a) -> Assignment {
        // 13 25 46
        return { a[2], a[4], a[0], a[5], a[1], a[3] };
    },
    [](const Assignment& a) -> Assignment {
        // 13 26 45
        return { a[2], a[5], a[0], a[4], a[3], a[1] };
    },
    [](const Assignment& a) -> Assignment {
        // 15 24 36
        return { a[4], a[3], a[5], a[1], a[0], a[2] };
    },
    [](const Assignment& a) -> Assignment {
        // 16 24 35
        return { a[5], a[3], a[4], a[1], a[2], a[0] };
    },
    [](const Assignment& a) -> Assignment {
        // 12 34 56
        return { a[1], a[0], a[3], a[2], a[5], a[4]  };
    },
    [](const Assignment& a) -> Assignment {
        // 14 23 56
        return { a[3], a[2], a[1], a[0], a[5], a[4] };
    }
};

bool is_equivalent_assignment(
    const Assignment& a,
    const Assignment& b
) {
    Assignment variant;
    // try c4 rotations
    for(const auto& c4_rotation_function: 
        EquivalentAssignmentGenerationInstructions::c4_rotations
    ) {
        // overwrite variant with a again
        variant = a;
        for(unsigned i = 0; i < 3; i++) {
            variant = c4_rotation_function(variant);
            if(std::equal(
                variant.begin(),
                variant.end(),
                b.begin()
            )) {
                return true;
            }
        }
    }

    // try c2 rotations
    for(const auto& c2_rotation_function:
        EquivalentAssignmentGenerationInstructions::c2_rotations
    ) {
        variant = c2_rotation_function(a);
        if(std::equal(
            variant.begin(),
            variant.end(),
            b.begin()
        )) {
            return true;
        }
    }

    return false;
}

int main() {
    std::vector<
        std::vector<char>
    > unique_assignments;

    std::vector<char> assignment = {'A', 'A', 'A', 'B', 'B', 'B'};

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
            std::cout << assignment << std::endl;
            unique_assignments.push_back(assignment);
        }
    } while(std::next_permutation(
        assignment.begin(),
        assignment.end()
    ));

    std::cout << " -> " << unique_assignments.size() << std::endl;

    return 0;
}
