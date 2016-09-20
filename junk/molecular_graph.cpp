#include "molecular_graph.hpp"
#include <deque>

namespace molecular_graph {
    // BondAdjacencyList
    uint8_t BondAdjacencyList::connected_components() {
        std::vector<uint8_t> component_count;

        /* here, component_count is used in the function by reference, but then
         * discarded
         */
        return _connected_components(component_count);
    }

    uint8_t BondAdjacencyList::_connected_components(
        std::vector<uint8_t>& component_count
    ) {
        component_count.resize(0);
        component_count.resize(_adjacencies.size(), 0);

        uint8_t group_number = 0;
        std::deque<size_t> visit_deque;
        for(unsigned i = 0; i < _adjacencies.size(); i++) {
            if(component_count[i] == group_number) {
                group_number += 1;
                // mark this vertex visited
                component_count[i] = group_number;
                for(const auto& index : _adjacencies[i]) {
                    visit_deque.push_back(index);
                }

                while(!visit_deque.empty()) {
                    // get the frontmost element
                    auto& element_index = visit_deque.front();

                    // mark the vertex visited
                    component_count[element_index] = group_number;

                    // add its connections to the visit deque
                    for(const auto& index : _adjacencies[i]) {
                        visit_deque.push_front(index);
                    }

                    // pop it from the deque
                    visit_deque.pop_front();
                }
            }
        }

        return group_number;
    }

    // MolecularGraph 
    auto MolecularGraph::num_atoms() {
        return _atoms.size();
    }

    auto MolecularGraph::num_groups() {
        return _adjacencies.connected_components();
    }

}
