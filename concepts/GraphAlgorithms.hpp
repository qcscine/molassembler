#include "AdjacencyList.hpp"

/* TODO
 * - Test
 * - Cycle detection -> see Tucker, Alan. "Applied Combinatorics p.49
 */

namespace MoleculeManip {

    /*!
     * Connected Components algorithm. Returns a vector of unsigned numbers 
     * that maps AtomIndexType -> Connected component group ID.
     * \param adjacencies The AdjacencyList instance to process
     * \returns A vector of unsigned numbers that maps AtomIndexType -> 
     * ConnectedComponentType.
     */
    std::vector<unsigned> connected_components(const AdjacencyList& adjacencies) {
        if(adjacencies.size() == 0) return {};

        std::vector<unsigned> visited (
            adjacencies.size(),
            0
        );
        std::vector<unsigned>::iterator vis_iter = visited.begin();
        std::vector<AtomIndexType> to_visit = {0};
        unsigned counter = 1;
        AtomIndexType current;

        while(to_visit.size() > 0) {
            // take an element
            current = to_visit[0];
            to_visit.pop_back();

            // mark it
            visited[current] = counter;

            // add it's connections to to_visit
            for(const auto& connected_index : adjacencies[current]) {
                if(visited[connected_index] != 0) {
                    to_visit.push_back(connected_index);
                }
            }

            // if to_visit is empty, increment the counter and add the next 
            // non-visited element to to_visit
            if(
                to_visit.size() == 0 &&
                vis_iter != visited.end()
            ) {
                // move vis_iter to the next unvisited element
                while(
                    vis_iter != visited.end() 
                    && *vis_iter != 0
                ) {
                    vis_iter++;
                }
                // if we're not at the end, then
                if(vis_iter != visited.end()) {
                    // increment the counter
                    counter += 1;
                    // and add it's index to to_visit
                    to_visit.push_back(
                        vis_iter - visited.begin()
                    );
                }
            }
        }

        // all atoms must be visited
        assert(std::accumulate(
            visited.begin(),
            visited.end(),
            true,
            [](const bool& carry, const unsigned& element) {
                return (
                    carry 
                    && element != 0
                );
            }
        ));

        return visited;
    }

    /*!
     * Returns the number of connected components in an AdjacencyList.
     * \param adjacencies The AdjacencyList instance to process
     * \returns The number of connected components in the AdjacencyList.
     */
    unsigned num_connected_components(const AdjacencyList& adjacencies) {
        auto visited = connected_components(adjacencies);
        if(visited.size() == 0) return 0;
        else return *max_element(visited.begin(), visited.end());
    }

    /*!
     * Constructs a list of AtomIndexType groups that are connected in an 
     * AdjacencyList.
     * \param adjacencies The AdjacencyList instance to process
     * \returns An unordered vector of vectors containing atom indices that are
     * connected.
     */
    std::vector<
        std::vector<AtomIndexType>
    > connected_component_groups(const AdjacencyList& adjacencies) {
        auto visited = connected_components(adjacencies);

        // guard against 0-length
        if(visited.size() == 0) return std::vector<
            std::vector<AtomIndexType>
        >();

        unsigned num_groups = *max_element(visited.begin(), visited.end());
        std::vector<
            std::vector<AtomIndexType>
        > groups (num_groups);

        // go through visited
        // it could be e.g. [1, 2, 1, 3, 4, 3]
        for(auto it = visited.begin(); it != visited.end(); it++) {
            // make space if the group number is bigger
            if(groups.size() < *it) {
                groups.resize(*it);
            }

            groups[(*it) - 1].push_back( 
                it - visited.begin() 
            );
        }

        return groups;
    } 
}
