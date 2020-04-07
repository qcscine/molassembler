/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Contains class to help with equal-size index partitioning
 *
 * Contains only the symmetry names and count for minimal header inclusion
 * necessities in dependencies
 */

#include <vector>

namespace Scine {
namespace shapes {

/**
 * @brief Given S * E distinguishable objects, this class helps enumerate all
 *   partitions into S groups of size E.
 *
 * The challenge with writing such a class is getting complexity less than
 * @math{\Theta((S\cdot E)!)} considering that the sets themselves are
 * indistinguishable. For instance, for the trivial example @math{S = 2, E = 2},
 * there are only three distinct solutions:
 *
 * - {1, 2}, {3, 4}
 * - {1, 3}, {2, 4}
 * - {1, 4}, {2, 3}
 *
 * There are two orders at play here, the sets are lexicographically ordered
 * and the elements within them are too.
 */
class Partitioner {
public:
  /*! @brief Constructor
   *
   * @param s Number of groups
   * @param e Number of elements per group
   *
   * @pre @p s is not zero and @p e is not zero
   *
   * @complexity{@math{\Theta(S\cdot E)}}
   */
  Partitioner(const unsigned s, const unsigned e);

  /*! @brief Advance underlying state to the next partition
   *
   * @complexity{@math{O(S\cdot E)}}
   *
   * @returns whether the underlying state reflects a new, valid partition
   */
  bool next_partition();

  /*! @brief Generate sets of element indices within their respective partitions
   *
   * @complexity{@math{\Theta(S\cdot E)}}
   */
  std::vector<
    std::vector<unsigned>
  > partitions() const;

  /*! @brief Check whether the sets as specified by the mapping are ordered
   *
   * @complexity{@math{O(S\cdot E)}}
   */
  static bool isOrderedMapping(const std::vector<unsigned>& mapping);

  //! Number of groups
  inline unsigned s() const {
    return S;
  }

  //! Number of elements per group
  inline unsigned e() const {
    return E;
  }

  //! Access to underlying flat map from element index to group index
  inline const std::vector<unsigned>& map() const {
    return mapping;
  }

private:
  //! Number of groups
  unsigned S;
  //! Number of elements per group
  unsigned E;
  //! Flat map from element index to group index
  std::vector<unsigned> mapping;

};

} // namespace shapes
} // namespace Scine
