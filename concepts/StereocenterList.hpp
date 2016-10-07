#include <vector>
#include <algorithm>
#include "common_typedefs.hpp"

namespace MoleculeManip {
    struct Stereocenter {
    /* Public members */
        std::vector<AtomIndexType> atoms;
        StereocenterType type;
        bool is_unspecified;

    /* Public member functions */
        /* Constructors */
        Stereocenter(
            const std::vector<AtomIndexType>& pass_atoms,
            const StereocenterType& pass_type,
            const bool& pass_is_unspecified
        ) :
            atoms (pass_atoms),
            type (pass_type),
            is_unspecified (pass_is_unspecified)
        {};
    };

    class StereocenterList {
    private:
        std::vector<Stereocenter> _stereocenters;

    public:
    /* Public member functions */
        /* Modification */
        void add(const Stereocenter& stereocenter) {
            _stereocenters.push_back(stereocenter);
        }

        void add(
            const std::vector<AtomIndexType> atoms,
            const StereocenterType& type,
            const bool& is_unspecified
        ) {
            _stereocenters.emplace_back(
                atoms,
                type,
                is_unspecified
            );
        }

        void remove_involving(const AtomIndexType& a);

        /* Information */
        std::vector<Stereocenter> get_all_matching_type(
            const StereocenterType& type
        ) const {
            std::vector<Stereocenter> matches;
            std::copy_if(
                _stereocenters.begin(),
                _stereocenters.end(),
                matches.begin(),
                [&type](const Stereocenter& stereocenter) {
                    return stereocenter.type == type;
                }
            );
            return matches;
        }

        /* Operators */
        std::vector<Stereocenter>::const_iterator begin() const {
            return _stereocenters.begin();
        }
        std::vector<Stereocenter>::const_iterator end() const {
            return _stereocenters.end();
        }

    };
}
