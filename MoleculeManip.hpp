#include "I_MolecularGraph.hpp"

namespace MoleculeManip {
    using MoleculeContainer = std::list<Molecule>;
    using unsigned_type = unsigned;

    // TODO: Implement
    struct SettingsCollection;

    /* Molecule boolean tests */
    // TODO: how to reference bond, is second argument for the following
    bool is_bond_rotatable(const Molecule& mol);
    bool is_bond_incrementable(const Molecule& mol);

    /* Map functions 
     *
     * These can be used thusly:
     *
     * MoleculeContainer initial_structures;
     * 
     * // MapFunctionContainer would be a STL Container of std::function
     * MapFunctionContainer functions = {
     *     permute_dummy_atoms,
     *     permute_metal_environment,
     *     ...
     * };
     *
     * SettingsCollection settings;
     *
     * unsigned_type stage = 0;
     *
     * // bind all functions to make them unary
     * for(auto& function: functions) {
     *     function = std::bind(
     *         function,
     *         std::placeholders::_1,
     *         std::cref(settings),
     *         std::cref(stage)
     *      );
     * }
     *
     * // fold the MapFunctionContainer on the set of initial structures
     * MoleculeContainer std::accumulate(
     *     functions.begin(),
     *     functions.end()
     *     initial_structures,
     *     [](
     *         const MoleculeContainer& carry,
     *         const std::function<
     *            MoleculeContainer(const MoleculeContainer&)
     *         > function
     *     ) {
     *         return function(carry);
     *     }
     *  );
     *
     */
    /*!
     * Permute all atoms with Delib::ElementType::none to metal molecules
     * persuant to the settings. This mapping function is only relevant to 
     * first shell generation. 
     * \param molecules A container with molecules to permute
     * \param settings The settings object
     * \param stage The current generation stage
     */
    MoleculeContainer map_permute_dummy_atoms(
        const MoleculeContainer& molecules,
        const SettingsCollection& settings,
        const unsigned& stage
    );

    /*!
     * Generate hydrogen environments around metal atoms, taking into account 
     * the existing ligands. Permutes over oxidation states, total charge, and
     * ... (update when implementated). This mapping function is only relevant 
     * to first shell generation.
     * \param molecules A container with molecules to permute
     * \param settings The settings object
     * \param stage The current generation stage
     */
    MoleculeContainer map_permute_metal_environment(
        const MoleculeContainer& molecules,
        const SettingsCollection& settings,
        const unsigned& stage
    );

    /*!
     * Permute and generate multidentate ligands.
     * \param molecules A container with molecules to permut
     * \param settings The settings object
     * \param stage The current generation stage
     */
    MoleculeContainer map_permute_denticity(
        const MoleculeContainer& molecules,
        const SettingsCollection& settings,
        const unsigned& stage
    );

    /*!
     * Permute and generate ring closures.
     * \param molecules A container with molecules to permute
     * \param settings The settings object
     * \param stage The current generation stage
     */
    MoleculeContainer map_permute_ring_closures(
        const MoleculeContainer& molecules,
        const SettingsCollection& settings,
        const unsigned& stage
    );

    /*!
     * Permute currently unspecified stereocenters.
     * \param molecules A container with molecules to permute
     * \param settings The settings object
     * \param stage The current generation stage
     */
    MoleculeContainer map_permute_unspecified_stereocenters(
        const MoleculeContainer& molecules,
        const SettingsCollection& settings,
        const unsigned& stage
    );

    /*!
     * Remove redundant structures from the MoleculeContainer.
     * \param molecules A container with molecules to permute
     * \param settings The settings object
     * \param stage The current generation stage
     */
    MoleculeContainer map_remove_redundants(
        const MoleculeContainer& molecules,
        const SettingsCollection& settings,
        const unsigned& stage
    );

    /*!
     * Sample the contained Molecules to a limited amount.
     * \param molecules A container with molecules to permute
     * \param settings The settings object
     * \param stage The current generation stage
     */
    MoleculeContainer map_sample(
        const MoleculeContainer& molecules,
        const SettingsCollection& settings,
        const unsigned& stage
    );

}
