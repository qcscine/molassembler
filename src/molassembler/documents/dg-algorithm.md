# notes

2.2 and 2.3 are independent, should be summarized as one step

# runDistanceGeometry algorithm steps

1. Figure out if the algorithm has to regenerate the distance bounds matrix at
   each step by looking through the list of stereocenters and checking if any
   are unassigned.

   - Yes: Leave DGData uninitialized
   - No: Initialize DGData immediately from gatherDGInformation()

2. As long as the number of failed refinements does not exceed the threshold to
   bailing out, for each structure:

   1. In case we need to assign stereocenters: Copy out Molecule, progressively
      assign stereocenters, and assign DGData from gatherDGInformation()
   2. Make a distances matrix from the distance bounds
   3. Propagate the chirality constraints using distance bounds and prototypes 
   4. Make a metric matrix from the distances matrix
   5. Embed the positions from the metric matrix
   6. Vectorized the positions into a dlib vector
   7. Invert y coordinates if desired and the proportion of correct chirality
      constraints is < 0.5
   8. If not all chirality constraints are correct, run a four-dimensional
      coordinate refinement without compression until there are no
      correctly-signed chirality constraints
   9. Finish refinement with a four-dimensional refinement including compression
   10. Check the postconditions for refinement failure
