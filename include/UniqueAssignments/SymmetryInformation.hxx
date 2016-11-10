namespace PermSymmetry {

template<typename T>
const std::vector<
  std::pair<
    std::function<
      std::vector<T>(
        const std::vector<T>&
      )
    >,
    uint8_t
  >
> Linear<T>::rotations = {
  std::make_pair(
    [](const std::vector<T>& a) -> std::vector<T> {
      // C2
      return {
        a[1],
        a[0]
      };
    },
    2
  )
};

template<typename T>
const std::vector<
  std::pair<
    std::function<
      std::vector<T>(
        const std::vector<T>&
      )
    >,
    uint8_t
  >
> TrigonalPlanar<T>::rotations = {
  std::make_pair(
    [](const std::vector<T>& a) -> std::vector<T> {
      // C3
      return {
        a[1],
        a[2],
        a[0]
      };
    },
    3
  ),
  std::make_pair(
    [](const std::vector<T>& a) -> std::vector<T> {
      // C2
      return {
        a[0],
        a[2],
        a[1]
      };
    },
    2
  )
};

template<typename T>
const std::vector<
  std::pair<
    std::function<
      std::vector<T>(
        const std::vector<T>&
      )
    >,
    uint8_t
  >
> Tetrahedral<T>::rotations = {
  std::make_pair(
    [](const std::vector<T>& a) -> std::vector<T> { 
      // this is for {2, 3, 4}
      return {
        a[0], 
        a[3], // 2 <- 4 = [3]
        a[1], // 3 <- 2 = [1]
        a[2]  // 4 <- 3 = [2]
      };
    },
    3
  ),
  std::make_pair(
    [](const std::vector<T>& a) -> std::vector<T> { 
      // this is for {1, 3, 4}
      return {
        a[2], // 1 <- 3 = [2]
        a[1], 
        a[3], // 3 <- 4 = [3]
        a[0]  // 4 <- 1 = [0]
      };
    },
    3
  ),
  std::make_pair(
    [](const std::vector<T>& a) -> std::vector<T> { 
      // this is for {1, 2, 4}
      return {
        a[3], // 1 <- 4 = [3]
        a[0], // 2 <- 1 = [0]
        a[2], 
        a[1]  // 4 <- 2 = [1]
      };
    },
    3
  ),
  std::make_pair(
    [](const std::vector<T>& a) -> std::vector<T> { 
      // this is for {1, 2, 3}
      return {
        a[1], // 1 <- 2 = [1]
        a[2], // 2 <- 3 = [2]
        a[0], // 3 <- 1 = [0]
        a[3]  
      };
    },
    3
  )
};

template<typename T>
const std::vector<
  std::pair<
    std::function<
      std::vector<T>(
        const std::vector<T>&
      )
    >,
    uint8_t
  >
> SquarePlanar<T>::rotations = {
  std::make_pair(
    [](const std::vector<T>& a) -> std::vector<T> { 
      // C4
      return {
        a[3], // 1 <- 4 = [3]
        a[0], // 2 <- 1 = [0]
        a[1], // 3 <- 2 = [1]
        a[2]  // 4 <- 3 = [2]
      };
    },
    4
  ),
  std::make_pair(
    [](const std::vector<T>& a) -> std::vector<T> { 
      // C2
      return {
        a[1], // 1 <- 2 = [1]
        a[0], // 2 <- 1 = [0]
        a[3], // 3 <- 4 = [3]
        a[2]  // 4 <- 3 = [2]
      };
    },
    2
  ),
  /* this might be redundant, I think one C2 and the C4 combine to make
   * anything you can rotate with this, but it doesn't really matter 
   * since the search algorithm fill figure this out.
   */
  std::make_pair(
    [](const std::vector<T>& a) -> std::vector<T> { 
      // C2
      return {
        a[3], // 1 <- 4 = [3]
        a[2], // 2 <- 3 = [2]
        a[1], // 3 <- 2 = [1]
        a[0]  // 4 <- 1 = [0]
      };
    },
    2
  )
};

template<typename T>
const std::vector<
  std::pair<
    std::function<
      std::vector<T>(
        const std::vector<T>&
      )
    >,
    uint8_t
  >
> SquarePyramidal<T>::rotations = {
  std::make_pair(
    [](const std::vector<T>& a) -> std::vector<T> { 
      // C4
      return {
        a[3], // 1 <- 4 = [3]
        a[0], // 2 <- 1 = [0]
        a[1], // 3 <- 2 = [1]
        a[2], // 4 <- 3 = [2]
        a[4]  
      };
    },
    4
  ),
};

template<typename T>
const std::vector<
  std::pair<
    std::function<
      std::vector<T>(
        const std::vector<T>&
      )
    >,
    uint8_t
  >
> TrigonalBiPyramidal<T>::rotations = {
  std::make_pair(
    [](const std::vector<T>& a) -> std::vector<T> { 
      // C3
      return {
        a[2], // 1 <- 3 = [2]
        a[0], // 2 <- 1 = [0]
        a[1], // 3 <- 2 = [1]
        a[3], 
        a[4] 
      };
    },
    3
  ),
  std::make_pair(
    [](const std::vector<T>& a) -> std::vector<T> { 
      // C2 on 1
      return {
        a[0], 
        a[2], // 2 <- 3 = [2]
        a[1], // 3 <- 2 = [1]
        a[4], // 4 <- 5 = [4]
        a[3]  // 5 <- 4 = [3]
      };
    },
    2
  ),
  std::make_pair(
    [](const std::vector<T>& a) -> std::vector<T> { 
      // C2 on 2
      return {
        a[2], // 1 <- 3 = [2]
        a[1], 
        a[0], // 3 <- 1 = [0]
        a[4], // 4 <- 5 = [4]
        a[3]  // 5 <- 4 = [3]
      };
    },
    2
  ),
  std::make_pair(
    [](const std::vector<T>& a) -> std::vector<T> { 
      // C2 on 3
      return {
        a[1], // 1 <- 2 = [1]
        a[0], // 2 <- 1 = [0]
        a[2], 
        a[4], // 4 <- 5 = [4]
        a[3]  // 5 <- 4 = [3]
      };
    },
    2
  )
};

template<typename T>
const std::vector<
  std::pair<
    std::function<
      std::vector<T>(
        const std::vector<T>&
      )
    >,
    uint8_t
  >
> TrigonalBiPyramidal<T>::pseudorotations = {
  std::make_pair(
    [](const std::vector<T>& a) -> std::vector<T> { 
      /* Pseudorotation with 1 fixed, combines with C3 to perform 
       * all pseudorotations, that's why not all are listed
       */
      return {
        a[0], 
        a[4], // 2 <- 5 = [4]
        a[3], // 3 <- 4 = [3]
        a[2], // 4 <- 3 = [2]
        a[1]  // 5 <- 2 = [1]
      };
    },
    2
  )
};

template<typename T>
const std::vector<
  std::pair<
    std::function<
      std::vector<T>(
        const std::vector<T>&
      )
    >,
    uint8_t
  >
> Octahedral<T>::rotations = {
  std::make_pair(
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
    4
  ),
  std::make_pair(
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
    4
  ),
  std::make_pair(
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
    4
  )
};

} // eo namespace PermSymmetry
