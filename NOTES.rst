Issues
------

- asin is discontinuous and inaccurate -> perhaps an implementation error! 
  Splitting into bounds where different approximations are used makes it
  unusable for integration. Must make sure each function is a single
  approximation, not multiple on different bounds


Incomplete
----------

- Rename Containers.h to Algorithms.h?
- Various TODOs littered throughout the code


Future
------

- Remove Set.h?
- More vector operations?
- Add more trig, inverse trig
- UpperTriangularMatrix could be way more useful -> add functionality
- multiple algorithms for every Math primitive -> investigate which are best in
  terms of numerical stability and speed
