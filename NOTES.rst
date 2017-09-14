Issues
------

- asin is discontinuous and inaccurate -> perhaps an implementation error! 
  Splitting into bounds where different approximations are used makes it
  unusable for integration. Must make sure each function is a single
  approximation, not multiple on different bounds


Incomplete
----------

- Finish Map
- Various TODOs littered throughout the code


Future
------

- current constructors of UpperTriangular force you to specify the size of the
  resulting matrix, but it could be inferred from the size of the array used for
  construction (and checked for correctness!)
- More vector operations?
- Add more trig, inverse trig
- UpperTriangularMatrix could be way more useful -> add functionality
- multiple algorithms for every Math primitive -> investigate which are best in
  terms of numerical stability and speed
