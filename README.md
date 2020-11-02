# Preon::Math

Preon::Math is a simple C++ math library that is used for all the basic math operations in [PreonLab](https://www.fifty2.eu/preonlab/). It includes vector, matrix, euler angles and quaternion classes and basic linear algebra operations. It also includes SIMD variants for the most common operations.

Preon::Math is tested to compile with Visual Studio 2015 and GCC 7.

## Changes
### 0.03 - 11.02.2020
- Added `SymmetricMatrix33`.
- Added `.clang-format` file and formatted all files with it.

### 0.02 - 02.20.2020
- Renamed `vec::zeroVector()` to `vec::zero()`
- Improved `double`-precision support (new cast methods in `Preon::Math::cast`, new templated operators).

## Roadmap
* [ ] Port remaining unit tests from QTest to catch2.
* [ ] Doxygen documentation.
* [ ] CI builds.
* [ ] Extend unit tests to cover all classes and methods.
* [x] Rename `vec::zeroVector()` to `vec::zero()` to match the `matrix::zero()` method.
* [ ] Add a `all.h` (or similar name) that includes all math classes/defines. Similar, `all_fwd.h` would be nice.
* [ ] Combine `math_misc.h` and `math_utils.h`.
* [ ] Compare `MatrixUtil::fromQuaternion` to `RotationUtils::matrixFromQuaternions` and decide on one variant. Delete the other.
* [ ] Remove `ErrorHandling.h` dependencies and exception throwing.

## License
Preon::Math is licensed under [MPL2](https://www.mozilla.org/en-US/MPL/2.0/), see `LICENSE`.

This means you can use Preon::Math in your application (even commercially) but you must include the license notice and at least provide a link to the source code. Furthermore, if you make adaptions to Preon::Math you have to share them under the same license.

This repo also contains the [catch2](https://github.com/catchorg/Catch2) unit test library which is licensed under the Boost Software License.
