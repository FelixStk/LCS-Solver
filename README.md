# LCS Solver

**LCS Solver** is a C++ library for solving various types of longest common subsequence (LCS) problems.

## Implemented Algorithms

The following table provides an overview of the implemented algorithms. In general:
- LCS algorithms return a `BaseCollector` Solution, which provides a type-erased iterator over the (one-based) longest common subsequence embeddings (`Points`)
- LLCS algorithms return a `UnsignedSolution`, which wrap the length of the longest common subsequence.

| Algorithm          | Description                                                                           | Method            | Time Complexity   |
|--------------------|---------------------------------------------------------------------------------------|-------------------|-------------------|
| `LLCS_STD_FL`      | The Classic LCS Problem                                                               | Folklore & DP     | O(N)              |
| `LLCS2_MC`         | Multiple Gap Constraint [1]                                                           | DP Matrices       | O(N * L)          |
| `LLCS2_MC_INC`     | Multiple Gap constraints where each is a superset of the previous [1]                 | Segment Tree & DP | O(N * log²(N))    |
| `LLCS2_MC_INC_E`   | Multiple Gap constraints where each is a superset of the previous [2]                 | Queues & DP       | O(N)              |
| `LLCS2_MC_O1_SYNC` | Multiple Gap Constraints where the number of gap bounds is known and synchronized [1] | 2D Buffer & DP    | O(N * h)          |
| `LLCS2_SR_MQ`      | Gap lengths depend on the symbol to the right of the gap [1]                          | 2D Buffer & DP    | O(N * σ)          |
| `LLCS2_SR_RMQ`     | Gap lengths depend on the symbol to the right of the gap [1]                          | RMQ & DP          | O(N * log(m))     |
| `LLCS2_SL_R`       | Gap lengths depend on the symbol to the left of the gap [1]                           | Wrapper/Mirror    | Same as SR        |
| `LLCS2_SA_MQ`      | Gap lengths depend on the surrounding symbols of the gap [1]                          | 2D Buffer & DP    | O(N * σ²)         |
| `LLCS2_SA_RMQ`     | Gap lengths depend on the surrounding symbols of the gap [1]                          | RMQ & DP          | O(N * σ * log(m)) |
| `LCS2_RT`          | LCS-Reconstruction based on marked points in LLCS2 Algorithm                          | L 2D-Rangetree    | O(N * log(N) * L) |
| `LCS2_STD_S`       | LCS-Reconstruction based on a stack for the Classic LCS Problem                       | Stack & DP        | O(N)              |
| `LCS2_STD_H`       | LCS-Reconstruction based on the hirschberg trick for the Classic LCS Problem          | DP                | O(N)              |

In this tabel for the time until the first solution is generated, the problem-dependent variables are:
- `N` is the product of the lengths of the two input strings,
- `L` is the length of the longest common subsequence,
- `h` is the number of distinct gap constraint in a gc tuple,
- `σ` is the size of the input alphabet, and
- `m` is the length of the longest input string.
- `t` is the number of (embedding-wise) different LCS solutions.

## Compiling the Project

1. **Requirements**. By default, CMake will be configured to automatically fetch the latest version of [Google Test](https://github.com/google/googletest/), [Google Benchmark](https://github.com/google/benchmark) and a [JSON library](https://github.com/nlohmann/json)from GitHub if no version is installed. You can disable this fetch behavior by setting the options beginning with FETCH to false. Ensure you have the following requirements:
   - C++23 compiler
   - [CMake](https://cmake.org/download/) (version >= 3.27)
   - [Doxygen](https://www.doxygen.nl/download.html) (version >= 1.9.7)
   - [Graphviz](https://graphviz.org/download/) (optional, only needed by Doxygen for generating class diagrams)

2. **Install missing packages**. If you are using a 64-bit Windows, you can install [MSYS2](https://www.msys2.org/) and execute the following command to meet the requirements:
    ```shell
    pacman -S mingw-w64-x86_64-doxygen mingw-w64-x86_64-cmake mingw-w64-x86_64-graphviz mingw-w64-x86_64-gtest mingw-w64-x86_64-benchmark mingw-w64-x86_64-doxygen mingw-w64-x86_64-nlohmann-json
    ```
    For Debian or Ubuntu the equivalent command would be:
    ```shell
    sudo apt install build-essential cmake doxygen graphviz libgtest-dev libbenchmark-dev nlohmann-json3-dev libgmock-dev
    ```
3. **Initial build**. To build the project in the `./cmake-build-debug/` directory, navigate to the main project directory and run:
    ```shell
    cmake -B cmake-build-debug -S .
    cmake --build cmake-build-debug/
    ```
   To build the project as release use the `-DCMAKE_BUILD_TYPE=Release` flag on the first cmake configure command. You can generate the documentation by executing either:
    ```shell
    doxygen ./docs/Doxyfile
    ```
   or
    ```shell
    cmake -B cmake-build-debug -S .
    cmake --build cmake-build-debug --target doxygen-docs
    ```
    Afterward, the documentation will be in [docs/html/index.html](docs/html/index.html)

## Demonstration Projects
- [`demo/solver`](demo/solver/README.md): Generation, processing and solving of problems defined in `lib/problems`
- [`demo/dna`](demo/dna/README.md): Run an algorithm on the [afproject genetree dataset](https://afproject.org/media/genetree/swisstree/dataset/swisstree.zip)
- [`demo/diffgc`](demo/diffgc/README.md): Run a lcs algorithm to generate the difference between two files. 

## Testing
You can run the unit tests using the GTest executable `LCS_Test` in the tests folder, or by using `ctest`:

- **Option 1**: Run the GTest executable:
   ```shell 
   cd cmake-build-debug/tests # necessary because of relative pathing in integration tests
   ./LCS_Test
   ```
  For more details, refer to the [GoogleTest documentation](https://google.github.io/googletest/advanced.html#running-test-programs-advanced-options), or call the test application in the `./build/src` directory with the `--help` argument.

- **Option 2**: Run `ctest` in the project's build directory:
   ```shell 
   cd cmake-build-debug
   ctest
   ```

- **Option 3**: Use an IDE to select and run specific tests.

A test report can be found in [`docs/Test_Results`](docs/Test_Results.htm)

## Installation
Currently, the lcs_solver lib depends on nlohmann_json. So this library will be also installed if it is not available. I strongly recommend installing it separately because it leads to a cleaner and quicker installation. Afterward install `lcs_solver` by running the following command:

```shell
  sudo cmake --install build/
```

To uninstall `lcs_solver`, use the information in `cmake-build-debug/install_manifest.txt`. For example, run the following command:

```shell
    sudo xargs rm -v < build/install_manifest.txt
    xargs -L1 dirname < build/install_manifest.txt | sort -r | sudo xargs -r rmdir -v 2>/dev/null
```

## References

- [1] D. Adamson, M. Kosche, T. Koß, F. Manea, and S. Siemer, "Longest common subsequence with gap constraints," 2023. https://arxiv.org/abs/2304.05270
- [2] D. Adamson, P. Sarnighausen-Cahn, M. Dumitran, M. Kosche, T. Koß, F. Manea, and S. Siemer, "Longest Common Subsequence with Gap Constraints,” Theory of Computing Systems, vol. 69,
  no. 2, p. 25, Jun. 2025. https://doi.org/10.1007/s00224-025-10223-0