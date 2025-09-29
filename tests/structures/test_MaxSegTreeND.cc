/*******************************************************************************
 * @file test_MaxSegTreeND.cc
 * @author Felix Steinkopp
 * @version 1.1
 * @brief Test the correctness of the MaxSegTreeND class
 ******************************************************************************/

#include <algorithm>
#include <random>
#include <string>
#include <vector>

#include "gtest/gtest.h"
#include "gmock/gmock-more-matchers.h"

#include "structures/Matrix.h"
#include "structures/MaxSegTreeND.h"
#include "util/CommonTypes.h"

namespace {

using ::lcs_solver::structures::MaxSegTreeND;

using ::testing::Lt;
using ::testing::Le;
using ::testing::Eq;
using ::testing::ContainerEq;
using ::testing::PrintToString;
using ::testing::ValuesIn;

using Pair = std::pair<size_t, size_t>;
using Index = std::vector<size_t>;
using Ranges = std::vector<Pair>;
using Matrix = lcs_solver::structures::Matrix<lcs_solver::util::uint>;

//== Parameter structure for the tests =========================================
using TestParam = std::tuple<
    std::string,                            // 0: name
    Matrix,                                 // 1: data
    std::vector<Index>,                     // 2: pointQueries
    std::vector<Ranges>,                    // 3: rangeQueries
    std::vector<std::tuple<Index, size_t>>, // 4: pointUpdates
    std::vector<std::tuple<Ranges, size_t>> // 5: rangeUpdates
>;

//== Definition of Testparameter ===============================================
std::vector<TestParam> getHardCodedParams() {
  std::vector<TestParam> params = {
      {
          "EmptyInitialValues",  // name of parameter struct
          {},  // Vector2D for Initialization
          {},  // Simple Updates (x, y, val}
          {},  // Ranged Updates (pair, pair, val)
          {},  // Simple Queries  (x,y )
          {}  // Ranged Queries (x-range, y-range)
      }
  };
  return params;
}

/*******************************************************************************
 * GenTestParams
 * @param initSeed unsigned int used for random number generation
 * @param nParams size_t for the number of parameters to generate
 * @param nDim std::vector for describing the dimensionality of the data matrix
 * @param nPointUpdates size_t for the number of point updates
 * @param nRangeUpdates size_t for the number of range updates
 * @param nPointQueries size_t for the number of point queries
 * @param nRangeQueries size_t for the number of range queries
 * @param Pair interval of possible values used data initialization
 * @return std::vector<TestParam>
 ******************************************************************************/
std::vector<TestParam> GenTestParams(const unsigned int initSeed,
                                     const size_t nParams = 5,
                                     const std::vector<size_t> nDim = {5, 10},
                                     const size_t nPointUpdates = 3,
                                     const size_t nRangeUpdates = 2,
                                     const size_t nPointQueries = 3,
                                     const size_t nRangeQueries = 3,
                                     const Pair numRange = {1, 100}) {
  std::vector<TestParam> params;
  for (size_t numParam = 0; numParam < nParams; ++numParam) {
    // Setup of generators, constants and empty TestParam x
    const size_t seed = initSeed + numParam;
    std::mt19937 gen(seed);
    std::uniform_int_distribution<size_t> dist(numRange.first, numRange.second);
    TestParam x;
    auto &[name, data, pointQueries, rangeQueries, pointUpdates, rangeUpdates] =
        x;

    // Generate Parameter: Name
    name = "Rnd_" + std::to_string(seed);

    // Generate Parameter: Data
    Matrix mat = Matrix(nDim);
    for (auto &val : mat.data)
      val = dist(gen);
    data = std::move(mat);

    // Generate Parameter: pointQueries
    auto genIdx = [&]() -> Index {
      std::vector<size_t> idx;
      idx.reserve(nDim.size());
      for (unsigned long long i : nDim)
        idx.emplace_back(dist(gen) % i);
      return idx;
    };
    std::generate_n(std::back_inserter(pointQueries), nPointQueries, genIdx);

    // Generate Parameter: rangeQueries
    auto genR = [&]() -> Ranges {
      Ranges r = Ranges();
      r.reserve(nDim.size());
      for (unsigned long long i : nDim) {
        size_t a = dist(gen) % i;
        size_t b = a + (dist(gen) % (i - a));
        r.emplace_back(a, b);
      }
      return r;
    };
    std::generate_n(std::back_inserter(rangeQueries), nRangeQueries, genR);

    // Generate Parameter: pointUpdates
    std::generate_n(std::back_inserter(pointUpdates), nPointUpdates,
                    [&]() -> std::tuple<Index, size_t> {
                      auto idx = genIdx();
                      auto val = dist(gen);
                      return {idx, val};
                    });

    // Generate Parameter: rangeUpdates
    std::generate_n(std::back_inserter(rangeUpdates), nRangeUpdates,
                    [&]() -> std::tuple<Ranges, size_t> {
                      auto r = genR();
                      auto val = dist(gen);
                      return {r, val};
                    });

    params.push_back(std::move(x));
  }
  return params;
}

/*******************************************************************************
 * Helper: applyRangeUpdate
 * @param mat Matrix to update
 * @param update tuple containing vector<Pair> for the range and a value
 ******************************************************************************/
void applyUpdate(Matrix& mat, const std::tuple<Ranges, size_t>& update){

  // Prepare looping through r[0] x r[1] x ... x r[r.size()-1]
  Index indices;
  auto &[r, val] = update;
  indices.reserve(r.size());
  for (size_t i = 0; i < mat.dim.size(); ++i)
    indices.emplace_back(r[i].first);
  bool error = mat.empty() || r.size() != mat.dim.size();

  // Apply the update to the matrix mat
  bool loop = true;
  while (!error && loop) {
    // Process the current combination
    for (size_t i = 0; i < mat.dim.size(); ++i) {
      error = error || r[i].first >= mat.dim[i];
      error = error || r[i].second >= mat.dim[i];
    }
    if (!error)
      mat[indices] = std::max<size_t>(mat[indices], val);

    // Increment the indices
    size_t dimension = r.size();
    while (dimension--) {
      if (++indices[dimension] <= r[dimension].second) {
        break;
      }
      if (dimension == 0) loop = false;
      indices[dimension] = r[dimension].first;
    }
  }
}


//==== Initialize Test Suites ==================================================

class MaxSegTreeNDTest : public testing::TestWithParam<TestParam> {};

const auto seed = std::random_device{}();
//const auto seed = 4025468157;
std::vector<TestParam> hardCodedParam = getHardCodedParams();
std::vector<TestParam> param1D = GenTestParams(seed, // initSeed
                                               6, // nParams
                                               {11}, // nDim
                                               4, // nPointUpdates
                                               20, // nRangeUpdates
                                               4, // nPointQueries
                                               4 //nRangeQueries
);
std::vector<TestParam> param2D = GenTestParams(seed, // initSeed
                                               5, // nParams
                                               {2, 11}, // nDim
                                               5, // nPointUpdates
                                               18, // nRangeUpdates
                                               5, // nPointQueries
                                               5 //nRangeQueries
);

INSTANTIATE_TEST_SUITE_P(
    HardCodedParameters, MaxSegTreeNDTest,
    ValuesIn(hardCodedParam),
    [](const testing::TestParamInfo<MaxSegTreeNDTest::ParamType> &info) {
      return std::get<0>(info.param);
    });

INSTANTIATE_TEST_SUITE_P(
    RandomParameters1D, MaxSegTreeNDTest,
    ValuesIn(param1D),
    [](const testing::TestParamInfo<MaxSegTreeNDTest::ParamType> &info) {
      return std::get<0>(info.param);
    });

INSTANTIATE_TEST_SUITE_P(
    RandomParameters2D, MaxSegTreeNDTest,
    ValuesIn(param2D),
    [](const testing::TestParamInfo<MaxSegTreeNDTest::ParamType> &info) {
      return std::get<0>(info.param);
    });

//==== Parametrized Test =======================================================
TEST_P(MaxSegTreeNDTest, SizeWithDimConstructor) {
  auto data = std::get<1>(GetParam());
  MaxSegTreeND st(data.dim);
  size_t expected = data.empty() ? 0 : 1;
  for (auto n : data.dim)
    expected *= 2 * n;
  EXPECT_THAT(st.size(), Le(expected));
}

TEST_P(MaxSegTreeNDTest, SizeWithMatrixConstructor) {
  auto data = std::get<1>(GetParam());
  MaxSegTreeND st(data);
  size_t expected = data.empty() ? 0 : 1;
  for (auto n : data.dim)
    expected *= 2 * n;
  size_t out = st.size();
  EXPECT_THAT(out, Le(expected));
}

TEST_P(MaxSegTreeNDTest, BuildQuery) {
  auto mat = std::get<1>(GetParam());
  MaxSegTreeND st(mat);
  auto leaf = Matrix(mat.dim);
  if (!mat.empty())
    for (size_t i = 0; i < mat.data.size(); ++i) {
      auto idx = mat.vectorizeIndex(i);
      leaf[i] = st.query(idx);

    }
  EXPECT_THAT(leaf.data, ContainerEq(mat.data))
            << "Exp M (as matrix): " << mat.DebugString() << "\n"
            << "Leafs (as matrix): " << leaf.DebugString() << "\n";
}

TEST_P(MaxSegTreeNDTest, RangeQuery) {
  auto mat = std::get<1>(GetParam());
  auto queries = std::get<3>(GetParam());
  MaxSegTreeND st(mat);

  for (const auto &r : queries) {

    auto result = st.query(r);
    size_t expected = 0;
    bool error = mat.empty() || r.size() != mat.dim.size();

    // Initialize indices
    Index indices = Index(r.size(), 0);
    for (size_t i = 0; i < mat.dim.size(); ++i) {
      indices[i] = r[i].first;
    }

    bool loop = true;
    while (!error && loop) {
      // Process the current combination
      for (size_t i = 0; i < mat.dim.size(); ++i) {
        error = error || r[i].first >= mat.dim[i];
        error = error || r[i].second >= mat.dim[i];
      }
      expected = error ? 0 : std::max<size_t>(expected, mat[indices]);

      // Increment the indices like in a mixed radix system
      size_t dimension = r.size();
      while (dimension--) {
        if (++indices[dimension] <= r[dimension].second) {
          break;
        }
        if (dimension == 0) loop = false;
        indices[dimension] = r[dimension].first;
      }
    }

    EXPECT_THAT(result, Eq(error ? 0 : expected))
              << "Query: " << PrintToString(r) << "\n"
              << "Initial Data: " << mat.DebugString();
  }
}

TEST_P(MaxSegTreeNDTest, PointUpdateTest) {
  auto mat = std::get<1>(GetParam());
  auto updates = std::get<4>(GetParam());
  MaxSegTreeND st(mat);

  auto out = Matrix(mat.dim, 0);      // space to store simple query
  bool error = mat.empty();
  for (const auto &update : updates) {
    auto [updateIdx, val] = update;
    st.update(updateIdx, val);
    for (size_t i = 0; i < mat.dim.size(); ++i) {
      error = error || updateIdx.size() != mat.dim.size();
      error = error || updateIdx[i] >= mat.dim[i];
    }
    if (!error) mat[updateIdx] = std::max<size_t>(mat[updateIdx], val);
    for (size_t pos = 0; pos < mat.data.size(); ++pos) {
      auto idx = mat.vectorizeIndex(pos);
      out[pos] = st.query(idx);
    }

    EXPECT_THAT(out.data, ContainerEq(mat.data))
              << "Update: " << PrintToString(update) << "\n"
              << "Initial Data: " << std::get<1>(GetParam()).DebugString() << "\n"
              << "Mat (as matrix): " << mat.DebugString() << "\n"
              << "Leafs (as matrix): " << out.DebugString() << "\n";
  }
}

TEST_P(MaxSegTreeNDTest, RangeUpdateTest1) {
  auto &mat = std::get<1>(GetParam());
  auto &updates = std::get<5>(GetParam());

  MaxSegTreeND st(mat.dim); // Initialized with MaxSegTreeND::neutralElement

  auto expected = Matrix(mat.dim, 0);
  auto out = Matrix(mat.dim, 0);
  auto prev = Matrix(mat.dim, 0);
  for (const auto &update : updates) {
    auto &[r, val] = update;
    // Apply Updates
    prev = expected;
    st.update(r, val);
    applyUpdate(expected, update);

    // Fill matrix with leaf data
    for (size_t pos = 0; pos < mat.data.size(); ++pos) {
      auto idx = mat.vectorizeIndex(pos);
      out[pos] = st.query(idx);
    }

    EXPECT_THAT(out.data, ContainerEq(expected.data))
              << "RangeUpdates: " << PrintToString(updates) << "\n"
              << "Update" << PrintToString(update) << "\n"
              << "Initial: (" << PrintToString(mat.dim) << ") Zero Matrix\n"
              << "Exp M (as matrix): " << expected.DebugString() << "\n"
              << "Leafs (as matrix): " << out.DebugString() << "\n"
              << "Prev  (as matrix): " << prev.DebugString() << "\n";
  } // end of foreach updates loop
}

TEST_P(MaxSegTreeNDTest, RangeUpdateTest2) {
  auto &mat = std::get<1>(GetParam());
  auto &updates = std::get<5>(GetParam());
  MaxSegTreeND st(mat.dim); // Initialized with MaxSegTreeND::neutralElement

  auto expected = Matrix(mat.dim, 0);
  for (const auto &update : updates) {
    auto &[r, val] = update;
    st.update(r, val);
    applyUpdate(expected, update);
  }

  auto out = Matrix(mat.dim, 0);
  for (size_t pos = 0; pos < mat.data.size(); ++pos) {
    auto idx = mat.vectorizeIndex(pos);
    out[pos] = st.query(idx);
  }

  EXPECT_THAT(out.data, ContainerEq(expected.data))
            << "RangeUpdates: " << PrintToString(updates) << "\n"
            << "Initial: (" << PrintToString(mat.dim) << ") Zero Matrix\n"
            << "Exp M (as matrix): " << expected.DebugString() << "\n"
            << "Leafs (as matrix): " << out.DebugString() << "\n";
}

}  // end of anonymous namespace