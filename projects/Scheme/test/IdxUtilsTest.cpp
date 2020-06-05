#include <Scheme/src/IdxUtils.h>
#include <gtest/gtest.h>

#include <array>


TEST(TestIdxUtils, TestConvolution) {
  const Grid1D x = {1, 2, 3, 4};
  const Grid1D y = {11, 12, 13, 14, 15, 16, 17, 18};
  const auto result = utils::conv_real<false>(x, y);
  const Grid1D expected_result = {11,  34,  70,  120, 130, 140,
                                              150, 160, 151, 122, 72};

  EXPECT_EQ(x.size() + y.size() - 1u, expected_result.size());
  EXPECT_EQ(result, expected_result);
}

TEST(TestIdxUtils, TestConvolution2) {
  const Grid1D x = {1, 3, 2};
  const Grid1D y = {4, 3, 2, 1};
  const auto result = utils::conv_real<false>(x, y);
  const Grid1D expected_result = {4, 15, 19, 13, 7, 2};

  EXPECT_EQ(x.size() + y.size() - 1u, expected_result.size());
  EXPECT_EQ(result, expected_result);
}

TEST(TestIdxUtils, TestCorrelation) {
  const Grid1D x = {1, 2, 3, 4};
  const Grid1D y = {11, 12, 13, 14, 15, 16, 17, 18};
  const auto result = utils::corr_real<false>(x, y);
  const Grid1D expected_result = {44,  81,  110, 130, 140, 150,
                                              160, 170, 104, 53,  18};

  EXPECT_EQ(result, expected_result);
}

TEST(TestIdxUtils, TestApplyOperationSumVectorized) {
  const Grid1D x = {1., 2., 3., 4.};
  const Grid1D y = {1., 2., 3., 4.};
  ASSERT_EQ(x.size(), y.size());
  const auto result = utils::apply_operation<true>(
      x, y, [](const auto& x, const auto& y) { return x + y; });

  EXPECT_EQ(x.size(), result.size());
  EXPECT_EQ(result, (Grid1D{2., 4., 6., 8.}));
}

TEST(TestIdxUtils, TestApplyOperationSum) {
  const Grid1D x = {1., 2., 3., 4.};
  const Grid1D y = {1., 2., 3., 4.};
  ASSERT_EQ(x.size(), y.size());
  const auto result = utils::apply_operation<false>(
      x, y, [](const auto& x, const auto& y) { return x + y; });

  EXPECT_EQ(x.size(), result.size());
  EXPECT_EQ(result, (Grid1D{2., 4., 6., 8.}));
}

TEST(TestIdxUtils, TestApplyOperationSumX2Vectorized) {
  const Grid1D x = {1., 2., 3., 4.};
  const Grid1D y = {1., 2., 3., 4.};
  ASSERT_EQ(x.size(), y.size());
  const auto result = utils::apply_operation<true>(
      x, y, [](const auto& x, const auto& y) { return x + 2. * y; });

  EXPECT_EQ(x.size(), result.size());
  EXPECT_EQ(result, (Grid1D{3., 6., 9., 12.}));
}

TEST(TestIdxUtils, TestApplyOperationSumX2) {
  const Grid1D x = {1., 2., 3., 4.};
  const Grid1D y = {1., 2., 3., 4.};
  ASSERT_EQ(x.size(), y.size());
  const auto result = utils::apply_operation<false>(
      x, y, [](const auto& x, const auto& y) { return x + 2. * y; });

  EXPECT_EQ(x.size(), result.size());
  EXPECT_EQ(result, (Grid1D{3., 6., 9., 12.}));
}

TEST(TestIdxUtils, TestApplyOperationDerivativeVectorized) {
  const Grid1D x = {1., 4., 6., 8.};
  const Grid1D y = {1., 2., 3., 4.};
  constexpr auto h = 2.;
  ASSERT_EQ(x.size(), y.size());
  const auto result = utils::apply_operation<true>(
      x, y, [](const auto& x, const auto& y) { return (x - y) / h; });

  EXPECT_EQ(x.size(), result.size());
  EXPECT_EQ(result, (Grid1D{0., 1., 1.5, 2.}));
}

TEST(TestIdxUtils, TestApplyOperationDerivative) {
  const Grid1D x = {1., 4., 6., 8.};
  const Grid1D y = {1., 2., 3., 4.};
  constexpr auto h = 2.;
  ASSERT_EQ(x.size(), y.size());
  const auto result = utils::apply_operation<false>(
      x, y, [](const auto& x, const auto& y) { return (x - y) / h; });

  EXPECT_EQ(x.size(), result.size());
  EXPECT_EQ(result, (Grid1D{0., 1., 1.5, 2.}));
}

TEST(TestIdxUtils, TestApplyConvolutionFactorVectorized) {
  const Grid1D x = {1., 2., 3., 4.};
  const auto result = utils::apply_conv_factor<true>(x, 1.);

  EXPECT_EQ(x.size(), result.size());
  EXPECT_EQ(result, (Grid1D{0., 2., 6., 12.}));
}

TEST(TestIdxUtils, TestApplyConvolutionFactor) {
  const Grid1D x = {1., 2., 3., 4.};
  const auto result = utils::apply_conv_factor<false>(x, 1.);

  EXPECT_EQ(x.size(), result.size());
  EXPECT_EQ(result, (Grid1D{0., 2., 6., 12.}));
}

TEST(TestIdxUtils, TestApplyCorrelationFactorVectorized) {
  const Grid1D x = {1., 2., 3., 4.};
  const auto result = utils::apply_corr_factor<true>(x, 1.);

  EXPECT_EQ(x.size(), result.size());
  EXPECT_EQ(result, (Grid1D{4., 6., 6., 4.}));
}

TEST(TestIdxUtils, TestApplyCorrelationFactor) {
  const Grid1D x = {1., 2., 3., 4.};
  const auto result = utils::apply_corr_factor<false>(x, 1.);

  EXPECT_EQ(x.size(), result.size());
  EXPECT_EQ(result, (Grid1D{4., 6., 6., 4.}));
}

TEST(TestIdxUtils, TestSummVectorized) {
  const Grid1D x = {11, 34, 70, 120, 130, 140, 150, 160, 151, 122, 72};
  const Grid1D y = {44, 81, 110, 130, 140, 150, 160, 170, 104, 53, 18};

  ASSERT_EQ(x.size(), y.size());

  const auto result = utils::summ_real<true>(x, y);
  const Grid1D expected_result = {55,  115, 180, 250, 270, 290,
                                  310, 330, 255, 175, 90};

  EXPECT_EQ(result, expected_result);
}

TEST(TestIdxUtils, TestSumm) {
  const Grid1D x = {11, 34, 70, 120, 130, 140, 150, 160, 151, 122, 72};
  const Grid1D y = {44, 81, 110, 130, 140, 150, 160, 170, 104, 53, 18};

  ASSERT_EQ(x.size(), y.size());

  const auto result = utils::summ_real<false>(x, y);
  const Grid1D expected_result = {55,  115, 180, 250, 270, 290,
                                  310, 330, 255, 175, 90};

  EXPECT_EQ(result, expected_result);
}
