#include <Core/MKL/Utils.h>
#include <gtest/gtest.h>

TEST(TestMKLCore, TestConvolution) {
  const fft::RealContainer x = {1, 2, 3, 4};
  const fft::RealContainer y = {11, 12, 13, 14, 15, 16, 17, 18};
  const auto result = fft::conv_real(x, y);
  const fft::RealContainer expected_result = {11,  34,  70,  120, 130, 140,
                                              150, 160, 151, 122, 72};

  EXPECT_EQ(x.size() + y.size() - 1u, expected_result.size());
  EXPECT_EQ(result, expected_result);
}

TEST(TestMKLCore, TestConvolution2) {
  const fft::RealContainer x = {1, 3, 2};
  const fft::RealContainer y = {4, 3, 2, 1};
  const auto result = fft::conv_real(x, y);
  const fft::RealContainer expected_result = {4, 15, 19, 13, 7, 2};

  EXPECT_EQ(x.size() + y.size() - 1u, expected_result.size());
  EXPECT_EQ(result, expected_result);
}

TEST(TestMKLCore, TestCorrelation) {
  const fft::RealContainer x = {1, 2, 3, 4};
  const fft::RealContainer y = {11, 12, 13, 14, 15, 16, 17, 18};
  const auto result = fft::corr_real(x, y);
  const fft::RealContainer expected_result = {44,  81,  110, 130, 140, 150,
                                              160, 170, 104, 53,  18};

  EXPECT_EQ(result, expected_result);
}
