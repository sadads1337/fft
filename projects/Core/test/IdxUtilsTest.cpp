#include <gtest/gtest.h>

#include <Core/IdxUtils.h>

#include <array>

TEST(TestIdxUtils, TestMiddleIdx)
{
	const auto idx = utils::get_idx(0., 1., 0.5, 10);
	EXPECT_EQ(idx, 5);
}

TEST(TestIdxUtils, TestApplyOperationSum)
{
	const Grid1D x = {1., 2., 3., 4.};
	const Grid1D y = {1., 2., 3., 4.};
	ASSERT_EQ(x.size(), y.size());
	const auto result = utils::apply_operation(x, y, [](const auto & x, const auto & y) { return x + y; });

	EXPECT_EQ(x.size(), result.size());
	EXPECT_EQ(result, (Grid1D{2., 4., 6., 8.}));
}

TEST(TestIdxUtils, TestApplyOperationSumX2)
{
	const Grid1D x = {1., 2., 3., 4.};
	const Grid1D y = {1., 2., 3., 4.};
	ASSERT_EQ(x.size(), y.size());
	const auto result = utils::apply_operation(x, y, [](const auto & x, const auto & y) { return x + 2. * y; });

	EXPECT_EQ(x.size(), result.size());
	EXPECT_EQ(result, (Grid1D{3., 6., 9., 12.}));
}

TEST(TestIdxUtils, TestApplyOperationDerivative)
{
	const Grid1D x = {1., 4., 6., 8.};
	const Grid1D y = {1., 2., 3., 4.};
	constexpr auto h = 2.;
	ASSERT_EQ(x.size(), y.size());
	const auto result = utils::apply_operation(x, y, [](const auto & x, const auto & y) { return (x - y) / h; });

	EXPECT_EQ(x.size(), result.size());
	EXPECT_EQ(result, (Grid1D{0., 1., 1.5, 2.}));
}

TEST(TestIdxUtils, TestApplyConvolutionFactor)
{
	const Grid1D x = {1., 2., 3., 4.};
	const auto result = utils::apply_conv_factor(x);

	EXPECT_EQ(x.size(), result.size());
	EXPECT_EQ(result, (Grid1D{0., 2., 6., 12.}));
}

TEST(TestIdxUtils, TestApplyCorrelationFactor)
{
	const Grid1D x = {1., 2., 3., 4.};
	const auto result = utils::apply_corr_factor(x);

	EXPECT_EQ(x.size(), result.size());
	EXPECT_EQ(result, (Grid1D{0., 2., 6., 12.}));
}
