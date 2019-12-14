#include <gtest/gtest.h>

#include <Core/IdxUtils.h>

TEST(TestIdxUtils, TestMiddleIdx)
{
	const auto idx = utils::get_idx(0., 1., 0.5, 10);
	EXPECT_EQ(idx, 5);
}
