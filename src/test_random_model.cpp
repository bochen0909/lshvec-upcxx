// Let Catch provide main():
#define CATCH_CONFIG_MAIN

#include <iostream>
#include "catch.hpp"
#include "model.h"
using namespace std;

template <typename VALUE_TYPE>
class RandomModel : public Model<VALUE_TYPE>
{
public:
	RandomModel(uint32_t dim, uint32_t neg_size, bool use_cbow) : Model<VALUE_TYPE>(dim, neg_size, use_cbow, 5)
	{
	}

	virtual void randomize_init() {}

protected:
	virtual void incr_word_count(uint32_t) {}
	virtual void add_vector_to_wi(const Vector<VALUE_TYPE> &, uint32_t, float)
	{
	}
	virtual void add_vector_to_wo(const Vector<VALUE_TYPE> &, uint32_t, float)
	{
	}
	virtual uint32_t getNegative(uint32_t)
	{
		return 1;
	}
	virtual Vector<VALUE_TYPE> get_vec_from_wi(uint32_t)
	{
		return Vector<VALUE_TYPE>(this->dim);
	}
	virtual Vector<VALUE_TYPE> get_vec_from_wo(uint32_t)
	{
		return Vector<VALUE_TYPE>(this->dim);
	}
};

TEST_CASE("random_model ", "[model]")
{

	RandomModel<float> model(10, 5, false);
	model.randomize_init();
	float loss = model.update({1, 2, 3, 4, 5}, 0.1f, true);
	std::cout << "loss=" << loss << "\n";
}

TEST_CASE("random_model_cbow ", "[model]")
{

	RandomModel<float> model(10, 5, true);
	model.randomize_init();
	float loss = model.update({1, 2, 3, 4, 5}, 0.1f, true);
	std::cout << "loss=" << loss << "\n";
}

TEST_CASE("single_node_model ", "[model]")
{

	SingleNodeModel<float> model(200, 10, 5, false, 5);
	model.randomize_init();
	float loss = model.update({1, 2, 3, 4, 5}, 0.1f, true);
	std::cout << "loss=" << loss << "\n";
}

TEST_CASE("single_node_model_cbow ", "[model]")
{

	SingleNodeModel<float> model(200, 10, 5, true, 5);
	model.randomize_init();
	float loss = model.update({1, 2, 3, 4, 5}, 0.1f, true);
	std::cout << "loss=" << loss << "\n";
}
