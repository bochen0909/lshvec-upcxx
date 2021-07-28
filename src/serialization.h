#ifndef _SRC_serialization_H_
#define _SRC_serialization_H_

#include <string>

template <typename T>
class SingleNodeModel;

void save_model(uint32_t this_epoch, const std::string &output_prefix, SingleNodeModel<float> &model, bool zip_output);

void save_vector(uint32_t this_epoch, const std::string &output_prefix, SingleNodeModel<float> &model, bool zip_output);

void save_vector_bin(uint32_t this_epoch, const std::string &output_prefix, SingleNodeModel<float> &model, bool zip_output);

#endif //_SRC_serialization_H_
