/*
 * MODEL.h
 *
 *  Created on: Jul 13, 2020
 *      Author: Bo Chen
 */

#ifndef SOURCE_DIRECTORY__SRC_SPARC_MODEL_UPCXX_H_
#define SOURCE_DIRECTORY__SRC_SPARC_MODEL_UPCXX_H_

#include <vector>
#include <cmath>
#include <iomanip>
#include "model.h"

template <>
struct upcxx::serialization<Vector<float>>
{
    template <typename Writer>
    static void serialize(Writer &writer, Vector<float> const &object)
    {
        auto &val = object.get_val();
        size_t n = val.size();

        writer.write(n);
        for (auto &x : val)
        {
            writer.write(x);
        }
    }
    template <typename Reader>
    static Vector<float> *deserialize(Reader &reader, void *storage)
    {

        size_t n_val = reader.template read<size_t>();
        Vector<float> *v = new (storage) Vector<float>();
        for (size_t n = 0; n < n_val; n++)
        {
            (*v).push_back(reader.template read<float>());
        }
        return v;
    }
};

template <typename VALUE_TYPE>
class UPCXXNodeModel;

template <typename VALUE_TYPE>
using UPCXXNodeModel_ptr = std::shared_ptr<UPCXXNodeModel<VALUE_TYPE>>;

template <typename VALUE_TYPE>
class UPCXXModel : public Model<VALUE_TYPE>
{
protected:
    uint32_t total_num_word;
    uint32_t num_rank;
    uint32_t this_rank;

    using dobj_node_model_t = upcxx::dist_object<UPCXXNodeModel_ptr<VALUE_TYPE>>;
    dobj_node_model_t *local_model;
    UPCXXNodeModel_ptr<VALUE_TYPE> model;

    uint32_t *bucket_start = 0;
    uint32_t *bucket_end = 0;
    uint32_t bucket = 0;

public:
    UPCXXModel(uint32_t total_num_word, uint32_t this_rank, uint32_t num_rank, uint32_t dim, uint32_t neg_size, bool use_cbow, uint32_t half_window)
        : Model<VALUE_TYPE>(dim, neg_size, use_cbow, half_window), total_num_word(total_num_word), num_rank(num_rank), this_rank(this_rank), local_model(0)
    {
        bucket_start = new uint32_t[num_rank];
        bucket_end = new uint32_t[num_rank];
        if (0)
        {
            bucket = total_num_word / num_rank;
            if (true)
            {
                uint32_t remainder = total_num_word - bucket * num_rank;
                uint32_t word_start = 0;
                for (uint32_t i = 0; i < num_rank; i++, remainder--)
                {

                    uint32_t word_end = word_start + bucket;
                    if (remainder > 0)
                    {
                        word_end++;
                    }
                    if (word_end > total_num_word)
                    {
                        word_end = total_num_word;
                    }
                    assert(word_start < word_end);
                    bucket_start[i] = word_start;
                    bucket_end[i] = word_end;

                    word_start += bucket;
                    if (remainder > 0)
                    {
                        word_start++;
                    }
                }
            }
        }
        else
        {
            bucket = total_num_word / num_rank;
            uint32_t rem = total_num_word % num_rank;
            if (rem > num_rank / 2)
            {
                bucket++;
            }
            if (bucket < 1)
            {
                bucket = 1;
            }

            for (uint32_t i = 0; i < num_rank; i++)
            {
                bucket_start[i] = i * bucket;
                bucket_end[i] = (i + 1) * bucket;
                if (bucket_end[i] > total_num_word)
                {
                    bucket_end[i] = total_num_word;
                }
                assert(bucket_start[i] <= bucket_end[i]);
            }
            bucket_end[num_rank - 1] = total_num_word;
        }
        uint32_t word_start = bucket_start[this_rank];
        uint32_t word_end = bucket_end[this_rank];

        myinfo("create model: total_num_word=%u, word_start=%u, word_end=%u, dim=%u, neg_size=%u, use_cbow=%u, half_window=%u",
               total_num_word, word_start, word_end, dim, neg_size, use_cbow ? 1 : 0, half_window);

        model = std::make_shared<UPCXXNodeModel<VALUE_TYPE>>(total_num_word, word_start, word_end, dim, neg_size, use_cbow, half_window);
        local_model = new dobj_node_model_t(model);
    }
    virtual ~UPCXXModel()
    {
        if (local_model)
        {
            delete local_model;
        }
        if (bucket_end)
        {
            delete[] bucket_end;
        }
        if (bucket_start)
        {
            delete[] bucket_start;
        }
    }
    void randomize_init()
    {
        model->randomize_init();
    }

    void uniform_init()
    {
        model->uniform_init();
    }

    void dump_vector_bin(const std::string &prefix, bool zip_output)
    {
        std::string filepath = prefix;
        if (zip_output)
        {
            filepath += ".gz";
        }

        if (zip_output)
        {
            ogzstream output(filepath.c_str());
            bitsery::OutputStreamAdapter bw(output);
            model->write_vec_bin(bw);
            output.flush();
        }
        else
        {
            std::ofstream output(filepath, std::ios::binary | std::ios::trunc);
            bitsery::OutputStreamAdapter bw(output);
            model->write_vec_bin(bw);
            output.flush();
        }
    }

    void dump_vector_text(const std::string &prefix, bool zip_output)
    {
        std::string filepath = prefix;
        if (zip_output)
        {
            filepath += ".gz";
        }

        if (zip_output)
        {
            ogzstream output(filepath.c_str());
            model->write_vec(output);
            output.flush();
        }
        else
        {
            std::ofstream output(filepath);
            model->write_vec(output);
            output.flush();
        }
    }

    void dump_word_counts(const std::string &prefix, bool zip_output)
    {
        std::string filepath = prefix;
        if (zip_output)
        {
            filepath += ".gz";
        }

        if (zip_output)
        {
            ogzstream output(filepath.c_str());
            model->write_word_count(output);
            output.flush();
        }
        else
        {
            std::ofstream output(filepath);
            model->write_word_count(output);
            output.flush();
        }
    }

protected:
    uint32_t _find_rank_by_word(uint32_t word, uint32_t start, uint32_t end)
    {
        assert(end >= start);
        if (start == end || start + 1 == end)
        {

            if (word >= bucket_start[start] && word < bucket_end[start])
            {
                return start;
            }
            if (word >= bucket_start[end] && word < bucket_end[end])
            {
                return start;
            }
            throw std::runtime_error("should not be here");
        }

        uint32_t mid = (start + end) / 2;
        if (word >= bucket_start[mid])
        {
            return _find_rank_by_word(word, mid, end);
        }
        else
        {
            return _find_rank_by_word(word, start, mid - 1);
        }
    }

    uint32_t get_target_rank(uint32_t word)
    {
        if (0)
        {
            return _find_rank_by_word(word, 0, num_rank);
        }
        else
        {
            return word / bucket;
        }
    }

    virtual void incr_word_count(uint32_t word)
    {
        upcxx_incr_word_count(word);
    }

    upcxx::future<> upcxx_incr_word_count(uint32_t word)
    {
        return upcxx::rpc(
            get_target_rank(word),
            [](dobj_node_model_t &lmodel, uint32_t word)
            {
                (*lmodel)->incr_word_count(word);
            },
            *local_model, word);
    }

    virtual void add_vector_to_wi(const Vector<VALUE_TYPE> &vec, uint32_t word, float a)
    {
        upcxx_add_vector_to_wi(vec, word, a);
    }

    upcxx::future<> upcxx_add_vector_to_wi(const Vector<VALUE_TYPE> &vec, uint32_t word, float a)
    {
        return upcxx::rpc(
            get_target_rank(word),
            [](dobj_node_model_t &lmodel, const Vector<VALUE_TYPE> &vec, uint32_t word, float a)
            {
                (*lmodel)->add_vector_to_wi(vec, word, a);
            },
            *local_model, vec, word, a);
    }

    virtual void add_vector_to_wo(const Vector<VALUE_TYPE> &vec, uint32_t word, float a)
    {
        upcxx_add_vector_to_wo(vec, word, a);
    }

    upcxx::future<> upcxx_add_vector_to_wo(const Vector<VALUE_TYPE> &vec, uint32_t word, float a)
    {
        return upcxx::rpc(
            get_target_rank(word),
            [](dobj_node_model_t &lmodel, const Vector<VALUE_TYPE> &vec, uint32_t word, float a)
            {
                (*lmodel)->add_vector_to_wo(vec, word, a);
            },
            *local_model, vec, word, a);
    }

    virtual uint32_t getNegative(uint32_t target)
    {
        while (true)
        {
            uint32_t i = (uint32_t)(sparc::myrand::uniform<double>() * total_num_word);
            if (i != target)
            {
                return i;
            }
        }
    }
    upcxx::future<Vector<VALUE_TYPE>> upcxx_get_vec_from_wi(uint32_t word)
    {
        return upcxx::rpc(
            get_target_rank(word),
            [](dobj_node_model_t &lmodel, uint32_t word)
            {
                return (*lmodel)->get_vec_from_wi(word);
            },
            *local_model, word);
    }
    upcxx::future<Vector<VALUE_TYPE>> upcxx_get_vec_from_wo(uint32_t word)
    {
        return upcxx::rpc(
            get_target_rank(word),
            [](dobj_node_model_t &lmodel, uint32_t word)
            {
                return (*lmodel)->get_vec_from_wo(word);
            },
            *local_model, word);
    }

    Vector<VALUE_TYPE> get_vec_from_wo(uint32_t word)
    {
        return upcxx_get_vec_from_wo(word).wait();
    }

    Vector<VALUE_TYPE> get_vec_from_wi(uint32_t word)
    {
        return upcxx_get_vec_from_wi(word).wait();
    }
};

template <typename VALUE_TYPE>
class UPCXXNodeModel : public Model<VALUE_TYPE>
{
protected:
    uint32_t total_num_word;
    uint32_t word_start; //include
    uint32_t word_end;   //exclude
    uint32_t num_word;
    std::vector<Vector<VALUE_TYPE>> wo_;
    std::vector<Vector<VALUE_TYPE>> wi_;
    std::vector<uint32_t> word_count;

public:
    void read_me(bitsery::InputStreamAdapter &br)
    {
        Model<VALUE_TYPE>::read_me(br);
        read(br, total_num_word);
        read(br, word_start);
        read(br, word_end);
        read(br, wo_);
        read(br, wi_);
        read(br, word_count);
    }

    void write_me(bitsery::OutputStreamAdapter &bw) const
    {
        Model<VALUE_TYPE>::write_me(bw);
        write(bw, total_num_word);
        write(bw, word_start);
        write(bw, word_end);
        write(bw, wo_);
        write(bw, wi_);
        write(bw, word_count);
    }

    template <typename OS>
    void write_vec(OS &output)
    {
        for (auto &x : wi_)
        {
            x.write_me(output);
        }
    }

    template <typename OS>
    void write_word_count(OS &output)
    {
        for (size_t i = 0; i < word_count.size(); i++)
        {
            output << i + word_start << " " << word_count.at(i) << "\n";
        }
    }

    void write_vec_bin(bitsery::OutputStreamAdapter &bw)
    {
        write(bw, wi_);
    }

public:
    UPCXXNodeModel() : Model<VALUE_TYPE>()
    {
    }
    UPCXXNodeModel(uint32_t total_num_word, uint32_t word_start, uint32_t word_end, uint32_t dim, uint32_t neg_size, bool use_cbow, uint32_t half_window)
        : Model<VALUE_TYPE>(dim, neg_size, use_cbow, half_window), total_num_word(total_num_word), word_start(word_start), word_end(word_end)
    {
        num_word = word_end - word_start;
        for (uint32_t i = 0; i < num_word; i++)
        {
            wo_.push_back(Vector<VALUE_TYPE>(dim));
            wi_.push_back(Vector<VALUE_TYPE>(dim));
        }
        word_count.resize(num_word);
    }

    void randomize_init()
    {
        for (uint32_t i = 0; i < num_word; i++)
        {
            wo_[i].randomize();
            wi_[i].randomize();
        }
    }

    void uniform_init()
    {
        for (uint32_t i = 0; i < num_word; i++)
        {
            wo_[i].uniform();
            wi_[i].uniform();
        }
    }

public:
    virtual void incr_word_count(uint32_t word)
    {
        word_count[word - word_start]++;
    }
    virtual void add_vector_to_wi(const Vector<VALUE_TYPE> &vec, uint32_t word, float a)
    {
        wi_[word - word_start].add(vec, a);
    }
    virtual void add_vector_to_wo(const Vector<VALUE_TYPE> &vec, uint32_t word, float a)
    {
        wo_[word - word_start].add(vec, a);
    }
    virtual uint32_t getNegative(uint32_t target)
    {
        while (true)
        {
            uint32_t i = (uint32_t)(sparc::myrand::uniform<double>() * total_num_word);
            if (i != target)
            {
                return i;
            }
        }
    }
    virtual Vector<VALUE_TYPE> get_vec_from_wi(uint32_t word)
    {
        return wi_.at(word - word_start);
    }
    virtual Vector<VALUE_TYPE> get_vec_from_wo(uint32_t word)
    {
        return wo_.at(word - word_start);
    }
};

template <typename VALUE_TYPE>
void write(bitsery::OutputStreamAdapter &bw, const UPCXXNodeModel<VALUE_TYPE> &obj)
{
    obj.write_me(bw);
}

#endif /* SOURCE_DIRECTORY__SRC_SPARC_MODEL_UPCXX_H_ */
