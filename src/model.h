/*
 * MODEL.h
 *
 *  Created on: Jul 13, 2020
 *      Author: Bo Chen
 */

#ifndef SOURCE_DIRECTORY__SRC_SPARC_MODEL_H_
#define SOURCE_DIRECTORY__SRC_SPARC_MODEL_H_

#include <vector>
#include <cmath>
#include <iomanip>
#include "static_block.hpp"
#include "utils.h"
#include "io.h"
template <typename VALUE_TYPE>
class Vector
{
    std::vector<VALUE_TYPE> val;

public:
    size_t get_dim()
    {
        return val.size();
    }
    void read_me(bitsery::InputStreamAdapter &br)
    {
        read(br, val);
    }

    void write_me(bitsery::OutputStreamAdapter &bw) const
    {
        write(bw, val);
    }
    template <typename OS>
    void write_me(OS &output)
    {
        for (size_t i = 0; i < val.size() - 1; i++)
        {
            output << std::fixed << std::setprecision(6) << val[i] << " ";
        }
        output << std::fixed << std::setprecision(6) << val[val.size() - 1] << "\n";
    }

    Vector()
    {
    }

    Vector(uint32_t dim)
    {
        val.resize(dim);
    }

    void randomize()
    {

        for (size_t i = 0; i < val.size(); i++)
        {
            val[i] = sparc::myrand::uniform<VALUE_TYPE>() * 2 - 1;
        }
    }

    void uniform()
    {
        float x = 1.0 / val.size();

        for (size_t i = 0; i < val.size(); i++)
        {
            val[i] = x;
        }
    }

    VALUE_TYPE at(size_t i)
    {
        return val.at(i);
    }

    void zero()
    {
        for (size_t i = 0; i < val.size(); i++)
        {
            val[i] = 0;
        }
    }
    void add(const Vector<VALUE_TYPE> &other)
    {
        for (size_t i = 0; i < val.size(); i++)
        {
            val[i] += other.val.at(i);
        }
    }

    void add(const Vector<VALUE_TYPE> &other, float alpha)
    {
        for (size_t i = 0; i < val.size(); i++)
        {
            val[i] += other.val.at(i) * alpha;
        }
    }

    void mul(float x)
    {
        for (size_t i = 0; i < val.size(); i++)
        {
            val[i] *= x;
        }
    }

    float dot(const Vector<VALUE_TYPE> &other) const
    {
        float r = 0;
        for (size_t i = 0; i < val.size(); i++)
        {
            r += val.at(i) * other.val.at(i);
        }
        return r;
    }
};

template <typename VALUE_TYPE>
inline void write(bitsery::OutputStreamAdapter &bw, const Vector<VALUE_TYPE> &val)
{
    val.write_me(bw);
}

template <typename VALUE_TYPE>
inline void read(bitsery::InputStreamAdapter &br, Vector<VALUE_TYPE> &val)
{
    val.read_me(br);
}

static int MAX_SIGMOID = 8;
static int SIGMOID_TABLE_SIZE = 512;
static int LOG_TABLE_SIZE = 512;
static float *t_sigmoid_;
static float *t_log_;

static_block
{

    if (true)
    {
        t_sigmoid_ = new float[SIGMOID_TABLE_SIZE + 1];
        for (int i = 0; i < SIGMOID_TABLE_SIZE + 1; i++)
        {
            float x = (1.0f * i * 2 * MAX_SIGMOID) / SIGMOID_TABLE_SIZE - MAX_SIGMOID;
            t_sigmoid_[i] = (float)(1.0 / (1.0 + exp(-x)));
        }
    }

    if (true)
    {
        t_log_ = new float[LOG_TABLE_SIZE + 1];
        for (int i = 0; i < LOG_TABLE_SIZE + 1; i++)
        {
            float x = (float)(((float)(i) + 1e-5) / LOG_TABLE_SIZE);
            t_log_[i] = (float)log(x);
        }
    }
};

template <typename VALUE_TYPE>
class OneSampleUpdator
{

protected:
    uint32_t dim;
    uint32_t neg_size;

public:
    uint32_t get_dim() const
    {
        return dim;
    }
    void write_me(bitsery::OutputStreamAdapter &bw) const
    {
        write(bw, dim);
        write(bw, neg_size);
    }
    void read_me(bitsery::InputStreamAdapter &br)
    {
        read(br, dim);
        read(br, neg_size);
    }

public:
    OneSampleUpdator()
    {
    }
    OneSampleUpdator(uint32_t dim, uint32_t neg_size) : dim(dim), neg_size(neg_size)
    {
    }
    virtual ~OneSampleUpdator()
    {
    }

protected:
    virtual void incr_word_count(uint32_t word) = 0;
    virtual void add_vector_to_wi(const Vector<VALUE_TYPE> &v, uint32_t word, float a) = 0;
    virtual void add_vector_to_wo(const Vector<VALUE_TYPE> &v, uint32_t word, float a) = 0;
    virtual uint32_t getNegative(uint32_t target) = 0;
    virtual Vector<VALUE_TYPE> get_vec_from_wi(uint32_t) = 0;
    virtual Vector<VALUE_TYPE> get_vec_from_wo(uint32_t) = 0;

    float update_one(uint32_t label, const std::vector<uint32_t> &words, float lr, bool update_wc)
    {
        if (words.empty())
            return 0;
        if (update_wc)
        {
            incr_word_count(label);
        }
        float loss = 0;
        if (true)
        {
            thread_local Vector<VALUE_TYPE> hidden_(dim);
            thread_local Vector<VALUE_TYPE> grad_(dim);
            computeHidden(words, hidden_);
            loss += negativeSampling(label, lr, hidden_, grad_);
            for (auto word : words)
            {
                add_vector_to_wi(grad_, word, 1.0f);
            }
        }
        return loss;
    }

private:
    void computeHidden(const std::vector<uint32_t> &words, Vector<VALUE_TYPE> &hidden_)
    {

        hidden_.zero();
        if (words.empty())
        {
            return;
        }

        for (auto word : words)
        {
            auto v = get_vec_from_wi(word);
            hidden_.add(v);
        }

        hidden_.mul(1.0f / words.size());
    }

    float negativeSampling(uint32_t target, float lr, Vector<VALUE_TYPE> &hidden_, Vector<VALUE_TYPE> &grad_)
    {
        float loss = 0.0f;
        grad_.zero();
        for (uint32_t n = 0; n <= neg_size; n++)
        {
            if (n == 0)
            {
                loss += binaryLogistic(target, true, lr, hidden_, grad_);
            }
            else
            {
                loss += binaryLogistic(getNegative(target), false, lr, hidden_, grad_);
            }
        }
        return loss / neg_size;
    }

    float binaryLogistic(uint32_t label, bool ytruth, float lr, Vector<VALUE_TYPE> &hidden_, Vector<VALUE_TYPE> &grad_)
    {
        auto v = get_vec_from_wo(label);
        float x = v.dot(hidden_);
        //fprintf(stderr, "x=%f %f %f \n", x, hidden_.at(0), v.at(0));
        float score = sigmoid(x);

        float alpha = lr * ((ytruth ? 1.0f : 0.0f) - score);
        grad_.add(v, alpha);
        add_vector_to_wo(hidden_, label, alpha); //wo[label]+=hidden_
        if (ytruth)
        {
            return -log(score);
        }
        else
        {
            return -log(1.0f - score);
        }
    }

    float log(float x)
    {
        if (x > 1.0)
        {
            return 0.0f;
        }
        int i = (int)(x * LOG_TABLE_SIZE);
        return t_log_[i];
    }

    float sigmoid(float x)
    {
        if (x < -MAX_SIGMOID)
        {
            return 0.0f;
        }
        else if (x > MAX_SIGMOID)
        {
            return 1.0f;
        }
        else
        {
            int i = (int)((x + MAX_SIGMOID) * SIGMOID_TABLE_SIZE / MAX_SIGMOID / 2);
            return t_sigmoid_[i];
        }
    }
};

template <typename VALUE_TYPE>
class Model : public OneSampleUpdator<VALUE_TYPE>
{
protected:
    bool use_cbow;
    uint32_t half_window;

public:
    void read_me(bitsery::InputStreamAdapter &br)
    {
        OneSampleUpdator<VALUE_TYPE>::read_me(br);
        read(br, use_cbow);
        read(br, half_window);
    }
    void write_me(bitsery::OutputStreamAdapter &bw) const
    {
        OneSampleUpdator<VALUE_TYPE>::write_me(bw);
        write(bw, use_cbow);
        write(bw, half_window);
    }

    void transform(const std::vector<uint32_t> &kmers, Vector<VALUE_TYPE> &vec)
    {
        vec.zero();
        if (kmers.empty())
        {
            return;
        }
        for (uint32_t kmer : kmers)
        {
            vec.add(this->get_vec_from_wi(kmer));
        }
        vec.mul(1.0f / kmers.size());
    }

public:
    Model() : OneSampleUpdator<VALUE_TYPE>()
    {
    }
    Model(uint32_t dim, uint32_t neg_size, bool use_cbow, uint32_t half_window) : OneSampleUpdator<VALUE_TYPE>(dim, neg_size), use_cbow(use_cbow), half_window(half_window)
    {
    }

    virtual void randomize_init() = 0;

    float update(const std::vector<uint32_t> &words, float lr, bool update_wc)
    {
        if (words.empty())
        {
            return 0;
        }
        float loss = 0;
        uint32_t count = 0;
        int ws = (int)words.size();
        int half_window = (int)this->half_window;
        std::vector<uint32_t> input_words;
        for (int i = 0; i < (int)words.size(); i++)
        {
            input_words.clear();
            uint32_t label = words.at(i);
            if (use_cbow)
            {

                for (int j = (int)(i - half_window); j < i + half_window; j++)
                {
                    if (j >= 0 && j != i && j < ws)
                    {
                        input_words.push_back(words.at(j));
                    }
                }
                loss += this->update_one(label, input_words, lr, update_wc);
                count++;
            }
            else
            {
                for (int j = (int)(i - half_window); j < i + half_window; j++)
                {
                    if (j >= 0 && j != i && j < ws)
                    {
                        input_words.clear();
                        input_words.push_back(words.at(j));
                        loss += this->update_one(label, input_words, lr, update_wc);
                        count++;
                    }
                }
            }
        }
        if (count == 0)
        {
            return 0;
        }
        else
        {
            return loss / (count);
        }
    }
};

template <typename VALUE_TYPE>
class SingleNodeModel : public Model<VALUE_TYPE>
{
protected:
    uint32_t num_word;
    std::vector<Vector<VALUE_TYPE>> wo_;
    std::vector<Vector<VALUE_TYPE>> wi_;
    std::vector<uint32_t> word_count;

public:
    void read_me(bitsery::InputStreamAdapter &br)
    {
        Model<VALUE_TYPE>::read_me(br);
        read(br, num_word);
        read(br, wo_);
        read(br, wi_);
        read(br, word_count);
    }

    void write_me(bitsery::OutputStreamAdapter &bw) const
    {
        Model<VALUE_TYPE>::write_me(bw);
        write(bw, num_word);
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
            output << i << " " << word_count.at(i) << "\n";
        }
    }

    void write_vec_bin(bitsery::OutputStreamAdapter &bw)
    {
        write(bw, wi_);
    }

public:
    SingleNodeModel() : Model<VALUE_TYPE>()
    {
    }
    SingleNodeModel(uint32_t num_word, uint32_t dim, uint32_t neg_size, bool use_cbow, uint32_t half_window) : Model<VALUE_TYPE>(dim, neg_size, use_cbow, half_window), num_word(num_word)
    {
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

protected:
    virtual void incr_word_count(uint32_t word)
    {
        word_count[word]++;
    }
    virtual void add_vector_to_wi(const Vector<VALUE_TYPE> &vec, uint32_t word, float a)
    {
        wi_[word].add(vec, a);
    }
    virtual void add_vector_to_wo(const Vector<VALUE_TYPE> &vec, uint32_t word, float a)
    {
        wo_[word].add(vec, a);
    }
    virtual uint32_t getNegative(uint32_t target)
    {
        while (true)
        {
            uint32_t i = (uint32_t)(sparc::myrand::uniform<double>() * num_word);
            if (i != target)
            {
                return i;
            }
        }
    }
    virtual Vector<VALUE_TYPE> get_vec_from_wi(uint32_t word)
    {
        return wi_.at(word);
    }
    virtual Vector<VALUE_TYPE> get_vec_from_wo(uint32_t word)
    {
        return wo_.at(word);
    }
};

template <typename VALUE_TYPE>
void write(bitsery::OutputStreamAdapter &bw, const SingleNodeModel<VALUE_TYPE> &obj)
{
    obj.write_me(bw);
}

#endif /* SOURCE_DIRECTORY__SRC_SPARC_MODEL_H_ */
