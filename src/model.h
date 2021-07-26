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
#include "static_block.hpp"
#include "utils.h"
template <typename VALUE_TYPE>
class Vector
{
    std::vector<VALUE_TYPE> val;

public:
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

private:
    Vector<VALUE_TYPE> hidden_;
    Vector<VALUE_TYPE> grad_;

public:
    OneSampleUpdator(uint32_t dim, uint32_t neg_size) : dim(dim), neg_size(neg_size), hidden_(dim), grad_(dim)
    {
    }
    virtual ~OneSampleUpdator()
    {
    }

protected:
    virtual void incr_word_count(uint32_t word) = 0;
    virtual void add_grad_to_wi(const Vector<VALUE_TYPE> &grad, uint32_t word, float a) = 0;
    virtual void add_grad_to_wo(const Vector<VALUE_TYPE> &grad, uint32_t word, float a) = 0;
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
            computeHidden(words);
            loss += negativeSampling(label, lr);
            for (auto word : words)
            {
                add_grad_to_wi(grad_, word, 1.0f);
            }
        }
        return loss;
    }

private:
    void computeHidden(const std::vector<uint32_t> &words)
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

    float negativeSampling(uint32_t target, float lr)
    {
        float loss = 0.0f;
        grad_.zero();
        for (uint32_t n = 0; n <= neg_size; n++)
        {
            if (n == 0)
            {
                loss += binaryLogistic(target, true, lr);
            }
            else
            {
                loss += binaryLogistic(getNegative(target), false, lr);
            }
        }
        return loss/neg_size;
    }

    float binaryLogistic(uint32_t label, bool ytruth, float lr)
    {
        auto v = get_vec_from_wo(label);
        float score = sigmoid(v.dot(hidden_));

        float alpha = lr * ((ytruth ? 1.0f : 0.0f) - score);
        grad_.add(v, alpha);
        add_grad_to_wo(hidden_, label, alpha); //wo[label]+=hidden_
        if (ytruth)
        {
            return -log(score);
        }
        else
        {
            return -log(1.0f - score);
        }
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
        for (int i = 0; i < (int)words.size(); i++)
        {
            std::vector<uint32_t> input_words;
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

protected:
    virtual void incr_word_count(uint32_t word)
    {
        word_count[word]++;
    }
    virtual void add_grad_to_wi(const Vector<VALUE_TYPE> &grad, uint32_t word, float a)
    {
        wi_[word].add(grad, a);
    }
    virtual void add_grad_to_wo(const Vector<VALUE_TYPE> &grad, uint32_t word, float a)
    {
        wo_[word].add(grad, a);
    }
    virtual uint32_t getNegative(uint32_t target)
    {
        return (uint32_t)(sparc::myrand::uniform<double>() * num_word);
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
#endif /* SOURCE_DIRECTORY__SRC_SPARC_MODEL_H_ */
