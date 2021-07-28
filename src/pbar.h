#ifndef _SRC_PBAR_H
#define _SRC_PBAR_H

#include <string>
#include "indicators.hpp"

template <int TYPE = 0>
class BaseBar
{
protected:
    std::string prefix;
    size_t item_size;
    indicators::BlockProgressBar *bar = 0;
    std::atomic<size_t> progress;

public:
    BaseBar(const std::string &prefix, size_t item_size) : prefix(prefix), item_size(item_size), progress(0)
    {
        if (TYPE == 0)
        {
            bar = new indicators::BlockProgressBar{
                indicators::option::BarWidth{50},
                indicators::option::PrefixText{prefix.c_str()},
                indicators::option::ForegroundColor{indicators::Color::blue},
                indicators::option::ShowElapsedTime{true},
                indicators::option::ShowRemainingTime{true},
                indicators::option::FontStyles{std::vector<indicators::FontStyle>{indicators::FontStyle::bold}},
                indicators::option::MaxProgress{item_size}};
        }
        else
        {
            bar = new indicators::BlockProgressBar{
                indicators::option::BarWidth{80},
                indicators::option::ShowElapsedTime{true},
                indicators::option::ShowRemainingTime{true},
                indicators::option::PrefixText{prefix.c_str()},
                indicators::option::ForegroundColor{indicators::Color::blue},
                indicators::option::FontStyles{
                    std::vector<indicators::FontStyle>{indicators::FontStyle::bold}},
                indicators::option::MaxProgress{item_size}};
        }
        indicators::show_console_cursor(false);
    }

    void end()
    {
        indicators::show_console_cursor(true);
        std::cout << std::endl;
    }
    virtual ~BaseBar()
    {

        if (bar)
            delete bar;
    }
    void tick()
    {
        if (TYPE == 0)
        {

            bar->tick();
        }
        else
        {

            bar->set_option(indicators::option::PostfixText{
                std::to_string(++progress) + "/" + std::to_string(item_size)});
            bar->tick();
        }
    }
};

template <int TYPE = 0>
class SimpleIndicator
{
protected:
    std::string prefix;
    size_t item_size;
    size_t percent;
    std::atomic<size_t> progress;

public:
    SimpleIndicator(const std::string &prefix, size_t item_size) : prefix(prefix), item_size(item_size), percent(0), progress{0}
    {
    }
    ~SimpleIndicator()
    {
    }
    void end()
    {
        std::cout << std::endl;
    }
    void tick(bool only_update = false)
    {
        size_t m;
        m = ++progress;
        if (TYPE == 0)
        {
            if (!only_update)
            {
                double this_percent = (m / double(item_size) * 100);
                std::cout << "\33[2K\r" << prefix << " " << std::fixed << std::setprecision(2) << this_percent << "%" << std::flush;
            }
        }
        else
        {
            if (!only_update)
            {
                size_t this_percent = (size_t)(m / double(item_size) * 1000);
                if (this_percent != percent)
                {
                    percent = this_percent;
                    std::cout << "\33[2K\r" << prefix << " " << std::fixed << std::setprecision(1) << this_percent / 10.0 << "%" << std::flush;
                }
            }
        }
    }
};

template <int TYPE = 0>
class SimpleUnknownIndicator
{
protected:
    std::string prefix;
    std::atomic<size_t> progress;

public:
    SimpleUnknownIndicator(const std::string &prefix) : prefix(prefix), progress{0}
    {
    }
    ~SimpleUnknownIndicator()
    {
        end();
    }
    void end()
    {
        std::cout << std::endl;
    }
    void tick(bool only_update = false)
    {

        ++progress;
        if (TYPE == 0)
        {
            if (!only_update)
            {
                std::cout << "\33[2K\r" << prefix << " " << std::fixed << std::setprecision(2) << progress << std::flush;
            }
        }
        else
        {
            if (!only_update)
            {
                std::cout << "\33[2K\r" << prefix << " " << std::fixed << std::setprecision(2) << progress << std::flush;
            }
        }
    }
};

using PBar = BaseBar<1>;
using PUnknownBar = SimpleUnknownIndicator<1>;

#endif