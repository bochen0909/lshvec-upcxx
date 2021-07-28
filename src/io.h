#ifndef __SRC_IO_H__
#define __SRC_IO_H__

#include <string>
#include <fstream>
#include <algorithm>
#include <vector>
#include <memory>
#include "bitsery/bitsery.h"
#include "bitsery/adapter/stream.h"
#include "bitsery/traits/vector.h"
#include "bitsery/traits/string.h"

#include "utils.h"
#include "gzstream.h"

inline void transform_seq(std::string &line)
{
    std::transform(line.begin(), line.end(), line.begin(), [](char c)
                   {
                       c = (char)::toupper(c);
                       if (c != 'A' && c != 'C' && c != 'G' && c != 'T')
                       {
                           c = 'N';
                       }
                       return c;
                   });
}
struct FastaRecord
{
    std::string id;
    std::string seq;
};

template <typename T>
T create_empty_record()
{
    return T{};
}

template <typename T>
bool record_is_valid(const T &)
{
    throw std::runtime_error("NA");
}

template <>
bool record_is_valid(const FastaRecord &t);
template <>
bool record_is_valid(const std::string &t);

std::vector<FastaRecord> read_fasta(const std::string &filepath);

template <typename RECORD_TYPE>
class ReaderBase
{
protected:
    std::string filepath;

public:
    ReaderBase(const std::string &filepath) : filepath(filepath)
    {
        if (!sparc::file_exists(filepath.c_str()))
        {
            throw std::runtime_error(std::string("file not exists: ") + filepath);
        }
    }

    virtual ~ReaderBase()
    {
        close();
    }

    virtual void close(){};

    virtual void open() = 0;

    virtual RECORD_TYPE next() = 0;

protected:
    virtual void read_one_record(RECORD_TYPE &ret) = 0;
};

template <typename T>
inline void bitsery_read(bitsery::InputStreamAdapter &br, T &val)
{
    return br.readBytes<sizeof(T)>(val);
}

template <typename RECORD_TYPE>
class BitseryReaderBase : public ReaderBase<RECORD_TYPE>
{
protected:
    std::ifstream *input = 0;
    bitsery::InputStreamAdapter *reader = 0;
    size_t n_readed = 0;
    size_t _size = 0;

public:
    BitseryReaderBase(const std::string &filepath) : ReaderBase<RECORD_TYPE>()
    {
    }

    virtual void close()
    {
        if (reader)
        {
            delete reader;
        }
        if (input)
        {
            delete input;
        }
        reader = 0;
        input = 0;
    }

    virtual void open()
    {
        input = new std::ifstream(this->filepath, std::ios::binary);
        reader = new bitsery::InputStreamAdapter{*input};
        bitsery_read(*reader, this->_size);
    }

    RECORD_TYPE next()
    {
        if (n_readed >= _size)
        {
            return create_empty_record<RECORD_TYPE>();
        }

        RECORD_TYPE ret;
        read_one_record(ret);
        n_readed++;
        return ret;
    }

protected:
    virtual void read_one_record(RECORD_TYPE &ret) = 0;
};

template <typename RECORD_TYPE>
class TextReaderBase : public ReaderBase<RECORD_TYPE>
{
protected:
    std::ifstream *input = 0;
    igzstream *gzinput = 0;
    bool is_stdin = false;
    bool is_gzip;

public:
    TextReaderBase(const std::string &filepath) : ReaderBase<RECORD_TYPE>(filepath)
    {
        is_gzip = sparc::endswith(filepath, ".gz");
        if (filepath == "-")
        {
            is_stdin = true;
        }
        open();
    }

    inline auto &read_line(std::string &line)
    {
        if (is_stdin)
        {
            return std::getline(std::cin, line);
        }
        else
        {
            if (is_gzip)
            {
                return std::getline(*gzinput, line);
            }
            else
            {
                return std::getline(*input, line);
            }
        }
    }

    virtual void open()
    {
        if (is_gzip)
        {
            gzinput = new igzstream(this->filepath.c_str());
        }
        else
        {
            input = new std::ifstream(this->filepath);
        }
    }

    virtual void close()
    {

        if (input)
        {
            delete input;
        }
        if (gzinput)
        {
            delete gzinput;
        }
        input = 0;
        gzinput = 0;
    }

protected:
    virtual void line_to_record(const std::string &, RECORD_TYPE &)
    {
        throw std::runtime_error("NA");
    }
};

class FastaTextReaderBase : public TextReaderBase<FastaRecord>
{
    FastaRecord record;

public:
    FastaTextReaderBase(const std::string &filepath) : TextReaderBase<FastaRecord>(filepath) {}
    FastaRecord next()
    {
        std::string line;
        while (this->read_line(line))
        {

            sparc::trim(line);
            if (line[0] == '>')
            {

                if (!record.id.empty())
                {
                    FastaRecord ret = record;
                    transform_seq(ret.seq);
                    record.id = line;
                    record.seq = "";
                    return ret;
                }
                else
                {
                    record.id = line;
                    record.seq = "";
                }
            }
            else
            {

                if (record.id.empty())
                {
                    throw std::runtime_error("failed to find header before sequence");
                }
                else
                {
                    record.seq += line;
                }
            }
        }

        if (record.id.empty())
        {
            return create_empty_record<FastaRecord>();
        }
        else
        {
            FastaRecord ret = record;
            transform_seq(ret.seq);
            record.id = "";
            record.seq = "";
            return ret;
        }
    }

protected:
    virtual void line_to_record(const std::string &, std::string &)
    {
        throw std::runtime_error("never be here");
    }
    virtual void read_one_record(FastaRecord &)
    {
        throw std::runtime_error("never be here");
    }
};

class FastqTextReaderBase : public TextReaderBase<FastaRecord>
{
    FastaRecord record;

public:
    FastqTextReaderBase(const std::string &filepath) : TextReaderBase<FastaRecord>(filepath) {}
    FastaRecord next()
    {
        std::string lines[4];
        std::string line;
        int i = 0;
        while (read_line(line))
        {

            sparc::trim(line);
            if (line.empty())
            {
                continue;
            }
            lines[i++] = line;
            if (i >= 4)
            {
                break;
            }
        }
        if (i == 0)
        {
            return create_empty_record<FastaRecord>();
        }
        if (i != 4)
        {
            throw std::runtime_error("failed to read fastq recrods, no enough rows");
        }
        if (lines[0][0] != '@' && lines[2][0] != '+')
        {
            throw std::runtime_error("failed to read fastq recrods");
        }
        FastaRecord ret;
        ret.id = lines[0];
        ret.seq = lines[1];
        transform_seq(ret.seq);
        return ret;
    }

protected:
    virtual void line_to_record(const std::string &, std::string &)
    {
        throw std::runtime_error("never be here");
    }
    virtual void read_one_record(FastaRecord &)
    {
        throw std::runtime_error("never be here");
    }
};

template <typename RECORD_TYPE>
class LineTextReaderBase : public TextReaderBase<RECORD_TYPE>
{
protected:
public:
    LineTextReaderBase(const std::string &filepath) : TextReaderBase<RECORD_TYPE>(filepath)
    {
    }

    RECORD_TYPE next()
    {
        std::string line;
        if (this->read_line(line))
        {
            RECORD_TYPE ret;
            line_to_record(line, ret);
            return ret;
        }
        else
        {
            return create_empty_record<RECORD_TYPE>();
        }
    }

protected:
    virtual void line_to_record(const std::string &line, RECORD_TYPE &ret) = 0;
    virtual void read_one_record(RECORD_TYPE &)
    {
        throw std::runtime_error("never be here");
    }
};

class SeqTextReaderBase : public LineTextReaderBase<FastaRecord>
{
public:
    SeqTextReaderBase(const std::string &filepath) : LineTextReaderBase<FastaRecord>(filepath) {}

protected:
    virtual void line_to_record(const std::string &line, FastaRecord &ret)
    {
        if (!line.empty())
        {
            ret.id = ">useless";
            ret.seq = line;
            transform_seq(ret.seq);
        }
    }
};

template <class BASE_READER, class RECORD>
class BatchReader
{
protected:
    std::vector<std::string> files;
    size_t curr_file_index = 0;
    std::shared_ptr<BASE_READER> reader = nullptr;

public:
    BatchReader(const std::vector<std::string> &files) : files(files)
    {
        for (auto &file : files)
        {
            if (!sparc::file_exists(file.c_str()))
            {
                throw std::runtime_error(
                    std::string("file not found: ") + file);
            }
        }
    }

    auto next(size_t batch_size)
    {
        std::vector<RECORD> ret;
        while (ret.size() < batch_size)
        {
            RECORD r;
            if (!_next(r))
            {
                break;
            }
            else
            {
                ret.push_back(r);
            }
        }
        return ret;
    }

protected:
    inline void check_reader()
    {
        if (reader == nullptr)
        {
            if (curr_file_index < files.size())
            {
                reader = std::make_shared<BASE_READER>(
                    files.at(curr_file_index));
                reader->open();
            }
        }
    }
    inline bool _next(RECORD &r)
    {
        check_reader();
        if (reader == nullptr)
        {
            return false;
        }

        r = reader->next();

        if (record_is_valid(r))
        {
            return true;
        }
        else
        {
            reader->close();
            reader = nullptr;
            curr_file_index++;
            return _next(r);
        }
    }
};

template <typename T>
inline void write(bitsery::OutputStreamAdapter &bw, T val)
{

    bw.writeBytes<sizeof(T)>(val);
}

template <>
inline void write(bitsery::OutputStreamAdapter &bw, double val)
{
    char *byteArray = reinterpret_cast<char *>(&val);
    for (uint32_t i = 0; i < sizeof(double); i++)
    {
        bw.writeBytes<sizeof(char)>(byteArray[i]);
    }
}

template <>
inline void write(bitsery::OutputStreamAdapter &bw, float val)
{
    char *byteArray = reinterpret_cast<char *>(&val);
    for (uint32_t i = 0; i < sizeof(float); i++)
    {
        bw.writeBytes<sizeof(char)>(byteArray[i]);
    }
}

template <>
inline void write(bitsery::OutputStreamAdapter &bw, bool val)
{
    int b = val ? 1 : 0;
    write(bw, b);
}

template <typename T>
inline void write(bitsery::OutputStreamAdapter &bw, std::vector<T> val)
{
    size_t n = val.size();
    write(bw, n);
    for (auto &x : val)
    {
        write(bw, x);
    }
}

#endif