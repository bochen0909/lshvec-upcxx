#include "model.h"
#include "serialization.h"

void save_model(uint32_t this_epoch, const std::string &output_prefix, SingleNodeModel<float> &model, bool zip_output)
{
    char txt[1024];
    sprintf(txt, "%s_%u.bin", output_prefix.c_str(), this_epoch);
    std::string filepath = txt;
    if (zip_output)
    {
        filepath += ".gz";
    }

    if (zip_output)
    {
        ogzstream output(filepath.c_str());
        bitsery::OutputStreamAdapter bw{output};
        write(bw, this_epoch);
        write(bw, model);
        output.flush();
    }
    else
    {
        std::ofstream output(txt, std::ios::binary | std::ios::trunc);
        bitsery::OutputStreamAdapter bw{output};
        write(bw, this_epoch);
        write(bw, model);
        output.flush();
    }
}

void save_vector(uint32_t this_epoch, const std::string &output_prefix, SingleNodeModel<float> &model, bool zip_output)
{
    char txt[1024];
    sprintf(txt, "%s_%u.vec", output_prefix.c_str(), this_epoch);
    std::string filepath = txt;
    if (zip_output)
    {
        filepath += ".gz";
    }

    if (zip_output)
    {
        ogzstream output(filepath.c_str());
        model.write_vec(output);
        output.flush();
    }
    else
    {
        std::ofstream output(txt);
        model.write_vec(output);
        output.flush();
    }
}
void save_vector_bin(uint32_t this_epoch, const std::string &output_prefix, SingleNodeModel<float> &model, bool zip_output)
{
    char txt[1024];
    sprintf(txt, "%s_%u.vec.bin", output_prefix.c_str(), this_epoch);
    std::string filepath = txt;
    if (zip_output)
    {
        filepath += ".gz";
    }

    if (zip_output)
    {
        ogzstream output(filepath.c_str());
        bitsery::OutputStreamAdapter bw(output);
        model.write_vec_bin(bw);
        output.flush();
    }
    else
    {
        std::ofstream output(txt, std::ios::binary | std::ios::trunc);
        bitsery::OutputStreamAdapter bw(output);
        model.write_vec_bin(bw);
        output.flush();
    }
}

void read_vec_bin(const std::string &binpath, std::vector<Vector<float>> &wi)
{
    if (sparc::endswith(binpath, ".gz"))
    {
        igzstream input(binpath.c_str());
        bitsery::InputStreamAdapter br(input);
        read(br, wi);
    }
    else
    {
        std::ifstream input(binpath, std::ios::binary);
        bitsery::InputStreamAdapter br(input);
        read(br, wi);
    }
}

void read_model(const std::string &modelpath, SingleNodeModel<float> &model)
{
    uint32_t this_epoch;
    if (sparc::endswith(modelpath, ".gz"))
    {
        igzstream input(modelpath.c_str());
        bitsery::InputStreamAdapter br(input);
        read(br, this_epoch);
        model.read_me(br);
    }
    else
    {
        std::ifstream input(modelpath, std::ios::binary);
        bitsery::InputStreamAdapter br(input);
        read(br, this_epoch);
        model.read_me(br);
    }
}