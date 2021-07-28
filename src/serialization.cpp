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
        model.write_vec_bin(output);
        output.flush();
    }
    else
    {
        std::ofstream output(txt, std::ios::binary | std::ios::trunc);
        model.write_vec_bin(output);
        output.flush();
    }

}