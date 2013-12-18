#ifndef DYNEARTHSOL3D_OUTPUT_HPP
#define DYNEARTHSOL3D_OUTPUT_HPP

class Output
{
private:
    const std::string &modelname;
    const double start_time;
    int frame;

public:
    Output(const Param& param, double start_time, int start_frame);
    ~Output();
    void write(const Variables& var);
};


void restart();

#endif
