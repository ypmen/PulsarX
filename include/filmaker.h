/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2021-05-28 11:43:41
 * @modify date 2021-05-28 11:43:41
 * @desc [description]
 */

#ifndef FILMAKER_H
#define FILMAKER_H

#include <vector>
#include <string>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp> 

#include "filterbank.h"
#include "databuffer.h"
#include "downsample.h"
#include "equalize.h"
#include "baseline.h"
#include "rfi.h"
#include "filterbankwriter.h"

using namespace boost::program_options;

class FilMaker
{
public:
    FilMaker();
    ~FilMaker();
    void prepare(DataBuffer<float> &databuffer);
    void run(DataBuffer<float> &databuffer);
public:
    //components
    Downsample downsample;
    Equalize equalize;
    BaseLine baseline;
    RFI rfi;
    FilterbankWriter filwriter;

    //downsample
    int td;
    int fd;

    //baseline
    float bswidth;

    //rfi
    std::vector<std::pair<double, double>> zaplist;
    std::vector<std::vector<std::string>> rfilist;
    double bandlimit;
    double widthlimit;
    double bandlimitKT;
    float threKadaneT;
    float threKadaneF;
    float threMask;
    string filltype;

    //filterbankwriter
    std::string source_name;
    int telescope_id;
    double src_raj;
    double src_dej;
    int ibeam;
    string rootname;
    int id;

    float outmean;
    float outstd;
    int outnbits;
};

void plan(variables_map &vm, std::vector<FilMaker> &filmakers);

#endif /* FILMAKER_H */
