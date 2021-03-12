/**
 * @author Yunpeng Men
 * @email ypmen@pku.edu.cn
 * @create date 2020-05-19 09:38:18
 * @modify date 2020-05-19 09:38:18
 * @desc [description]
 */

#ifndef PULSARSEARCH
#define PULSARSEARCH

#include <utility>
#include <string.h>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp> 

#include "filterbank.h"
#include "subdedispersion.h"
#include "databuffer.h"
#include "downsample.h"
#include "equalize.h"
#include "rfi.h"

using namespace std;
using namespace boost::program_options;

class PulsarSearch
{
public:
    PulsarSearch();
    ~PulsarSearch();
    void prepare(DataBuffer<float> &databuffer);
    void run(DataBuffer<float> &databuffer);
public:
    //components
    Downsample downsample;
    Equalize equalize;
    RFI rfi;
    RealTime::SubbandDedispersion dedisp;

    //downsample
    int td;
    int fd;

    //rfi
    vector<pair<double, double>> zaplist;
    vector<vector<string>> rfilist;
    double bandlimit;
    double widthlimit;
    double bandlimitKT;
    float threKadaneT;
    float threKadaneF;
    float threMask;

    //dedispere
    double dms;
    double ddm;
    long int ndm;

    int ibeam;
    string rootname;
    int id;

    float outmean;
    float outstd;
    int outnbits;

    Filterbank fildedisp;
};

void plan(variables_map &vm, vector<PulsarSearch> &search);

#endif /* PULSARSEARCH */
