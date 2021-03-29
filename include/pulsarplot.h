/**
 * @author Yunpeng Men
 * @email ypmen@pku.edu.cn
 * @create date 2020-06-14 15:54:14
 * @modify date 2020-06-14 15:54:14
 * @desc [description]
 */

#ifndef PULSARPLOT
#define PULSARPLOT

#include <map>

#include "utils.h"
#include "gridsearch.h"

namespace Pulsar
{
    class PulsarPlot
    {
    public:
        PulsarPlot();
        ~PulsarPlot();
        void plot(const DedispersionLite &dedisp, const ArchiveLite &archive, GridSearch &gridsearch, std::map<std::string, std::string> &obsinfo, int id, const string &rootname, bool plotx=false);        
    };
}

#endif /* PULSARPLOT */
