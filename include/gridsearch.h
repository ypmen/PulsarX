/**
 * @author Yunpeng Men
 * @email ypmen@pku.edu.cn
 * @create date 2020-06-09 17:27:20
 * @modify date 2020-06-09 17:27:20
 * @desc [description]
 */

#ifndef GRIDSEARCH
#define GRIDSEARCH

#include <vector>

#include "archivelite.h"

using namespace std;

namespace Pulsar
{
    class GridSearch
    {
    public:
        GridSearch();
        GridSearch(const GridSearch &gridsearch);
        GridSearch & operator=(const GridSearch &gridsearch);
        ~GridSearch();
        void prepare(ArchiveLite &arch);
        void runFFdot();
        void runDM();
        bool bestprofiles();
        void dumpFFdot(const string &rootname)
        {
            ofstream outfile;
            outfile.open(rootname+".ffdot", ofstream::binary);

            outfile.write((char *)(&mxsnr_ffdot[0]), sizeof(float)*ndf1*ndf0);

            outfile.close();
        }
        void dumpDM(const string &rootname)
        {
            ofstream outfile;
            outfile.open(rootname+".dm", ofstream::binary);

            outfile.write((char *)(&vsnr_dm[0]), sizeof(float)*nddm);

            outfile.close();
        }
        void dumpProfiles(const string &rootname)
        {
            ofstream outfile;
            outfile.open(rootname+".sft", ofstream::binary);

            outfile.write((char *)(&profiles[0]), sizeof(float)*nsubint*nchan*nbin);

            outfile.close();
        }
    public:
        float get_chisq(vector<float> &pro);
        void get_snr_width(double &snr, double &width);
        void get_rms();
        void subints_normalize();
    public:
        double df0start;
        double df0step;
        int ndf0;
        double df1start;
        double df1step;
        int ndf1;
        double ddmstart;
        double ddmstep;
        int nddm;
    public:
        double f0;
        double f1;
        double dm;
        double bestdf0;
        double bestdf1;
        double bestddm;
        bool ffdotsearch;
        bool dmsearch;
        int nsubint;
        int nchan;
        int nbin;
        double mean;
        double var;
        vector<double> ffold;
        vector<double> tsuboff;
        vector<double> frequencies;
        vector<float> profiles;
        vector<double> profile;
        vector<float> mxsnr_ffdot; 
        vector<float> vsnr_dm;
    };
}

#endif /* GRIDSEARCH */
