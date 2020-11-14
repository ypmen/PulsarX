/**
 * @author Yunpeng Men
 * @email ypmen@pku.edu.cn
 * @create date 2020-06-09 18:29:34
 * @modify date 2020-06-09 18:29:34
 * @desc [description]
 */

#include <algorithm>

#include "dedisperse.h"
#include "gridsearch.h"
#include "dedispersionlite.h"

using namespace std;
using namespace Pulsar;

GridSearch::GridSearch()
{
    ffdotsearch = false;
    dmsearch = false;
    f0 = 0.;
    f1 = 0.;
    dm = 0.;
    bestdf0 = 0.;
    bestdf1 = 0.;
    bestddm = 0.;
    df0start = 0.;
    df0step = 0.;
    ndf0 = 0;
    df1start = 0.;
    df1step = 0.;
    ndf1 = 0;
    ddmstart = 0.;
    ddmstep = 0.;
    nddm = 0;
    nsubint = 0;
    nchan = 0;
    nbin = 0;
    mean = 0.;
    var = 0.;
}

GridSearch::GridSearch(const GridSearch &gridsearch)
{
    ffdotsearch = gridsearch.ffdotsearch;
    dmsearch = gridsearch.dmsearch;
    f0 = gridsearch.f0;
    f1 = gridsearch.f1;
    dm = gridsearch.dm;
    bestdf0 = gridsearch.bestdf0;
    bestdf1 = gridsearch.bestdf1;
    bestddm = gridsearch.bestddm;
    df0start = gridsearch.df0start;
    df0step = gridsearch.df0step;
    ndf0 = gridsearch.ndf0;
    df1start = gridsearch.df1start;
    df1step = gridsearch.df1step;
    ndf1 = gridsearch.ndf1;
    ddmstart = gridsearch.ddmstart;
    ddmstep = gridsearch.ddmstep;
    nddm = gridsearch.nddm;
    nsubint = gridsearch.nsubint;
    nchan = gridsearch.nchan;
    nbin = gridsearch.nbin;
    mean = gridsearch.mean;
    var = gridsearch.var;
    ffold = gridsearch.ffold;
    tsuboff = gridsearch.tsuboff;
    frequencies = gridsearch.frequencies;
    profiles = gridsearch.profiles;
    profile = gridsearch.profile;
    mxsnr_ffdot = gridsearch.mxsnr_ffdot;
    vsnr_dm = gridsearch.vsnr_dm;
}

GridSearch & GridSearch::operator=(const GridSearch &gridsearch)
{
    ffdotsearch = gridsearch.ffdotsearch;
    dmsearch = gridsearch.dmsearch;
    f0 = gridsearch.f0;
    f1 = gridsearch.f1;
    dm = gridsearch.dm;
    bestdf0 = gridsearch.bestdf0;
    bestdf1 = gridsearch.bestdf1;
    bestddm = gridsearch.bestddm;
    df0start = gridsearch.df0start;
    df0step = gridsearch.df0step;
    ndf0 = gridsearch.ndf0;
    df1start = gridsearch.df1start;
    df1step = gridsearch.df1step;
    ndf1 = gridsearch.ndf1;
    ddmstart = gridsearch.ddmstart;
    ddmstep = gridsearch.ddmstep;
    nddm = gridsearch.nddm;
    nsubint = gridsearch.nsubint;
    nchan = gridsearch.nchan;
    nbin = gridsearch.nbin;
    mean = gridsearch.mean;
    var = gridsearch.var;
    ffold = gridsearch.ffold;
    tsuboff = gridsearch.tsuboff;
    frequencies = gridsearch.frequencies;
    profiles = gridsearch.profiles;
    profile = gridsearch.profile;
    mxsnr_ffdot = gridsearch.mxsnr_ffdot;
    vsnr_dm = gridsearch.vsnr_dm;

    return *this;
}

GridSearch::~GridSearch(){}

void GridSearch::prepare(ArchiveLite &arch)
{
    f0 = arch.f0;
    f1 = arch.f1;
    dm = arch.dm;

    nsubint = arch.profiles.size();
    nchan = arch.nchan;
    nbin = arch.nbin;

    frequencies = arch.frequencies;
    ffold.resize(nsubint, 0.);
    tsuboff.resize(nsubint, 0.);
    profiles.resize(nsubint*nchan*nbin, 0.);

    double tcentre = 0.;
    for (long int l=0; l<nsubint; l++)
    {
        tcentre += arch.profiles[l].offs_sub;
    }
    tcentre /= nsubint;

    int npol = arch.npol;
    for (long int l=0; l<nsubint; l++)
    {
        tsuboff[l] = arch.profiles[l].offs_sub - tcentre;
        ffold[l] = arch.profiles[l].ffold;

        for (long int k=0; k<npol; k++)
        {
            for (long int j=0; j<nchan; j++)
            {
                for (long int i=0; i<nbin; i++)
                {
                   profiles[l*nchan*nbin+j*nbin+i] += arch.profiles[l].data[k*nchan*nbin+j*nbin+i];
                }
            }
        }
    }

    subints_normalize();
    get_rms();
}

void GridSearch::runFFdot()
{
    ffdotsearch = true;

    vector<float> mxtph(nsubint*nbin, 0.);

    for (long int k=0; k<nsubint; k++)
    {
        for (long int j=0; j<nchan; j++)
        {
            for (long int i=0; i<nbin; i++)
            {
                mxtph[k*nbin+i] += profiles[k*nchan*nbin+j*nbin+i];
            }
        }
    }

    mxsnr_ffdot.resize(ndf1*ndf0, 0.);

    int if1 = -1;
    int if0 = -1;
    float maxsnr = 0.;
    for (long int k1=0; k1<ndf1; k1++)
    {
        double df1 = df1start + k1*df1step;
        for (long int k0=0; k0<ndf0; k0++)
        {

            double df0 = df0start + k0*df0step;

            vector<float> pro(nbin, 0.);
            float *pmxtph = &mxtph[0];
            for (long int j=0; j<nsubint; j++)
            {
                int delayn = round((df0*tsuboff[j]+0.5*df1*tsuboff[j]*tsuboff[j])*nbin);
                delayn %= nbin;
                if (delayn<0) delayn += nbin;
                float *pd = pmxtph-delayn+nbin;
                for (long int i=0; i<delayn; i++)
                {
                    pro[i] += pd[i];
                }
                pd = pmxtph-delayn;
                for (long int i=delayn; i<nbin; i++)
                {
                    pro[i] += pd[i];
                }
                
                pmxtph += nbin;
            }

            mxsnr_ffdot[k1*ndf0+k0] = get_chisq(pro);

            if (mxsnr_ffdot[k1*ndf0+k0] > maxsnr)
            {
                maxsnr = mxsnr_ffdot[k1*ndf0+k0];
                if1 = k1;
                if0 = k0;
            }
        }
    }

    bestdf1 = df1start + if1*df1step;
    bestdf0 = df0start + if0*df0step;
}

void GridSearch::runDM()
{
    dmsearch = true;

    vector<float> mxfph(nchan*nbin, 0.); 

    double ffoldmean = 0.;

    long int m = 0;
    for (long int k=0; k<nsubint; k++)
    {
        ffoldmean += ffold[k];
        long int l = 0;
        for (long int j=0; j<nchan; j++)
        {
            for (long int i=0; i<nbin; i++)
            {
                mxfph[l++] += profiles[m++];
            }
        }
    }
    ffoldmean /= nsubint;

    double fmin = 1e6;
    double fmax = 0.;
    for (long int j=0; j<nchan; j++)
    {
        fmax = frequencies[j]>fmax? frequencies[j]:fmax;
        fmin = frequencies[j]<fmin? frequencies[j]:fmin;
    }

    vsnr_dm.resize(nddm, 0.);

    int idm = -1;
    float maxsnr = 0.;
    for (long int k=0; k<nddm; k++)
    {
        vector<float> pro(nbin, 0.);
        
        double ddm = ddmstart + k*ddmstep;
        
        float *pmxfph = &mxfph[0];
        for (long int j=0; j<nchan; j++)
        {
            int delayn = round(DedispersionLite::dmdelay(ddm, fmax, frequencies[j])*ffoldmean*nbin);
            delayn %= nbin;
            if (delayn<0) delayn += nbin;
            int nleft = nbin - delayn;
            float *pd = pmxfph + delayn;
            for (long int i=0; i<nleft; i++)
            {
                pro[i] += pd[i];
            }
            pd = pmxfph + delayn - nbin;
            for (long int i=nleft; i<nbin; i++)
            {
                pro[i] += pd[i];
            }
            pmxfph += nbin;
        }

        vsnr_dm[k] = get_chisq(pro);

        if (vsnr_dm[k] > maxsnr)
        {
            maxsnr = vsnr_dm[k];
            idm = k;
        }
    }

    bestddm = ddmstart + idm*ddmstep;
}

bool GridSearch::bestprofiles()
{
    if (!dmsearch and !ffdotsearch) return false;

    double fmin = 1e6;
    double fmax = 0.;
    for (long int j=0; j<nchan; j++)
    {
        fmax = frequencies[j]>fmax? frequencies[j]:fmax;
        fmin = frequencies[j]<fmin? frequencies[j]:fmin;
    }

    int maxtdelayn = round(2*(bestdf0*tsuboff[nsubint-1])*nbin);
    int maxfdelayn = round(DedispersionLite::dmdelay(bestddm, fmax, fmin)*(f0+bestdf0)*nbin);

    if (abs(maxtdelayn)==0 and abs(maxfdelayn)==0)
    {
        bestddm = 0.;
        bestdf0 = 0;
        bestdf1 = 0.;

        return true;
    }

    vector<float> profiles2(nsubint*nchan*nbin, 0.);

    for (long int k=0; k<nsubint; k++)
    {
        int tdelayn = round((bestdf0*tsuboff[k]+0.5*bestdf1*tsuboff[k]*tsuboff[k])*nbin);
        for (long int j=0; j<nchan; j++)
        {
            int fdelayn = round(DedispersionLite::dmdelay(bestddm, fmax, frequencies[j])*(f0+bestdf0+(f1+bestdf1)*tsuboff[k])*nbin);
            
            int delayn = (tdelayn-fdelayn)%nbin;
            if (delayn<0) delayn += nbin;
            for (long int i=0; i<delayn; i++)
            {
                profiles2[k*nchan*nbin+j*nbin+i] = profiles[k*nchan*nbin+j*nbin+(i-delayn+nbin)];
            }
            for (long int i=delayn; i<nbin; i++)
            {
                profiles2[k*nchan*nbin+j*nbin+i] = profiles[k*nchan*nbin+j*nbin+(i-delayn)];
            }
        }
    }

    profiles = profiles2;

    dm += bestddm;
    f0 += bestdf0;
    f1 += bestdf1;

    bestddm = 0.;
    bestdf0 = 0;
    bestdf1 = 0.;

    return true;
}

float GridSearch::get_chisq(vector<float> &pro)
{  
    if (pro.empty()) return 0.;

    float temp_mean = 0.;
    float temp_var = 0.;

    int size = pro.size();
    for (long int i=0; i<size; i++)
    {
        temp_mean += pro[i];
        temp_var += pro[i]*pro[i]; 
    }

    temp_mean /= size;
    temp_var /= size;
    temp_var -= temp_mean*temp_mean;

    return temp_var/var;
}

/**
 * @brief get the snr after tscrunch and fscrunch
 * 
 * @return snr
 */
void GridSearch::get_snr_width(double &snr, double &width)
{
    profile.resize(nbin, 0.);
    for (long int k=0; k<nsubint; k++)
    {
        for (long int j=0; j<nchan; j++)
        {
            for (long int i=0; i<nbin; i++)
            {
                profile[i] += profiles[k*nchan*nbin+j*nbin+i];
            }
        }
    }

    double boxsum = 0.;
    for (long int i=0; i<nbin/2; i++)
    {
        boxsum += profile[i];
    }
    double min = boxsum;
    long int istart = 0;
    long int iend = nbin/2;
    for (long int i=0; i<nbin; i++)
    {
        boxsum -= profile[i];
        boxsum += profile[(i+nbin/2)%nbin];
        if (boxsum < min)
        {
            min = boxsum;
            istart = i+1;
            iend = nbin/2+i+1;
        }
    }

    double tmp_mean = min/(nbin/2);
    // double tmp_var = 0.;
    // for (long int i=istart; i<iend; i++)
    // {
    //     double tmp = profile[i%nbin];
    //     tmp_var += (tmp-tmp_mean)*(tmp-tmp_mean);
    // }
    // tmp_var /= (nbin/2);

    for (long int i=0; i<nbin; i++)
    {
        profile[i] -= tmp_mean;
        //profile[i] /= sqrt(tmp_var);
        profile[i] /= sqrt(var);
    }

    snr = 0.;
    width = 0;
    for (long int w=1; w<=nbin/2; w++)
    {
        double boxsum = 0.;
        for (long int i=0; i<w; i++)
        {
            boxsum += profile[i];
        }
        double max = boxsum;
        for (long int i=0; i<nbin; i++)
        {
            boxsum -= profile[i];
            boxsum += profile[(i+w)%nbin];
            if (boxsum > max)
            {
                max = boxsum;
                width = w;
            }
        }

        snr = (max/sqrt(w))>snr ? (max/sqrt(w)):snr;
    }

    width /= f0*nbin;
}

/**
 * @brief get the mean,rms,profile after tscrunch and fscrunch
 * 
 */
void GridSearch::get_rms()
{
    var = 0.;
    mean = 0.;
    for (long int k=0; k<nsubint; k++)
    {
        for (long int j=0; j<nchan; j++)
        {
            double boxsum = 0.;
            for (long int i=0; i<nbin/2; i++)
            {
                boxsum += profiles[k*nchan*nbin+j*nbin+i];
            }
            double min = boxsum;
            long int istart = 0;
            long int iend = nbin/2;
            for (long int i=0; i<nbin; i++)
            {
                boxsum -= profiles[k*nchan*nbin+j*nbin+i];
                boxsum += profiles[k*nchan*nbin+j*nbin+(i+nbin/2)%nbin];
                if (boxsum < min)
                {
                    min = boxsum;
                    istart = i+1;
                    iend = nbin/2+i+1;
                }
            }

            double tmp_mean = min/(nbin/2);
            double tmp_var = 0.;
            for (long int i=istart; i<iend; i++)
            {
                double tmp = profiles[k*nchan*nbin+j*nbin+i%nbin];
                tmp_var += (tmp-tmp_mean)*(tmp-tmp_mean);
            }
            tmp_var /= (nbin/2);

            mean += tmp_mean;
            var += tmp_var;
        }
    }
}

/**
 * @brief normalize the subints
 * 
 */
void GridSearch::subints_normalize()
{
    for (long int k=0; k<nsubint; k++)
    {
        double boxsum = 0.;
        for (long int j=0; j<nchan; j++)
        {
            for (long int i=0; i<nbin/2; i++)
            {
                boxsum += profiles[k*nchan*nbin+j*nbin+i];
            }
        }
        double min = boxsum;
        long int istart = 0;
        long int iend = nbin/2;

        for (long int i=0; i<nbin; i++)
        {
            for (long int j=0; j<nchan; j++)
            {
                boxsum -= profiles[k*nchan*nbin+j*nbin+i];
                boxsum += profiles[k*nchan*nbin+j*nbin+(i+nbin/2)%nbin];
            }

            if (boxsum < min)
            {
                min = boxsum;
                istart = i+1;
                iend = nbin/2+i+1;
            }
        }

        double tmp_mean = min/((nbin/2)*nchan);
        double tmp_var = 0.;
        for (long int j=0; j<nchan; j++)
        {
            for (long int i=istart; i<iend; i++)
            {
                double tmp = profiles[k*nchan*nbin+j*nbin+i%nbin];
                tmp_var += (tmp-tmp_mean)*(tmp-tmp_mean);
            }
        }
        tmp_var /= (nbin/2)*nchan;

        for (long int j=0; j<nchan; j++)
        {
            for (long int i=0; i<nbin; i++)
            {
                profiles[k*nchan*nbin+j*nbin+i] -= tmp_mean;
                profiles[k*nchan*nbin+j*nbin+i] /= sqrt(tmp_var);
            }
        }
    }
}