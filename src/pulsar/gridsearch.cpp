/**
 * @author Yunpeng Men
 * @email ypmen@pku.edu.cn
 * @create date 2020-06-09 18:29:34
 * @modify date 2020-06-09 18:29:34
 * @desc [description]
 */

#include <algorithm>

#include "constants.h"
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
    clfd_q = -1.;

    nsubint = 0;
    nchan = 0;
    nbin = 0;
    mean = 0.;
    var = 0.;

    snr = 0.;
    width = 0.;
    p0 = 0.;
    p1 = 0.;
    acc = 0.;
    err_f0 = 0.;
    err_f1 = 0.;
    err_p0 = 0.;
    err_p1 = 0.;
    err_dm = 0.;
    err_acc = 0.;
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
    clfd_q = gridsearch.clfd_q;
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

    snr = gridsearch.snr;
    width = gridsearch.width;
    p0 = gridsearch.p0;
    p1 = gridsearch.p1;
    acc = gridsearch.acc;
    err_f0 = gridsearch.err_f0;
    err_f1 = gridsearch.err_f1;
    err_p0 = gridsearch.err_p0;
    err_p1 = gridsearch.err_p1;
    err_dm = gridsearch.err_dm;
    err_acc = gridsearch.err_dm;
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
    clfd_q = gridsearch.clfd_q;
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

    snr = gridsearch.snr;
    width = gridsearch.width;
    p0 = gridsearch.p0;
    p1 = gridsearch.p1;
    acc = gridsearch.acc;
    err_f0 = gridsearch.err_f0;
    err_f1 = gridsearch.err_f1;
    err_p0 = gridsearch.err_p0;
    err_p1 = gridsearch.err_p1;
    err_dm = gridsearch.err_dm;
    err_acc = gridsearch.err_dm;

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

    if (clfd_q > 0)
        clfd();
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

    int maxtdelayn = round(2*(bestdf0*tsuboff[nsubint-1]+0.5*bestdf1*tsuboff[nsubint-1]*tsuboff[nsubint-1])*nbin);
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
void GridSearch::get_snr_width()
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
    width = 0.;
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
 * @brief Calculate the derived parameters and errors
 * 
 * @param obsinfo 
 */
void GridSearch::get_error(std::map<std::string, std::string> &obsinfo)
{
    /** get error f0,f1,p0,p1,a */
    double obslen = stod(obsinfo["Obslen"]);
    double toaerr = width/(pow(M_PI, 0.25)*snr);
    err_f0 = sqrt(12)*toaerr*f0/obslen;
    err_f1 = sqrt(320)*toaerr*f0/(obslen*obslen);
    err_dm = 1./4.148741601e3*(frequencies[0]*frequencies.back())/sqrt((frequencies[0]/frequencies.back()+frequencies.back()/frequencies[0]+1)/3.-1)*toaerr;
    
    if (abs(f1)<err_f1/1000) f1 = 0.;
    
    p0 = 1/f0;
    p1 = -f1/(f0*f0);
    err_p0 = abs(err_f0/(f0*f0));
    err_p1 = abs(err_f1/(f0*f0));
    acc = f1/f0*CONST_C;
    err_acc = abs((err_f1*f0-err_f0*f1)/(f0*f0)*CONST_C);
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

    var = var==0? 1:var;
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
        tmp_var = tmp_var==0? 1:tmp_var;

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

/**
 * @brief remove rfi based on clfd
 * 
 */
void GridSearch::clfd()
{
    /**
     * @brief apply to time-phase image
     * 
     */
    vector<float> tfstd(nsubint*nchan, 0.);
    vector<float> tfstd_sort(nsubint*nchan, 0.);
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

            tfstd_sort[k*nchan+j] = tfstd[k*nchan+j] = sqrt(tmp_var);
        }
    }

    std::nth_element(tfstd_sort.begin(), tfstd_sort.begin()+nsubint*nchan/4, tfstd_sort.end(), std::less<float>());
    float Q1 = tfstd_sort[nsubint*nchan/4];
    std::nth_element(tfstd_sort.begin(), tfstd_sort.begin()+nsubint*nchan/4, tfstd_sort.end(), std::greater<float>());
    float Q3 = tfstd_sort[nsubint*nchan/4];

    float R = Q3-Q1;
    float vmin = Q1-clfd_q*R;
    float vmax = Q3+clfd_q*R;

    for (long int k=0; k<nsubint; k++)
    {
        for (long int j=0; j<nchan; j++)
        {
            if (tfstd[k*nchan+j]<vmin or tfstd[k*nchan+j]>vmax)
            {
                for (long int i=0; i<nbin; i++)
                {
                    profiles[k*nchan*nbin+j*nbin+i] = 0.;
                }
            }
        }
    }

    /**
     * @brief apply to time
     * 
     */

    vector<float> tstd(nsubint, 0.);
    vector<float> tstd_sort(nsubint, 0.);

    for (long int k=0; k<nsubint; k++)
    {
        vector<float> tmppro(nbin, 0.);
        for (long int j=0; j<nchan; j++)
        {
            for (long int i=0; i<nbin; i++)
            {
                tmppro[i] += profiles[k*nchan*nbin+j*nbin+i];
            }
        }
       
        double boxsum = 0.;
        for (long int i=0; i<nbin/2; i++)
        {
            boxsum += tmppro[i];
        }
        double min = boxsum;
        long int istart = 0;
        long int iend = nbin/2;
        for (long int i=0; i<nbin; i++)
        {
            boxsum -= tmppro[i];
            boxsum += tmppro[(i+nbin/2)%nbin];
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
            double tmp = tmppro[i%nbin];
            tmp_var += (tmp-tmp_mean)*(tmp-tmp_mean);
        }
        tmp_var /= (nbin/2);

        tstd_sort[k] = tstd[k] = sqrt(tmp_var);
    }

    std::nth_element(tstd_sort.begin(), tstd_sort.begin()+nsubint/4, tstd_sort.end(), std::less<float>());
    Q1 = tstd_sort[nsubint/4];
    std::nth_element(tstd_sort.begin(), tstd_sort.begin()+nsubint/4, tstd_sort.end(), std::greater<float>());
    Q3 = tstd_sort[nsubint/4];

    R = Q3-Q1;
    vmin = Q1-clfd_q*R;
    vmax = Q3+clfd_q*R;

    for (long int k=0; k<nsubint; k++)
    {
        if (tstd[k]<vmin or tstd[k]>vmax)
        {
            for (long int j=0; j<nchan; j++)
            {
                for (long int i=0; i<nbin; i++)
                {
                    profiles[k*nchan*nbin+j*nbin+i] = 0.;
                }
            }
        }
    }

    /**
     * @brief apply to frequency
     * 
     */

    vector<float> fstd(nchan, 0.);
    vector<float> fstd_sort(nchan, 0.);

    vector<float> fph(nchan*nbin, 0.);
    for (long int k=0; k<nsubint; k++)
    {
        for (long int j=0; j<nchan; j++)
        {
            for (long int i=0; i<nbin; i++)
            {
                fph[j*nbin+i] += profiles[k*nchan*nbin+j*nbin+i];
            }
        }
    }

    for (long int j=0; j<nchan; j++)
    {
        double boxsum = 0.;
        for (long int i=0; i<nbin/2; i++)
        {
            boxsum += fph[j*nbin+i];
        }
        double min = boxsum;
        long int istart = 0;
        long int iend = nbin/2;
        for (long int i=0; i<nbin; i++)
        {
            boxsum -= fph[j*nbin+i];
            boxsum += fph[(i+nbin/2)%nbin];
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
            double tmp = fph[j*nbin+i%nbin];
            tmp_var += (tmp-tmp_mean)*(tmp-tmp_mean);
        }
        tmp_var /= (nbin/2);

        fstd_sort[j] = fstd[j] = sqrt(tmp_var);
    }

    std::nth_element(fstd_sort.begin(), fstd_sort.begin()+nchan/4, fstd_sort.end(), std::less<float>());
    Q1 = fstd_sort[nchan/4];
    std::nth_element(fstd_sort.begin(), fstd_sort.begin()+nchan/4, fstd_sort.end(), std::greater<float>());
    Q3 = fstd_sort[nchan/4];

    R = Q3-Q1;
    vmin = Q1-clfd_q*R;
    vmax = Q3+clfd_q*R;

    for (long int k=0; k<nsubint; k++)
    {
        for (long int j=0; j<nchan; j++)
        {
            if (fstd[j]<vmin or fstd[j]>vmax)
            {
                for (long int i=0; i<nbin; i++)
                {
                    profiles[k*nchan*nbin+j*nbin+i] = 0.;
                }
            }
        }
    }
}