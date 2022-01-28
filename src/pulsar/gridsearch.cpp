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

    double tcentre = (arch.ref_epoch-arch.start_mjd).to_second();

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
    
    if (!zaplist.empty())
        zap();
    
    //subints_normalize();
    if (f0*arch.profiles[0].tsubint >= 1)
        normalize();

    if (clfd_q > 0)
        clfd3();

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

    long int m = 0;
    for (long int k=0; k<nsubint; k++)
    {
        long int l = 0;
        for (long int j=0; j<nchan; j++)
        {
            for (long int i=0; i<nbin; i++)
            {
                mxfph[l++] += profiles[m++];
            }
        }
    }

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
        
        if (dm+ddm >= 0)
        {
            float *pmxfph = &mxfph[0];
            for (long int j=0; j<nchan; j++)
            {
                int delayn = round(DedispersionLite::dmdelay(ddm, fmax, frequencies[j])*f0*nbin);
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
        }
        else
        {
            vsnr_dm[k] = 0.;
        }

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

    int maxdelayn = -1;
    for (long int k=0; k<nsubint; k++)
    {
        int tdelayn = round((bestdf0*tsuboff[k]+0.5*bestdf1*tsuboff[k]*tsuboff[k])*nbin);
        for (long int j=0; j<nchan; j++)
        {
            int fdelayn = round(DedispersionLite::dmdelay(bestddm, fmax, frequencies[j])*(f0+bestdf0+(f1+bestdf1)*tsuboff[k])*nbin);
            int delayn = (tdelayn-fdelayn)%nbin;
            maxdelayn = abs(delayn) > maxdelayn ? abs(delayn) : maxdelayn;
        }
    }

    if (maxdelayn < 2)
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

    double tmp_mean = 0.;
    double tmp_var = 0.;

    get_mean_var2<std::vector<double>::iterator>(profile.begin(), nbin, tmp_mean, tmp_var);

    tmp_var = tmp_var==0. ? 1.:tmp_var;

    for (long int i=0; i<nbin; i++)
    {
        profile[i] -= tmp_mean;
        profile[i] /= sqrt(tmp_var);
    }

    snr = 0.;
    width = 0.;
    pulsespan.resize(2, 0);
    for (long int w=1; w<=nbin/2; w++)
    {
        double boxsum = 0.;
        for (long int i=0; i<w; i++)
        {
            boxsum += profile[i];
        }
        int istart = 0, iend=0;
        double max = boxsum;
        for (long int i=0; i<nbin; i++)
        {
            boxsum -= profile[i];
            boxsum += profile[(i+w)%nbin];
            if (boxsum > max)
            {
                max = boxsum;
                istart = i+1;
                iend = ((i+1+w)%nbin);
            }
        }

        if ((max/sqrt(w))>snr)
        {
            snr = (max/sqrt(w));
            width = w;
            pulsespan[0] = istart;
            pulsespan[1] = iend;
        }
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
    acc = -f1/f0*CONST_C;
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
            double tmp_mean = 0.;
            double tmp_var = 0.;

            get_mean_var<std::vector<float>::iterator>(profiles.begin()+k*nchan*nbin+j*nbin, nbin, tmp_mean, tmp_var);

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
        double tmp_mean = 0.;
        double tmp_var = 0.;

        get_mean_var<std::vector<float>::iterator>(profiles.begin()+k*nchan*nbin, nchan, nbin, tmp_mean, tmp_var);

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
 * @brief normalize all profiles
 * 
 */
void GridSearch::normalize()
{
    for (long int k=0; k<nsubint; k++)
    {
        for (long int j=0; j<nchan; j++)
        {
            double tmp_mean = 0.;
            double tmp_var = 0.;

            get_mean_var<std::vector<float>::iterator>(profiles.begin()+k*nchan*nbin+j*nbin, nbin, tmp_mean, tmp_var);

            tmp_var = tmp_var==0? 1:tmp_var;
            tmp_var = std::sqrt(tmp_var);

            for (long int i=0; i<nbin; i++)
            {
                profiles[k*nchan*nbin+j*nbin+i] -= tmp_mean;
                profiles[k*nchan*nbin+j*nbin+i] /= tmp_var;
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
    vector<float> tfskewness(nsubint*nchan, 0.);
    vector<float> tfskewness_sort(nsubint*nchan, 0.);
    vector<float> tfkurtosis(nsubint*nchan, 0.);
    vector<float> tfkurtosis_sort(nsubint*nchan, 0.);
    for (long int k=0; k<nsubint; k++)
    {
        for (long int j=0; j<nchan; j++)
        {        
            double tmp_skewness = 0.;
            double tmp_kurtosis = 0.;

            get_skewness_kurtosis<std::vector<float>::iterator>(profiles.begin()+k*nchan*nbin+j*nbin, nbin, tmp_skewness, tmp_kurtosis);

            tfskewness_sort[k*nchan+j] = tfskewness[k*nchan+j] = tmp_skewness;
            tfkurtosis_sort[k*nchan+j] = tfkurtosis[k*nchan+j] = tmp_kurtosis;
        }
    }

    std::nth_element(tfskewness_sort.begin(), tfskewness_sort.begin()+nsubint*nchan/4, tfskewness_sort.end(), std::less<float>());
    float Q1_skewness = tfskewness_sort[nsubint*nchan/4];
    std::nth_element(tfskewness_sort.begin(), tfskewness_sort.begin()+nsubint*nchan/4, tfskewness_sort.end(), std::greater<float>());
    float Q3_skewness = tfskewness_sort[nsubint*nchan/4];

    float R_skewness = Q3_skewness-Q1_skewness;
    float vmin_skewness = Q1_skewness-clfd_q*R_skewness;
    float vmax_skewness = Q3_skewness+clfd_q*R_skewness;

    std::nth_element(tfkurtosis_sort.begin(), tfkurtosis_sort.begin()+nsubint*nchan/4, tfkurtosis_sort.end(), std::less<float>());
    float Q1_kurtosis = tfkurtosis_sort[nsubint*nchan/4];
    std::nth_element(tfkurtosis_sort.begin(), tfkurtosis_sort.begin()+nsubint*nchan/4, tfkurtosis_sort.end(), std::greater<float>());
    float Q3_kurtosis = tfkurtosis_sort[nsubint*nchan/4];

    float R_kurtosis = Q3_kurtosis-Q1_kurtosis;
    float vmin_kurtosis = Q1_kurtosis-clfd_q*R_kurtosis;
    float vmax_kurtosis = Q3_kurtosis+clfd_q*R_kurtosis;

    for (long int k=0; k<nsubint; k++)
    {
        for (long int j=0; j<nchan; j++)
        {
            if (tfskewness[k*nchan+j]<vmin_skewness or tfskewness[k*nchan+j]>vmax_skewness or tfkurtosis[k*nchan+j]<vmin_kurtosis or tfkurtosis[k*nchan+j]>vmax_kurtosis)
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

    vector<float> tskewness(nsubint, 0.);
    vector<float> tskewness_sort(nsubint, 0.);
    vector<float> tkurtosis(nsubint, 0.);
    vector<float> tkurtosis_sort(nsubint, 0.);

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
       
        double tmp_skewness = 0.;
        double tmp_kurtosis = 0.;

        get_skewness_kurtosis<std::vector<float>::iterator>(tmppro.begin(), nbin, tmp_skewness, tmp_kurtosis);

        tskewness_sort[k] = tskewness[k] = tmp_skewness;
        tkurtosis_sort[k] = tkurtosis[k] = tmp_kurtosis;
    }

    std::nth_element(tskewness_sort.begin(), tskewness_sort.begin()+nsubint/4, tskewness_sort.end(), std::less<float>());
    Q1_skewness = tskewness_sort[nsubint/4];
    std::nth_element(tskewness_sort.begin(), tskewness_sort.begin()+nsubint/4, tskewness_sort.end(), std::greater<float>());
    Q3_skewness = tskewness_sort[nsubint/4];

    R_skewness = Q3_skewness-Q1_skewness;
    vmin_skewness = Q1_skewness-clfd_q*R_skewness;
    vmax_skewness = Q3_skewness+clfd_q*R_skewness;

    std::nth_element(tkurtosis_sort.begin(), tkurtosis_sort.begin()+nsubint/4, tkurtosis_sort.end(), std::less<float>());
    Q1_kurtosis = tkurtosis_sort[nsubint/4];
    std::nth_element(tkurtosis_sort.begin(), tkurtosis_sort.begin()+nsubint/4, tkurtosis_sort.end(), std::greater<float>());
    Q3_kurtosis = tkurtosis_sort[nsubint/4];

    R_kurtosis = Q3_kurtosis-Q1_kurtosis;
    vmin_kurtosis = Q1_kurtosis-clfd_q*R_kurtosis;
    vmax_kurtosis = Q3_kurtosis+clfd_q*R_kurtosis;

    for (long int k=0; k<nsubint; k++)
    {
        if (tskewness[k]<vmin_skewness or tskewness[k]>vmax_skewness or tkurtosis[k]<vmin_kurtosis or tkurtosis[k]>vmax_kurtosis)
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

    vector<float> fskewness(nchan, 0.);
    vector<float> fskewness_sort(nchan, 0.);
    vector<float> fkurtosis(nchan, 0.);
    vector<float> fkurtosis_sort(nchan, 0.);

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
        double tmp_skewness = 0.;
        double tmp_kurtosis = 0.;

        get_skewness_kurtosis<std::vector<float>::iterator>(fph.begin()+j*nbin, nbin, tmp_skewness, tmp_kurtosis);

        fskewness_sort[j] = fskewness[j] = tmp_skewness;
        fkurtosis_sort[j] = fkurtosis[j] = tmp_kurtosis;
    }

    std::nth_element(fskewness_sort.begin(), fskewness_sort.begin()+nchan/4, fskewness_sort.end(), std::less<float>());
    Q1_skewness = fskewness_sort[nchan/4];
    std::nth_element(fskewness_sort.begin(), fskewness_sort.begin()+nchan/4, fskewness_sort.end(), std::greater<float>());
    Q3_skewness = fskewness_sort[nchan/4];

    R_skewness = Q3_skewness-Q1_skewness;
    vmin_skewness = Q1_skewness-clfd_q*R_skewness;
    vmax_skewness = Q3_skewness+clfd_q*R_skewness;

    std::nth_element(fkurtosis_sort.begin(), fkurtosis_sort.begin()+nchan/4, fkurtosis_sort.end(), std::less<float>());
    Q1_kurtosis = fkurtosis_sort[nchan/4];
    std::nth_element(fkurtosis_sort.begin(), fkurtosis_sort.begin()+nchan/4, fkurtosis_sort.end(), std::greater<float>());
    Q3_kurtosis = fkurtosis_sort[nchan/4];

    R_kurtosis = Q3_kurtosis-Q1_kurtosis;
    vmin_kurtosis = Q1_kurtosis-clfd_q*R_kurtosis;
    vmax_kurtosis = Q3_kurtosis+clfd_q*R_kurtosis;

    for (long int k=0; k<nsubint; k++)
    {
        for (long int j=0; j<nchan; j++)
        {
            if (fskewness[j]<vmin_skewness or fskewness[j]>vmax_skewness or fkurtosis[j]<vmin_kurtosis or fkurtosis[j]>vmax_kurtosis)
            {
                for (long int i=0; i<nbin; i++)
                {
                    profiles[k*nchan*nbin+j*nbin+i] = 0.;
                }
            }
        }
    }
}

void GridSearch::clfd2()
{
    /**
     * @brief apply to frequency
     * 
     */

    vector<float> fvar(nchan, 0.);
    vector<float> fvar_sort(nchan, 0.);
    vector<float> fmean(nchan, 0.);
    vector<float> fmean_sort(nchan, 0.);

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
        double tmp_var = 0.;
        double tmp_mean = 0.;

        for (long int i=0; i<nbin; i++)
        {
            tmp_mean += fph[j*nbin+i];
            tmp_var += fph[j*nbin+i]*fph[j*nbin+i];
        }

        tmp_mean /= nbin;
        tmp_var /= nbin;
        
        tmp_var -= tmp_mean*tmp_mean;

        fvar_sort[j] = fvar[j] = tmp_var;
        fmean_sort[j] = fmean[j] = tmp_mean;
    }

    std::nth_element(fvar_sort.begin(), fvar_sort.begin()+nchan/4, fvar_sort.end(), std::less<float>());
    float Q1_var = fvar_sort[nchan/4];
    std::nth_element(fvar_sort.begin(), fvar_sort.begin()+nchan/4, fvar_sort.end(), std::greater<float>());
    float Q3_var = fvar_sort[nchan/4];

    float R_var = Q3_var-Q1_var;
    float vmin_var = Q1_var-clfd_q*R_var;
    float vmax_var = Q3_var+clfd_q*R_var;

    std::nth_element(fmean_sort.begin(), fmean_sort.begin()+nchan/4, fmean_sort.end(), std::less<float>());
    float Q1_mean = fmean_sort[nchan/4];
    std::nth_element(fmean_sort.begin(), fmean_sort.begin()+nchan/4, fmean_sort.end(), std::greater<float>());
    float Q3_mean = fmean_sort[nchan/4];

    float R_mean = Q3_mean-Q1_mean;
    float vmin_mean = Q1_mean-clfd_q*R_mean;
    float vmax_mean = Q3_mean+clfd_q*R_mean;

    for (long int k=0; k<nsubint; k++)
    {
        for (long int j=0; j<nchan; j++)
        {
            if (fvar[j]<vmin_var or fvar[j]>vmax_var or fmean[j]<vmin_mean or fmean[j]>vmax_mean)
            {
                for (long int i=0; i<nbin; i++)
                {
                    profiles[k*nchan*nbin+j*nbin+i] = 0.;
                }
            }
        }
    }
}

void GridSearch::clfd3()
{
    /**
     * @brief apply to frequency
     * 
     */

    vector<float> fvar(nchan, 0.);
    vector<float> fvar_sort(nchan, 0.);
    vector<float> fmean(nchan, 0.);
    vector<float> fmean_sort(nchan, 0.);

    vector<float> fvardiff(nchan, 0.);
    vector<float> fvardiff_sort(nchan, 0.);
    vector<float> fmeandiff(nchan, 0.);
    vector<float> fmeandiff_sort(nchan, 0.);

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

    int nchan_corr = std::ceil(std::abs(bandcorr/(frequencies[1]-frequencies[0]))/2);

    vector<float> fvar_ext(nchan+2*nchan_corr, -std::numeric_limits<float>::infinity());
    vector<float> fmean_ext(nchan+2*nchan_corr, -std::numeric_limits<float>::infinity());

    for (long int j=0; j<nchan; j++)
    {
        double tmp_var = 0.;
        double tmp_mean = 0.;

        for (long int i=0; i<nbin; i++)
        {
            tmp_mean += fph[j*nbin+i];
            tmp_var += fph[j*nbin+i]*fph[j*nbin+i];
        }

        tmp_mean /= nbin;
        tmp_var /= nbin;
        
        tmp_var -= tmp_mean*tmp_mean;

        fvar_sort[j] = fvar[j] = tmp_var;
        fmean_sort[j] = fmean[j] = tmp_mean;

        fvar_ext[j+nchan_corr] = tmp_var;
        fmean_ext[j+nchan_corr] = tmp_mean;
    }

    for (long int j=0; j<nchan; j++)
    {
        fvardiff_sort[j] = fvardiff[j] = fvar_ext[j+nchan_corr] - std::max(fvar_ext[j], fvar_ext[j+2*nchan_corr]);
        fmeandiff_sort[j] = fmeandiff[j] = fmean_ext[j+nchan_corr] - std::max(fmean_ext[j], fmean_ext[j+2*nchan_corr]);
    }

    std::nth_element(fvardiff_sort.begin(), fvardiff_sort.begin()+nchan/4, fvardiff_sort.end(), std::less<float>());
    float Q1_vardiff = fvardiff_sort[nchan/4];
    std::nth_element(fvardiff_sort.begin(), fvardiff_sort.begin()+nchan/4, fvardiff_sort.end(), std::greater<float>());
    float Q3_vardiff = fvardiff_sort[nchan/4];

    float R_vardiff = Q3_vardiff-Q1_vardiff;
    float vmin_vardiff = Q1_vardiff-clfd_q*R_vardiff;
    float vmax_vardiff = Q3_vardiff+clfd_q*R_vardiff;

    std::nth_element(fmeandiff_sort.begin(), fmeandiff_sort.begin()+nchan/4, fmeandiff_sort.end(), std::less<float>());
    float Q1_meandiff = fmeandiff_sort[nchan/4];
    std::nth_element(fmeandiff_sort.begin(), fmeandiff_sort.begin()+nchan/4, fmeandiff_sort.end(), std::greater<float>());
    float Q3_meandiff = fmeandiff_sort[nchan/4];

    float R_meandiff = Q3_meandiff-Q1_meandiff;
    float vmin_meandiff = Q1_meandiff-clfd_q*R_meandiff;
    float vmax_meandiff = Q3_meandiff+clfd_q*R_meandiff;

    std::nth_element(fvar_sort.begin(), fvar_sort.begin()+nchan/4, fvar_sort.end(), std::less<float>());
    float Q1_var = fvar_sort[nchan/4];
    std::nth_element(fvar_sort.begin(), fvar_sort.begin()+nchan/4, fvar_sort.end(), std::greater<float>());
    float Q3_var = fvar_sort[nchan/4];

    float R_var = Q3_var-Q1_var;
    float vmin_var = Q1_var-clfd_q*R_var;
    float vmax_var = Q3_var+clfd_q*R_var;

    std::nth_element(fmean_sort.begin(), fmean_sort.begin()+nchan/4, fmean_sort.end(), std::less<float>());
    float Q1_mean = fmean_sort[nchan/4];
    std::nth_element(fmean_sort.begin(), fmean_sort.begin()+nchan/4, fmean_sort.end(), std::greater<float>());
    float Q3_mean = fmean_sort[nchan/4];

    float R_mean = Q3_mean-Q1_mean;
    float vmin_mean = Q1_mean-clfd_q*R_mean;
    float vmax_mean = Q3_mean+clfd_q*R_mean;

    std::vector<float> weights(nchan, 1.);
    for (long int j=0; j<nchan; j++)
    {
        if ((fvardiff[j]>vmax_vardiff or fmeandiff[j]>vmax_meandiff) and
                (fvar[j]<vmin_var or fvar[j]>vmax_var or fmean[j]<vmin_mean or fmean[j]>vmax_mean)
            )
        {
            weights[j] = 0.;
        }
    }

    for (long int k=0; k<nsubint; k++)
    {
        for (long int j=0; j<nchan; j++)
        {
            for (long int i=0; i<nbin; i++)
            {
                profiles[k*nchan*nbin+j*nbin+i] *= weights[j];
            }
        }
    }
}

void GridSearch::zap()
{
    if (zaplist.empty()) return;

    std::vector<float> weights(nchan, 1.);

    for (long int j=0; j<nchan; j++)
    {
        for (auto k=zaplist.begin(); k!=zaplist.end(); ++k)
        {
            if (frequencies[j]>=(*k).first and frequencies[j]<=(*k).second)
            {
                weights[j] = 0.;
            }
        }
    }

    for (long int k=0; k<nsubint; k++)
    {
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
        for (long int j=0; j<nchan; j++)
        {
            for (long int i=0; i<nbin; i++)
            {
                profiles[k*nchan*nbin+j*nbin+i] *= weights[j];
            }
        }
    }
}
