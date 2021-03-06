/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2021-06-01 18:08:09
 * @modify date 2021-06-01 18:08:09
 * @desc [description]
 */

#include <limits>
#include "preprocesslite.h"
#include "utils.h"
#include "dedisperse.h"

void PreprocessLite::prepare(DataBuffer<float> &databuffer)
{
    nsamples = databuffer.nsamples/td;
    nchans = databuffer.nchans/fd;

    resize(nsamples, nchans);
    tsamp = databuffer.tsamp*td;

    fill(frequencies.begin(), frequencies.end(), 0.);
    for (long int j=0; j<nchans; j++)
    {
        for (long int k=0; k<fd; k++)
        {
            frequencies[j] += databuffer.frequencies[j*fd+k];
        }
        frequencies[j] /= fd;
    }

    closable = false;
}

DataBuffer<float> * PreprocessLite::run(DataBuffer<float> &databuffer)
{
    if (closable) open();

    std::vector<float> chkurtosis(databuffer.nchans, 0.), chskewness(databuffer.nchans, 0.), chmean(databuffer.nchans, 0.), chstd(databuffer.nchans, 0.);
    std::vector<double> chmean1(databuffer.nchans, 0.), chmean2(databuffer.nchans, 0.), chmean3(databuffer.nchans, 0.), chmean4(databuffer.nchans, 0.);

    for (long int i=0; i<databuffer.nsamples; i++)
    {        
        for (long int j=0; j<databuffer.nchans; j++)
        {
            double tmp1 = databuffer.buffer[i*databuffer.nchans+j];
            double tmp2 = tmp1*tmp1;
            double tmp3 = tmp2*tmp1;
            double tmp4 = tmp2*tmp2;
            chmean1[j] += tmp1;
            chmean2[j] += tmp2;
            chmean3[j] += tmp3;
            chmean4[j] += tmp4;
        }
    }

    for (long int j=0; j<databuffer.nchans; j++)
    {
        chmean1[j] /= databuffer.nsamples;
        chmean2[j] /= databuffer.nsamples;
        chmean3[j] /= databuffer.nsamples;
        chmean4[j] /= databuffer.nsamples;

        double tmp = chmean1[j]*chmean1[j];

        chmean[j] = chmean1[j];
        chstd[j] = chmean2[j]-tmp;
        
        if (chstd[j] != 0)
        {
            chskewness[j] = chmean3[j]-3.*chmean2[j]*chmean1[j]+2.*tmp*chmean1[j];
            chkurtosis[j] = chmean4[j]-4.*chmean3[j]*chmean1[j]+6.*chmean2[j]*tmp-3.*tmp*tmp;

            chkurtosis[j] /= chstd[j]*chstd[j];
            chkurtosis[j] -= 3.;

            chskewness[j] /= chstd[j]*std::sqrt(chstd[j]);
        }
        else
        {
            chstd[j] = 1.;
            chkurtosis[j] = std::numeric_limits<float>::max();
            chskewness[j] = std::numeric_limits<float>::max();
        }

        chstd[j] = std::sqrt(chstd[j]);
    }

    /* calculate mean and std of chkurtosis and chskewness */
    std::vector<float> kurtosis_sort = chkurtosis;
    std::nth_element(kurtosis_sort.begin(), kurtosis_sort.begin()+kurtosis_sort.size()/4, kurtosis_sort.end(), std::less<float>());
    float kurtosis_q1 = kurtosis_sort[kurtosis_sort.size()/4];
    std::nth_element(kurtosis_sort.begin(), kurtosis_sort.begin()+kurtosis_sort.size()/4, kurtosis_sort.end(), std::greater<float>());
    float kurtosis_q3 =kurtosis_sort[kurtosis_sort.size()/4];
    float kurtosis_R = kurtosis_q3-kurtosis_q1;

    std::vector<float> skewness_sort = chskewness;
    std::nth_element(skewness_sort.begin(), skewness_sort.begin()+skewness_sort.size()/4, skewness_sort.end(), std::less<float>());
    float skewness_q1 = skewness_sort[skewness_sort.size()/4];
    std::nth_element(skewness_sort.begin(), skewness_sort.begin()+skewness_sort.size()/4, skewness_sort.end(), std::greater<float>());
    float skewness_q3 =skewness_sort[skewness_sort.size()/4];
    float skewness_R = skewness_q3-skewness_q1;


    std::vector<float> weights(databuffer.nchans, 0.);
    for (long int j=0; j<databuffer.nchans; j++)
    {
        if (chkurtosis[j]>=kurtosis_q1-thresig*kurtosis_R && chkurtosis[j]<=kurtosis_q3+thresig*kurtosis_R && chskewness[j]>=skewness_q1-thresig*skewness_R && chskewness[j]<=skewness_q3+thresig*skewness_R)
        {
            weights[j] = 1.;
        }
    }

    if (td == 1 && fd == 1)
    {
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
        for (long int i=0; i<nsamples; i++)
        {
            for (long int j=0; j<nchans; j++)
            {
                buffer[i*nchans+j] = weights[j]*(databuffer.buffer[i*databuffer.nchans+j]-chmean[j])/chstd[j];
            }
        }
    }
    else if (fd == 1)
    {
        std::fill(buffer.begin(), buffer.end(), 0);
        for (long int i=0; i<nsamples; i++)
        {
            for (long int m=0; m<td; m++)
            {
                for (long int j=0; j<nchans; j++)
                {
                    buffer[i*nchans+j] += weights[j]*(databuffer.buffer[(i*td+m)*databuffer.nchans+j]-chmean[j])/chstd[j];
                }
            }
        }
    }
    else
    {
        std::fill(buffer.begin(), buffer.end(), 0);
        for (long int i=0; i<nsamples; i++)
        {
            for (long int m=0; m<td; m++)
            {
                for (long int l=0; l<fd; l++)
                {
                    for (long int j=0; j<nchans; j++)
                    {
                        buffer[i*nchans+j] += weights[j*fd+l]*(databuffer.buffer[(i*td+m)*databuffer.nchans+(j*fd+l)]-chmean[j*fd+l])/chstd[j*fd+l];
                    }
                }
            }
        }
    }

    if (td == 1 && fd == 1)
        equalized = true;
    else
        equalized = false;

    counter += nsamples;

    databuffer.isbusy = false;
    isbusy = true;

    if (databuffer.closable) databuffer.close();

    return this;
}