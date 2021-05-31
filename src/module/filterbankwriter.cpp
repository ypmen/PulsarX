/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2021-05-30 10:54:12
 * @modify date 2021-05-30 10:54:12
 * @desc [description]
 */

#include <iostream>

#include "filterbankwriter.h"

FilterbankWriter::FilterbankWriter()
{
    outstd = 3.;
    outmean = 10.;
}

FilterbankWriter::FilterbankWriter(const Filterbank &temp)
{
    outstd = 3.;
    outmean = 10.;
    
    fil = temp;
}

FilterbankWriter::FilterbankWriter(const FilterbankWriter &filwriter) : fil(filwriter.fil)
{
    outstd = filwriter.outstd;
    outmean = filwriter.outmean;
}

FilterbankWriter & FilterbankWriter::operator=(const FilterbankWriter &filwriter)
{
    outstd = filwriter.outstd;
    outmean = filwriter.outmean;

    fil = filwriter.fil;

    return *this;
}

FilterbankWriter::~FilterbankWriter()
{
    outfile.close();
}

void FilterbankWriter::prepare(DataBuffer<float> &databuffer)
{
    fil.nifs = 1;
    fil.nchans = databuffer.nchans;
    fil.tsamp = databuffer.tsamp;
    fil.fch1 = databuffer.frequencies.front();
    fil.foff = (databuffer.frequencies.back()-databuffer.frequencies.front())/(databuffer.nchans-1);

    if (!fil.write_header())
    {
        std::cerr<<"Error: Can not write filterbank header"<<std::endl;
        exit(-1);
    }
    fil.close();

    outfile.open(fil.filename, std::ios::binary|std::ios::app);
}

void FilterbankWriter::run(DataBuffer<float> &databuffer)
{
    // write data to filterbank
    switch (fil.nbits)
    {
    case 8:
    {
        double tmpmean=0., tmpstd=0.;
        for (long int i=0; i<databuffer.nsamples; i++)
        {
            for (long int j=0; j<databuffer.nchans; j++)
            {
                tmpmean += databuffer.buffer[i*databuffer.nchans+j];
                tmpstd += databuffer.buffer[i*databuffer.nchans+j]*databuffer.buffer[i*databuffer.nchans+j];
            }
        }
        tmpmean /= databuffer.nsamples*databuffer.nchans;
        tmpstd /= databuffer.nsamples*databuffer.nchans;
        tmpstd -= tmpmean*tmpmean;
        tmpstd = std::sqrt(tmpstd);

        float scl = outstd/tmpstd;
        float offs = outmean-scl*tmpmean;

        std::vector<unsigned char> data8bit(databuffer.nsamples*databuffer.nchans, 0);
        for (long int i=0; i<databuffer.nsamples; i++)
        {
            for (long int j=0; j<databuffer.nchans; j++)
            {
                float tmp = scl*databuffer.buffer[i*databuffer.nchans+j]+offs;
                tmp = std::max(0.f, tmp);
                tmp = std::min(255.f, tmp);
                data8bit[i*databuffer.nchans+j] = tmp;
            }
        }

        outfile.write((char *)(data8bit.data()), sizeof(unsigned char)*databuffer.nsamples*databuffer.nchans);

    };break;
    case 32:
    {
        outfile.write((char *)(databuffer.buffer.data()), sizeof(float)*databuffer.nsamples*databuffer.nchans);
    };break;
    default:
    {
        std::cerr<<"Error: data type unsupported"<<endl;
        exit(-1);
    }; break;
    }

    databuffer.isbusy = false;
    if (databuffer.closable) databuffer.close();
}