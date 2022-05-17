/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2021-05-30 10:54:12
 * @modify date 2021-05-30 10:54:12
 * @desc [description]
 */

#ifdef __AVX2__
#include "avx2.h"
#endif

#include <iostream>

#include "filterbankwriter.h"
#include <assert.h>

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
	if (fil.nbits == 32)
	{
		outfile.write((char *)(databuffer.buffer.data()), sizeof(float)*databuffer.nsamples*databuffer.nchans);
		return;
	}

	double tmpmean=0., tmpstd=0.;
#ifndef __AVX2__
	for (long int i=0; i<databuffer.nsamples; i++)
	{
		for (long int j=0; j<databuffer.nchans; j++)
		{
			tmpmean += databuffer.buffer[i*databuffer.nchans+j];
			tmpstd += databuffer.buffer[i*databuffer.nchans+j]*databuffer.buffer[i*databuffer.nchans+j];
		}
	}
#else
	if (databuffer.nsamples*databuffer.nchans % 4 == 0)
	{
		PulsarX::accumulate_mean_var2(tmpmean, tmpstd, databuffer.buffer.data(), databuffer.nsamples*databuffer.nchans);
	}
	else
	{
		for (long int i=0; i<databuffer.nsamples; i++)
		{
			for (long int j=0; j<databuffer.nchans; j++)
			{
				tmpmean += databuffer.buffer[i*databuffer.nchans+j];
				tmpstd += databuffer.buffer[i*databuffer.nchans+j]*databuffer.buffer[i*databuffer.nchans+j];
			}
		}
	}
#endif
	tmpmean /= databuffer.nsamples*databuffer.nchans;
	tmpstd /= databuffer.nsamples*databuffer.nchans;
	tmpstd -= tmpmean*tmpmean;
	tmpstd = std::sqrt(tmpstd);

	float scl = outstd/tmpstd;
	float offs = outmean-scl*tmpmean;

	// write data to filterbank
	switch (fil.nbits)
	{
	case 8:
	{
		std::vector<unsigned char> data8bit(databuffer.nsamples*databuffer.nchans, 0);
#ifndef __AVX2__
		for (long int i=0; i<databuffer.nsamples; i++)
		{
			for (long int j=0; j<databuffer.nchans; j++)
			{
				float tmp = std::round(scl*databuffer.buffer[i*databuffer.nchans+j]+offs);
				tmp = std::max(0.f, tmp);
				tmp = std::min(255.f, tmp);
				data8bit[i*databuffer.nchans+j] = tmp;
			}
		}
#else
		if (databuffer.nsamples*databuffer.nchans % 8 == 0)
		{
			PulsarX::scale(data8bit.data(), databuffer.buffer.data(), scl, offs, 8, databuffer.nsamples*databuffer.nchans);
		}
		else
		{
			for (long int i=0; i<databuffer.nsamples; i++)
			{
				for (long int j=0; j<databuffer.nchans; j++)
				{
					float tmp = std::round(scl*databuffer.buffer[i*databuffer.nchans+j]+offs);
					tmp = std::max(0.f, tmp);
					tmp = std::min(255.f, tmp);
					data8bit[i*databuffer.nchans+j] = tmp;
				}
			}
		}
#endif

		outfile.write((char *)(data8bit.data()), sizeof(unsigned char)*databuffer.nsamples*databuffer.nchans);
	};break;
	case 4:
	{
		assert(databuffer.nchans % 2 == 0);

		std::vector<unsigned char> data4bit(databuffer.nsamples*databuffer.nchans/2, 0);
#ifndef __AVX2__
		for (long int i=0; i<databuffer.nsamples; i++)
		{
			for (long int j=0; j<databuffer.nchans/2; j++)
			{
				float tmp0 = std::round(scl*databuffer.buffer[i*databuffer.nchans+j*2+0]+offs);
				float tmp1 = std::round(scl*databuffer.buffer[i*databuffer.nchans+j*2+1]+offs);
				tmp0 = std::max(0.f, tmp0);
				tmp0 = std::min(15.f, tmp0);
				tmp1 = std::max(0.f, tmp1);
				tmp1 = std::min(15.f, tmp1);

				data4bit[i*databuffer.nchans/2+j] = ((unsigned char)tmp0) + (((unsigned char)tmp1) << 4);
			}
		}
#else
		if (databuffer.nsamples*databuffer.nchans % 8 == 0)
		{
			PulsarX::scale(data4bit.data(), databuffer.buffer.data(), scl, offs, 4, databuffer.nsamples*databuffer.nchans);
		}
		else
		{
			for (long int i=0; i<databuffer.nsamples; i++)
			{
				for (long int j=0; j<databuffer.nchans/2; j++)
				{
					float tmp0 = std::round(scl*databuffer.buffer[i*databuffer.nchans+j*2+0]+offs);
					float tmp1 = std::round(scl*databuffer.buffer[i*databuffer.nchans+j*2+1]+offs);
					tmp0 = std::max(0.f, tmp0);
					tmp0 = std::min(15.f, tmp0);
					tmp1 = std::max(0.f, tmp1);
					tmp1 = std::min(15.f, tmp1);

					data4bit[i*databuffer.nchans/2+j] = ((unsigned char)tmp0) + (((unsigned char)tmp1) << 4);
				}
			}
		}
#endif

		outfile.write((char *)(data4bit.data()), sizeof(unsigned char)*databuffer.nsamples*databuffer.nchans/2);

	};break;
	case 2:
	{
		assert(databuffer.nchans % 4 == 0);

		std::vector<unsigned char> data2bit(databuffer.nsamples*databuffer.nchans/4, 0);
#ifndef __AVX2__
		for (long int i=0; i<databuffer.nsamples; i++)
		{
			for (long int j=0; j<databuffer.nchans/4; j++)
			{
				float tmp0 = std::round(scl*databuffer.buffer[i*databuffer.nchans+j*4+0]+offs);
				float tmp1 = std::round(scl*databuffer.buffer[i*databuffer.nchans+j*4+1]+offs);
				float tmp2 = std::round(scl*databuffer.buffer[i*databuffer.nchans+j*4+2]+offs);
				float tmp3 = std::round(scl*databuffer.buffer[i*databuffer.nchans+j*4+3]+offs);
				tmp0 = std::max(0.f, tmp0);
				tmp0 = std::min(3.f, tmp0);
				tmp1 = std::max(0.f, tmp1);
				tmp1 = std::min(3.f, tmp1);
				tmp2 = std::max(0.f, tmp2);
				tmp2 = std::min(3.f, tmp2);
				tmp3 = std::max(0.f, tmp3);
				tmp3 = std::min(3.f, tmp3);

				data2bit[i*databuffer.nchans/4+j] = ((unsigned char)tmp0) + (((unsigned char)tmp1) << 2) + (((unsigned char)tmp2) << 4) + (((unsigned char)tmp3) << 6);
			}
		}
#else
		if (databuffer.nsamples*databuffer.nchans % 8 == 0)
		{
			PulsarX::scale(data2bit.data(), databuffer.buffer.data(), scl, offs, 2, databuffer.nsamples*databuffer.nchans);
		}
		else
		{
			for (long int i=0; i<databuffer.nsamples; i++)
			{
				for (long int j=0; j<databuffer.nchans/4; j++)
				{
					float tmp0 = std::round(scl*databuffer.buffer[i*databuffer.nchans+j*4+0]+offs);
					float tmp1 = std::round(scl*databuffer.buffer[i*databuffer.nchans+j*4+1]+offs);
					float tmp2 = std::round(scl*databuffer.buffer[i*databuffer.nchans+j*4+2]+offs);
					float tmp3 = std::round(scl*databuffer.buffer[i*databuffer.nchans+j*4+3]+offs);
					tmp0 = std::max(0.f, tmp0);
					tmp0 = std::min(3.f, tmp0);
					tmp1 = std::max(0.f, tmp1);
					tmp1 = std::min(3.f, tmp1);
					tmp2 = std::max(0.f, tmp2);
					tmp2 = std::min(3.f, tmp2);
					tmp3 = std::max(0.f, tmp3);
					tmp3 = std::min(3.f, tmp3);

					data2bit[i*databuffer.nchans/4+j] = ((unsigned char)tmp0) + (((unsigned char)tmp1) << 2) + (((unsigned char)tmp2) << 4) + (((unsigned char)tmp3) << 6);
				}
			}
		}
#endif

		outfile.write((char *)(data2bit.data()), sizeof(unsigned char)*databuffer.nsamples*databuffer.nchans/4);

	};break;
	case 1:
	{
		assert(databuffer.nsamples*databuffer.nchans % 8 == 0);

		std::vector<unsigned char> data1bit(databuffer.nsamples*databuffer.nchans/8, 0);
#ifndef __AVX2__
		for (long int i=0; i<databuffer.nsamples; i++)
		{
			for (long int j=0; j<databuffer.nchans/8; j++)
			{
				float tmp0 = std::round(scl*databuffer.buffer[i*databuffer.nchans+j*8+0]+offs);
				float tmp1 = std::round(scl*databuffer.buffer[i*databuffer.nchans+j*8+1]+offs);
				float tmp2 = std::round(scl*databuffer.buffer[i*databuffer.nchans+j*8+2]+offs);
				float tmp3 = std::round(scl*databuffer.buffer[i*databuffer.nchans+j*8+3]+offs);
				float tmp4 = std::round(scl*databuffer.buffer[i*databuffer.nchans+j*8+4]+offs);
				float tmp5 = std::round(scl*databuffer.buffer[i*databuffer.nchans+j*8+5]+offs);
				float tmp6 = std::round(scl*databuffer.buffer[i*databuffer.nchans+j*8+6]+offs);
				float tmp7 = std::round(scl*databuffer.buffer[i*databuffer.nchans+j*8+7]+offs);
				tmp0 = std::max(0.f, tmp0);
				tmp0 = std::min(1.f, tmp0);
				tmp1 = std::max(0.f, tmp1);
				tmp1 = std::min(1.f, tmp1);
				tmp2 = std::max(0.f, tmp2);
				tmp2 = std::min(1.f, tmp2);
				tmp3 = std::max(0.f, tmp3);
				tmp3 = std::min(1.f, tmp3);
				tmp4 = std::max(0.f, tmp4);
				tmp4 = std::min(1.f, tmp4);
				tmp5 = std::max(0.f, tmp5);
				tmp5 = std::min(1.f, tmp5);
				tmp6 = std::max(0.f, tmp6);
				tmp6 = std::min(1.f, tmp6);
				tmp7 = std::max(0.f, tmp7);
				tmp7 = std::min(1.f, tmp7);

				data1bit[i*databuffer.nchans/8+j] = ((unsigned char)tmp0) + (((unsigned char)tmp1) << 1) + (((unsigned char)tmp2) << 2) + (((unsigned char)tmp3) << 3) + (((unsigned char)tmp4) << 4) + (((unsigned char)tmp5) << 5) + (((unsigned char)tmp6) << 6) + (((unsigned char)tmp7) << 7);
			}
		}
#else
		if (databuffer.nchans % 8 == 0)
		{
			PulsarX::scale(data1bit.data(), databuffer.buffer.data(), scl, offs, 1, databuffer.nsamples*databuffer.nchans);
		}
		else
		{
			for (long int i=0; i<databuffer.nsamples; i++)
			{
				for (long int j=0; j<databuffer.nchans/8; j++)
				{
					float tmp0 = std::round(scl*databuffer.buffer[i*databuffer.nchans+j*8+0]+offs);
					float tmp1 = std::round(scl*databuffer.buffer[i*databuffer.nchans+j*8+1]+offs);
					float tmp2 = std::round(scl*databuffer.buffer[i*databuffer.nchans+j*8+2]+offs);
					float tmp3 = std::round(scl*databuffer.buffer[i*databuffer.nchans+j*8+3]+offs);
					float tmp4 = std::round(scl*databuffer.buffer[i*databuffer.nchans+j*8+4]+offs);
					float tmp5 = std::round(scl*databuffer.buffer[i*databuffer.nchans+j*8+5]+offs);
					float tmp6 = std::round(scl*databuffer.buffer[i*databuffer.nchans+j*8+6]+offs);
					float tmp7 = std::round(scl*databuffer.buffer[i*databuffer.nchans+j*8+7]+offs);
					tmp0 = std::max(0.f, tmp0);
					tmp0 = std::min(1.f, tmp0);
					tmp1 = std::max(0.f, tmp1);
					tmp1 = std::min(1.f, tmp1);
					tmp2 = std::max(0.f, tmp2);
					tmp2 = std::min(1.f, tmp2);
					tmp3 = std::max(0.f, tmp3);
					tmp3 = std::min(1.f, tmp3);
					tmp4 = std::max(0.f, tmp4);
					tmp4 = std::min(1.f, tmp4);
					tmp5 = std::max(0.f, tmp5);
					tmp5 = std::min(1.f, tmp5);
					tmp6 = std::max(0.f, tmp6);
					tmp6 = std::min(1.f, tmp6);
					tmp7 = std::max(0.f, tmp7);
					tmp7 = std::min(1.f, tmp7);

					data1bit[i*databuffer.nchans/8+j] = ((unsigned char)tmp0) + (((unsigned char)tmp1) << 1) + (((unsigned char)tmp2) << 2) + (((unsigned char)tmp3) << 3) + (((unsigned char)tmp4) << 4) + (((unsigned char)tmp5) << 5) + (((unsigned char)tmp6) << 6) + (((unsigned char)tmp7) << 7);
				}
			}
		}
#endif

		outfile.write((char *)(data1bit.data()), sizeof(unsigned char)*databuffer.nsamples*databuffer.nchans/8);

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