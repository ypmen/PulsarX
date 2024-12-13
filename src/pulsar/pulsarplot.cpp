/**
 * @author Yunpeng Men
 * @email ypmen@pku.edu.cn
 * @create date 2020-06-14 15:59:59
 * @modify date 2020-06-14 15:59:59
 * @desc [description]
 */

#include <fstream>
#include <string.h>
#include <sstream>
#include <cmath>

#include "config.h"

#include "constants.h"
#include "pulsarplot.h"
#ifdef __AVX2__
#include "avx_mathfun.h"
#endif

#include "plotx.h"

using namespace std;
using namespace Pulsar;

PulsarPlot::PulsarPlot(){}

PulsarPlot::~PulsarPlot(){}

void PulsarPlot::plot(const ArchiveLite &archive, GridSearch &gridsearch, std::map<std::string, std::string> &obsinfo, int id, const string &rootname, bool plotx, bool save2fits)
{
	string basename;
	if (obsinfo.find("Basename") == obsinfo.end())
	{
		stringstream ss_id;
		ss_id << setw(5) << setfill('0') << id;
		string s_id = ss_id.str();

		basename = rootname + "_" + obsinfo["Date"] + "_" + obsinfo["Beam"] + "_" + s_id;
	}
	else
	{
		basename = obsinfo["Basename"];
	}
	string figname = basename + ".png";

	long int nsubint = gridsearch.nsubint;
	long int nchan = gridsearch.nchan;
	long int nbin = gridsearch.nbin;

	long int ndm = gridsearch.nddm;
	long int nf1 = gridsearch.ndf1;
	long int nf0 = gridsearch.ndf0;

	vector<float> vp, vph, vt, vf, vdm, vsnr_dm, vdf0, vsnr_f0, vdf1, vsnr_f1;
	vp.resize(nbin*2, 0.);
	vph.resize(nbin*2, 0.);
	vt.resize(nsubint);
	vf.resize(nchan);
	vdm.resize(ndm);
	vsnr_dm.resize(ndm);
	vdf0.resize(nf0);
	vsnr_f0.resize(nf0);
	vdf1.resize(nf1);
	vsnr_f1.resize(nf1);

	vector<float> mxtph, mxfph, mxsnr_ffdot;
	mxtph.resize(nsubint*nbin*2, 0.);
	mxfph.resize(nchan*nbin*2, 0.);
	mxsnr_ffdot.resize(nf1*nf0, 0.);

	double snr = gridsearch.snr;
	double width = gridsearch.width;
	double f0 = gridsearch.f0;
	double f1 = gridsearch.f1;
	double dm = gridsearch.dm;
	double err_f0 = gridsearch.err_f0;
	double err_f1 = gridsearch.err_f1;
	double err_dm = gridsearch.err_dm;
	double p0 = gridsearch.p0;
	double p1 = gridsearch.p1;
	double err_p0 = gridsearch.err_p0;
	double err_p1 = gridsearch.err_p1;
	double acc = gridsearch.acc;
	double err_acc = gridsearch.err_acc;
	std::vector<double> pulsespan(2, 0.);
	if (gridsearch.pulsespan.back() > gridsearch.pulsespan.front())
	{
		pulsespan[0] = gridsearch.pulsespan.front()*1./nbin;
		pulsespan[1] = gridsearch.pulsespan.back()*1./nbin;
	}
	else
	{
		pulsespan[0] = gridsearch.pulsespan.front()*1./nbin;
		pulsespan[1] = (gridsearch.pulsespan.back()+nbin)*1./nbin;
	}

	for (long int i=0; i<2*nbin; i++)
	{
		vp[i] = gridsearch.profile[i%nbin];
		vph[i] = i*1./nbin;
	}

	for (long int k=0; k<nsubint; k++)
	{
		vt[k] = gridsearch.tsuboff[k]-gridsearch.tsuboff[0];
		for (long int i=0; i<nbin*2; i++)
		{
			mxtph[k*nbin*2+i] = 0.;
		}
	}

	for (long int j=0; j<nchan; j++)
	{
		vf[j] = gridsearch.frequencies[j];
		for (long int i=0; i<nbin*2; i++)
		{
			mxfph[j*nbin*2+i] = 0.;
		}
	}

	for (long int k=0; k<nsubint; k++)
	{
		for (long int j=0; j<nchan; j++)
		{
			for (long int i=0; i<nbin; i++)
			{
				mxfph[j*nbin*2+i] += gridsearch.profiles[k*nchan*nbin+j*nbin+i];
				mxfph[j*nbin*2+i+nbin] += gridsearch.profiles[k*nchan*nbin+j*nbin+i];

				mxtph[k*nbin*2+i] += gridsearch.profiles[k*nchan*nbin+j*nbin+i];
				mxtph[k*nbin*2+i+nbin] += gridsearch.profiles[k*nchan*nbin+j*nbin+i];
			}
		}
	}

	std::vector<float> vddm(ndm, 0.);
	for (long int k=0; k<ndm; k++)
	{
		vddm[k] = gridsearch.ddmstart + k*gridsearch.ddmstep;
		vdm[k] = gridsearch.dm + gridsearch.ddmstart + k*gridsearch.ddmstep;
		if (gridsearch.dmsearch)
			vsnr_dm[k] = gridsearch.vsnr_dm[k];
		else
			vsnr_dm[k] = 0.;
	}

	std::vector<float> vchisq_dm(ndm, 0.);
	get_dm_chisq_curve(vchisq_dm, vddm, gridsearch.frequencies, f0, width*f0, nbin);
	float a=0., b=0.;
	get_bestfit(a, b, vsnr_dm, vchisq_dm);
	for (long int k=0; k<ndm; k++)
	{
		vchisq_dm[k] = a*vchisq_dm[k]+b;
	}

	int if1 = nf1/2;
	int if0 = nf0/2;
	double maxsnr = -1.;
	for (long int k0=0; k0<nf0; k0++)
	{
		vdf0[k0] = gridsearch.df0start + k0*gridsearch.df0step;
		if (std::abs(vdf0[k0])<std::abs(0.5*gridsearch.df0step))
			if0 = k0;
		vsnr_f0[k0] = 0.;
	}

	for (long int k1=0; k1<nf1; k1++)
	{
		vdf1[k1] = gridsearch.df1start + k1*gridsearch.df1step;
		if (std::abs(vdf1[k1])<std::abs(0.5*gridsearch.df1step))
			if1 = k1;
		vsnr_f1[k1] = 0.;
	}

	for (long int k1=0; k1<nf1; k1++)
	{
		for (long int k0=0; k0<nf0; k0++)
		{
			if (gridsearch.ffdotsearch)
				mxsnr_ffdot[k1*nf0+k0] = gridsearch.mxsnr_ffdot[k1*nf0+k0];
			else
				mxsnr_ffdot[k1*nf0+k0] = 0.;
			
			//vsnr_f0[k0] += mxsnr_ffdot[k1*nf0+k0];
			//vsnr_f1[k1] += mxsnr_ffdot[k1*nf0+k0];

			// if (mxsnr_ffdot[k1*nf0+k0] > maxsnr)
			// {
			//     maxsnr = mxsnr_ffdot[k1*nf0+k0];
			//     if1 = k1;
			//     if0 = k0;
			// }
		}
	}

	maxsnr = mxsnr_ffdot[if1*nf0+if0];

	for (long int k0=0; k0<nf0; k0++)
	{
		//vsnr_f0[k0] /= nf1;
		vsnr_f0[k0] = mxsnr_ffdot[if1*nf0+k0];
	}

	std::vector<float> vchisq_f0(nf0, 0.);
	get_f0_chisq_curve(vchisq_f0, vdf0, gridsearch.tsuboff, f0, width*f0, nbin);
	get_bestfit(a, b, vsnr_f0, vchisq_f0);
	for (long int k=0; k<nf0; k++)
	{
		vchisq_f0[k] = a*vchisq_f0[k]+b;
	}

	for (long int k1=0; k1<nf1; k1++)
	{
		//vsnr_f1[k1] /= nf0;
		vsnr_f1[k1] = mxsnr_ffdot[k1*nf0+if0];
	}

	std::vector<float> vchisq_f1(nf1, 0.);
	get_f1_chisq_curve(vchisq_f1, vdf1, gridsearch.tsuboff, f0, width*f0, nbin);
	get_bestfit(a, b, vsnr_f1, vchisq_f1);
	for (long int k=0; k<nf1; k++)
	{
		vchisq_f1[k] = a*vchisq_f1[k]+b;
	}

	/**
	 * @brief format print the parameters
	 * 
	 */
	std::string s_f0, s_f1, s_p0, s_p1, s_acc, s_dm, s_snr;
	format_val_err(s_f0, f0, err_f0);
	format_val_err(s_f1, f1, err_f1, "sci");
	format_val_err(s_p0, p0, err_p0);
	format_val_err(s_p1, p1, err_p1, "sci");
	format_val_err(s_acc, acc, err_acc);
	format_val_err(s_dm, dm, err_dm);

	std:stringstream ss_snr;
	ss_snr<<fixed<<setprecision(2)<<snr;
	s_snr = ss_snr.str();

	std::string s_gl, s_gb;
#ifdef HAVE_SOFA
	s_gl = obsinfo["GL"];
	s_gb = obsinfo["GB"];
#endif

	std::string s_ymw16_maxdm, s_ymw16_dist;

	std::stringstream ss_ymw16_maxdm;
	ss_ymw16_maxdm<<fixed<<setprecision(1)<<stod(obsinfo["MaxDM_YMW16"]);
	s_ymw16_maxdm = ss_ymw16_maxdm.str();

	std::stringstream ss_ymw16_dist;
	ss_ymw16_dist<<fixed<<setprecision(1)<<stod(obsinfo["Dist_YMW16"]);
	s_ymw16_dist = ss_ymw16_dist.str();

	//DM smearing
	double fcentre = std::stod(obsinfo["Fcentre"]);
	double bandwidth = std::stod(obsinfo["Bandwidth"]);
	int nchans = std::stoi(obsinfo["Nchan"]);
	double dmsmear_phase = abs(DedispersionLite::dmdelay(dm, fcentre-0.5*bandwidth, fcentre+0.5*bandwidth))/nchans*f0;

	std::string fontsize_label = "0.7";
	std::string fontsize_ticklabel = "0.7";

	/** plot */
	PlotX::Figure fig(11.75, 1.);
	fig.set_background_color("black");
	fig.set_default_color("white");

	//profile
	float smear_left = 1.-0.5*dmsmear_phase;
	float smear_right = 1.+0.5*dmsmear_phase;
	if (smear_left < 0.) smear_left = 0;
	if (smear_right > 2.) smear_right = 2.;
	PlotX::Axes ax1(0.08, 0.38, 0.76, 0.96);
	ax1.axvspan(smear_left, smear_right, 0., 1., {{"color", "lightgray"}, {"filled", "true"}});
	ax1.axhline(-0.5, 0., 1., {{"color", "red"},{"linestyle", "--"}});
	ax1.axhline(0.5, 0., 1., {{"color", "red"},{"linestyle", "--"}});
	ax1.axvline(pulsespan[0], 0., 1., {{"color", "red"},{"linestyle", "--"}});
	ax1.axvline(pulsespan[1], 0., 1., {{"color", "red"},{"linestyle", "--"}});
	ax1.plot(vph, vp);
	ax1.autoscale(true, "x", true);
	ax1.set_ylabel("Flux");
	ax1.set_fontsize_label(fontsize_label);
	ax1.set_fontsize_ticklabel(fontsize_ticklabel);
	ax1.label(true, false, false, false);
	fig.push(ax1);

	//dynamic image
	PlotX::Axes ax2_twin(0.08, 0.38, 0.42, 0.76);
	ax2_twin.pcolor(vph.size(), vf.size(), mxfph, "viridis");
	ax2_twin.label(false, true, false, false);
	ax2_twin.set_fontsize_label(fontsize_label);
	ax2_twin.set_fontsize_ticklabel(fontsize_ticklabel);
	fig.push(ax2_twin);

	PlotX::Axes ax2(0.08, 0.38, 0.42, 0.76);
	ax2.pcolor(vph, vf, mxfph, "viridis");
	ax2.set_ylabel("Frequency (MHz)");
	ax2.set_fontsize_label(fontsize_label);
	ax2.set_fontsize_ticklabel(fontsize_ticklabel);
	ax2.label(true, false, false, false);
	fig.push(ax2);

	//subint image
	PlotX::Axes ax3_twin(0.08, 0.38, 0.08, 0.42);
	ax3_twin.pcolor(vph.size(), vt.size(), mxtph, "viridis");
	ax3_twin.label(false, true, false, false);
	ax3_twin.set_fontsize_label(fontsize_label);
	ax3_twin.set_fontsize_ticklabel(fontsize_ticklabel);
	fig.push(ax3_twin);

	PlotX::Axes ax3(0.08, 0.38, 0.08, 0.42);
	ax3.pcolor(vph, vt, mxtph, "viridis");
	ax3.set_xlabel("Phase");
	ax3.set_ylabel("Tint (s)");
	ax3.label(true, false, true, false);
	ax3.set_fontsize_label(fontsize_label);
	ax3.set_fontsize_ticklabel(fontsize_ticklabel);
	fig.push(ax3);

	double xpos = archive.f0-gridsearch.f0;
	double ypos = archive.f1-gridsearch.f1;
	double dmpos = archive.dm;
	xpos = xpos <vdf0[0]? vdf0[0]: xpos;
	xpos = xpos >vdf0.back()? vdf0.back(): xpos;
	ypos = ypos <vdf1[0]? vdf1[0]: ypos;
	ypos = ypos >vdf1.back()? vdf1.back(): ypos;
	dmpos = dmpos <vdm[0]? vdm[0]:dmpos;
	dmpos = dmpos >vdm.back()? vdm.back():dmpos;
	//chi2-f0
	PlotX::Axes ax4(0.48, 0.78, 0.6, 0.72);
	ax4.plot(vdf0, vsnr_f0);
	ax4.plot(vdf0, vchisq_f0, {{"color", "yellow"}});
	ax4.axvline(xpos, 0, 1, {{"color", "red"}});
	ax4.annotate("P0 (s) = "+s_p0, 0.1, 1.1, {{"xycoords","axes fraction"}, {"fontsize", "0.7"}});
	ax4.set_xlabel("F0 - " + to_string(gridsearch.f0) + " (Hz)");
	ax4.set_ylabel("\\gx\\u2");
	ax4.set_fontsize_label(fontsize_label);
	ax4.set_fontsize_ticklabel(fontsize_ticklabel);
	ax4.autoscale(true, "x", true);
	fig.push(ax4);

	//chi2-f0-f1
	PlotX::Axes ax5(0.48, 0.78, 0.24, 0.52);
	ax5.pcolor(vdf0, vdf1, mxsnr_ffdot, "viridis");
	ax5.cross(xpos, ypos, 4.);
	ax5.label(false, false, false, false);
	fig.push(ax5);

	//chi2-f1
	std::stringstream ss_f1;
	ss_f1<<fixed<<setprecision(5)<<scientific<<f1;

	PlotX::Axes ax6(0.86, 0.98, 0.24, 0.52);
	ax6.plot(vsnr_f1, vdf1);
	ax6.plot(vchisq_f1, vdf1, {{"color", "yellow"}});
	ax6.axhline(ypos, 0, 1, {{"color", "red"}});
	ax6.annotate("P1 (s/s) = "+s_p1, -0.8, 1.1, {{"xycoords","axes fraction"}, {"fontsize", "0.7"}});
	ax6.set_xlabel("\\gx\\u2");
	ax6.set_ylabel("F1 -" + ss_f1.str() + " (Hz/s)");
	ax6.set_fontsize_label(fontsize_label);
	ax6.set_fontsize_ticklabel(fontsize_ticklabel);
	ax6.autoscale(true, "y", true);
	ax6.label(true, false, true, false);
	fig.push(ax6);

	//chi2-dm
	PlotX::Axes ax7(0.48, 0.78, 0.08, 0.2);
	ax7.plot(vdm, vsnr_dm);
	ax7.plot(vdm, vchisq_dm, {{"color", "yellow"}});
	ax7.axvline(dmpos, 0, 1, {{"color", "red"}});
	ax7.annotate("DM (pc/cc) = "+s_dm, 0.1, 1.1, {{"xycoords","axes fraction"}, {"fontsize", "0.7"}});
	ax7.set_xlabel("DM (pc/cc)");
	ax7.set_ylabel("\\gx\\u2");
	ax7.set_fontsize_label(fontsize_label);
	ax7.set_fontsize_ticklabel(fontsize_ticklabel);
	ax7.autoscale(true, "x", true);
	fig.push(ax7);

	std::string fontsize = "0.6"; 
	//metadata
	PlotX::Axes ax8(0.48, 0.98, 0.76, 0.96);
	ax8.label(false, false, false, false);
	ax8.majorticks_off();
	ax8.majorticks_off();

	ax8.annotate("Telescope = " + obsinfo["Telescope"], 0.02, 0.91, {{"fontsize", fontsize}});
	ax8.annotate("Beam = " + obsinfo["Beam"], 0.02, 0.80, {{"fontsize", fontsize}});
	ax8.annotate("RA = " + obsinfo["RA"], 0.02, 0.69, {{"fontsize", fontsize}});
	ax8.annotate("DEC = " + obsinfo["DEC"], 0.02, 0.58, {{"fontsize", fontsize}});
	ax8.annotate("P0 (s) = " + s_p0, 0.02, 0.47, {{"fontsize", fontsize}});
	ax8.annotate("P1 (s/s) = " + s_p1, 0.02, 0.36, {{"fontsize", fontsize}});
	ax8.annotate("DM (pc/cc) = " + s_dm, 0.02, 0.25, {{"fontsize", fontsize}});
	ax8.annotate("acc (m/s/s) = " + s_acc, 0.02, 0.14, {{"fontsize", fontsize}});
	ax8.annotate("S/N = " + s_snr, 0.02, 0.03, {{"fontsize", fontsize}});
	ax8.annotate("Source name = " + obsinfo["Source_name"], 0.46, 0.91, {{"fontsize", fontsize}});
	ax8.annotate("Date (MJD) = " + obsinfo["Date"], 0.46, 0.80, {{"fontsize", fontsize}});
	ax8.annotate("GL (deg) = " + s_gl, 0.46, 0.69, {{"fontsize", fontsize}});
	ax8.annotate("GB (deg) = " + s_gb, 0.46, 0.58, {{"fontsize", fontsize}});
	ax8.annotate("F0 (Hz) = " + s_f0, 0.46, 0.47, {{"fontsize", fontsize}});
	ax8.annotate("F1 (Hz/s) = " + s_f1, 0.46, 0.36, {{"fontsize", fontsize}});
	ax8.annotate("MaxDM YMW16  (pc/cc) = " + s_ymw16_maxdm, 0.46, 0.25, {{"fontsize", fontsize}});
	ax8.annotate("Distance YMW16 (pc) = " + s_ymw16_dist, 0.46, 0.14, {{"fontsize", fontsize}});
	ax8.annotate("Pepoch (MJD) = " + obsinfo["Pepoch"], 0.46, 0.03, {{"fontsize", fontsize}});

	ax8.annotate(obsinfo["Filename"], 0.07, 0.97, {{"xycoords","figure fraction"}, {"fontsize", "0.45"}});
	fig.push(ax8);

	fig.save(figname+"/PNG");

	if (save2fits)
	{
		fig.savepx(basename + ".px");
	}
}

void PulsarPlot::get_dm_chisq_curve(std::vector<float> &vchisq, const std::vector<float> &vddm, const std::vector<double> &frequencies, double f0, double boxphw, int nbin)
{
	int ndm = vddm.size();
	int nchan = frequencies.size();
	double fc = 0.5*(*std::min_element(frequencies.begin(), frequencies.end())+*std::max_element(frequencies.begin(), frequencies.end()));

	vchisq.resize(ndm, 0.);
	for (long int k=0; k<ndm; k++)
	{
		float ddm = vddm[k];
#ifndef __AVX2__
		std::vector<float> mprofile(nbin, 0.);
#else
		vector<float, boost::alignment::aligned_allocator<float, 32>> mprofile(nbin, 0.);
#endif
		for (long int j=0; j<nchan; j++)
		{
			float delay = f0*CONST_DM*ddm*(1./(frequencies[j]*frequencies[j])-1./(fc*fc));

#ifndef __AVX2__
			for (long int i=0; i<nbin; i++)
			{
				float phi = i*1./nbin+delay;
				phi -= std::floor(phi);
				phi -= 0.5;
				phi /= 0.5*boxphw;
				mprofile[i] += std::exp(-phi*phi);
			}
#else
			if (nbin % 8 == 0)
			{
				__m256 avx_delay = _mm256_set_ps(delay,delay,delay,delay,delay,delay,delay,delay);
				__m256 avx_0 = _mm256_set_ps(0,0,0,0,0,0,0,0);
				__m256 avx_nbin = _mm256_set_ps(nbin,nbin,nbin,nbin,nbin,nbin,nbin,nbin);
				__m256 avx_1_2 = _mm256_set_ps(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5);
				__m256 avx_boxphw = _mm256_set_ps(boxphw,boxphw,boxphw,boxphw,boxphw,boxphw,boxphw,boxphw);
				avx_boxphw = _mm256_mul_ps(avx_boxphw, avx_1_2);
				for (long int i=0; i<nbin; i+=8)
				{
					__m256 avx_phi = _mm256_set_ps(i,i+1,i+2,i+3,i+4,i+5,i+6,i+7);
					avx_phi = _mm256_div_ps(avx_phi, avx_nbin);
					avx_phi = _mm256_add_ps(avx_phi, avx_delay);
					avx_phi = _mm256_sub_ps(avx_phi, _mm256_floor_ps(avx_phi));
					avx_phi = _mm256_sub_ps(avx_phi, avx_1_2);
					avx_phi = _mm256_div_ps(avx_phi, avx_boxphw);
					__m256 avx_x = exp256_ps(_mm256_sub_ps(avx_0, _mm256_mul_ps(avx_phi, avx_phi)));
					__m256 avx_pro = _mm256_load_ps(mprofile.data()+i);
					avx_pro = _mm256_add_ps(avx_pro, avx_x);
					_mm256_store_ps(mprofile.data()+i, avx_pro);
				}
			}
			else
			{
				for (long int i=0; i<nbin; i++)
				{
					float phi = i*1./nbin+delay;
					phi -= std::floor(phi);
					phi -= 0.5;
					phi /= 0.5*boxphw;
					mprofile[i] += std::exp(-phi*phi);
				}
			}
#endif
		}

		for (long int i=0; i<nbin; i++)
		{
			vchisq[k] += mprofile[i]*mprofile[i];
		}
	}
}

void PulsarPlot::get_f0_chisq_curve(std::vector<float> &vchisq, const std::vector<float> &vdf0, const std::vector<double> &tsuboff, double f0, double boxphw, int nbin)
{
	int nf0 = vdf0.size();
	int nsubint = tsuboff.size();

	vchisq.resize(nf0, 0.);
	for (long int k=0; k<nf0; k++)
	{
		float df0 = vdf0[k];
#ifndef __AVX2__
		std::vector<float> mprofile(nbin, 0.);
#else
		vector<float, boost::alignment::aligned_allocator<float, 32>> mprofile(nbin, 0.);
#endif
		for (long int j=0; j<nsubint; j++)
		{
			float delay = df0*tsuboff[j];

#ifndef __AVX2__
			for (long int i=0; i<nbin; i++)
			{
				float phi = i*1./nbin+delay;
				phi -= std::floor(phi);
				phi -= 0.5;
				phi /= 0.5*boxphw;
				mprofile[i] += std::exp(-phi*phi);
			}
#else
			if (nbin % 8 == 0)
			{
				__m256 avx_delay = _mm256_set_ps(delay,delay,delay,delay,delay,delay,delay,delay);
				__m256 avx_0 = _mm256_set_ps(0,0,0,0,0,0,0,0);
				__m256 avx_nbin = _mm256_set_ps(nbin,nbin,nbin,nbin,nbin,nbin,nbin,nbin);
				__m256 avx_1_2 = _mm256_set_ps(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5);
				__m256 avx_boxphw = _mm256_set_ps(boxphw,boxphw,boxphw,boxphw,boxphw,boxphw,boxphw,boxphw);
				avx_boxphw = _mm256_mul_ps(avx_boxphw, avx_1_2);
				for (long int i=0; i<nbin; i+=8)
				{
					__m256 avx_phi = _mm256_set_ps(i,i+1,i+2,i+3,i+4,i+5,i+6,i+7);
					avx_phi = _mm256_div_ps(avx_phi, avx_nbin);
					avx_phi = _mm256_add_ps(avx_phi, avx_delay);
					avx_phi = _mm256_sub_ps(avx_phi, _mm256_floor_ps(avx_phi));
					avx_phi = _mm256_sub_ps(avx_phi, avx_1_2);
					avx_phi = _mm256_div_ps(avx_phi, avx_boxphw);
					__m256 avx_x = exp256_ps(_mm256_sub_ps(avx_0, _mm256_mul_ps(avx_phi, avx_phi)));
					__m256 avx_pro = _mm256_load_ps(mprofile.data()+i);
					avx_pro = _mm256_add_ps(avx_pro, avx_x);
					_mm256_store_ps(mprofile.data()+i, avx_pro);
				}
			}
			else
			{
				for (long int i=0; i<nbin; i++)
				{
					float phi = i*1./nbin+delay;
					phi -= std::floor(phi);
					phi -= 0.5;
					phi /= 0.5*boxphw;
					mprofile[i] += std::exp(-phi*phi);
				}
			}
#endif
		}

		for (long int i=0; i<nbin; i++)
		{
			vchisq[k] += mprofile[i]*mprofile[i];
		}
	}
}

void PulsarPlot::get_f1_chisq_curve(std::vector<float> &vchisq, const std::vector<float> &vdf1, const std::vector<double> &tsuboff, double f0, double boxphw, int nbin)
{
	int nf1 = vdf1.size();
	int nsubint = tsuboff.size();

	vchisq.resize(nf1, 0.);
	for (long int k=0; k<nf1; k++)
	{
		float df1 = vdf1[k];
#ifndef __AVX2__
		std::vector<float> mprofile(nbin, 0.);
#else
		vector<float, boost::alignment::aligned_allocator<float, 32>> mprofile(nbin, 0.);
#endif
		for (long int j=0; j<nsubint; j++)
		{
			float delay = 0.5*df1*tsuboff[j]*tsuboff[j];

#ifndef __AVX2__
			for (long int i=0; i<nbin; i++)
			{
				float phi = i*1./nbin+delay;
				phi -= std::floor(phi);
				phi -= 0.5;
				phi /= 0.5*boxphw;
				mprofile[i] += std::exp(-phi*phi);
			}
#else
			if (nbin % 8 == 0)
			{
				__m256 avx_delay = _mm256_set_ps(delay,delay,delay,delay,delay,delay,delay,delay);
				__m256 avx_0 = _mm256_set_ps(0,0,0,0,0,0,0,0);
				__m256 avx_nbin = _mm256_set_ps(nbin,nbin,nbin,nbin,nbin,nbin,nbin,nbin);
				__m256 avx_1_2 = _mm256_set_ps(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5);
				__m256 avx_boxphw = _mm256_set_ps(boxphw,boxphw,boxphw,boxphw,boxphw,boxphw,boxphw,boxphw);
				avx_boxphw = _mm256_mul_ps(avx_boxphw, avx_1_2);
				for (long int i=0; i<nbin; i+=8)
				{
					__m256 avx_phi = _mm256_set_ps(i,i+1,i+2,i+3,i+4,i+5,i+6,i+7);
					avx_phi = _mm256_div_ps(avx_phi, avx_nbin);
					avx_phi = _mm256_add_ps(avx_phi, avx_delay);
					avx_phi = _mm256_sub_ps(avx_phi, _mm256_floor_ps(avx_phi));
					avx_phi = _mm256_sub_ps(avx_phi, avx_1_2);
					avx_phi = _mm256_div_ps(avx_phi, avx_boxphw);
					__m256 avx_x = exp256_ps(_mm256_sub_ps(avx_0, _mm256_mul_ps(avx_phi, avx_phi)));
					__m256 avx_pro = _mm256_load_ps(mprofile.data()+i);
					avx_pro = _mm256_add_ps(avx_pro, avx_x);
					_mm256_store_ps(mprofile.data()+i, avx_pro);
				}
			}
			else
			{
				for (long int i=0; i<nbin; i++)
				{
					float phi = i*1./nbin+delay;
					phi -= std::floor(phi);
					phi -= 0.5;
					phi /= 0.5*boxphw;
					mprofile[i] += std::exp(-phi*phi);
				}
			}
#endif
		}

		for (long int i=0; i<nbin; i++)
		{
			vchisq[k] += mprofile[i]*mprofile[i];
		}
	}
}