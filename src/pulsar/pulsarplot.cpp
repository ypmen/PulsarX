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

#include "config.h"

#include "constants.h"
#include "pulsarplot.h"

#ifdef HAVE_PYTHON
    #include "matplotlibcpp.h"
    namespace plt = matplotlibcpp;
#endif

#ifdef HAVE_PLOTX
    #include "plotx.h"
#endif

using namespace std;
using namespace Pulsar;

PulsarPlot::PulsarPlot(){}

PulsarPlot::~PulsarPlot(){}

void PulsarPlot::plot(const DedispersionLite &dedisp, const ArchiveLite &archive, GridSearch &gridsearch, std::map<std::string, std::string> &obsinfo, int id, const string &rootname, bool plotx)
{
    stringstream ss_id;
    ss_id << setw(5) << setfill('0') << id;
    string s_id = ss_id.str();

    string basename = rootname + "_" + obsinfo["Date"] + "_" + obsinfo["Beam"] + "_" + s_id;
    string figname = basename + ".png";

    long int nsubint = gridsearch.nsubint;
    long int nchan = gridsearch.nchan;
    long int nbin = gridsearch.nbin;

    long int ndm = gridsearch.nddm;
    long int nf1 = gridsearch.ndf1;
    long int nf0 = gridsearch.ndf0;

    vector<float> vp, vph, vt, vf, vdm, vsnr_dm, vf0, vdf0, vsnr_f0, vf1, vdf1, vsnr_f1;
    vp.resize(nbin*2, 0.);
    vph.resize(nbin*2, 0.);
    vt.resize(nsubint);
    vf.resize(nchan);
    vdm.resize(ndm);
    vsnr_dm.resize(ndm);
    vf0.resize(nf0);
    vdf0.resize(nf0);
    vsnr_f0.resize(nf0);
    vf1.resize(nf1);
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

    for (long int k=0; k<ndm; k++)
    {
        vdm[k] = gridsearch.dm + gridsearch.ddmstart + k*gridsearch.ddmstep;
        if (gridsearch.dmsearch)
            vsnr_dm[k] = gridsearch.vsnr_dm[k];
        else
            vsnr_dm[k] = 0.;
    }

    for (long int k0=0; k0<nf0; k0++)
    {
        vf0[k0] = (gridsearch.f0 + gridsearch.df0start + k0*gridsearch.df0step);
        vdf0[k0] = vf0[k0]-gridsearch.f0;
        vsnr_f0[k0] = 0.;
    }

    for (long int k1=0; k1<nf1; k1++)
    {
        vf1[k1] = (gridsearch.f1 + gridsearch.df1start + k1*gridsearch.df1step);
        vdf1[k1] = vf1[k1]-gridsearch.f1;
        vsnr_f1[k1] = 0.;
    }

    int if1 = -1;
    int if0 = -1;
    double maxsnr = -1.;
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

            if (mxsnr_ffdot[k1*nf0+k0] > maxsnr)
            {
                maxsnr = mxsnr_ffdot[k1*nf0+k0];
                if1 = k1;
                if0 = k0;
            }
        }
    }

    for (long int k0=0; k0<nf0; k0++)
    {
        //vsnr_f0[k0] /= nf1;
        vsnr_f0[k0] = mxsnr_ffdot[if1*nf0+k0];
    }

    for (long int k1=0; k1<nf1; k1++)
    {
        //vsnr_f1[k1] /= nf0;
        vsnr_f1[k1] = mxsnr_ffdot[k1*nf0+if0];
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
    double dmsmear_phase = abs(DedispersionLite::dmdelay(dm, dedisp.frequencies[0], dedisp.frequencies.back()))/dedisp.nchans*f0;

    if (!plotx)
    {
#ifdef HAVE_PYTHON
        /** plot */
        long int nrows = 12;
        long int ncols = 10;

        plt::figure_size(1000, 1000, 100);

        //profile
        plt::subplot2grid(nrows, ncols, 0, 0, 2, 4);

        plt::plot(vph, vp);
        plt::axvspan(1.-0.5*dmsmear_phase, 1.+0.5*dmsmear_phase, 0., 1., {{"color", "lightgrey"}});
        plt::autoscale(true, "x", true);
        plt::xlabel("Phase");
        plt::ylabel("Flux");
        plt::tick_params({{"which", "both"}, {"direction", "in"}, 
                        {"bottom", "on"}, {"top", "on"}, {"left", "on"}, {"right", "on"}, 
                        {"labelbottom", "on"}, {"labeltop", ""}, {"labelleft", "on"}, {"labelright", ""}});

        //dynamic image
        plt::subplot2grid(nrows, ncols, 2, 0, 5, 4);
        plt::pcolormesh(vph, vf, mxfph);
        plt::xlabel("Phase");
        plt::ylabel("Frequency (MHz)");
        plt::axis("auto");
        plt::tick_params({{"which", "both"}, {"direction", "in"}, 
                        {"bottom", "on"}, {"top", "on"}, {"left", "on"}, {"right", "on"}, 
                        {"labelbottom", "on"}, {"labeltop", ""}, {"labelleft", "on"}, {"labelright", ""}});
        plt::twinx();
        vector<float> vfid(nchan);
        if (vf[1]<vf[0])
        {
            for (long int i=0; i<nchan; i++)
            {
                vfid[i] = nchan-1-i;
            }
        }
        else
        {
            for (long int i=0; i<nchan; i++)
            {
                vfid[i] = i;
            }        
        }
        plt::pcolormesh(vph, vfid, mxfph);
        //plt::ylabel("Index");
        plt::axis("auto");
        plt::tick_params({{"which", "both"}, {"direction", "in"}, 
                        {"bottom", "on"}, {"top", "on"}, {"left", "on"}, {"right", "on"}, 
                        {"labelbottom", ""}, {"labeltop", ""}, {"labelleft", ""}, {"labelright", "on"}});    

        //subint image
        plt::subplot2grid(nrows, ncols, 7, 0, 5, 4);
        plt::pcolormesh(vph, vt, mxtph);
        plt::xlabel("Phase");
        plt::ylabel("Tint (s)");
        plt::axis("auto");
        plt::tick_params({{"which", "both"}, {"direction", "in"}, 
                        {"bottom", "on"}, {"top", "on"}, {"left", "on"}, {"right", "on"}, 
                        {"labelbottom", "on"}, {"labeltop", ""}, {"labelleft", "on"}, {"labelright", ""}});
        plt::twinx();
        vector<float> vtid(nsubint);
        if (vt[1] > vt[0])
        {
            for (long int i=0; i<nsubint; i++)
            {
                vtid[i] = i;
            }    
        }
        else
        {
            for (long int i=0; i<nsubint; i++)
            {
                vtid[i] = nsubint-1-i;
            }
        }
        plt::pcolormesh(vph, vtid, mxtph);
        //plt::ylabel("Index");
        plt::axis("auto");
        plt::tick_params({{"which", "both"}, {"direction", "in"}, 
                        {"bottom", "on"}, {"top", "on"}, {"left", "on"}, {"right", "on"}, 
                        {"labelbottom", ""}, {"labeltop", ""}, {"labelleft", ""}, {"labelright", "on"}});

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
        plt::subplot2grid(nrows, ncols, 4, 4, 2, 4);
        plt::plot(vdf0, vsnr_f0);
        plt::axvline(xpos, 0, 1, {{"color", "red"}});
        plt::annotate("P0 (s) = "+s_p0, 0.25, 1.1, {{"xycoords","axes fraction"}, {"annotation_clip", ""}, {"fontsize", "11"}});
        plt::xlabel("F0 - " + to_string(gridsearch.f0) + " (Hz)");
        plt::ylabel("$\\chi^2$");
        plt::autoscale(true, "x", true);
        plt::ticklabel_format("x", "sci", 0, 0);
        plt::tick_params({{"which", "both"}, {"direction", "in"}, 
                        {"bottom", "on"}, {"top", "on"}, {"left", "on"}, {"right", "on"}, 
                        {"labelbottom", "on"}, {"labeltop", ""}, {"labelleft", ""}, {"labelright", "on"}});

        //chi2-f0-f1
        plt::subplot2grid(nrows, ncols, 6, 4, 4, 4);
        plt::pcolormesh(vdf0, vdf1, mxsnr_ffdot);
        plt::axis("auto");

        plt::plot(vector<double>{xpos}, vector<double>{ypos}, {{"markersize", "40"}, {"marker", "+"}, {"color", "red"}});
        plt::tick_params({{"which", "both"}, {"direction", "in"}, 
                        {"bottom", "on"}, {"top", "on"}, {"left", "on"}, {"right", "on"}, 
                        {"labelbottom", ""}, {"labeltop", ""}, {"labelleft", ""}, {"labelright", ""}});

        //chi2-f1
        std::stringstream ss_f1;
        ss_f1<<fixed<<setprecision(5)<<scientific<<f1;
        plt::subplot2grid(nrows, ncols, 6, 8, 4, 4);
        plt::plot(vsnr_f1, vdf1);
        plt::axhline(ypos, 0, 1, {{"color", "red"}});
        plt::annotate("P1 (s/s) = "+s_p1, -0.5, 1.1, {{"xycoords","axes fraction"}, {"annotation_clip", ""}, {"fontsize", "11"}});
        plt::xlabel("$\\chi^2$");
        plt::ylabel("F1 -" + ss_f1.str() + " (Hz/s)");
        plt::autoscale(true, "y", true);
        plt::ticklabel_format("y", "sci", 0, 0);
        plt::tick_params({{"which", "both"}, {"direction", "in"}, 
                        {"bottom", "on"}, {"top", "on"}, {"left", "on"}, {"right", "on"}, 
                        {"labelbottom", "on"}, {"labeltop", ""}, {"labelleft", "on"}, {"labelright", ""}});

        //chi2-dm
        plt::subplot2grid(nrows, ncols, 10, 4, 2, 4);
        plt::plot(vdm, vsnr_dm);
        plt::axvline(dmpos, 0, 1, {{"color", "red"}});
        plt::annotate("DM (pc/cc) = "+s_dm, 0.25, 1.1, {{"xycoords","axes fraction"}, {"annotation_clip", ""}, {"fontsize", "11"}});
        plt::xlabel("DM (pc/cc)");
        plt::ylabel("$\\chi^2$");
        plt::autoscale(true, "x", true);
        plt::tick_params({{"which", "both"}, {"direction", "in"}, 
                        {"bottom", "on"}, {"top", "on"}, {"left", "on"}, {"right", "on"}, 
                        {"labelbottom", "on"}, {"labeltop", ""}, {"labelleft", ""}, {"labelright", "on"}});

        std::string fontsize = "11"; 
        //metadata
        plt::subplot2grid(nrows, ncols, 0, 4, 4, 6);

        plt::text(0.02, 0.91, "Telescope = " + obsinfo["Telescope"], {{"fontsize", fontsize}});
        plt::text(0.02, 0.80, "Beam = " + obsinfo["Beam"], {{"fontsize", fontsize}});
        plt::text(0.02, 0.69, "RA = " + obsinfo["RA"], {{"fontsize", fontsize}});
        plt::text(0.02, 0.58, "DEC = " + obsinfo["DEC"], {{"fontsize", fontsize}});
        plt::text(0.02, 0.47, "P0 (s) = " + s_p0, {{"fontsize", fontsize}});
        plt::text(0.02, 0.36, "P1 (s/s) = " + s_p1, {{"fontsize", fontsize}});
        plt::text(0.02, 0.25, "DM (pc/cc) = " + s_dm, {{"fontsize", fontsize}});
        plt::text(0.02, 0.14, "acc (m/s/s) = " + s_acc, {{"fontsize", fontsize}});
        plt::text(0.02, 0.03, "S/N = " + s_snr, {{"fontsize", fontsize}});

        plt::text(0.46, 0.91, "Source name = " + obsinfo["Source_name"], {{"fontsize", fontsize}});
        plt::text(0.46, 0.80, "Date (MJD) = " + obsinfo["Date"], {{"fontsize", fontsize}});
        plt::text(0.46, 0.69, "GL (deg) = " + s_gl, {{"fontsize", fontsize}});
        plt::text(0.46, 0.58, "GB (deg) = " + s_gb, {{"fontsize", fontsize}});
        plt::text(0.46, 0.47, "F0 (Hz) = " + s_f0, {{"fontsize", fontsize}});
        plt::text(0.46, 0.36, "F1 (Hz/s) = " + s_f1, {{"fontsize", fontsize}});
        plt::text(0.46, 0.25, "MaxDM YMW16  (pc/cc) = " + s_ymw16_maxdm, {{"fontsize", fontsize}});
        plt::text(0.46, 0.14, "Distance YMW16 (pc) = " + s_ymw16_dist, {{"fontsize", fontsize}});
        plt::text(0.46, 0.03, "Pepoch (MJD) = " + obsinfo["Pepoch"], {{"fontsize", fontsize}});

        plt::tick_params({{"which", "both"}, {"direction", "in"}, 
                        {"bottom", ""}, {"top", ""}, {"left", ""}, {"right", ""}, 
                        {"labelbottom", ""}, {"labeltop", ""}, {"labelleft", ""}, {"labelright", ""}});

        plt::tight_layout(0, 0, 1, 0.97);
        
        plt::annotate(obsinfo["Filename"], 0.07, 0.97, {{"xycoords","figure fraction"}, {"annotation_clip", ""}, {"fontsize", "8"}});

        plt::save(figname);

        plt::close();
#endif
    }
    else
    {
#ifdef HAVE_PLOTX
        std::string fontsize_label = "0.7";
        std::string fontsize_ticklabel = "0.7";

        /** plot */
        PlotX::Figure fig(11.75, 1.);
        fig.set_background_color("black");
        fig.set_defalut_color("white");

        //profile
        PlotX::Axes ax1(0.08, 0.38, 0.76, 0.96);
        ax1.plot(vph, vp);
        ax1.axvspan(1.-0.5*dmsmear_phase, 1.+0.5*dmsmear_phase, 0., 1., {{"color", "lightgray"}, {"filled", "true"}});
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
#endif
    }
}
