#!/usr/bin/env python3
"""
	clean archive file from RFI using statistical methods
"""

import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.gridspec import GridSpec
import astropy.io.fits as fits
from scipy.stats import skew, kurtosis, shapiro
from scipy.optimize import minimize, curve_fit
from pybaselines.whittaker import airpls
from scipy.signal import savgol_filter
import argparse


def plot_channel_statistics(vf, spectra, statfile):
	# check channel statistics
	fig = plt.figure()
	gs = GridSpec(6, 2, figure=fig, width_ratios=[4, 1])

	ax1 = fig.add_subplot(gs[0, 0])
	x = np.mean(spectra, 1)
	ax1.plot(vf, x)
	ax1.set_title('mean')
	ax1.set_xlabel('Frequency (MHz)')
	ax1_2 = fig.add_subplot(gs[0, 1])
	_ = ax1_2.hist(x, bins=100)

	ax2 = fig.add_subplot(gs[1, 0])
	x = peak2peak(spectra)
	ax2.plot(vf, x)
	ax2.set_title('std')
	ax2.set_xlabel('Frequency (MHz)')
	ax2_2 = fig.add_subplot(gs[1, 1])
	_ = ax2_2.hist(x, bins=100)

	ax3 = fig.add_subplot(gs[2, 0])
	x = skew(spectra, 1)
	ax3.plot(vf, x)
	ax3.set_title('skew')
	ax3.set_xlabel('Frequency (MHz)')
	ax3_2 = fig.add_subplot(gs[2, 1])
	_ = ax3_2.hist(x, bins=100)

	ax4 = fig.add_subplot(gs[3, 0])
	x = kurtosis(spectra, 1)
	ax4.plot(vf, x)
	ax4.set_title('kurtosis')
	ax4.set_xlabel('Frequency (MHz)')
	ax4_2 = fig.add_subplot(gs[3, 1])
	_ = ax4_2.hist(x, bins=100)

	#plt.tight_layout()
	plt.savefig(statfile)
	plt.show()

def smooth_profile(x, lam=100.):
	res = minimize(lambda s: np.mean((x - s) * (x - s)) + lam * np.mean((np.roll(s, 1) - s) * (np.roll(s, 1) - s)), x, method='BFGS')
	return res.x

def compute_scloffs(dat, s):
	ex = np.sum(dat, 1)
	se = np.sum(s)
	sx = np.sum(s*dat, 1)
	ee = len(s)
	ss = np.sum(s*s)

	tmp = se*se-ss*ee

	va = np.zeros(dat.shape[0])
	vb = np.zeros(dat.shape[0])

	idx = tmp != 0

	va[idx] = ((ex*se-sx*ee))[idx]/tmp[idx]
	vb[idx] = ((sx*se-ex*ss))[idx]/tmp[idx]

	return va, vb

def remove_signal(spectra, threshold=0.5):
	x = skew(spectra, 1)
	x_baseline, _ = airpls(x, lam=1e6, diff_order=2, max_iter=50, tol=1e-3)
	idx_skew = (x - x_baseline) < threshold

	x = kurtosis(spectra, 1)
	x_baseline, _ = airpls(x, lam=1e6, diff_order=2, max_iter=50, tol=1e-3)
	idx_kurt = (x - x_baseline) < threshold

	idx = np.bitwise_and(idx_skew, idx_kurt)

	profile = np.sum(spectra[idx], 0)
	profile_smooth = smooth_profile(profile)

	va, vb = compute_scloffs(spectra, profile_smooth)

	residuals = spectra - (va[:, None] * profile_smooth[None, :] + vb[:, None])

	return residuals	

def autocorrelation(data, lag=1):
	data1 = np.roll(data, lag, axis=1)
	vmean = np.mean(data, 1)[:, None]
	return np.sum((data1 - vmean) * (data - vmean), 1) / np.sum((data - vmean) * (data - vmean), 1)

def peak2peak(data):
	return np.max(data, 1) - np.min(data, 1)

def indentify_rfi_channels(spectra, threshold=1.5):
	x = np.mean(spectra, 1)
	x_baseline, _ = airpls(x, lam=1e6, diff_order=2, max_iter=50, tol=1e-3)
	x -= x_baseline
	q75, q25 = np.percentile(x, [75, 25])
	iqr = q75 - q25
	idx_mean = np.bitwise_or(x < q25 - threshold * iqr, x > q75 + threshold * iqr)

	x = np.std(spectra, 1)
	x_baseline, _ = airpls(x, lam=1e6, diff_order=2, max_iter=50, tol=1e-3)
	x -= x_baseline
	q75, q25 = np.percentile(x, [75, 25])
	iqr = q75 - q25
	idx_std = np.bitwise_or(x < q25 - threshold * iqr, x > q75 + threshold * iqr)
    
	x = skew(spectra, 1)
	q75, q25 = np.percentile(x, [75, 25])
	iqr = q75 - q25
	idx_skew = np.bitwise_or(x < q25 - threshold * iqr, x > q75 + threshold * iqr)

	x = kurtosis(spectra, 1)
	q75, q25 = np.percentile(x, [75, 25])
	iqr = q75 - q25
	idx_kurt = np.bitwise_or(x < q25 - threshold * iqr, x > q75 + threshold * iqr)

	x = autocorrelation(spectra)
	q75, q25 = np.percentile(x, [75, 25])
	iqr = q75 - q25
	idx_ac = np.bitwise_or(x < q25 - threshold * iqr, x > q75 + threshold * iqr)

	x = peak2peak(spectra)
	q75, q25 = np.percentile(x, [75, 25])
	iqr = q75 - q25
	idx_p2p = np.bitwise_or(x < q25 - threshold * iqr, x > q75 + threshold * iqr)

	idx = idx_mean + idx_std + idx_skew + idx_kurt + idx_ac + idx_p2p
	idx = idx > 0

	return idx

def parse_options():
	parser = argparse.ArgumentParser(description="Clean archive file from RFI using statistical methods")
	parser.add_argument("infile", type=str, help="Input archive file (psrchive format)")
	parser.add_argument("--plot", action='store_true', help="Output statistics plot (optional)")
	parser.add_argument("--threshold", type=float, default=1.5, help="Threshold for RFI detection (in sigma)")
	args = parser.parse_args()
	return args

def main():
	options = parse_options()
	infile = options.infile
	outfile = infile + '.clean'
	statfile = outfile.replace('.clean', '.png')

	os.system(f"cp {infile} {outfile}")

	# read data
	f = fits.open(outfile, mode='update')

	subint_hdu = f['SUBINT']
	nsub, npol, nchan, nbin = subint_hdu.data['DATA'].shape
 
	vf = subint_hdu.data['DAT_FREQ'][0]
	vph = np.arange(nbin)/nbin

	data = subint_hdu.data['DATA'] * subint_hdu.data['DAT_SCL'].reshape(nsub, npol, nchan, 1) + subint_hdu.data['DAT_OFFS'].reshape(nsub, npol, nchan, 1)
	if subint_hdu.header['POL_TYPE'] == 'AABBCRCI':
		XX, YY = np.sum(data, axis=0)[0:2]
		spectra = XX + YY
	else:
		spectra = np.sum(subint_hdu.data['DATA'][:, 0, :, :], axis=0)

	idx_rfi = indentify_rfi_channels(spectra, threshold=options.threshold)

	subint_hdu.data['DAT_WTS'][:, idx_rfi] = 0.
	f.close()

	if options.plot:
		plot_channel_statistics(vf[np.bitwise_not(idx_rfi)], spectra[np.bitwise_not(idx_rfi)], statfile)

if __name__ == "__main__":
	main()