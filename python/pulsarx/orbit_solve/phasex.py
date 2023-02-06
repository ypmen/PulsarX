#!/usr/bin/env python

__author__ = "Yunpeng Men"
__date__ = "2022/12/14"
__license__ = "GPLv3"
__version__ = "0.0.0"
__maintainer__ = "Yunpeng Men"
__email__ = "ypmen@mpifr-bonn.mpg.de"

import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
import tqdm
import orbit_utils
import logging
import argparse
import yaml, re

loader = yaml.SafeLoader
loader.add_implicit_resolver(
	u'tag:yaml.org,2002:float',
	re.compile(u'''^(?:
	 [-+]?(?:[0-9][0-9_]*)\\.[0-9_]*(?:[eE][-+]?[0-9]+)?
	|[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)
	|\\.[0-9_]+(?:[eE][-+][0-9]+)?
	|[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\\.[0-9_]*
	|[-+]?\\.(?:inf|Inf|INF)
	|\\.(?:nan|NaN|NAN))$''', re.X),
	list(u'-+0123456789.'))

log = logging.getLogger("orbit_solver")
FORMAT = "[%(levelname)s - %(asctime)s - %(filename)s:%(lineno)s] %(message)s"
logging.basicConfig(format=FORMAT)
log.setLevel("INFO")

def parse_timfile(timfile):
	dat = np.loadtxt(timfile, dtype='float')

	vTbmean = np.mean(dat[:, 0])

	idx = dat[:, -1]

	vTb_arr_test = dat[:, 0] - vTbmean
	vTberr_arr_test = dat[:, 1]
	vobs_ssb_arr_test = dat[:, 2:5]

	vTb_arr_fit = vTb_arr_test[idx != 0]
	vTberr_arr_fit = vTberr_arr_test[idx != 0]
	vobs_ssb_arr_fit = vobs_ssb_arr_test[idx != 0, :]

	return vTb_arr_fit, vTberr_arr_fit, vobs_ssb_arr_fit, vTb_arr_test, vTberr_arr_test, vobs_ssb_arr_test, vTbmean

def get_opts():
	parser = argparse.ArgumentParser()
	parser.add_argument("-y", "--yaml", dest="yaml", type=str, help="configuration file", default="yaml.conf", metavar="YAML")
	options = parser.parse_args()
	with open(options.yaml, 'r') as f:
		config = yaml.safe_load(f)
		log.info(config)
	return config

config = get_opts()

vTb_arr_fit, vTberr_arr_fit, vobs_ssb_arr_fit, vTb_arr_test, vTberr_arr_test, vobs_ssb_arr_test, vTbmean = parse_timfile(config['data']['path'])

F0_min = config['parameter_range']['F0_min']
F0_max = config['parameter_range']['F0_max']
F1_min = config['parameter_range']['F1_min']
F1_max = config['parameter_range']['F1_max']
DRA_min = config['parameter_range']['DRA_min']
DRA_max = config['parameter_range']['DRA_max']
DDEC_min = config['parameter_range']['DDEC_min']
DDEC_max = config['parameter_range']['DDEC_max']
A1_min = config['parameter_range']['A1_min']
A1_max = config['parameter_range']['A1_max']
PB_min = config['parameter_range']['PB_min']
PB_max = config['parameter_range']['PB_max']
T0_min = config['parameter_range']['T0_min']-vTbmean
T0_max = config['parameter_range']['T0_max']-vTbmean
ECC_min = config['parameter_range']['ECC_min']
ECC_max = config['parameter_range']['ECC_max']
OM_min = config['parameter_range']['OM_min'] / 360. * 2. * np.pi
OM_max = config['parameter_range']['OM_max'] / 360. * 2. * np.pi
OMDOT_min = config['parameter_range']['OMDOT_min'] / 360. * 2. * np.pi
OMDOT_max = config['parameter_range']['OMDOT_max'] / 360. * 2. * np.pi

f0 = config['pulsar']['F0'][0]
f1 = config['pulsar']['F1'][0]
dα = config['pulsar']['DRA'][0]
dδ = config['pulsar']['DDEC'][0]
asini_c = config['pulsar']['A1'][0]
Pb = config['pulsar']['PB'][0]
T0 = config['pulsar']['T0'][0]-vTbmean
e = config['pulsar']['ECC'][0]
omega = config['pulsar']['OM'][0] / 360. * 2. * np.pi
omega_dot = config['pulsar']['OMDOT'][0] / 360. * 2. * np.pi

# calculate search grid
vpherr_test = vTberr_arr_test * 1e-6 * 0.5 * (F0_max + F0_min)
vpherr_fit = vTberr_arr_fit * 1e-6 * 0.5 * (F0_max + F0_min)

Tspan = vTb_arr_fit.max() - vTb_arr_fit.min()

f0_step = 1. / (Tspan * 86400.)
f1_step = 2. / (Tspan * 86400.)**2
dα_step = 1. / (500. * f0)
dδ_step = 1. / (500. * f0)
asini_c_step = 1. / f0
Pb_step = 1.
T0_step = 1.
if asini_c != 0. :
	Pb_step = Pb * Pb / (2. * np.pi * f0 * asini_c * Tspan)
	T0_step = Pb / (2. * np.pi * f0 * asini_c)
e_step = 0. # 1. / (asini_c * f0)
omega_step = 0. # 1. / (asini_c * f0 * e)
omega_dot_step = 0.

ngrid_f0 = int((F0_max - F0_min) / f0_step) * 3
ngrid_f1 = int((F1_max - F1_min) / f1_step) * 3
ngrid_dα = int((DRA_max - DRA_min) / dα_step) * 1
ngrid_dδ = int((DDEC_max - DDEC_min) / dδ_step) * 1
ngrid_asini_c = int((A1_max - A1_min) / asini_c_step) * 3
ngrid_Pb = int((PB_max - PB_min) / Pb_step) * 6
ngrid_T0 = int((T0_max - T0_min) / T0_step) * 3
ngrid_e = 1
ngrid_omega = 1
ngrid_omega_dot = 1

ngrid_f0 = np.max([ngrid_f0, 1])
ngrid_f1 = np.max([ngrid_f1, 1])
ngrid_dα = np.max([ngrid_dα, 1])
ngrid_dδ = np.max([ngrid_dδ, 1])
ngrid_asini_c = np.max([ngrid_asini_c, 1])
ngrid_Pb = np.max([ngrid_Pb, 1])
ngrid_T0 = np.max([ngrid_T0, 1])
ngrid_e = np.max([ngrid_e, 1])
ngrid_omega = np.max([ngrid_omega, 1])
ngrid_omega_dot = np.max([ngrid_omega_dot, 1])

vf0 = np.linspace(F0_min, F0_max, ngrid_f0)
vf1 = np.linspace(F1_min, F1_max, ngrid_f1)
vdα = np.linspace(DRA_min, DRA_max, ngrid_dα)
vdδ = np.linspace(DDEC_min, DDEC_max, ngrid_dδ)
vasini_c = np.linspace(A1_min, A1_max, ngrid_asini_c)
vPb = np.linspace(PB_min, PB_max, ngrid_Pb)
vT0 = np.linspace(T0_min, T0_max, ngrid_T0)
ve = np.linspace(ECC_min, ECC_max, ngrid_e)
vomega = np.linspace(OM_min, OM_max, ngrid_omega)
vomega_dot = np.linspace(OMDOT_min, OMDOT_max, ngrid_omega_dot)

ra = config['pulsar']['ra']
dec = config['pulsar']['dec']

def ra2α(ra):
	h,m,s = ra.split(':')
	return (float(h)/24.+float(m)/1440.+float(s)/86400.)*360.*np.pi/180.

def dec2δ(dec):
	if dec[0] == '-':
		sn = -1.
		d,m,s = dec[1:].split(':')
	else:
		sn = 1.
		d,m,s = dec.split(':')
	return sn*(float(d)/360.+float(m)/21600.+float(s)/1296000.)*360.*np.pi/180.

α = ra2α(ra)
δ = dec2δ(dec)

fit_f0 = config['pulsar']['F0'][1] != 0
fit_f = config['pulsar']['F1'][1] != 0
fit_ra = config['pulsar']['DRA'][1] != 0
fit_dec = config['pulsar']['DDEC'][1] != 0
fit_asini_c = config['pulsar']['A1'][1] != 0
fit_Pb = config['pulsar']['PB'][1] != 0
fit_T0 = config['pulsar']['T0'][1] != 0
fit_e = config['pulsar']['ECC'][1] != 0
fit_omega = config['pulsar']['OM'][1] != 0
fit_omega_dot = config['pulsar']['OMDOT'][1] != 0

if_fit = np.array([fit_f0, fit_f, fit_ra, fit_dec, fit_asini_c, fit_Pb, fit_T0, fit_e, fit_omega, fit_omega_dot])

fitpar = [vf0, vf1, vdα, vdδ, vasini_c, vPb, vT0, ve, vomega, vomega_dot]

dpar = np.mean(vpherr_fit) * np.array([(par.max()-par.min())/len(par) for par in fitpar]) * 10

def get_phase(vTb_arr, f0, f1, dα, dδ, asini_c, Pb, T0, e, omega, omega_dot, vobs_ssb_arr, α, δ):
	vph_arr = orbit_utils.timing2(vTb_arr, dα, dδ, f0, f1, asini_c, Pb, T0, e, omega, omega_dot, vobs_ssb_arr, α, δ)
	return vph_arr

def get_chisq(vTb_arr, p0, vpherr, vobs_ssb_arr, α, δ):
	vph = get_phase(vTb_arr, *p0, vobs_ssb_arr, α, δ)
	vchi = vph / vpherr
	chisq = np.sum(vchi*vchi)
	return chisq

def plot(vTb_arr, p0, vpherr, vobs_ssb_arr, α, δ):
	vph = get_phase(vTb_arr, *p0, vobs_ssb_arr, α, δ)
	fig = plt.figure()
	plt.errorbar(vTb_arr, vph, yerr=vpherr, linestyle='', marker='o')
	plt.xlabel(f'MJD-{vTbmean}')
	plt.ylabel(f'Residual in pulse periods (days)')

def optimize_orbit(f0):
	p_best_test = None
	chisq_best_test = np.sum(1./(vpherr_test*vpherr_test))

	p_best_fit = None
	chisq_best_fit = np.sum(1./(vpherr_fit*vpherr_fit))

	vidx = []
	vlen = []
	vpar = []

	p0 = np.array([f0, f1, dα, dδ, asini_c, Pb, T0, e, omega, omega_dot])

	for k in np.arange(1, len(if_fit)):
		if if_fit[k]:
			vidx.append(k)
			vlen.append(len(fitpar[k]))
			vpar.append(fitpar[k])

	vlen = np.array(vlen)
	ndim = len(vlen)
	N = vlen.prod()
	for k in np.arange(N):
		cnt = 1
		n = k
		while cnt <= ndim:
			idx, n = n % vlen[ndim-cnt], n // vlen[ndim-cnt]
			p0[vidx[ndim-cnt]] = vpar[ndim-cnt][idx]
			cnt += 1
			
			p = orbit_utils.lmfit2(vTb_arr_fit, np.zeros(len(vTb_arr_fit)), p0, vpherr_fit, dpar, if_fit, vobs_ssb_arr_fit, α, δ)

			chisq_fit = get_chisq(vTb_arr_fit, p, vpherr_fit, vobs_ssb_arr_fit, α, δ)
			if chisq_fit < chisq_best_fit:
				p_best_fit = p
				chisq_best_fit = chisq_fit

			p = orbit_utils.lmfit2(vTb_arr_test, np.zeros(len(vTb_arr_test)), p, vpherr_test, dpar, if_fit, vobs_ssb_arr_test, α, δ)
			
			chisq_test = get_chisq(vTb_arr_test, p, vpherr_test, vobs_ssb_arr_test, α, δ)
			if chisq_test < chisq_best_test:
				p_best_test = p
				chisq_best_test = chisq_test

	return p_best_test, chisq_best_test, p_best_fit, chisq_best_fit

def main():
	log.info(f"f0 = {f0}\t{if_fit[0]}")
	log.info(f"f1 = {f1}\t{if_fit[1]}")
	log.info(f"dα = {dα}\t{if_fit[2]}")
	log.info(f"dδ = {dδ}\t{if_fit[3]}")
	log.info(f"asini_c = {asini_c}\t{if_fit[4]}")
	log.info(f"Pb = {Pb}\t{if_fit[5]}")
	log.info(f"T0 = {T0+vTbmean}\t{if_fit[6]}")
	log.info(f"e = {e}\t{if_fit[7]}")
	log.info(f"omega = {omega}\t{if_fit[8]}")
	log.info(f"omega_dot = {omega_dot}\t{if_fit[9]}")
	
	log.info(f"time span = {Tspan} days")

	log.info(f"reference epoch = {vTbmean:.15f}")

	if if_fit[0]:
		log.info(f"f0 grid: {F0_min} {F0_max} {ngrid_f0}")
	if if_fit[1]:
		log.info(f"f1 grid: {F1_min} {F1_max} {ngrid_f1}")
	if if_fit[2]:
		log.info(f"dα grid: {DRA_min} {DRA_max} {ngrid_dα}")
	if if_fit[3]:
		log.info(f"dδ grid: {DDEC_min} {DDEC_max} {ngrid_dδ}")
	if if_fit[4]:
		log.info(f"asini_c grid: {A1_min} {A1_max} {ngrid_asini_c}")
	if if_fit[5]:
		log.info(f"Pb grid: {PB_min} {PB_max} {ngrid_Pb}")
	if if_fit[6]:
		log.info(f"T0 grid: {T0_min} {T0_max} {ngrid_T0}")
	if if_fit[7]:
		log.info(f"e grid: {ECC_min} {ECC_max} {ngrid_e}")
	if if_fit[8]:
		log.info(f"omega grid: {OM_min} {OM_max} {ngrid_omega}")
	if if_fit[9]:
		log.info(f"omega_dot grid: {OMDOT_min} {OMDOT_max} {ngrid_omega_dot}")
 
	# search
	p_best = None
	chisq_best_test = np.sum(1./(vpherr_test*vpherr_test))
	vchisq_test = np.ones(ngrid_f0) * np.sum(1./(vpherr_test*vpherr_test))
	vchisq_fit = np.ones(ngrid_f0) * np.sum(1./(vpherr_fit*vpherr_fit))
	vp = np.empty((len(vf0), len(fitpar)))
	if not config['plot']['plot_only']:
		if fit_f0:
			with Pool(config['run']['cpus']) as pool:
				res = list(tqdm.tqdm(pool.imap(optimize_orbit, vf0), total=len(vf0)))
				for k,r in enumerate(res):
					p_best_test_local, chisq_best_test_local, p_best_fit_local, chisq_best_fit_local = r
					vp[k] = p_best_test_local

					if chisq_best_test_local < vchisq_test[k]:
						vchisq_test[k] = chisq_best_test_local
					if chisq_best_test_local < chisq_best_test:
						p_best = p_best_test_local
						chisq_best_test = chisq_best_test_local

					if chisq_best_fit_local < vchisq_fit[k]:
						vchisq_fit[k] = chisq_best_fit_local
		else:
			p_best, chisq_best_test, p_best_fit, chisq_best_fit = optimize_orbit(f0)

	np.savetxt('best.txt', p_best)
	plot(vTb_arr_test, p_best, vpherr_test, vobs_ssb_arr_test, α, δ)
	plt.savefig('best.png')

if __name__ == "__main__":
	main()

