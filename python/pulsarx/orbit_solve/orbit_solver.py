#!/usr/bin/env python
"""
solve the orbit based on the epoch-F0-F1
to-do: 1) the error of parameter is not checked; 2) add F2 search; 3) add OMDOT search
"""
from email.policy import default
import numpy as np
import matplotlib.pyplot as plt
import logging
import argparse
import pymultinest
from scipy.optimize import fsolve
import orbit_utils
import yaml
import os

log = logging.getLogger("orbit_solver")
FORMAT = "[%(levelname)s - %(asctime)s - %(filename)s:%(lineno)s] %(message)s"
logging.basicConfig(format=FORMAT)
log.setLevel("INFO")

const_G = 6.6743e-11
const_c = 299792458
const_solarmass = 1.9891e30

def get_opts():
	parser = argparse.ArgumentParser()
	parser.add_argument("-y", "--yaml", dest="yaml", type=str, help="configuration file", metavar="YAML")
	options = parser.parse_args()
	with open(options.yaml, 'r') as f:
		config = yaml.safe_load(f)
		log.info(config)
	return config

config = get_opts()

F0_min = config['parameter_range']['F0_min']
F0_max = config['parameter_range']['F0_max']
A1_min = np.log10(config['parameter_range']['A1_min'])
A1_max = np.log10(config['parameter_range']['A1_max'])
PB_min = np.log10(config['parameter_range']['PB_min'])
PB_max = np.log10(config['parameter_range']['PB_max'])
ECC_min = config['parameter_range']['ECC_min']
ECC_max = config['parameter_range']['ECC_max']

dat = np.loadtxt(config['data']['path'])

vT = dat[:, 0]
vf0 = dat[:, 1]
vf0err = dat[:, 2]
if config['data']['fit_F1']:
	vf1 = dat[:, 3]
	vf1err = dat[:, 4]

def myprior(cube, ndim, nparams):
	# f0
	cube[0] = F0_min + cube[0] * (F0_max-F0_min)
	# log10(a*sin(i)/c)
	cube[1] = A1_min + cube[1] * (A1_max-A1_min)
	# log10(Pb)
	cube[2] = PB_min + cube[2] * (PB_max-PB_min)
	# T0 / Pb
	cube[3] = 0. + cube[3] * 1.
	if not config['binary']['circular']:
		# e
		cube[4] = ECC_min + cube[4] * (ECC_max-ECC_min)
		# omega
		cube[5] = 0 + cube[5] * 2. * np.pi

def myloglike(cube, ndim, nparams):
	return likli(cube)

def likli(vpar):
	f0 = vpar[0]
	asini_c = 10**vpar[1]
	Pb = 10**vpar[2]
	T0 = vpar[3] * Pb
	e = 0.
	omega = 0.
	if not config['binary']['circular']:
		e = vpar[4]
		omega = vpar[5]
	
	vffdot = orbit_utils.compute_f0_f1(vT, f0, asini_c, Pb, T0, e, omega)
 
	vf0_s, vf1_s = vffdot[:, 0], vffdot[:, 1]
 
	chi2 = orbit_utils.get_chi2(np.concatenate([vf0, vf0_s, vf0err]).reshape(3, -1))

	if config['data']['fit_F1']:
		chi2 += orbit_utils.get_chi2(np.concatenate([vf1, vf1_s, vf1err]).reshape(3, -1))
  
	return chi2

def get_mass_function(asini_c, Pb, e):
	return Pb*86400.*(2.*np.pi/Pb/86400.*asini_c*const_c)**3/(2.*np.pi*const_G)*(1.-e*e)**1.5/const_solarmass

def get_minimum_mass(f):
	m1 = 1.4
	def func(m2):
		return f * (m1 + m2)**2 - m2**3
	root = fsolve(func, 10.)
	return root[0]

def plot(vpar):
	f0 = vpar[0]
	asini_c = 10**vpar[1]
	Pb = 10**vpar[2]
	T0 = vpar[3] * Pb
	e = 0.
	omega = 0.
	if not config['binary']['circular']:
		e = vpar[4]
		omega = vpar[5]

	fmass = get_mass_function(asini_c, Pb, e)

	mass_min = get_minimum_mass(fmass)

	log.info(f"binary parameters:\nF0\t{f0}\nA1\t{asini_c}\nPb\t{Pb}\nT0\t{T0}\nECC\t{e}\nOM\t{omega*180./np.pi}")
	log.info(f"mass_function = {fmass}, mass_min={mass_min}")

	vT_s = T0 + np.arange(1000)/1000. * Pb
	vffdot_s = orbit_utils.compute_f0_f1(vT_s, f0, asini_c, Pb, T0, e, omega)
	vf0_s, vf1_s = vffdot_s[:, 0], vffdot_s[:, 1]
 
	vffdot_t = orbit_utils.compute_f0_f1(vT, f0, asini_c, Pb, T0, e, omega)
	vf0_t, vf1_t = vffdot_t[:, 0], vffdot_t[:, 1]

	deltaf0 = vf0 - vf0_t		
 
	fig = plt.figure(figsize=(12, 9))
	fig.suptitle(rf"$\nu = {f0:.6f} (Hz), \frac{{a \sin i}}{{c}} = {asini_c:.2f} (ls), P_b = {Pb:.5f} (day), e = {e:2f}, f = {fmass:.8f} (M_\odot), M_{{min}} = {mass_min:.2f} (M_\odot)$", fontsize=14)

	if config['data']['fit_F1']:
		log.info(f"plot F0 and F1 fitting")
  
		deltaf1 = vf1 - vf1_t

		ax1 = plt.subplot(221)
		ax2 = plt.subplot(223, sharex = ax1)
	
		ax3 = plt.subplot(222)
		ax4 = plt.subplot(224, sharex = ax3)
	
		ax1.plot((vT_s-T0)/Pb % 1., vf0_s, 'r')
		ax1.errorbar((vT-T0)/Pb % 1., vf0, yerr=vf0err, fmt='bo')
		ax1.set_ylabel(r'$f0 (Hz)$')
		
		ax2.errorbar((vT-T0)/Pb % 1., deltaf0, yerr=vf0err, fmt='bo')
		ax2.set_xlabel(r'Mean anomaly (phase)')
		ax2.set_ylabel(r'$\Delta f0 (Hz)$')
	
		ax3.plot((vT_s-T0)/Pb % 1., vf1_s, 'r')
		ax3.errorbar((vT-T0)/Pb % 1., vf1, yerr=vf1err, fmt='bo')
		ax3.set_ylabel(r'$f1 (Hz)$')
		
		ax4.errorbar((vT-T0)/Pb % 1., deltaf1, yerr=vf1err, fmt='bo')
		ax4.set_xlabel(r'Mean anomaly (phase)')
		ax4.set_ylabel(r'$\Delta f1 (Hz)$')
	
		ax2.set_title(rf"$\chi^2 = {np.mean(deltaf0*deltaf0/(vf0err*vf0err)):.2f}$")
		ax4.set_title(rf"$\chi^2 = {np.mean(deltaf1*deltaf1/(vf1err*vf1err)):.2f}$")
	else:
		log.info(f"plot F0 fitting")
  
		ax1 = plt.subplot(211)
		ax2 = plt.subplot(212, sharex = ax1)
	
		ax1.plot((vT_s-T0)/Pb % 1., vf0_s, 'r')
		ax1.errorbar((vT-T0)/Pb % 1., vf0, yerr=vf0err, fmt='bo')
		ax1.set_ylabel(r'$f0 (Hz)$')
		
		ax2.errorbar((vT-T0)/Pb % 1., deltaf0, yerr=vf0err, fmt='bo')
		ax2.set_xlabel(r'Mean anomaly (phase)')
		ax2.set_ylabel(r'$\Delta f0 (Hz)$')
	
		ax2.set_title(rf"$\chi^2 = {np.mean(deltaf0*deltaf0/(vf0err*vf0err)):.2f}$")

	log.info(f"save fitting plot orbit_fitting.png")
	plt.savefig('orbit_fitting.png')

def main():
	
	out_dir = config['multinest']['out_dir']
	rootname = config['multinest']['rootname']
 
	try:
		os.mkdir(out_dir)
	except OSError as error:
		log.warning(f"{out_dir} already exists")
 
	if not config['plot']['plot_only']:
		ndim = 4
		if not config['binary']['circular']:
			ndim = 6
		pymultinest.run(myloglike, myprior, ndim,
						n_params=None,
						n_clustering_params=None,
						wrapped_params=None,
						importance_nested_sampling=False,
						multimodal=False,
						const_efficiency_mode=False,
						n_live_points=config['multinest']['n_live_points'],
						evidence_tolerance=config['multinest']['evidence_tolerance'],
						sampling_efficiency=config['multinest']['sampling_efficiency'],
						n_iter_before_update=100,
						null_log_evidence=-1e+90,
						max_modes=100,
						mode_tolerance=-1e+90,
						outputfiles_basename=out_dir+'/'+rootname,
						seed=-1,
						verbose=True,
						resume=False,
						context=0,
						write_output=True,
						log_zero=-1e+100,
						max_iter=0,
						init_MPI=False,
						dump_callback=None)

	samples = np.loadtxt(out_dir+'/'+rootname+'.txt')
	id = np.argmax(-samples[:, 1])
	vpar0 = samples[id, 2:]
	plot(vpar0)

if __name__ == "__main__":
	main()