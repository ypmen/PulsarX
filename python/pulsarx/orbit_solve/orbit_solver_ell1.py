#!/usr/bin/env python
"""
solve the orbit based on the epoch-F0-F1 (ELL1 model)
"""
from email.policy import default
import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 15
mpl.rcParams['ytick.labelsize'] = 15
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['axes.labelsize'] = 15
import numpy as np
import matplotlib.pyplot as plt
import logging
import argparse
import pymultinest
from scipy.optimize import fsolve
import orbit_utils
import yaml
import os, re

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
EPS1_min = (config['parameter_range']['EPS1_min'])
EPS1_max = (config['parameter_range']['EPS1_max'])
EPS2_min = (config['parameter_range']['EPS2_min'])
EPS2_max = (config['parameter_range']['EPS2_max'])

dat = np.loadtxt(config['data']['path'])

vT = dat[:, 0]
T_c = np.mean(vT)
vT -= np.mean(vT)
vf0 = dat[:, 1]
vf0err = dat[:, 2]
if config['data']['fit_F1']:
	vf1 = dat[:, 3]
	vf1err = dat[:, 4]

def get_ffdot(T, f0, asini_c, Pb, Tasc, eps1, eps2):
	vffdot = orbit_utils.compute_f0_f1_ell1(T, f0, asini_c, Pb, Tasc, eps1, eps2)
	return vffdot[:, 0], vffdot[:, 1]

def myprior(cube, ndim, nparams):
	# f0
	cube[0] = F0_min + cube[0] * (F0_max-F0_min)
	# log10(a*sin(i)/c)
	cube[1] = A1_min + cube[1] * (A1_max-A1_min)
	# log10(Pb)
	cube[2] = PB_min + cube[2] * (PB_max-PB_min)
	# Tasc / Pb
	cube[3] = -0.5 + cube[3] * 1.
	# eps1
	cube[4] = EPS1_min + cube[4] * (EPS1_max-EPS1_min)
	# eps2
	cube[5] = EPS2_min + cube[5] * (EPS2_max-EPS2_min)

def myloglike(cube, ndim, nparams):
	return likli(cube)

def likli(vpar):
	f0 = vpar[0]
	asini_c = 10**vpar[1]
	Pb = 10**vpar[2]
	Tasc = vpar[3] * Pb
	eps1 = vpar[4]
	eps2 = vpar[5]
	
	vf0_s, vf1_s = get_ffdot(vT, f0, asini_c, Pb, Tasc, eps1, eps2)
  
	chi2 = orbit_utils.get_chi2(np.concatenate([vf0, vf0_s, vf0err]).reshape(3, -1))

	if config['data']['fit_F1']:
		chi2 += orbit_utils.get_chi2(np.concatenate([vf1, vf1_s, vf1err]).reshape(3, -1))

	if config['parameter_range']['EFAC']:
		return -0.5 * len(vT) * np.log(chi2)
	else:
		return -0.5 * chi2  

def get_mass_function(asini_c, Pb, e):
	return (4.*np.pi**2)*(asini_c*const_c)**3/(Pb*86400.)**2/const_G/const_solarmass

def get_minimum_mass(f):
	m1 = 1.4
	def func(m2):
		return f * (m1 + m2)**2 - m2**3
	root = fsolve(func, 10.)
	return root[0]

def format_val_err(val, err, style="plain", low=-5., high=6.):
	if style == "sci" and (np.abs(val) < 10**low or np.abs(val) > 10**high):
		n = 0
		if val != 0.:
			n = np.int_(np.floor(np.log10(np.abs(val))))
			val *= 10**(-n)
			err *= 10**(-n)
			
			if err > 1.:
				return f"{val:.0f}({err:.0f})e{n}"
			else:
				n = np.int_(-np.floor(np.log10(err))+1)
				return f"{val:.{n}f}({err*10**n:.0f})e{n}"
		else:
			n = np.int_(np.floor(np.log10(np.abs(err))))
			val *= 10**(-n)
			err *= 10**(-n)
			return f"{val:.0f}({err:.0f})e{n}"
	else:
		if err >= 1. or (val == 0. and err == 0.):
			return f"{val:.0f}({err:.0f})"
		else:
			if (err < np.abs(val) and err != 0.) or val == 0.:
				n = np.int_(-np.floor(np.log10(err))+1)
				return f"{val:.{n}f}({err*10**n:.0f})"
			else:
				n = np.int_(-np.floor(np.log10(np.abs(err)))+1)
				return f"{val:.{n}f}({err*10**n:.0f})"

def circular_transform(samples):
	phi_mean = np.arctan2(np.mean(np.sin(samples*2*np.pi)), np.mean(np.cos(samples*2*np.pi)))/(2*np.pi)
	phimod = (samples-phi_mean)%1.
	ind = phimod>0.5
	phimod[ind] = phimod[ind]-1.
	return phimod+phi_mean

def plot(vpar, vpar_err):
	f0 = vpar[0]
	f0_err = vpar_err[0]
	asini_c = vpar[1]
	asini_c_err = vpar_err[1]
	Pb = vpar[2]
	Pb_err = vpar_err[2]
	Tasc = vpar[3]
	Tasc_err = vpar_err[3]
	eps1 = vpar[4]
	eps1_err = vpar_err[4]
	eps2 = vpar[5]
	eps2_err = vpar_err[5]

	s_f0 = format_val_err(f0, f0_err, style='sci')
	s_asini_c = format_val_err(asini_c, asini_c_err, style='sci')
	s_Pb = format_val_err(Pb, Pb_err, style='sci')
	s_Tasc = format_val_err(Tasc + T_c, Tasc_err, style='sci')
	s_eps1 = format_val_err(eps1, eps1_err, style='sci')
	s_eps2 = format_val_err(eps2, eps2_err, style='sci')

	fmass = get_mass_function(asini_c, Pb, 0.)

	mass_min = get_minimum_mass(fmass)

	log.info(f"binary parameters:\n" f"F0\t{f0}\t{f0_err}\n" f"A1\t{asini_c}\t{asini_c_err}\n" f"PB\t{Pb}\t{Pb_err}\n" f"Tasc\t{Tasc + T_c}\t{Tasc_err}\n" f"EPS1\t{eps1}\t{eps1_err}\n" f"EPS2\t{eps2}\t{eps2_err}\n")
	log.info(f"mass_function = {fmass}, mass_min={mass_min}")

	if not config['plot']['mean_anomaly']:
		vT_s = np.linspace(vT.min()-1., vT.max()+1., 1000)
	else:
		vT_s = Tasc + np.arange(1000)/1000. * Pb

	vf0_s, vf1_s = get_ffdot(vT_s, f0, asini_c, Pb, Tasc, eps1, eps2)
 
	vf0_t, vf1_t = get_ffdot(vT, f0, asini_c, Pb, Tasc, eps1, eps2)

	deltaf0 = vf0 - vf0_t
 
	fig = plt.figure(figsize=(12, 9))
	fig.suptitle(f"A1={s_asini_c}, PB={s_Pb}, Tasc={s_Tasc}, EPS1={s_eps1}, EPS2={s_eps2}\n"f"F0={s_f0} (Hz), f={fmass:2f}, Massmin={mass_min:2f}", fontsize=14)

	if config['data']['fit_F1']:
		log.info(f"plot F0 and F1 fitting")
  
		deltaf1 = vf1 - vf1_t

		ax1 = plt.subplot(221)
		ax2 = plt.subplot(223, sharex = ax1)
	
		ax3 = plt.subplot(222)
		ax4 = plt.subplot(224, sharex = ax3)
	
		if config['plot']['mean_anomaly']:
			ax1.plot((vT_s-Tasc)/Pb % 1., vf0_s, 'r')
			ax1.errorbar((vT-Tasc)/Pb % 1., vf0, yerr=vf0err, fmt='bo')
			ax1.set_ylabel(r'$f0 (Hz)$')
		
			ax2.errorbar((vT-Tasc)/Pb % 1., deltaf0, yerr=vf0err, fmt='bo')
			ax2.set_xlabel(r'Mean anomaly (phase)')
			ax2.set_ylabel(r'$\Delta f0 (Hz)$')
	
			ax3.plot((vT_s-Tasc)/Pb % 1., vf1_s, 'r')
			ax3.errorbar((vT-Tasc)/Pb % 1., vf1, yerr=vf1err, fmt='bo')
			ax3.set_ylabel(r'$f1 (Hz)$')
		
			ax4.errorbar((vT-Tasc)/Pb % 1., deltaf1, yerr=vf1err, fmt='bo')
			ax4.set_xlabel(r'Mean anomaly (phase)')
			ax4.set_ylabel(r'$\Delta f1 (Hz)$')
		else:
			ax1.plot(vT_s, vf0_s, 'r')
			ax1.errorbar(vT, vf0, yerr=vf0err, fmt='bo')
			ax1.set_ylabel(r'$f0 (Hz)$')
		
			ax2.errorbar(vT, deltaf0, yerr=vf0err, fmt='bo')
			ax2.set_xlabel(r'T (day)')
			ax2.set_ylabel(r'$\Delta f0 (Hz)$')
	
			ax3.plot(vT_s, vf1_s, 'r')
			ax3.errorbar(vT, vf1, yerr=vf1err, fmt='bo')
			ax3.set_ylabel(r'$f1 (Hz)$')
		
			ax4.errorbar(vT, deltaf1, yerr=vf1err, fmt='bo')
			ax4.set_xlabel(r'T (day)')
			ax4.set_ylabel(r'$\Delta f1 (Hz)$')
	
		ax2.set_title(rf"$\chi^2 = {np.mean(deltaf0*deltaf0/(vf0err*vf0err)):.2f}$", fontsize=15)
		ax4.set_title(rf"$\chi^2 = {np.mean(deltaf1*deltaf1/(vf1err*vf1err)):.2f}$", fontsize=15)
	else:
		log.info(f"plot F0 fitting")
  
		ax1 = plt.subplot(211)
		ax2 = plt.subplot(212, sharex = ax1)
	
		if config['plot']['mean_anomaly']:
			ax1.plot((vT_s-Tasc)/Pb % 1., vf0_s, 'r')
			ax1.errorbar((vT-Tasc)/Pb % 1., vf0, yerr=vf0err, fmt='bo')
			ax1.set_ylabel(r'$f0 (Hz)$')
			
			ax2.errorbar((vT-Tasc)/Pb % 1., deltaf0, yerr=vf0err, fmt='bo')
			ax2.set_xlabel(r'Mean anomaly (phase)')
			ax2.set_ylabel(r'$\Delta f0 (Hz)$')
		else:
			ax1.plot(vT_s, vf0_s, 'r')
			ax1.errorbar(vT, vf0, yerr=vf0err, fmt='bo')
			ax1.set_ylabel(r'$f0 (Hz)$')
			
			ax2.errorbar(vT, deltaf0, yerr=vf0err, fmt='bo')
			ax2.set_xlabel(r'T (day)')
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

	samples = np.loadtxt(out_dir+'/'+rootname+'post_equal_weights.dat')
	samples[:, 1] = 10**samples[:, 1]
	samples[:, 2] = 10**samples[:, 2]
	samples[:, 3] = circular_transform(samples[:, 3])
	samples[:, 4] = samples[:, 4]
	samples[:, 5] = samples[:, 5]
	id = np.argmax(samples[:, 6])
	vpar0 = samples[id, :-1]
	samples[:, 3] *= samples[:, 2]
	vpar0 = samples[id, :-1]
	vpar0_err = np.std(samples[:, :-1], 0)

	samples_full = np.loadtxt(out_dir+'/'+rootname+'.txt')
	samples_full[:, 2+1] = 10**samples_full[:, 2+1]
	samples_full[:, 2+2] = 10**samples_full[:, 2+2]
	samples_full[:, 2+3] = circular_transform(samples_full[:, 2+3])
	samples_full[:, 2+4] = samples_full[:, 2+4]
	samples_full[:, 2+5] = samples_full[:, 2+5]
	id = np.argmax(-0.5*samples_full[:, 1])
	vpar0 = samples_full[id, 2]

	plot(vpar0, vpar0_err)

if __name__ == "__main__":
	main()
