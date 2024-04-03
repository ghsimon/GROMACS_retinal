import numpy
cimport numpy
from scipy import interpolate

ctypedef numpy.float_t DTYPE_t

class diffusionC():
	'''
	Diffusion from umbrella sampling simulation.

	Reference:
	Gerhard Hummer 2005 New J. Phys.7 34
	https://iopscience.iop.org/article/10.1088/1367-2630/7/1/034/pdf
	'''

	def __init__(self, int Nsteps, double dt, int stride):

		# Number of timesteps per simulation
		self.Nsteps = Nsteps

		# Timestep
		self.dt     = dt

		# Number of timesteps used by PLUMED to print out the data
		self.stride = stride

	def compute(self, numpy.ndarray[DTYPE_t, ndim=1] x):

		cdef int t = 0
		cdef int i

		cdef DTYPE_t avgx  = numpy.mean(x)
		cdef numpy.ndarray[DTYPE_t, ndim=1] corr  = numpy.zeros(self.Nsteps)


		while corr[t-1] >= 0:
			for i in numpy.arange(self.Nsteps-(t)):
				corr[t] = corr[t] + (x[i] - avgx) * (x[i + t] - avgx)

			corr[t] = corr[t] / (self.Nsteps - 1)

			t += 1

		return corr[0] ** 2 / numpy.trapz(corr * self.stride * self.dt)

	def spline(xj, xhalf, diff_arr):
		tck      = interpolate.splrep(xj, diff_arr, s=0.0000018)
		diff_spl = interpolate.splev(xhalf, tck)

		return diff_spl
