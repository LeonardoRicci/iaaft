import numpy
import numpy.random
import numpy.fft

def generate_surrogate_iaaft(original_data, M = 1, detrend = False, verbose = True):
	"""
	surrogate_data = generate_surrogate_iaaft(original_data, M = 1, detrend = False, verbose = True)

		Generates, by means of the Iterative Amplitude Adjusted Fourier Transform
		(IAAFT) algorithm, surrogate time series corresponding to the original
		time series original_data.

	Options:

	'M'		Number of surrogate time series to be generated. By
			default, M = 1. Each of the M surrogate time series
			is a row in the returned array, which therefore
			has size MxL, where L is the original sequence
			length.

	'detrend'	Specifies whether the time series has to be detrended
			prior to surrogate generation. Default is False.

	'verbose'	Sets the verbosity of the function. If True (default),
			all messages are displayed. If False, only critical errors
			are displayed.

	NSE Laboratory, Dept. of Physics, University of Trento

	REFERENCE:
	The IAAFT (Iterative Amplitude-Adjusted Fourier Transform) algorithm
	for surrogate generation was originally proposed by T. Schreiber
	and A. Schmitz in Phys. Rev. Lett. 77 (1996), 635,
	doi:10.1103/PhysRevLett.77.635
	"""


	# --- Input validation
	if (not (isinstance(original_data,numpy.ndarray))):
			print('ERROR (in generate_surrogate_iaaft): function argument "original_data" must be a one-dimensional numpy array.')
			return False
	if (original_data.ndim != 1):
			print('ERROR (in generate_surrogate_iaaft): function argument "original_data" must be a one-dimensional numpy array.')
			return False
	L = original_data.size
	try:
		M = int(M)
	except ValueError:
		print('ERROR (in generate_surrogate_iaaft): function argument "M" must be a positive integer or convertible to a positive integer.')
		return False
	if M < 1:
		print('ERROR (in generate_surrogate_iaaft): function argument "M" must be a positive integer or convertible to a positive integer.')
		return False
	if (not (isinstance(detrend,(bool)))):
		print('ERROR (in generate_surrogate_iaaft): function option "detrend" must be a boolean.')
		return False
	if (not (isinstance(verbose,(bool)))):
		print('ERROR (in generate_surrogate_iaaft): function option "verbose" must be a boolean.')
		return False

	# --- Provide some information
	if verbose is True:
		print('\n### Starting IAAFT routine ###')

	# --- Detrend, if requested
	if (detrend):
		x0 = original_data[0]
		xN = original_data[L-1]
		if (verbose):
			print('\n### Detrending time series ###')
		time_idx = numpy.arange(L)
		original_temp = numpy.copy(original_data);
		original_data = original_temp - time_idx*(xN-x0)/(L-1.0);

	# --- Assess spectrum & distribution of original sequence
	spectrum_magnitude_original = numpy.absolute(numpy.fft.rfft(original_data))
	distribution_original = numpy.sort(original_data)

 	# --- Initialize output array (each of the M generated surrogates is a row)
	surrogate_data = numpy.empty_like(original_data);

	# --- Define check of convergence (sequence did not change with respect to last iteration)
	def check_convergence(sequence_this, sequence_prev):
		z = numpy.sum(numpy.square(sequence_prev - sequence_this))
		I = numpy.sum(numpy.square(sequence_this))
		if (z / I <= 1e-6):
			return True
		else:
			return False

	# --- Generate surrogates: iterate M times
	for iter in range(0, M):

		if verbose is True:
			print('# Surrogate number ', iter + 1, ' out of ', M, '.')

		# --- Starting conditions
		numpy.random.seed()
		run_iterations = True
		surrogate_timeseries = numpy.copy(original_data)
		numpy.random.shuffle(surrogate_timeseries)
		surrogate_timeseries_previous = numpy.copy(surrogate_timeseries)

		# --- Iterative algorithm
		n_iter = 0
		while (run_iterations == True):

			n_iter = n_iter + 1;

			spectrum_surrogate = numpy.fft.rfft(surrogate_timeseries);
			phases_surrogate = numpy.angle(spectrum_surrogate);
			spectrum_surrogate = spectrum_magnitude_original * (numpy.cos(phases_surrogate) + 1j*numpy.sin(phases_surrogate));
			surrogate_timeseries = numpy.fft.irfft(spectrum_surrogate, L);

			order_elements = numpy.argsort(surrogate_timeseries)
			surrogate_timeseries[order_elements] = distribution_original

			if (check_convergence(surrogate_timeseries, surrogate_timeseries_previous)):
				break

			surrogate_timeseries_previous = numpy.copy(surrogate_timeseries)

		if verbose is True:
			print('Process converged after ', n_iter, ' iterations.\n')

		# --- Append sequence to the set of those already generated
		if (iter == 0):
			surrogate_data = surrogate_timeseries
		else:
			surrogate_data = numpy.vstack((surrogate_data, surrogate_timeseries))

	return surrogate_data
