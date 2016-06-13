from math import sqrt
from math import factorial
from operator import *
from numpy import *

### Defining the DeltaJ function that will be used in Wigner6J

def DeltaJ(a, b, c):

	Total = 0
		
	while True:
		if (a+b-c) < 0:
			break
		elif (a-b+c) < 0:
			break
		elif (-a+b+c) < 0:
			break
		elif (a+b+c+1) < 0:
			break
	
		Total = sqrt( float(factorial(a+b-c)*factorial(a-b+c)*factorial(-a+b+c)) / float(factorial(a+b+c+1)) )
		
	#print "Total: ", Total
		
		break

	return Total
	


	
### Defining Wigner 6J function	
	
def Wigner6J(j1, j2, j3, j4, j5, j6):
	
	# Wigner 6J angular momentum coupling

	# {j1 j2 j3} =	{J_upper F_upper I}
	# {j4 j5 j6}	{F_lower J_lower 1}
	
	Delta_Total = DeltaJ(j1, j2, j3)*DeltaJ(j1, j5, j6)*DeltaJ(j4, j2, j6)*DeltaJ(j4, j5, j3)
		
	Wigner_Total = 0
	z = 0
	
	while True:
		
		if (j1+j2+j4+j5-z) < 0:
			break
		elif (j2+j3+j5+j6-z) < 0:
			break
		elif (j3+j1+j6+j4-z) < 0:
			break
	
		while True:
							
			if (z-j1-j2-j3) < 0:
				break
			elif (z-j1-j5-j6) < 0:
				break
			elif (z-j4-j2-j6) < 0:
				break
			elif (z-j4-j5-j3) < 0:
				break
		
			Wigner1 = float(factorial(z-j1-j2-j3))
			Wigner2 = float(factorial(z-j1-j5-j6))
			Wigner3 = float(factorial(z-j4-j2-j6))
			Wigner4 = float(factorial(z-j4-j5-j3))
			Wigner5 = float(factorial(j1+j2+j4+j5-z))
			Wigner6 = float(factorial(j2+j3+j5+j6-z))
			Wigner7 = float(factorial(j3+j1+j6+j4-z))
				
			Wigner_Denominator = Wigner1*Wigner2*Wigner3*Wigner4*Wigner5*Wigner6*Wigner7
			Wigner_Total = Wigner_Total + (float(-1)**z*factorial(float(z+1)))/Wigner_Denominator
						
			z = z+1
			
			if (j1+j2+j4+j5-z) <= 0:
				break
			elif (j2+j3+j5+j6-z) <= 0:
				break
			elif (j3+j1+j6+j4-z) <= 0:
				break
		z = z+1

		Total = float(Delta_Total*Wigner_Total)	
	return Total
	



	
### Defining the Gaussian function for the HF spectrum
	
def Gaussian(x, HFS_frequency, FWHM, intensity):		
		return float(intensity)*exp(- 0.5*((HFS_frequency-x)/(FWHM/2.355))**2) # Gaussian function
		
		
	
### Defining the Lorentzian function for the HF spectrum

def Lorentzian(x, HFS_frequency, gamma, intensity):		
	return intensity*(gamma**2/((x-HFS_frequency)**2 + gamma**2)) # Lorentzian function	
	
	

### Defining the Voigt function for the HF spectrum

def pseudoVoigt(x, HFS_frequency, FWHM, intensity, eta):
	Gauss = exp(-0.6931*((x-HFS_frequency)/(FWHM/2))**2)
	Lorentz = 1/(1+((x-HFS_frequency)/(FWHM/2))**2)
	Voigt = eta*Lorentz + (1-eta)*Gauss
	return intensity*Voigt # Voigt function	
	
	

	
### Defining the Crystalball function

def Crystalball(x_array, x0, N, sigma, alpha, n):
	y_array = []
	for i in range(len(x_array)):
		x = x_array[i]
		t = (x-x0)/sigma
		if (alpha < 0):
			t = -t
		if (t >= -abs(alpha)):
			y =  exp(-0.5*t*t)
		else:
			a =  ((n/abs(alpha))**n)*exp(-0.5*abs(alpha)*abs(alpha))
			b = n/abs(alpha) - abs(alpha)
			y = a/(b - t)**n
		y_array.append(N*y)
	return array(y_array)
	

### Defining an exponential pseudoVoigt function

def expoVoigt(x_array, x0, intensity, FWHM, alpha, eta):
	y_array = []
	for i in range(len(x_array)):
		x = x_array[i]
		t = (x-x0)/FWHM
		if (alpha < 0):
			t = -t
		if (t >= -abs(alpha)):
			y =  pseudoVoigt(x, x0, FWHM, intensity, eta)*exp(-0.5*t*t)
		else:
			y = pseudoVoigt(x, x0, FWHM, intensity, eta)
		y_array.append(y)
	return array(y_array)
	
	
	
	
### Defining the HF function which simulates the HF spectrum		
						
def HF_function(I, J_lower, J_upper, centroid_frequency, A_lower, A_upper, B_lower, B_upper):
	
	# Calculates the F values for a J_lower to J_upper transition
	
	HFS_frequency = []; HF_intensity = []
	F_lower_min = pos(I - J_lower)
	F_lower_max = pos(I + J_lower)
	F_upper_min = pos(I - J_upper)
	F_upper_max = pos(I + J_upper)
	
	while F_lower_min < (F_lower_max +1) :
		F_upper_min = pos(I - J_upper)
	
		while F_upper_min < (F_upper_max +1) :
			F_lower = F_lower_min
			F_upper = F_upper_min 
			F_delta = F_upper - F_lower
			
			if (-1 <= F_delta <= 1):
			
				K_lower = F_lower*(F_lower+1)-I*(I+1)-J_lower*(J_lower+1)
				alpha_lower	 = K_lower/2
				if I <= 0.5 :
					beta_lower = 0
				elif J_lower <= 0.5 :
					beta_lower = 0
				else:
					beta_lower = (3*K_lower*(K_lower+1)-4*I*(I+1)*J_lower*(J_lower+1))/(8*I*(2*I-1)*J_lower*(2*J_lower-1))
						
				K_upper = F_upper*(F_upper+1)-I*(I+1)-J_upper*(J_upper+1)
				alpha_upper = K_upper/2
				if I <= 0.5 :
					beta_upper = 0
				elif J_upper <= 0.5 :
					beta_upper = 0
				else:
					beta_upper = (3*K_upper*(K_upper+1)-4*I*(I+1)*J_upper*(J_upper+1))/(8*I*(2*I-1)*J_upper*(2*J_upper-1))
		
				HFS_frequency.append(centroid_frequency + alpha_upper*A_upper + beta_upper*B_upper - alpha_lower*A_lower - beta_lower*B_lower)
				HF_intensity.append((2*F_lower+1)*(2*F_upper+1)*Wigner6J(F_lower, F_upper, 1, J_upper, J_lower, I)**2)
				
			F_upper_min = F_upper_min +1
				
		F_lower_min = F_lower_min +1
			
	return HFS_frequency, HF_intensity



### Defining the intensity for each of the HF peaks

def HF_intensity(I, J_lower, J_upper, F_lower, F_upper):

	# Intensity ratio = (2F_lower+1)(2F_upper+1){F_lower F_upper 1}
	#											{J_upper J_lower I}
	
	Intensity = (2*F_lower+1)*(2*F_upper+1)*Wigner6J(F_lower, F_upper, 1, J_upper, J_lower, I)**2
	
	if Intensity == 0:
		print("Intensity = 0")
	
	return Intensity



### Defining Doppler correction from lab frame to rest frame

def Doppler_correction(freq_range_lab, mass, iscool_voltage):		
	alpha	= iscool_voltage/(mass*931.494061*10**6)
	freq_range_rest = freq_range_lab*( 1 + alpha - sqrt(2*alpha + alpha*alpha))
	return freq_range_rest



### Defining wavenumber to frequency conversion

def Frequency_conversion(wavenumber, mass, iscool_voltage, harmonic, frequency_correction):
	
	frequency_rest_frame = array([]); frequency_lab_frame =  array([])
	
	c		= 299792458.0								# Speed of light [ms-1] in a vacuum
	e		= 2.7182818284								# Maths constant	
	alpha	= iscool_voltage/(mass*931.494061*10**6)	# alpha = eV/mc*c - in units of e/c**2
	
	# Convert from wavenumber to frequency in lab frame
	frequency_lab_frame = harmonic*wavenumber*c/10**4	# Wavenumber doubled as reading taken at fundamental frequency (calculate in MHz)
	
	# Convert frequency from lab frame to rest frame
	frequency_rest_frame = frequency_lab_frame*( 1 + alpha - sqrt(2*alpha + alpha*alpha))
	
	# Convert to relative frequency
	frequency_relative = frequency_rest_frame - frequency_correction
	
	return frequency_relative


### Defining the HF structure

def HFS(I, J_lower, J_upper, CF, A_lower, A_upper, B_lower, B_upper, FWHM, Int, Bkgnd, x):

	for i in range(len(HF_function(I, J_lower, J_upper, CF, A_lower, A_upper, B_lower, B_upper)[0])):
		HFS_frequency = HF_function(I, J_lower, J_upper, CF, A_lower, A_upper, B_lower, B_upper)[0][i]
		intensity_1 = HF_function(I, J_lower, J_upper, CF, A_lower, A_upper, B_lower, B_upper)[1][i]
		intensity_1 = intensity_1*Int
		Bkgnd = Bkgnd + Gaussian(x, HFS_frequency, FWHM, intensity_1)
	
	return Bkgnd

	