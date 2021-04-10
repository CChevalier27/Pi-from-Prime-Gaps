#import libraries
import sympy
import numpy as np
import scipy.integrate as integrate

def approximate_prime(primes): #find approximate value for largest prime in primes (same method as in video)
    primegaps_dict = {}
    for i in range(len(primes)-1):
        gap = primes[i+1] - primes[i]
        if gap in primegaps_dict: 
            primegaps_dict[gap] += 1
        else:
            primegaps_dict[gap] = 1

    primegaps = np.asarray(list(primegaps_dict.keys()))
    frequency = np.asarray(list(primegaps_dict.values()))

    m, c = np.polyfit(primegaps, np.log(frequency), 1)

    approx_prime = integrate.quad(lambda x,m,c: x*np.exp(m*x+c), 0, primegaps[-1], args=(m,c))[0] + 3
    return approx_prime

prime_approximations = [2,3]
primes = [3]
max=1000

#populate prime_approximations
for i in range(3,max):
    primes.append(sympy.prime(i))
    prime_approximations.append(approximate_prime(primes))
    if(i%100 == 0):
        print(str(int(100*i/max)) + "% complete.")

#calculate zeta(2) from the euler product using our prime approximations
zeta_2 = 1
for prime in prime_approximations:
    zeta_2 *= pow((1-pow(prime,-2)),-1)

approx_pi = np.sqrt(6*zeta_2)
print("\npi = " + str(approx_pi))
print("%.2f percent error." % (100*abs(approx_pi/np.pi - 1)))