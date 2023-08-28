import numpy as np
from numpy.random import standard_normal
import matplotlib.pyplot as plt

def qam_symbol_generator(M):
  '''
  Generates QAM constellation symbols.
  Input --> M: QAM constellation order
  Output ---> M-QAM alphabet
  If M=4 || M=8 M-Psk constellation
  '''
  N = np.log2(M)
  if N != np.round(N):
    raise ValueError("M must be 2^n for n=0,1,2...")
  m = np.arange(M)
  c = np.sqrt(M)
  b = -2*(np.array(m)%c) + c-1
  a = 2*np.floor(np.array(m) / c ) - c+1
  s = list((a+1j*b)) #list of QAM symbols
  return s
def H_channel(Nr,Nt):
  '''
  Generates the H channel matrix having:
  Inputs
  Nt --> number of Tx antennas
  Nr ---> number of Rx antennas
  Output: H --> channel matrix of complex coefficients
  '''
  return (1/np.sqrt(2))*(standard_normal((Nr,Nt))+1j*(standard_normal((Nr,Nt))))
def awgn_noise(No, Nr):
  '''
  Generates AWGN additive noise:
  Inputs:
  Nr ---> Number of Rx antennas
  No --> Noise spectrum density No = Es / ((10**(EbNo/10))*np.log2(M))
  Output:
  vector of AWGN noise
  '''
  return np.sqrt(No/2)*(standard_normal((Nr,1))+1j*standard_normal((Nr,1)))
def biterr_calculation(mod_idx, demod_idx, M):
  '''
  Calculates the number of error bits (or different) having 2 lists with the
  indices of the transmitted and detected QAM symbols.
  Inputs:
  mod_bits --> list with the indices of the transmitted symbols
  demod_bits --> list with the indices of the decoded symbols
  Outputs:
  bit_err --> the number of bits which are different from the indices lists
  pe --> the probability of bit errors
  '''
  bits_per_symbol = int(np.log2(M))
  if isinstance(mod_idx,list) and isinstance(demod_idx,list):
    if len(mod_idx) != len(demod_idx):
       raise ValueError("The lists of indices must have the same lenght")
    mod_bits = ''.join([np.binary_repr(dec, width=bits_per_symbol) for dec in mod_idx])
    demod_bits = ''.join([np.binary_repr(dec, width=bits_per_symbol) for dec in demod_idx])
  else:
     mod_bits = np.binary_repr(mod_idx, width=bits_per_symbol)
     demod_bits = np.binary_repr(demod_idx, width=bits_per_symbol)
  bit_err = sum(1 for a, b in zip(mod_bits, demod_bits) if a != b)
  return bit_err

Nt = 1
Nr = 1
SNR_dB = np.arange(start=0, stop=32, step=2)
SNR_l = 10**(SNR_dB/10)
be = 0
BER = []
n_iter = 10**4

M = 4
bpsym = np.log2(M)
y = np.array(qam_symbol_generator(M))
FN = 1/np.sqrt((2/3)*(M-1))
y = FN*y # normalización de referencia

suma = 0
for q in range(M):
    pot1 = np.sqrt((y[q].real)**2+(y[q].imag)**2)
    suma += pot1
pot = suma/M
y = y/pot

x = y[0] # tomar un símbolo con norma unitaria
g = []
for snr in SNR_l:
    for k in range(n_iter):
        H1 = H_channel(Nr,Nt)
        dato1 = np.sqrt(snr)*(x*H1.T)
        n = awgn_noise(1,Nr)
        #n = 0
        r1 = dato1+n
        s1 = np.abs(r1-np.sqrt(snr)*y*H1.T)**2
        s = s1
        index = np.argmin(s)
        if index != 0:
            a = biterr_calculation(0,index,M)
            be += a
    BER.append(be/(n_iter*bpsym))
    be = 0
    b = snr
    g.append(0.5*(1-np.sqrt((b/2)/(1+(b/2)))))

print(f'BER = {BER}')

fig,ax = plt.subplots(nrows=1, ncols=1)
ax.semilogy(SNR_dB, BER, 'o-r', label='Simulated')
ax.semilogy(SNR_dB, g, '--*k', label='Theoretical')
ax.set_xlabel('SNR (dB)')
ax.set_ylabel('ABEP')
ax.set_title('4-QAM Modulation SISO System $N_t=1$, $N_r=1$')
ax.legend()
plt.grid()
plt.show()