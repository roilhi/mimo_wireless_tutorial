import numpy as np
from numpy.random import standard_normal
import matplotlib.pyplot as plt
import math

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
def biterr_calculation(mod_idx, demod_idx, bit_width):
  '''
  Calculates the number of error bits (or different) having 2 lists with the
  indices of the transmitted and detected QAM symbols.
  Inputs:
  mod_bits --> list with the indices of the transmitted symbols
  demod_bits --> list with the indices of the decoded symbols
  bit_width --> Number of bits to be used for conversion
  Outputs:
  bit_err --> the number of bits which are different from the indices lists
  pe --> the probability of bit errors
  '''
  if isinstance(mod_idx,list) and isinstance(demod_idx,list):
    if len(mod_idx) != len(demod_idx):
       raise ValueError("The lists of indices must have the same lenght")
    mod_bits = ''.join([np.binary_repr(dec, width=bit_width) for dec in mod_idx])
    demod_bits = ''.join([np.binary_repr(dec, width=bit_width) for dec in demod_idx])
  else:
     mod_bits = np.binary_repr(mod_idx, width=bit_width)
     demod_bits = np.binary_repr(demod_idx, width=bit_width)
  bit_err = sum(1 for a, b in zip(mod_bits, demod_bits) if a != b)
  return bit_err

Nt = 1
Nr = 2
SNR_dB = np.arange(start=0, stop=20, step=2)
SNR_l = 10**(SNR_dB/10)
be = 0
BER = []
#g = []
Pe = []
n_iter = 10**4

M = 4
bpsym = int(np.log2(M))
y = np.array(qam_symbol_generator(M))
FN = 1/np.sqrt((2/3)*(M-1))
y = FN*y # normalización de referencia

suma = 0
for q in range(M):
    pot1 = np.sqrt((y[q].real)**2+(y[q].imag)**2)
    suma += pot1
pot = suma/M
y = y/pot

x = np.array([[y[1]]]) # tomar un símbolo con norma unitaria

for snr in SNR_l:
    for k in range(n_iter):
        H = H_channel(Nr,Nt)
        datos = np.sqrt(snr)*np.matmul(H,x)
        n = awgn_noise(1,Nr)
        #n = 0
        r = datos+n
        #s = np.abs(r-np.sqrt(snr)*y*H.T)**2
        s1 = np.abs(r[0]-np.sqrt(snr)*y*H[0])**2
        s2 = np.abs(r[1]-np.sqrt(snr)*y*H[1])**2
        s = s1+s2
        index = np.argmin(s)
        if index != 1:
            a = biterr_calculation(1,index,bpsym)
            be += a
    BER.append(be/(n_iter*bpsym*Nt))
    be = 0
    b = snr
    z = 0.5*(1-np.sqrt((b/2)/(1+(b/2))))
    #g.append(z)
    ss = 0
    for nn in range(Nr):
       ss+= math.comb((Nr-1+nn),nn)*((1-z)**nn)
    Pe.append((z**Nr)*ss)
       

fig,ax = plt.subplots(nrows=1, ncols=1)
ax.semilogy(SNR_dB, BER, 'o-r', label='Simulated')
ax.semilogy(SNR_dB, Pe, '--*k', label='Theoretical')
ax.set_xlabel('SNR (dB)')
ax.set_ylabel('ABEP')
ax.set_title('4-QAM SIMO system ($N_t=1$, $N_r=2$)')
ax.legend()
plt.grid()
plt.show()