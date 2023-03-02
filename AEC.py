

import numpy as np
import scipy as sci
from scipy import linalg as lin
from scipy import fftpack as fft
from scipy import io as sio
from matplotlib import pyplot as plt
import scipy.io.wavfile as wav
import IPython.display as ipd

wave_file2 = 'female.wav'
wave_file = 'male.wav'
h_file = 'h.wav'


# read and rescale input signals
freq,   x_far  = wav.read(wave_file);
freq2,  x_near = wav.read(wave_file2);
f_freq, h      = wav.read(h_file);
x_far = x_far/32768.;
x_near = x_near/32768.;
h = h/32768.;
N = len(x_far);

# near-end speaker active in [140000;180000] only
x_near[1:140000-1] = 0;
x_near[180000+1:] = 0;




# convolve far-end with RIR
y_far = np.convolve(x_far,h,mode='full');
y_near = x_near;
y = y_far[0:len(x_far)] + y_near;

# plot signals
plt.figure(figsize=(15,10))
plt.subplot(411)
plt.plot(x_far)
plt.xticks([])
plt.title('Clean far-end signal')
plt.subplot(412)
plt.plot(x_near)
plt.xticks([])
plt.title('Near-end signal')
plt.subplot(413)
plt.plot(y_far)
plt.xticks([])
plt.title('Recorded far-end signal (''echo'')')
plt.subplot(414)
plt.plot(y)
plt.xticks([])
plt.title('Recorded signal (''echo'' + ''near-end'')')
plt.show()




ipd.Audio(y.reshape(-1), rate=freq)




ipd.Audio(x_near.reshape(-1), rate=freq)




ipd.Audio(x_far.reshape(-1), rate=freq)




L = 2048;
NFFT = 2*L;

buf = np.zeros(NFFT); # input signal buffer (contains 2L samples)
W = np.zeros(NFFT); # frequency-domain filter coefficients
P = np.zeros(NFFT); # estimated PDS of X
y_hat = np.zeros(N); # signal buffer for adaptive filter output (= estimated echo)
e = np.zeros(N); # error signal (= "microphone signal" - "estimated echo")

alpha = 0.5
delta = 0.000001
mu = 0.2
gamma = 0.9

n_blocks = int(N/L)
x_old = np.zeros(L)

for k_block in range(n_blocks):
    idx = np.arange(L) + k_block*L
 
    x_b = np.concatenate([x_old, x_far[idx]]); #size NFFT
    d_b = y[idx] #distorted signal 
    x_old = x_far[idx] #replace
    
    #Q5 should the FFT be of the Rxx or the X?
    X_b = np.fft.fft(x_b, NFFT)

    Y_hat_b = W * X_b ; 
    y_hat_b = np.real(np.fft.ifft(Y_hat_b)[L:])
    y_hat[idx] = y_hat_b;  
  
    # * estimate PSD for far-end signal (for normalization in normalized LMS)
    P =  gamma*P+(1-gamma)*np.abs(X_b)**2

    e_b = d_b - y_hat_b;
    e[idx] = e_b;
    e_b_fft = np.concatenate([np.zeros(L),e_b]) #pad with L zeros
    E_b = np.fft.fft(e_b_fft, NFFT)

    w1 = np.fft.ifft(E_b*X_b.conj()/(P+1e-06))
    w1[L:] = 0 
    W = W + 2 * mu * np.fft.fft(w1, NFFT) 




ipd.Audio(e.reshape(-1), rate=freq)




ipd.Audio(y_hat.reshape(-1), rate=freq)




plt.figure(figsize=(15,8))
plt.subplot(411)
plt.plot(y_hat)
plt.xticks([])
plt.title('Adaptive filter output')
plt.subplot(412)
plt.plot(e)
plt.title('Error signal')
plt.xticks([])
plt.subplot(413)
plt.plot(e-y_near)
plt.title('Residual echo signal (= "error signal" - "near-end speech")')
plt.xticks([])
plt.show()



def calc_rho(x, d, L, gamma, c_11, c_12, c_22):
    
    X = np.fft.fft(x, 2*L)[L:]
    D = np.fft.fft(d, 2*L)[L:]
    
    c_11 = gamma*c_11 + (1-gamma)*np.abs(X)**2
    c_22 = gamma*c_22 + (1-gamma)*np.abs(D)**2
    c_12 = gamma*c_12 + (1-gamma)*(X*D.conj())
    
    coh_b = np.abs(c_12)**2/(c_11*c_22)
    wts = (1.0/np.sum(c_11))*c_11
    rho = np.sum(wts*coh_b)
    
    return rho, c_11, c_12, c_22
    


def aec_double_talk(x, y, N, L, gamma1=0.9, gamma2=0.6,stepsize=0.2, threshold=0.6):
    
    n_blocks = int(N/L)
    x_old = np.zeros(L)
    NFFT = 2*L
    
    buf = np.zeros(NFFT); # input signal buffer (contains 2L samples)
    W = np.zeros(NFFT); # frequency-domain filter coefficients
    P = np.zeros(NFFT); # estimated PDS of X
    y_hat = np.zeros(N); # signal buffer for adaptive filter output (= estimated echo)
    e = np.zeros(N); # error signal (= "microphone signal" - "estimated echo")

    # PSDs for open-loop detector
    col_12 = np.zeros(L); # cross-PSD of X and D
    col_11 = np.ones(L); # auto-PSD of X
    col_22 = np.ones(L); # auto-PSD of D

    # PSDs for closed-loop detector
    ccl_12 = np.zeros(L); # cross-PSD of Y_hat and D
    ccl_11 = np.ones(L); # auto-PSD of Y_hat
    ccl_22 = np.ones(L); # auto-PSD of D (can be estimated only once and shared with open loop detection)

    n_blocks = int(N/L);

    rho_ol = np.zeros(n_blocks);
    rho_cl = np.zeros(n_blocks);
    flagAdapt = np.zeros(n_blocks);
    mse = np.zeros(n_blocks)
    
    for k_block in range(n_blocks):
        idx = np.arange(L) + k_block*L
 
        x_b = np.concatenate([x_old, x_far[idx]]); 
        d_b = y[idx]
        x_old = x_far[idx]
    
        X_b = np.fft.fft(x_b, NFFT)
    
        Y_hat_b = W * X_b ; 
        y_hat_b = np.real(np.fft.ifft(Y_hat_b)[L:])
        y_hat[idx] = y_hat_b;  #store the output
    
        # * estimate PSD for far-end signal (for normalization in normalized LMS)
        P =  gamma1*P + (1-gamma1)*np.abs(X_b)**2
    
        rho_ol[k_block], col_11, col_12, col_22 = calc_rho(x_far[idx], d_b,L, gamma2, col_11, col_12, col_22) 
        rho_cl[k_block], ccl_11, ccl_12, ccl_22 = calc_rho(y_hat_b, d_b, L, gamma2, ccl_11, ccl_12, ccl_22)

        e_b = d_b - y_hat_b;
        e[idx] = e_b; #store the error
        e_b_fft = np.concatenate([np.zeros(L),e_b]) #pad with L zeros
        E_b = np.fft.fft(e_b_fft, NFFT)
    
        if ((rho_ol[k_block] > threshold and rho_cl[k_block] > threshold) or k_block < 30):
            flagAdapt[k_block] = 1
            w1 = (np.fft.ifft(E_b*X_b.conj()/(P+np.finfo(float).eps)))
            w1[L:] = 0 #gradient constraint
            W = W + 2 * stepsize * np.fft.fft(w1, NFFT)  



        mse[k_block] = 10 * np.log10(np.mean(np.square(y-y_hat))/np.mean(np.square(y)))
    
    return y_hat, e, mse,rho_ol, rho_cl, flagAdapt


L=2048 
mse = np.zeros((int(N/L)))
g_list = np.arange(0.5, 0.8, 0.1)
mu = np.arange(0.1, 0.6, 0.1)

plt.figure(figsize=(15,10))
for i, g in enumerate(g_list):
    for j, m in enumerate(mu):
        y_hat, e, mse, rho_ol, rho_cl, flagAdapt = aec_double_talk(x_far, y, N, L, gamma1=0.9, gamma2=g, stepsize=m, threshold=0.6)
        plt.plot(mse,label='g2 = {0:.2f}, stepsize={1:0.2f}'.format(g, m))

plt.legend()
plt.show()    




ipd.Audio(e.reshape(-1), rate=freq)




ipd.Audio(y_hat.reshape(-1), rate=freq)




y_hat, e, mse, rho_ol, rho_cl, flagAdapt = aec_double_talk(x_far, y, N, L=512, gamma1=0.9, gamma2= 0.6,
                                                            stepsize=0.2, threshold=0.2)

plt.figure(figsize=(15,10))
plt.subplot(411)
plt.plot(rho_ol,'r-*')
plt.plot(rho_cl,'b-*')
plt.title('Magnitude square coherence')
plt.xticks([])
plt.subplot(412)
plt.plot(flagAdapt)
plt.title('Adaptation flag')
plt.xticks([])
plt.subplot(413)
plt.plot(y_near)
plt.title('Near-end signal')
plt.xticks([])
plt.subplot(414)
plt.plot(y,label='y')
plt.plot(y_hat,label='y_hat')
plt.title('Near-end signal')
plt.legend()
plt.xticks([])
plt.show()




ipd.Audio(y_hat.reshape(-1), rate=freq)

