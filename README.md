# Frequency-domain-Adaptive-Filtering-for-Auto-Echo-Cancellation
# Introduction
I implemented the frequency domain adaptation echo cancellation based on the NLMS algorithm.
The main difference to the NLMS algorithm is that all the operations are done in
the frequency domain. Frequency domain adaptive filters are both computationally efficient
and also generate signals that are uncorrelated(orthogonal). Below are the description of the
algorithm used and some experimental results.
# Algorithm
Two tasks were completed. The first one was a straightforward block based extension of the
NLMS algorithm to the frequency domain. As seen in Fig 5 of the reference [1], the following
equations were implemented after we concatenate the x_far block of previous loop to the
current x_far block. Here x_far is the uncorrupted microphone signal of first speaker. x_near
is the second speaker and y is the near speaker added to far speaker x_far convolved with
Room Impulse Response(RIR). The echo is estimated by the set of coefficients W and d is the
distorted signal that is measured which is equal to y. The estimated output of the adaptive
filter is represented as y_hat. Variables, matrices or vectors with small letters represent time
domain, big letters represent frequency domain. The time index k or block index nblock is
omitted for brevity.
![image](https://user-images.githubusercontent.com/126625677/222429087-11722df5-c1e6-4001-b17c-ac1374c82f99.png)

The last block of y is saved which is calculated as the inverse FFT of the signal multiplied with
weights in freqency domain.

![image](https://user-images.githubusercontent.com/126625677/222429275-ae52e734-a49b-4a8a-96d8-4a570ba69e28.png)

The variables Γ(gamma) is a forgetting factor and P is the PSD of far end signal x_far.

![image](https://user-images.githubusercontent.com/126625677/222429386-d7120035-4870-42b7-8f37-74c9f8cbc6bb.png)

The weights are computed as the inverse FFT of the convolved error and concatenated x signal.
A gradient constraint is applied by setting the second half of the weights to zero and the adaptive
weights are updated. The main update equation is seen in eqn 4.

![image](https://user-images.githubusercontent.com/126625677/222429482-ce830e07-9b81-44ad-9bb7-3a6223058bf4.png)

In the second task, we implemented the double talk detector with open-loop(using auto and
cross PSDs of x and distorted signal d) and closed-loop(using auto and cross PSDs of y_hat and
distorted signal d) detectors. The eqns 5 show the computation of the open loop detector where
we replace x with y_hat for the closed loop detector calculation. Hence this was implemented
as small function that is called twice. The forgetting factor gamma mentioned here is different
from the gamma used to calculate P.

![image](https://user-images.githubusercontent.com/126625677/222429587-2b860b81-3154-403f-9df6-c2455a4c3208.png)

C11 is Auto-PSD of x, C22 is Auto-PSD of d, C12 is cross-PSD of x and d. We additionally
used a weighting factor of the C11 divided by the sum of all C11 to calculate a pseudo-coherence
[2] that showed better coherence results than without using weights.

![image](https://user-images.githubusercontent.com/126625677/222429663-e5605ee3-1cdc-49f5-9f75-cc969a3c31ac.png)


Simple thresholding using the vales of ρ are used to check if the filter weights should be
adapted or not. We adapt the weights only if both the ρ values are above a threshold where
the appropriate threshold is derived by plotting the ρ for all the blocks.
Lastly, to check the quality of our output, we used an MSE criterion between y and y_hat,
normalized with total energy of target y.

# Implementation

The first task was directly implemented using the reference Fig 5 in [1]. We had to find an
appropriate γ value, the forgetting factor for the recursive computation of the PSD of the
undistorted first speaker. A γ value of 0.9 showed good results.
For the second task, the provided reference was not very clear on how to implement the open
and closed loop detectors in the frequency domain. We were able to implement the algorithm
with the help of the tutor and hints from additional material online. As mentioned previously
we experimented with weighting the Magnitude Square Coherence(MSC) with weights based on
the input signal power and found the pseudo-coherence [2] to be a better indicated than simply
averaging the MSCs over the frequency bins. The results of both tasks are shown below.
The main challenges were in resolving all the notations used between the different reference
materials and understanding the implementation of the open and closed loop detectors. Hints
from the tutor were helpful in finding the right hyper parameters for better double talk detection.

# Evaluation Results

The input near and far signals are shown below.

![image](https://user-images.githubusercontent.com/126625677/222430009-493e4bf8-95b8-42bf-92f3-aec17b287b42.png)
Figure 1: Adaptive Filter Output, Error Signal, Residual Error Signal

Time domain results of the first task with stepsize μ=0.2 and PSD forgetting factor Γ=0.9(see
Figure 2). We can see that initially the echo cancellation results are not very good, but the
error signal has very small magnitude after approximately 4sec of adapting the weights. This is
due to the initialization of the weights , where the weights updates take time to converge due
to the update step being not ideal.

![image](https://user-images.githubusercontent.com/126625677/222430053-e0e7e9e5-93c2-4c2e-982b-e1e534677c71.png)
Figure 2: Adaptive Filter Output, Error Signal, Residual Error Signal

Here we see the MSE output results with 3 different block sizes. We can see good results
with appropriate step sizes and forgetting factors. So we can conclude that we do not always
need large filters or block sizes to perform good echo cancellation. (alteast for this audio input
data). The MSE plots are shown with different μ, Γ and L(see Figure 3) and the xaxis shows
the block number.

![image](https://user-images.githubusercontent.com/126625677/222430147-b23d6a73-ab34-45f5-88b0-36b959a426ec.png)

![image](https://user-images.githubusercontent.com/126625677/222430175-48b91fed-536d-4b41-b026-721d9b2618fc.png)

![image](https://user-images.githubusercontent.com/126625677/222430195-d815ffab-b903-42f8-88df-f8f4cb31d7c9.png)

Figure 3: MSE with L=512, L=1024 and L=2048
The last figure shows the plot of the MSCs and the adaptation flag based on thresholding of
the MSCs. The last plot is the time domain difference between distorted and the estimated signal
which shows we can stop the updation of the weights when the second speaker starts talking
and restart again when only the far speaker talks, hence effectively implementing frequency
domain double talk detection and echo cancellation. Magnitude Square Coherence, Adaptation
Flag, Near-end Signal with L=2048, γ1=0.9, γ2=0.6, μ=0.2(see Figure 4)
![image](https://user-images.githubusercontent.com/126625677/222430252-fda833c0-461b-4c1e-b809-13233ad4a3f7.png)
Figure 4: Magnitude Square Coherence, Adaptation Flag, Near-end Signal, Y - estimated Y

For reference we can also see the output results with L=512 and γ2=0.2. As we see, it gets
harder to pick a good threshold when the block sizes are smaller, as the open loop MSC values
are noisier and lower overall. The closed loop MSCs are still a good indicator of double talk
here, so we could further adapt the threshold or use only closed loop MSC. The audio of the
larger block sizes sounded a bit better.

![image](https://user-images.githubusercontent.com/126625677/222430385-a8a91f37-50c8-4ba9-9965-28ac868b3322.png)
Figure 5: Magnitude Square Coherence, Adaptation Flag, Near-end Signal, Y - estimated Y

# Bibliography
[1] J. J. Shynk et al., “Frequency-domain and multirate adaptive filtering,” IEEE Signal processing
magazine, vol. 9, no. 1, pp. 14–37, 1992.

[2] T. Gänsler and J. Benesty, “A frequency-domain double-talk detector based on a normalized
cross-correlation vector,” Signal Processing, vol. 81, no. 8, pp. 1783–1787, 2001.
