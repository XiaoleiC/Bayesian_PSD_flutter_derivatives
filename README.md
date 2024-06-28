# Bayesian_PSD_flutter_derivatives
Matlab code for Bayesian PSD identification of flutter derivatives operated in turbulent flow

## Installation
Please copy the whole package under MATLAB platform. "Bayesian_PSD.m" is the main function. For the purpose of illustration, people can run the main function directly.

## Theoretical Background
Please refer to our papar [Chu et al. (2022)](https://doi.org/10.1016/j.ymssp.2021.108782)
![paper](/readmeFigures/paper.png)

## Usage

### Input structural parameters
![structure](/readmeFigures/structure.png)

### Input wind field parameters (Only mean wind speed currently)
![windField](/readmeFigures/windField.png)

### Upload vibration signals from the sectional test
![buffetingSignal](/readmeFigures/buffetingSignal.png)

### FFT
![fft](/readmeFigures/fft.png)


### Optimize probability density function using Simulated Annealing. The initial value and band can be referred to Theodorsen values or your best guess
![optimize](/readmeFigures/optimize.png)

### Compare the reconstructed PSD from identified parameters against actual ones
![validate](/readmeFigures/validate.png)

### Validate the identified parameters through fitness parameter lambda of the reconstructed PSD
![fitness](/readmeFigures/fitness.png)

### Uncertainty quantification
Only the most probable values identification currently, if you need to quantify the uncertainty, please use MCMC.

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## For MCMC
We uploaded the file "buffeting_disp2.rar". You should be able to run the main function "UQ_FlutterDerivatives.m". Here I employ the AIES-based MCMC, but it lacks stability and often does not converge. I really recommend the sequential Monte Carlo to generate samples adaptively, which is promising.

## License
[MIT](https://choosealicense.com/licenses/mit/)
