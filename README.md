# Bayesian_PSD_flutter_derivatives
Matlab code for Bayesian PSD identification of flutter derivatives operated in turbulent flow

## Installation
Please copy the whole package under MATLAB platform. 

## Theoretical Background
Please refer to our papar [Chu et al. (2022)](https://doi.org/10.1016/j.ymssp.2021.108782))
[paper](readmeFigures/paper.bmp)

## Usage

### Input structural parameters
[structure](readmeFigures/structure.bmp)

### Input wind field parameters (Only mean wind speed currently)
[windField](readmeFigures/windField.bmp)

### Upload vibration signals from sectional test
[buffetingSignal](readmeFigures/buffetingSignal.bmp)

### FFT
[fft](readmeFigures/fft.bmp)


### Optimize probability density function using Simulated Annealing. The initial value and band can be refereed to Theodorsen values or your best guess
[optimize](readmeFigures/optimize.bmp)


### Compare the reconstructed PSD from identified parameters against actual ones
[validate](readmeFigures/validate.bmp)

### Validate the identified parameters through fitness parameter lambda of the reconstructed PSD
[fitness](readmeFigures/fitness.bmp)

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.


## License
[MIT](https://choosealicense.com/licenses/mit/)
