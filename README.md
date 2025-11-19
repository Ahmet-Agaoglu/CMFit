# CMFit: A Tool for Line Fitting Under Bounded Errors 

## Authors
[Ahmet Agaoglu](https://github.com/Ahmet-Agaoglu)

## Description
This repository provides a MATLAB implementation of a fast and accurate method for line fitting under bounded measurement errors. The approach is based on the Corridor Method, which computes the Feasible Parameter Set (FPS) exactly and enables two types of parameter estimation:

1. **Unknown-error case:** the centroid of the FPS  
2. **Normally distributed errors:** the likelihood-maximizing point inside the FPS

Both estimators rely on the exact FPS geometry, ensuring reliable regression when only bounded-error information is available.

## Files Included
This repository contains three MATLAB functions:

- **CorridorMethod.m**  
  Computes the FPS as a convex polygon in slope-intercept space.

- **UnknownErrorEstimate.m**  
  Computes the centroid of the FPS for the unknown-error case.

- **NormalErrorEstimate.m**  
  Computes the optimal likelihood-maximizing estimate for normally distributed errors.

All functions support optional FPS visualization and plotting of estimation results.

### Inputs/Outputs
- **x** — Independent variable  
- **y** — Dependent variable  
- **lb** — Lower bounds on y 
- **ub** — Upper bounds on y

- **V** - Vertex set of the FPS
- **e_1** — FPS centroid  
- **e_2** — Likelihood-optimal solution

CorridorMethod.m accepts x,lb,ub and returns V.
The following functions call CorridorMethod.m to get V.
UnknownErrorEstimate.m accepts x,lb,ub and returns e_1.
NormalErrorEstimate.m accepts x,y,lb,ub and returns e_2.

## Installation
1. Clone or download this repository.
2. Set the folder as your MATLAB working directory.

### Dependencies
- MATLAB R2014 or later  

## Usage
To compute e_1 and e_2, call UnknownErrorEstimate.m and NormalErrorEstimate.m, respectively.

### Example

```matlab
x = 1:10;                         % define independent variable vector
slope = 2.1; 
intercept = -0.45;               % define true parameters
y_t = slope*x + intercept;       % calculate true response values

maxErr = 0.5;                    % define error bound
y = y_t + (rand(size(x))*2 - 1)*maxErr;   % define noisy response

lb = y - maxErr;                 % measurement lower bound
ub = y + maxErr;                 % measurement upper bound

e_1 = UnknownErrorEstimate(x, y, lb, ub);   % estimate for unknown error case
e_2 = NormalErrorEstimate(x, y, lb, ub);    % estimate for normal error case
```


## Reproducibility

To reproduce the experimental results presented in the accompanying study, follow the steps below:

1. Prepare synthetic or real-world datasets where measurement errors are bounded by known lower and upper limits.

2. Run UnknownErrorEstimate.m to get e_1 and NormalErrorEstimate.m to get e_2.

3. Compare the resulting estimates with baseline methods: Ordinary Least Squares (OLS): polyfit, Constrained Linear Least Squares (CLLS): lsqlin, Regularized Chebyshev Center (RCC). For validation, reproduce the accuracy (MAE) and runtime results by looping over multiple dataset sizes (N = 20, 50, 100, 200, 500) and multiple repetitions (e.g., 100 runs per configuration).

## Citation
If you use this code in your research, please cite:
@article{Agaoglu2025, title={CMFit: A Fast and Accurate Tool for Line Fitting Under Bounded Errors}, author={Ahmet Agaoglu and Namik Ciblak}, journal={Software Impacts}, year={2025}, doi={}}


## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

