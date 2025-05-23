# Iron-Overload-Quant

````markdown
# Quantification of Hepatic Iron Overload in Thalassemia Patients Using MRI

This repository contains code and resources for implementing and evaluating two MRI T2* mapping techniques—Alternating Direction Method of Multipliers (ADMM) and Levenberg–Marquardt (LM)—to quantify hepatic iron overload in transfusion-dependent thalassemia patients. The project was carried out as part of the Computational Modelling for Biomedical Imaging (QHLX2) module at University College London.

---

## Table of Contents

- [Background](#background)  
- [Features](#features)  
- [Requirements](#requirements)  
- [Installation](#installation)  
- [Usage](#usage)  
- [Data](#data)  
- [Methods](#methods)  
  - [ADMM-Based Fitting](#admm-based-fitting)  
  - [Levenberg–Marquardt Fitting](#levenberg–marquardt-fitting)  
- [Results](#results)  
- [Contributing](#contributing)  
- [License](#license)  
- [References](#references)  

---

## Background

Regular blood transfusions in β-thalassemia major lead to progressive iron overload in the liver. MRI-based T2* relaxometry offers a non-invasive alternative to biopsy for monitoring liver iron concentration (LIC). This project compares two fitting techniques—ADMM and LM—to assess their accuracy, robustness to noise, and clinical relevance in both synthetic phantoms and patient data.

## Features

- Synthetic dataset generation with known ground-truth parameters  
- Spatially-regularized ADMM solver with total variation (TV) penalty  
- Pixel-wise Levenberg–Marquardt mono-exponential fitting via MATLAB’s `lsqcurvefit`  
- Quantitative and visual evaluation on synthetic phantoms and four patient datasets  
- ROI-based analysis and comparison against scanner-console estimates  

## Requirements

- MATLAB R2020a or later  
- Optimization Toolbox  
- Image Processing Toolbox  
- (Optional) Python with NumPy and Matplotlib for auxiliary scripts  

## Installation

1. **Clone this repository**  
   ```bash
   git clone https://github.com/angusgcunningham/Iron-Overload-Quant.git
   cd Iron-Overload-Quant
````

2. **Add project folders to MATLAB path**

   ```matlab
   addpath(genpath('src'))
   ```
3. Ensure required MATLAB toolboxes are installed and licensed.

## Usage

### Synthetic Phantom Evaluation

1. **Generate synthetic phantom**

   ```matlab
   createPhantoms;  % Produces a 4-quadrant phantom with known T2* values
   ```
2. **Run ADMM fitting on noise-free data**

   ```matlab
   results_ADMM = relaxationEst(syntheticData);
   ```
3. **Run LM fitting on noisy data**

   ```matlab
   results_LM = task_11_synthetic_LM(noisyData);
   ```
4. **Inspect results**
   Error tables and comparison figures are saved in `results/phantom/`.

### Patient Data Analysis

1. **Load & preprocess DICOM volumes**

   ```matlab
   [vol, info] = dicomreadPatient('data/patientX');
   ```
2. **Select liver ROIs**

   * Use interactive GUI or provided masks in `data/masks/`.
3. **Run ADMM fitting**

   ```matlab
   maps_ADMM = relaxationEst(vol);
   ```
4. **Run LM fitting**

   ```matlab
   maps_LM = task_11_LM_model(vol);
   ```
5. **Review outputs**
   Spatial maps: `results/patients/`
   ROI summary tables: `results/tables/`

## Data

* **Synthetic phantoms**: `data/phantoms/`
* **Patient DICOMs**: `data/patient1/`, `data/patient2/`, etc.
* **Preprocessed MAT files**: `data/preprocessed/`

## Methods

### ADMM-Based Fitting

* **Model**: $S(t) = a \, e^{-r t}$
* **Regularization**: TV penalties on amplitude $a$ and rate $r$
* **Solver**: Augmented Lagrangian with variable splitting
* **Iterations**: Up to 5,000 ADMM cycles or convergence

### Levenberg–Marquardt Fitting

* **Algorithm**: MATLAB’s `lsqcurvefit` (LM)
* **Constraints**: $a \ge 0$, $0 \le r \le 0.5$
* **Post-processing**: Median filtering of T2\* maps

## Results

* **Synthetic data**:

  * Both methods achieve <0.4% error under noise-free and noisy conditions.
  * ADMM shows marginally greater robustness at edges and low-SNR voxels.
* **Patient studies**:

  * Mean T2\* ranged from \~5 ms (severe overload) to \~15 ms (near-normal).
  * LM estimates matched console values in 3/4 cases; ADMM produced smoother spatial maps.

All figures and quantitative tables are in the `results/` directory.

## References

1. Hernando, D. et al., “Quantification of liver iron with MRI: state of the art and remaining challenges,” *Journal of Magnetic Resonance Imaging* (JMRI), 2014.
2. Bondestam, S. et al., “Magnetic resonance imaging of transfusional hepatic iron overload,” *British Journal of Radiology* (BJR), 1994.
3. Tipirneni-Sajja, A. et al., “Robust R2\* mapping in the presence of high iron concentrations,” *Magnetic Resonance in Medicine* (MRM), 2020.

