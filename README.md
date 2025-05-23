# Iron-Overload-Quant

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

- Synthetic dataset generation with known ground truth parameters  
- Implementation of a spatially-regularized ADMM solver with total variation (TV) penalty  
- Pixel-wise Levenberg–Marquardt mono-exponential fitting via MATLAB's `lsqcurvefit`  
- Quantitative and visual evaluation on phantoms and four patient datasets  
- ROI-based analysis and comparison against scanner console estimates  

## Requirements

- MATLAB R2020a or later  
- Optimization Toolbox  
- Image Processing Toolbox  
- (Optional) Python with NumPy and Matplotlib for auxiliary scripts  

## Installation

1. Clone this repository:  
   ```bash
   git clone https://github.com/yourusername/QHLX2-Iron-Quantification.git
   cd QHLX2-Iron-Quantification
