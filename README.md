# admr: Aggregate Data Modeling in R

admr (Aggregate Data Modeling in R) is an open-source R package designed to facilitate pharmacometric modeling using summary-level data. 
It enables users to work with aggregate data, such as mean observations and variance-covariance matrices, to fit pharmacokinetic and pharmacodynamic (PK/PD) models efficiently. 
This package implements a newly developed Expectation-Maximization (EM) algorithm to enhance computational performance and provides tools for advanced modeling applications.

Features
- Expectation-Maximization Algorithm: Efficiently fits pharmacometric models to aggregate data, improving speed and scalability compared to Monte Carlo methods.
- Aggregate Data Flexibility: Allows for the integration of summary-level data from diverse sources, including published literature and simulated models.
- Meta-Analysis Support: Facilitates model-based meta-analyses by enabling the combination of summary data across studies.
- R Integration: Fully compatible with R, leveraging popular pharmacometric modeling libraries like rxode2.
- Open-Source: Developed for accessibility and ease of use by the pharmacometric community.

## Installation

This is an *R* package. [*R*](https://www.r-project.org/) is required,
[*RStudio*](https://www.rstudio.com/) is recommended.

As *admr* is still hosted on Github, we need to install it via 
the [`devtools`](https://devtools.r-lib.org/) official CRAN package.

```r
devtools::install_github("hiddevandebeek/admr")
```
After installing, we have to attach the package as usual:

```r
library(admr)
```

## Basic Usage
for the future

## License
This project is licensed under the MIT License. See the LICENSE file for details.
