# Monte Carlo Simulations for Insurance Portfolio

## Overview
This repository contains the code for a final assignment project in the "Simulation Techniques" course. The project involves applying **Monte Carlo simulations** to model risks in an **insurance portfolio**. The main objective is to simulate the probability of accidents and calculate insurance reimbursements based on various parameters, with a focus on the effect of weather conditions.

The simulations are implemented using **Markov chains** and Monte Carlo methods, where two hypotheses were tested to explore the impact of weather on the accident probabilities and insurance reimbursements.

## Key Sections of the Report:
- **Introduction**: Overview of the project and its objectives.
- **Background and Context**: Explanation of the modeling approach using Markov chains and the insurance risk analysis.
- **Parameters Setup**: Definition and interpretation of the key parameters used in the simulation.
- **Distribution of the Law of Pr**: Methods used to simulate the law of insurance claims and optimizations.
- **Hypothesis A & B**: Simulation of weather and accident models under different assumptions.
- **Stationary Measure and Validity Domain**: Analysis of the stationary distribution and validation using the Gelman-Rubin criteria.
- **Hypothesis B with C++**: Performance improvements using C++ integration for key functions.
- **Suggestions**: Ideas for further improvements and optimizations of the model.
- **Conclusion**: Summary of findings and reflections on the project’s performance.
- **Final Code**: The final version of all the implemented functions.

## Code Implementation

The code includes the following functions and features:
- **R Functions**: Implementations of the weather and accident models, both for **Hypothesis A** (same weather for all motorcyclists) and **Hypothesis B** (independent weather for each motorcyclist).
- **Performance Optimizations**: Integrations with C++ (via Rcpp) to speed up time-consuming operations like `rbinom` and `sample` functions.
- **Stationary Measure**: Implementation of a stationary measure based on Markov chain theory to optimize weather generation.
- **Gelman-Rubin Test**: Used to verify convergence of the Markov chain to the stationary distribution.

## Libraries Used
- **Rcpp**: For C++ integration within R to improve performance in generating weather and accident data.
- **dplyr**: For data manipulation and handling of results (if needed).
- **ggplot2**: For visualization of the results (e.g., probability distributions and cumulative distribution functions).

## Performance Results
- The optimizations made using **C++** integration significantly improved the performance of key functions. The final tests show a considerable decrease in simulation time.

## Final Code
All the final versions of the functions, optimized and selected for the final report, are included in the repository. The code is structured in an easy-to-understand way for future improvements or adaptations.

## License
This project is open-source and can be freely used and modified under the [MIT License](https://opensource.org/licenses/MIT).

---

### Note:
You can find the complete, final version of the code for each function at the end of the provided report.
