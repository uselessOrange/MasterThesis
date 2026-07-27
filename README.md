# MasterThesis

MATLAB research project on estimating refrigeration appliance energy performance from minimal sensor data by reconstructing compressor operation. The research was conducted for AMC Vibro in collaboration with ES System K.

## What Is In This Repo

- `Compressor_Time_Prediction/`: peak-valley compressor-state prediction from internal temperature.
- `calib_peak_valley/`: calibration workflow for cases where raw peak-valley timing does not align with real compressor switching.
- `GrayBox_Compressor_Time_Prediction/`: gray-box digital twin for compressor prediction, simulation, and analysis.
- `Functions/`: shared MATLAB functions used across the project.
- `Data/`: experimental datasets used in the thesis work.
- `RelevantPapers/`: background literature.
- `thesis.pdf`: full thesis with methodology, experiments, and conclusions.

## Project Summary

The repo compares two ways of estimating compressor operation, which is then used as the basis for simplified energy estimation.

- The peak-valley approach is lightweight and aimed at low-cost deployment with minimal sensing.
- The calibration workflow improves peak-valley on systems where temperature extrema are shifted relative to true compressor transitions.
- The gray-box approach uses a simplified physical thermal model and is better suited for simulation, diagnostics, and model-based analysis.

## Key Results

- Peak-valley achieved strong accuracy on favorable datasets using mainly internal temperature.
- Calibration reduced the difficult Dataset 3 case from `12.90%` MAPE to `3.36%`.
- The gray-box method achieved competitive prediction accuracy while also providing physically interpretable model parameters.

## Tech Used

MATLAB, signal processing, anomaly handling, time-series preprocessing, genetic-algorithm optimization, state-space modeling, and model validation.

## Where To Start

- Read `thesis.pdf` for the full research context.
- Start in `Compressor_Time_Prediction/` for the lightweight method.
- Go to `calib_peak_valley/` for the calibration protocol.
- Go to `GrayBox_Compressor_Time_Prediction/` for the digital twin work.
