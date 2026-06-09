#pragma once
#include "common.h"

namespace PyMKF {

// Utility functions for processing and calculations
double resolve_dimension_with_tolerance(json dimensionWithToleranceJson, std::string preferredValue = "Nominal");
json calculate_basic_processed_data(json waveformJson);
json calculate_harmonics(json waveformJson, double frequency);
json calculate_sampled_waveform(json waveformJson, double frequency);
json calculate_processed_data(json signalDescriptorJson, json sampledWaveformJson, bool includeDcComponent);

// Power calculation utilities  
double calculate_instantaneous_power(json excitationJson);
double calculate_rms_power(json excitationJson);

// Reflection utilities
json calculate_reflected_secondary(json primaryExcitationJson, double turnRatio);
json calculate_reflected_primary(json secondaryExcitationJson, double turnRatio);

// Array conversion utilities
py::list list_of_list_to_python_list(std::vector<std::vector<double>> arrayOfArrays);
std::vector<double> python_list_to_vector(py::list pythonList);

// Signal processing utilities
json standardize_signal_descriptor(json signalDescriptorJson, double frequency);
json create_waveform(json processedJson, double frequency);
json calculate_processed(json harmonicsJson, json waveformJson);
json scale_waveform_time_to_frequency(json waveformJson, double newFrequency);
json scale_excitation_time_to_frequency(json excitationJson, double newFrequency);
json calculate_induced_voltage(json excitationJson, double magnetizingInductance);
json calculate_induced_current(json excitationJson, double magnetizingInductance);
bool check_requirement(json requirementJson, double value);
json get_main_harmonic_indexes(json harmonicsJson, double threshold, int mainHarmonicIndex);
json get_excitation_harmonic_indexes(json excitationJson, double threshold);

// Register all utility bindings
void register_utils_bindings(py::module& m);

} // namespace PyMKF