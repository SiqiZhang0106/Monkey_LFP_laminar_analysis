# Monkey_LFP_laminar_analysis
We analyzed laminar LFP data recorded from bi-hemispheric premotor and primary motor cortex of monkeys during a unimanual reaching task.

Task-related epochs in the following analysis: Baseline, Target onset, Go cue, Reach onset, Enter target, Return on, Return peak

#Epoched_PSD_fooof: aims to calculate the laminar PSD and also remove the aperiodic components by FOOOF. The oscillatory power and aperiodic parameters (offset and slope could be found changing in different electrode depths (layers), and also changing in different epochs.

#Epoched_beta_average_trials: provides the time-frequency analysis by Superlets algorithm, then tries to extract the power in beta frequency band (13-30)Hz, and finally plots the distribution of beta power from the superficial to deep layers by a gradient figure.

#Reach_onset_LFP_CSD: applies Current Source Density (CSD) analysis to the LFP signals in the reach onset epoch when the movement is on.

#csd_burst_lM1/ csd_burst_rM1: firstly detects the beta bursts in each trials among all sessions (Note here: Target onset epoch and Enter target epoch are used to detect beta bursts here, longer than the previous epochs), the electrode used to detect the bursts on are based on the flipped layer of common averaging LFP. Then the burst windows are applied to all electrodes to calculate the CSD among layers.
