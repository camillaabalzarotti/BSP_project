# Biomedical Signal Processing project
This project, developed by Camilla Balzarotti (me), Emma Battaglia, Ester Penati and Sara Rebuzzi, as part of the Biomedical Signal Processing and Medical Images course at Politecnico di Milano, investigates how brain activity changes between resting state and cognitive engagement during mental arithmetic tasks.

It focuses on identifying variations in EEG spectral power and functional connectivity (coherence) across different frequency bands and brain regions.

🔗 A more detailed explanation is available in the pdf located in the presentation folder. 

## Goals
* Understand how brain activity shifts from rest to task-based engagement.
* Identify patterns of power distribution across EEG frequency bands.
* Compare connectivity between different brain regions under varying cognitive loads.

## Materials and methods

**Subjects** : 6 participants

**Recordings**: EEG signals during

**Resting state**: 3 min

**Arithmetic tasks**: 1 min

**Frequency bands:**
* Theta (4–7 Hz)
* Alpha (8–12 Hz)
* Beta (13–30 Hz)

**Channels**: 19, grouped into 6 brain regions

## Analyses
### Power Spectral Density (PSD)
* Welch’s method, 50% overlap
* Cumulative Spectral Analysis (CSA)
* Visualization: waterfall plots, pie charts

### Magnitude-Squared Coherence (MSC)
* *mscohere* function with Hamming windows (50% overlap)
* Visualizations: coherence matrices (*imagesc*) and topographical maps (*topoplot*)
* Implemented in MATLAB (EEGLAB toolbox)
