## Method for cluster correcting hierarchical/linear mixed effects models in intracranial EEG timeseries

Our research aims to uncover mechanistic explanations of the neural basis of human behavior, that is, move from where to how. Our goals are multifaceted: (1) advance fundamental science by discovering new knowledge using rigorous, reproducible methods; and (2) advance translational applications in neurotechnology, precision medicine, and product development that are grounded in rigorous science. Intracranial EEG (iEEG) is an invaluable tool to study the human brain because it provides data with the high spatiotemporal resolution and signal quality traditionally limited to animal neurophysiology. Analyzing iEEG data on the group level is complex because every neurosurgical patient dataset is different and sample sizes tend to be small. By modeling random effects - that is, iEEG channels (sub-samples of a brain region) and/or experimental trials (sub-samples of a task) nested in subjects, hierarchical models can harmonize iEEG data while maximizing statistical power.  

Because the data being modeled are often timeseries (or frequency series or time-frequency representations), the model outputs must be corrected for multiple comparisons. Nonparametric clustering is recommended for multiple comparison correction in neural data because the test statistic incorporates biophysically motivated constraints (i.e., adjacent datapoints are not independent) and captures any biases in the data. These scripts provide a clustering method for hierarchical models in iEEG timeseries (or frequency series) that preserves the random effects structure. They may be adjusted to work with time-frequency representations.

The method is described in:
- Dede, AJO, Cross, ZR, Gray, SM, Kelly, JP, Yin, Q, Vahidi, P, Asano, E, Schuele, SU, Rosenow, JM, Wu, JY, Lam, SK, Raskin, JS, Lin, JJ, Kim McManus, O, Sattar, S, Shaikhouni, A, King-Stephens, D, Weber, PB, Laxer, KD, Brunner, P, Roland, JL, Saez, I, Girgis, F, Knight, RT, Ofen, N, Johnson, EL. Declarative memory through the lens of single-trial peaks in high frequency power. _bioRxiv_ (2025). [DOI](https://doi.org/10.1101/2025.01.02.631123)

Publications or other papers using these scripts should cite the publication above.

### Run using Python or MATLAB:  

Python:
- Software: Python 3.12.12
- Environment: Google Colab
- Package versions:
  - NumPy 2.0.2
  - Pandas 2.2.2
  - SciPy 1.16.2
  - Statsmodels 0.14.5
- Notebooks:
  - cluster_1way.ipynb
  - cluster_2way.ipynb
- Subfunctions: cluster_utils.py
  - Upload to Colab File panel

MATLAB:
- Software: MATLAB 9.7 (last tested with R2024b)
- Scripts:
  - cluster_1way.m
  - cluster_2way.m
- Subfunctions: subfunctions.zip
  - Unzip folder

Notes:
- Use cluster_1way for models with one fixed effect with two levels. It was designed for the fixed effect of hit vs. miss on a memory task, and may be used as a template for other models with an equivalent design.
- Use cluster_2way for models with two fixed effects and their interaction. It was designed for the fixed effects of (1) hit vs. miss on a memory task and (2) one brain region vs. another, and may be used as a template for models with more than two fixed effects and/or fixed effects with more than two levels.
- Running 1000 permutations requires fitting 1000 models per timepoint (sample data file: 139 timepoints × 1000 permutations = 139,000 models), which can take several hours or more depending on hardware allocation. Set nperm to a lower number for prototyping. 
