## Method for cluster correcting hierarchical/linear mixed effects models in intracranial EEG timeseries

Our research aims to uncover mechanistic explanations of the neural basis of human behavior, that is, move from where to how. Our goals are multifaceted: (1) advance fundamental science by discovering new knowledge using rigorous, reproducible methods; and (2) advance translational applications in neurotechnology, precision medicine, and product development that are grounded in rigorous science. Intracranial EEG (iEEG) is an invaluable tool to study the human brain because it provides data with the high spatiotemporal resolution and signal quality traditionally limited to animal neurophysiology. Analyzing iEEG data on the group level is complex because every neurosurgical patient dataset is different and sample sizes tend to be small. By modeling random effects - that is, iEEG channels (sub-samples of a brain region) and/or experimental trials (sub-samples of a task) nested in subjects, hierarchical models can harmonize iEEG data while maximizing statistical power. This particularly relevant for studies of inter-individual variability, as described in:
- Johnson, EL, Knight, RT. How can iEEG be used to study inter-individual and developmental differences? in _Intracranial EEG: A Guide for Cognitive Neuroscientists_ (ed N Axmacher) (Springer, Cham, 2023). [DOI](https://doi.org/10.1007/978-3-031-20910-9_10)

Because the data being modeled are often timeseries (or frequency series or time-frequency representations), the model outputs must be corrected for multiple comparisons. Nonparametric clustering is recommended for multiple comparison correction in neural data because the test statistic incorporates biophysically motivated constraints (i.e., adjacent datapoints are not independent) and captures any biases in the data. These scripts provide a clustering method for hierarchical models in iEEG timeseries (or frequency series) that preserves the random effects structure. They may be adjusted to work with time-frequency representations.

The method is described in:
- Dede, AJO, Cross, ZR, Gray, SM, Kelly, JP, Yin, Q, Vahidi, P, Asano, E, Schuele, SU, Rosenow, JM, Wu, JY, Lam, SK, Raskin, JS, Lin, JJ, Kim McManus, O, Sattar, S, Shaikhouni, A, King-Stephens, D, Weber, PB, Laxer, KD, Brunner, P, Roland, JL, Saez, I, Girgis, F, Knight, RT, Ofen, N, Johnson, EL. Declarative memory through the lens of single-trial peaks in high frequency power. _bioRxiv_ (2025). [DOI](https://doi.org/10.1101/2025.01.02.631123)

Publications or other papers using these scripts should cite the publication(s) above.

Software:
- MATLAB 9.7 (last tested with R2024b)

Notes:
- Use the function cluster_1way for models with one fixed effect. It was designed for the fixed effect of hit vs. miss on a memory task, and may be used as a template for fixed effects with more than two levels.
- Use the function cluster_2way for models with two fixed effects and their interaction. It was designed for the fixed effects of (1) hit vs. miss on a memory task and (2) brain region, and may be used as a template for models with more than two fixed effects and/or fixed effects with more than two levels.
- Both functions call the subfunction cluster_test.
