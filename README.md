# Epidemic-strategies-amid-vaccine-uncertainty

Code for reproducing the results presented in 'When should lockdown be implemented? Devising cost-effective strategies for managing epidemics amid vaccine uncertainty' (N.J. Doyle, F. Cumming, R.N. Thompson, M.J. Tildesley).

All code written in MATLAB, compatible with version R2022b.

Order to run for reproducibility:
1. define_params.m to save variables for simulating model
2. Simulations.m to define the four control strategies and save thresholds
3. dists_joint.m to compute strategy costs for different vaccine outcomes
4. FigS1S2_jointdists.m to generate vaccine joint probability distributions
5. Remaining figure scripts
