import pymc as pm
import numpy as np
import arviz as az
import pytensor.tensor as pt
import matplotlib.pyplot as plt

print(f"Running on PyMC v{pm.__version__}")

# ----------------------------------------------------------------------
# Part 1: System Representation & Synthetic Data Generation
# (Includes Steps 1, 2, 3, 4 from your prompt)
# ----------------------------------------------------------------------

# --- Ground Truth Parameters (THE "ANSWER KEY" FOR THE INFERENCE) ---
n_particles = 10000 # Your suggestion: 10,000 particles
particle_box_size = 100.0 # e.g., particles in a 100x100x100 Å cube

n_total_crosslinks = 200 # Total crosslinks to generate
true_fp_rate = 0.30     # 30% of crosslinks will be false positives
true_overlength_rate = 0.10 # 10% of *true* crosslinks will be "over-length"
true_sigma = 10.0       # The "fuzziness" (std dev) of over-length links (in Å)

l_tri = 25.0            # The "true" max length (in Å) for all 3 arms of the linker

# --- Calculate data composition ---
n_fp = int(n_total_crosslinks * true_fp_rate)
n_tp = n_total_crosslinks - n_fp
n_tp_overlength = int(n_tp * true_overlength_rate)
n_tp_good = n_tp - n_tp_overlength

print(f"Generating synthetic trifunctional data:")
print(f"  {n_total_crosslinks} total crosslinks")
print(f"    - {n_fp} False Positives (FPs)")
print(f"    - {n_tp} True Positives (TPs), which are:")
print(f"        - {n_tp_good} 'Good' TPs")
print(f"        - {n_tp_overlength} 'Over-length' TPs")

# --- 1. System Representation (Ground Truth "Universe") ---
particles_xyz = np.random.rand(n_particles, 3) * particle_box_size

def get_triplet_distances(p_idx, particles):
    """Helper function to get all 3 pairwise distances of a triplet."""
    p_a = particles[p_idx[0]]
    p_b = particles[p_idx[1]]
    p_c = particles[p_idx[2]]
    d_ab = np.linalg.norm(p_a - p_b)
    d_ac = np.linalg.norm(p_a - p_c)
    d_bc = np.linalg.norm(p_b - p_c)
    return d_ab, d_ac, d_bc

# --- 2, 3, 4. Generate & Corrupt Synthetic Data, Simulate ID Scores ---
observed_distances_ab = []
observed_distances_ac = []
observed_distances_bc = []
observed_id_scores = []
data_types = [] # To keep track for plotting/validation

# 2a. Generate "Good" True Triplets
print(f"  Generating {n_tp_good} 'Good' True Positives...")
for _ in range(n_tp_good):
    while True:
        p_idx = np.random.choice(n_particles, 3, replace=False)
        d_ab, d_ac, d_bc = get_triplet_distances(p_idx, particles_xyz)
        # Condition for a "good" trifunctional crosslink
        if d_ab < l_tri and d_ac < l_tri and d_bc < l_tri:
            observed_distances_ab.append(d_ab)
            observed_distances_ac.append(d_ac)
            observed_distances_bc.append(d_bc)
            observed_id_scores.append(np.clip(np.random.normal(0.85, 0.1), 0, 1)) # High ID score
            data_types.append("TP_Good")
            break

# 2b & 3b. Generate "Over-length" True Triplets
print(f"  Generating {n_tp_overlength} 'Over-length' True Positives...")
for _ in range(n_tp_overlength):
    while True:
        p_idx = np.random.choice(n_particles, 3, replace=False)
        d_ab, d_ac, d_bc = get_triplet_distances(p_idx, particles_xyz)
        # Initially satisfy the geometry
        if d_ab < l_tri and d_ac < l_tri and d_bc < l_tri:
            # Corrupt one distance by adding noise to simulate "over-length"
            # Your prompt: "move particle C so that d(B,C) is now 25 + max_over-length_deviation"
            # We'll simulate this by adding absolute normal noise to d_bc
            d_bc_corrupt = d_bc + np.abs(np.random.normal(0, true_sigma))

            observed_distances_ab.append(d_ab)
            observed_distances_ac.append(d_ac)
            observed_distances_bc.append(d_bc_corrupt) # This one is corrupted
            observed_id_scores.append(np.clip(np.random.normal(0.80, 0.1), 0, 1)) # High ID score (slightly lower on average)
            data_types.append("TP_Overlength")
            break

# 3a. Generate False Positives
print(f"  Generating {n_fp} False Positives...")
for _ in range(n_fp):
    p_idx = np.random.choice(n_particles, 3, replace=False)
    d_ab, d_ac, d_bc = get_triplet_distances(p_idx, particles_xyz)
    observed_distances_ab.append(d_ab)
    observed_distances_ac.append(d_ac)
    observed_distances_bc.append(d_bc)
    observed_id_scores.append(np.clip(np.random.normal(0.2, 0.1), 0, 1)) # Low ID score
    data_types.append("FP")

d_ab_obs = np.array(observed_distances_ab)
d_ac_obs = np.array(observed_distances_ac)
d_bc_obs = np.array(observed_distances_bc)
id_scores_obs = np.array(observed_id_scores)
print("Synthetic data generation complete.")

# --- 1c. Empirically Calculate True Beta (Random-match probability) ---
# This is crucial for fixing beta_tri in the model
print("Empirically calculating true beta_tri (random-match probability)...")
n_test_triplets = 100_000
satisfied_count = 0
for _ in range(n_test_triplets):
    p_idx = np.random.choice(n_particles, 3, replace=False)
    d_ab, d_ac, d_bc = get_triplet_distances(p_idx, particles_xyz)
    if d_ab < l_tri and d_ac < l_tri and d_bc < l_tri:
        satisfied_count += 1
true_beta_empirical = satisfied_count / n_test_triplets
print(f"  Empirical Beta_tri = {true_beta_empirical:.8f}") # Beta_tri will be very small

# ----------------------------------------------------------------------
# Part 2: Bayesian Model Definition (PyMC) - Step 5
# ----------------------------------------------------------------------

with pm.Model() as trifunctional_model:
    # --- Nuisance Parameters (Priors) ---
    
    # Sigma for the "Flexible" / "Over-length" True Positives
    # (True value = true_sigma = 10.0)
    sigma_flexible = pm.HalfNormal("sigma_flexible", sigma=20.0, initval=true_sigma) # A slightly wider prior for more exploration
    
    # Sigma for the "Good" / "Rigid" True Positives
    # Fixed to a small value, as per discussion, representing a sharp cutoff
    # A small HalfNormal prior allows for slight variation but keeps it constrained
    sigma_good = pm.HalfNormal("sigma_good", sigma=2.0, initval=2.0)

    # Probability that a TP is "Flexible" (our true_overlength_rate)
    # (True value = true_overlength_rate = 0.10)
    p_flexible = pm.Beta("p_flexible", alpha=2.0, beta=10.0, initval=true_overlength_rate) # Prior centered around 0.2

    # Beta_tri (FP random match probability) is a fixed value
    beta_tri = true_beta_empirical

    # --- Isotonic Regression (Logistic Model for Calibrated Alpha) ---
    # The 'id_score' provides an empirical estimate of precision.
    # We calibrate this with a logistic function.
    logistic_intercept = pm.Normal("logistic_intercept", mu=-3.0, sigma=2.0) # Adjusted prior for potentially wider range
    logistic_slope = pm.HalfNormal("logistic_slope", sigma=3.0, initval=5.0) # Adjusted prior
    
    calibrated_alpha = pm.math.sigmoid(
        logistic_intercept + logistic_slope * id_scores_obs
    )
    pm.Deterministic("calibrated_alpha", calibrated_alpha)

    # --- The Likelihood Function: Mixture Model for True Positives ---
    
    # P(Data | FP): The probability of a False Positive satisfying the geometry by chance
    p_satisfied_fp = beta_tri

    # P(Data | TP): This is a mixture of "good" and "flexible" True Positives
    
    # Stack all observed distances for convenience
    all_distances = pt.stack([d_ab_obs, d_ac_obs, d_bc_obs], axis=1)
    
    # Find the maximum distance for each crosslink triplet
    # This is critical for the "trifunctional geometry" constraint
    d_max_obs = pt.max(all_distances, axis=1)
    
    # Probability for the "Good" TP component (sharp cutoff)
    # invlogit ensures probability between 0 and 1. l_tri - d_max > 0 means within bounds.
    p_tp_good_component = pm.math.invlogit((l_tri - d_max_obs) / sigma_good)

    # Probability for the "Flexible" TP component (wider cutoff)
    p_tp_flexible_component = pm.math.invlogit((l_tri - d_max_obs) / sigma_flexible)

    # Combine the two TP components using p_flexible
    p_satisfied_tp = (1.0 - p_flexible) * p_tp_good_component + p_flexible * p_tp_flexible_component
    
    # --- The Final Combined Likelihood ---
    # P(Observed '1' | alpha, p_satisfied_tp, p_satisfied_fp)
    # = alpha * P(Observed '1' | TP) + (1 - alpha) * P(Observed '1' | FP)
    likelihood_per_xl = (
        calibrated_alpha * p_satisfied_tp +
        (1.0 - calibrated_alpha) * p_satisfied_fp
    )

    # --- The "Observation" ---
    # We observe all crosslinks, so the observed value for each Bernoulli trial is 1.
    obs = pm.Bernoulli(
        "obs",
        p=likelihood_per_xl,
        observed=np.ones(n_total_crosslinks)
    )

print("Trifunctional Model definition complete.")

# ----------------------------------------------------------------------
# Part 3: Sampling (MCMC) - Step 5 continued
# ----------------------------------------------------------------------
print("Starting MCMC sampling...")
with trifunctional_model:
    # Use higher target_accept to handle potentially complex posterior landscape
    # Increased tune and draws for better convergence of a more complex model
    idata = pm.sample(
        3000, 
        tune=2000, 
        cores=4, 
        target_accept=0.95, 
        init="adapt_diag"
    )
print("Sampling complete.")

# ----------------------------------------------------------------------
# Part 4: Analysis & Validation - Step 6
# ----------------------------------------------------------------------
print("Analyzing results...")

# --- 6a. Check Parameter Recovery ---
summary = az.summary(
    idata, 
    var_names=[
        "sigma_flexible", 
        "sigma_good", 
        "p_flexible", 
        "logistic_slope", 
        "logistic_intercept"
    ]
)
print(summary)

print("\n--- Ground Truth vs. Inferred ---")
print(f"sigma_flexible: True={true_sigma:.2f}, Inferred={summary.loc['sigma_flexible', 'mean']:.2f} (HDI_94%=[{summary.loc['sigma_flexible', 'hdi_3%']:.2f}, {summary.loc['sigma_flexible', 'hdi_97%']:.2f}])")
print(f"p_flexible: True={true_overlength_rate:.2f}, Inferred={summary.loc['p_flexible', 'mean']:.2f} (HDI_94%=[{summary.loc['p_flexible', 'hdi_3%']:.2f}, {summary.loc['p_flexible', 'hdi_97%']:.2f}])")
print(f"beta_tri (fixed): {true_beta_empirical:.8f}")

# Trace plots for visual convergence check
az.plot_trace(idata, var_names=["sigma_flexible", "p_flexible", "logistic_slope", "logistic_intercept"])
plt.suptitle("Parameter Trace Plots (Check for Convergence)", y=1.02)
plt.tight_layout()
plt.show()

# --- 6b. The "Isotonic Regression" Plot ---
mean_calibrated_alpha = idata.posterior["calibrated_alpha"].mean(
    dim=("chain", "draw")
)
# Define color map for different data types
color_map = {"TP_Good": "green", "TP_Overlength": "blue", "FP": "red"}
colors = [color_map[dt] for dt in data_types]

plt.figure(figsize=(10, 6))
plt.scatter(
    id_scores_obs,
    mean_calibrated_alpha,
    c=colors,
    alpha=0.7,
    label="Inferred Precision (per XL)"
)
plt.axvline(0.5, color='gray', linestyle='--', label='Approx. TP/FP cutoff')
plt.axhline(1.0 - true_fp_rate, color='black', linestyle=':', label=f'True Avg. Precision ({1.0 - true_fp_rate:.2f})')

plt.title("Model Validation: Calibrated Precision (α) vs. ID Score for Trifunctional XLs")
plt.xlabel("Synthetic ID Score (Input)")
plt.ylabel("Inferred Calibrated Precision (α)")
plt.legend(handles=[
    plt.Line2D([0], [0], marker='o', color='w', mfc='green', label='TP (Good)'),
    plt.Line2D([0], [0], marker='o', color='w', mfc='blue', label='TP (Over-length)'),
    plt.Line2D([0], [0], marker='o', color='w', mfc='red', label='FP (Noise)')
])
plt.grid(True, linestyle='--', alpha=0.6)
plt.show()

# --- 6c. Optional: Histogram of d_max_obs ---
# This helps visualize the data and confirm the two TP populations
plt.figure(figsize=(10, 6))
d_max_all = np.array([max(d) for d in zip(d_ab_obs, d_ac_obs, d_bc_obs)])
plt.hist(d_max_all[np.array(data_types) == "TP_Good"], bins=30, alpha=0.6, label="TP (Good)", color="green", density=True)
plt.hist(d_max_all[np.array(data_types) == "TP_Overlength"], bins=30, alpha=0.6, label="TP (Over-length)", color="blue", density=True)
plt.hist(d_max_all[np.array(data_types) == "FP"], bins=30, alpha=0.6, label="FP (Noise)", color="red", density=True)
plt.axvline(l_tri, color='purple', linestyle='--', label=f'L_TRI = {l_tri:.1f} Å')
plt.title("Distribution of Max Pairwise Distances (d_max) for Triplet Crosslinks")
plt.xlabel("Max Pairwise Distance (Å)")
plt.ylabel("Density")
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)
plt.show()