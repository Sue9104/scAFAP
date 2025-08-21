import os
import uuid
import socket
import pytensor
unique_id = uuid.uuid4()
hostname = socket.gethostname()
cache_dir = os.path.join(os.path.expanduser('~'), f'.pymc_cache/{hostname}_{unique_id}')
pytensor.config.compiledir = cache_dir
os.environ['NUMBA_CACHE_DIR'] = cache_dir
os.environ['PYTENSOR_CACHE_DIR'] = cache_dir


import warnings
warnings.filterwarnings("ignore", category=UserWarning, message="The figure layout has changed to tight")

import argparse
import arviz as az
import numpy as np
import pymc as pm
import seaborn as sns
import xarray as xr
    
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import norm
from xarray_einstats.stats import XrContinuousRV
import pytensor.tensor as pt
import nutpie

print(f"Running on PyMC v{pm.__version__}")
print(f'PYMC cache dir: {cache_dir}')
# Set random seed for reproducibility
RANDOM_SEED = 8927
rng = np.random.default_rng(RANDOM_SEED)
az.style.use("arviz-darkgrid")
def generate_gmm_data_file(
    filename: str = "simulated_gmm_data.txt",
    n_total_samples: int = 1000,
    component_weights: list[float] = None,
    component_means: list[float] = None,
    component_sds: list[float] = None,
    delimiter: str = ',',
    random_seed: int = 8927
) -> None:
    """
    Generates data from a Gaussian Mixture Model (GMM) and saves it to a file
    in 'value<tab>count' format.

    Args:
        filename (str): The name of the output file. Defaults to "simulated_gmm_data.txt".
        n_total_samples (int): The total number of data points to generate. Defaults to 1000.
        component_weights (list[float]): List of weights for each Gaussian component.
                                         Must sum to 1. If None, defaults to [0.33, 0.67].
        component_means (list[float]): List of means for each Gaussian component.
                                       If None, defaults to [300, 400].
        component_sds (list[float]): List of standard deviations for each Gaussian component.
                                     If None, defaults to [60, 50].
        delimiter (str): The delimiter to use between value and count in the output file.
                         Defaults to a tab ('\t').
        random_seed (int): Seed for the random number generator for reproducibility. Defaults to 8927.
    """
    
    rng = np.random.default_rng(random_seed)

    # Set default parameters if not provided
    if component_weights is None:
        component_weights = [0.33, 0.67]
    if component_means is None:
        component_means = [300, 500]
    if component_sds is None:
        component_sds = [50, 50]

    # Validate inputs
    num_components = len(component_weights)
    if not (num_components == len(component_means) == len(component_sds)):
        raise ValueError("Number of component weights, means, and standard deviations must be equal.")
    if not np.isclose(sum(component_weights), 1.0):
        raise ValueError("Component weights must sum to 1.")
    if n_total_samples <= 0:
        raise ValueError("n_total_samples must be a positive integer.")

    # --- Generate the Data ---
    # Step 1: Choose a component for each data point based on weights
    component_indices = rng.choice(num_components, size=n_total_samples, p=component_weights)

    # Step 2: Generate data points from the chosen components
    # For each data point, sample from the normal distribution corresponding to its chosen component
    simulated_data = rng.normal(
        loc=np.array(component_means)[component_indices],
        scale=np.array(component_sds)[component_indices],
        size=n_total_samples
    )

    print(f"Generated {n_total_samples} raw data points from GMM with {num_components} components.")
    print(f"Overall mean of simulated data: {np.mean(simulated_data):.2f}")
    print(f"Overall std dev of simulated data: {np.std(simulated_data):.2f}")

    # --- Process data into 'value count' format ---
    # Round to create discrete "positions" for the file output
    rounded_data = np.round(simulated_data, 0)
    
    # Get unique positions and their counts
    unique_positions, counts = np.unique(rounded_data, return_counts=True)
    
    # --- Save to File ---
    try:
        with open(filename, "w") as f:
            for pos, count in zip(unique_positions, counts):
                f.write(f"{pos}{delimiter}{count}\n")
        print(f"\nData saved to '{filename}' with '{delimiter}' delimiter.")
        print(f"File contains {len(unique_positions)} unique (rounded) positions.")
    except IOError as e:
        print(f"Error writing to file '{filename}': {e}")
    except Exception as e:
        print(f"An unexpected error occurred during file writing: {e}")
# generate_gmm_data_file(
#     filename="simulated_data.txt",
#     n_total_samples=1000, # A small number for demonstration
#     component_weights=[0.5, 0.5],
#     component_means=[300, 500],
#     component_sds=[50, 50]
# )

# --- 1. Function to Read and Expand Data ---
def read_and_expand_data(file, delimiter=','):
    """
    Reads a file with 'value count' per line and expands it into a 1D NumPy array.

    Args:
        file (str): The path to the input file.
        delimiter (str): The character(s) separating value and count. Defaults to space.

    Returns:
        numpy.ndarray: A 1D NumPy array with values repeated by their counts.
                       Returns an empty array if the file is not found or no valid data.
    """
    expanded_list = []
    data = np.genfromtxt(
        file, delimiter=delimiter, encoding=None, invalid_raise=False,
        dtype=[('value', int), ('count', int)]
    )
    expanded = np.repeat(data['value'], data['count'])
    return(expanded)

# --- 2. Main Script Logic ---
def run_gmm_analysis(input_file_path, clusters, sigma = 50, mu_known = [], output_prefix="gmm", infer = False):
    """
    Performs GMM analysis using PyMC, generates plots, and saves results.

    Args:
        input_file_path (str): Path to the input text file (position, counts).
        clusters (int): Number of Gaussian components (K > 1).
        mu_known (int): mu init values.
        output_prefix (str): Prefix for output files (e.g., gmm_results.txt, gmm_plots.pdf).
    """

    # Read and expand data
    x = read_and_expand_data(input_file_path, delimiter=',')
    if x.size == 0:
        print("Error: No valid data found in the input file. Exiting.")
        return
    if clusters <= 1:
        print("Error: K (number of clusters) must be greater than 1.")
        return
    n_known = len(mu_known)  
    # total_sd = x.std() / clusters
    sigma_init = sigma
    total_mean = x.mean()
    # --- 3. Define the GMM Model ---
    # Using the 'GMM marginal' structure from your notebook, which uses tau
    with pm.Model(coords={"cluster": np.arange(clusters), "obs_id": np.arange(len(x))}) as model:
        # Weights for the mixture components
        w = pm.Dirichlet("w", np.ones(clusters), dims="cluster")
        
        sigma = pm.HalfNormal("sigma", sigma = sigma_init * np.sqrt(np.pi/2), dims="cluster") 
        # Calculate penalty for sigma outside [30, 100]
        penalty_low = pm.math.maximum(0, 30 - sigma)**2
        penalty_high = pm.math.maximum(0, sigma - 100)**2
        penalty = penalty_low + penalty_high
        strength = 5
        pm.Potential("sigma_penalty", - strength * pm.math.sum(penalty) )
    
        # Means of the mixture components
        if (n_known == 0):
            mu = pm.Normal(
                "mu", mu = total_mean, sigma = sigma_init, dims="cluster", 
                transform=pm.distributions.transforms.ordered, 
                initval = np.linspace(x.min(), x.max(), clusters)
            )
        elif (n_known == clusters):
            if infer:
                mu = pm.Normal(
                    "mu", mu = mu_known, sigma = 5, dims="cluster",
                    transform=pm.distributions.transforms.ordered,
                    initval = mu_known
                )
            else:
                mu = pm.Data("mu", mu_known, dims="cluster")
        elif (n_known > 0):
            mu_unknown = pm.Normal(
                "mu_unknown", mu = total_mean, sigma = sigma_init, shape = clusters - n_known
            )
            mu_components = []
            for i in range(n_known):
                mu_components.append(mu_known[i]) 
            for i in range(clusters - n_known):
                mu_components.append(mu_unknown[i]) 
            mu = pm.Deterministic("mu", pt.stack(mu_components), dims="cluster")
        obs = pm.NormalMixture("obs", w=w, mu=mu, sigma = sigma, observed=x, dims="obs_id")
        # components = [pm.Normal.dist(mu=mu[i], tau=tau[i]) for i in range(clusters)]
        # likelihood = pm.Mixture("likelihood", w=w, comp_dists=components, observed=x, dims="obs_id")
        
    # pm.model_to_graphviz(model) # Commented out for script, but useful for debugging

    # --- 4. Sample from the Posterior ---
    print("\nStarting MCMC sampling...")
    with model:
        try:
            trace = pm.sample(
                draws=2000, tune=1000, random_seed=RANDOM_SEED,
                nuts_sampler='nutpie', 
                n_init = 3000, init='auto', cores = 2
            )
        except RuntimeError as e:
            if "All initialization points failed" in str(e):
                print("Default initialization failed. Retrying with 'advi'...")
                trace = pm.sample(
                    draws=2000, tune=1000, random_seed=RANDOM_SEED,
                    nuts_sampler='nutpie', 
                    n_init = 3000, init='advi', cores = 2
                )
    #     # Sample posterior predictive samples, slow, could be hided if not plot_ppc
    #     # ppc_trace = pm.sample_posterior_predictive(trace, extend_inferencedata=True, random_seed=RANDOM_SEED)

    az.summary(trace).to_csv(f"{output_prefix}.model.csv")
    print("MCMC sampling complete.")

    # --- 5. Generate Posterior Probability File ---
    print(f"\nGenerating posterior probabilities and saving to '{output_prefix}.probs.txt'...")
    xi_coord_name = "x_pdf"
    # Create xi_array covering data range with a buffer, for smoother PDFs
    linspace_values = np.linspace(x.min(), x.max(), len(x)).round()
    xi_array = xr.DataArray(linspace_values, coords={xi_coord_name: linspace_values}, dims=[xi_coord_name])

    post = trace.posterior
    # Calculate PDF components (weighted individual Gaussians)
    if (n_known == clusters):
        mu = mu_known
    else:
        mu = post["mu"]
    pdf_components = XrContinuousRV(norm, mu, post["sigma"]).pdf(xi_array) * post["w"]
    pdf_components = pdf_components.mean(dim=["chain", "draw"])
    # Calculate overall PDF (sum of components)
    pdf = pdf_components.sum("cluster")
    member_probs = (pdf_components / pdf)
    
    output_data_path = f"{output_prefix}.probs.txt"
    with open(output_data_path, 'w') as f:
        header = f"xi_array\t" + "\t".join([f"Cluster_{i}_Prob" for i in range(clusters)])
        f.write(header + "\n")
        # Write data rows
        for i, xi_val in enumerate(xi_array.values):
            current_row_probs = member_probs.isel({xi_coord_name: i}).values
            line_data = [f"{int(xi_val)}"] + [f"{prob:.4f}" for prob in current_row_probs]
            f.write("\t".join(line_data) + "\n")
    print(f"Posterior probabilities saved to '{output_data_path}'.")


    # --- 6. Generate All Plots into a Single PDF ---
    print(f"Generating plots and saving to '{output_prefix}.pdf'...")
    pdf_plots_path = f"{output_prefix}.pdf"

    with PdfPages(pdf_plots_path) as pdf_saver:
        # --- Plot 2: Arviz Trace Plot ---
        fig2 = az.plot_trace(trace, var_names=["w", "sigma"]).ravel()[0].figure
        fig2.set_size_inches(9, 6)
        fig2.suptitle("MCMC Trace Plots", y=1.02) # Add a title to the figure
        plt.tight_layout()
        pdf_saver.savefig(fig2)
        plt.close(fig2)

        # --- Plot 3: Arviz Posterior Plot ---
        plt.rcParams['font.size'] = 8
        axes3 = az.plot_posterior(trace, var_names=["w", "sigma"], grid=(2, clusters))
        fig3 = axes3.ravel()[0].figure
        fig3.set_size_inches(9, 6)
        fig3.suptitle("Posterior Distributions", y=1.02)
        for ax_row in axes3:
            for ax in ax_row:
                # 调整子图标题（即变量名 'w', 'mu', 'tau'）
                ax.set_title(ax.get_title(), fontsize=8)
                # 调整轴标签（例如 'Sample', 'Value', 'Density'）
                ax.set_xlabel(ax.get_xlabel(), fontsize=8)
                ax.set_ylabel(ax.get_ylabel(), fontsize=8)
                # 调整刻度标签（轴上的数字）
                ax.tick_params(axis='x', labelsize=6)
                ax.tick_params(axis='y', labelsize=6)
        plt.tight_layout()
        pdf_saver.savefig(fig3)
        plt.close(fig3)

        # --- Plot 4: Arviz PPC Plot ---
        # slow, could be hided if not plot_ppc
        # axes_ppc = az.plot_ppc(trace)
        # if isinstance(axes_ppc, np.ndarray):
        #     fig4 = axes_ppc.flatten()[0].figure
        # else: # It's a single Axes object
        #     fig4 = axes_ppc.figure
        # fig4.set_size_inches(9, 6)
        # fig4.suptitle("Posterior Predictive Check", y=1.02)
        # plt.tight_layout()
        # pdf_saver.savefig(fig4)
        # plt.close(fig4)

        # --- Plot 5: Combined PDF and Membership Probabilities ---
        fig5, ax = plt.subplots(4, 1, figsize=(9, 10), sharex=True)

        # 0: Empirical histogram (re-doing for the combined plot)
        ax[0].hist(x, bins=20, density=True, alpha=0.7, color='skyblue')
        ax[0].set(title="Raw Data", ylabel="Density")
        for mu_i in mu_known:
            ax[0].axvline(mu_i, color='red', linestyle='--', label=f'Known μ: {mu_i:.2f}')
        # 1: Individual PDF Components
        pdf_components.plot.line(x=xi_coord_name, hue="cluster", ax=ax[1])
        ax[1].set(title="Posterior Mean of Individual PDF", ylabel="Probability")
        # 2: Overall PDF (Sum of components)
        pdf.plot.line(x=xi_coord_name, ax=ax[2], color='darkorange')
        ax[2].set(title="Posterior Mean of Overall PDF", ylabel="Probability")
        for mu_i in mu_known:
            ax[2].axvline(mu_i, color='red', linestyle='--', label=f'Known μ: {mu_i:.2f}')
        # 3: Group membership probabilities
        # mean_group_membership_probs already calculated above
        member_probs.plot.line(x=xi_coord_name, hue="cluster", ax=ax[3])
        ax[3].set(title="Posterior Mean of Group Membership Probabilities", xlabel="x", ylabel="Probability")
        # Ensure x-axis ticks are visible on all plots, but labels only on the bottom
        for i in range(len(ax)):
            ax[i].set_xlabel("")
            ax[i].tick_params(axis='x', which='both', bottom=True, top=True, labelbottom=False)
        ax[-1].tick_params(axis='x', labelbottom=True)
        ax[-1].set_xlabel("x")
        plt.tight_layout()
        pdf_saver.savefig(fig5)
        plt.close(fig5)

    print("\nAnalysis complete.")

# --- 7. Command Line Argument Parsing ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Perform Bayesian Gaussian Mixture Model (GMM) analysis on 2-column text data."
    )
    parser.add_argument(
        "infile",
        type=str,
        help="Input text file path (two columns: position<tab>counts)."
    )
    parser.add_argument(
        "K",
        type=int,
        help="Number of Gaussian components (K > 1)."
    )
    parser.add_argument(
        "--sigma",
        type=int,
        default=50,
        help="Evaluated sd (e.g., --sigma 50)."
    )
    parser.add_argument(
        "--init",
        nargs='*',  # Accepts more than 1 arguments
        type=int,
        default=[],
        help="Required list of integers (e.g., --init 1 2)."
    )
    parser.add_argument(
        "--infer",
        action='store_true',
        help="Enable inference mode. Defaults to False."
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default="gmm_analysis",
        help="Prefix for output files (e.g., gmm_analysis.probs.txt, gmm_analysis.pdf)."
    )
    
    args = parser.parse_args()
    run_gmm_analysis(args.infile, args.K, args.sigma, np.array(args.init), args.prefix, args.infer)
