import streamlit as st
import numpy as np
import plotly.graph_objects as go
from scipy.special import gammaln

def calculate_poisson_samples(lambda_val, sample_size, total_species=350000, mean_copies=30):
    """Calculate number of samples needed for 95% probability of capturing species."""
    max_k = min(1000, int(lambda_val * 5))
    k_range = np.arange(max_k + 1)
    
    # Calculate log poisson probabilities
    log_poisson = k_range * np.log(lambda_val) - lambda_val - gammaln(k_range + 1)
    
    # Calculate probability of missing in one sample
    prob_miss = np.power(1 - k_range/(total_species * mean_copies), sample_size)
    log_prob_miss = np.log(np.maximum(prob_miss, 1e-308))
    
    # Sum the exponentials using log-sum-exp trick
    sum_exp = -np.inf
    for log_term in log_poisson + log_prob_miss:
        sum_exp = np.log(np.exp(sum_exp) + np.exp(log_term))
    
    return int(np.ceil(np.log(0.05) / sum_exp))

def calculate_capture_probability(copy_num, sample_size, num_samples, total_species, mean_copies):
    """Calculate probability of capturing a species with given copy number in any sample."""
    # Probability of missing in a single sample
    prob_miss_single = np.power(1 - copy_num/(total_species * mean_copies), sample_size)
    # Probability of missing in all samples
    prob_miss_all = np.power(prob_miss_single, num_samples)
    # Return probability of capturing in at least one sample
    return (1 - prob_miss_all) * 100

def main():
    st.title("Poisson Sample Size Calculator")
    st.write("""
    Calculate the number of samples needed to capture species with a given copy number,
    and visualize capture probabilities across different copy numbers.
    """)
    
    # Input parameters
    col1, col2 = st.columns(2)
    
    with col1:
        copy_number = st.number_input(
            "Minimum Copy Number per Species (λ)",
            min_value=1,
            value=30
        )
        sample_size = st.number_input(
            "Phage per Sample",
            min_value=1,
            value=8000
        )
        
    with col2:
        total_species = st.number_input(
            "Total Number of Species",
            min_value=1,
            value=350000
        )
        mean_copies = st.number_input(
            "Mean Copies per Species",
            min_value=1,
            value=30
        )
    
    # Calculate results when button is pressed
    if st.button("Calculate"):
        # Calculate required samples
        samples_needed = calculate_poisson_samples(
            copy_number,
            sample_size,
            total_species,
            mean_copies
        )
        
        st.write(f"## Results")
        st.write(f"### {samples_needed:,} samples needed")
        st.write(
            "This will give you a 95% probability of capturing species "
            f"with {copy_number} or more copies."
        )
        
        # Generate probability curve data
        max_copy_num = copy_number * 2
        copy_numbers = np.arange(1, max_copy_num + 1)
        probabilities = [
            calculate_capture_probability(
                cn, sample_size, samples_needed, total_species, mean_copies
            )
            for cn in copy_numbers
        ]
        
        # Create plot
        fig = go.Figure()
        
        # Add probability curve
        fig.add_trace(go.Scatter(
            x=copy_numbers,
            y=probabilities,
            mode='lines',
            name='Detection Probability',
            line=dict(color='#2563eb', width=2)
        ))
        
        # Add 95% reference line
        fig.add_trace(go.Scatter(
            x=[1, max_copy_num],
            y=[95, 95],
            mode='lines',
            name='95% Threshold',
            line=dict(color='#dc2626', width=1, dash='dash')
        ))
        
        # Update layout
        fig.update_layout(
            title='Detection Probability by Copy Number',
            xaxis_title='Copy Number',
            yaxis_title='Detection Probability (%)',
            yaxis_range=[0, 100],
            height=500,
            showlegend=True,
            hovermode='x unified'
        )
        
        st.plotly_chart(fig, use_container_width=True)
        
        # Add formula explanation
        st.write("### Formula Details")
        st.latex(r"n = \left\lceil\frac{\ln(0.05)}{\ln\left(\sum_k \frac{\lambda^k e^{-\lambda}}{k!} \times (1 - \frac{k}{N \times M})^S\right)}\right\rceil")
        st.latex(r"P(\text{detection in n samples}) = 1 - \left(1 - \frac{k}{N \times M}\right)^{S \times n}")
        
        st.write("Where:")
        st.write("""
        - λ (lambda): Minimum copy number per species
        - N: Total number of species
        - M: Mean copies per species
        - S: Sample size (phage per sample)
        - n: Number of samples needed
        - k: Copy number (for probability curve)
        """)

if __name__ == "__main__":
    main()