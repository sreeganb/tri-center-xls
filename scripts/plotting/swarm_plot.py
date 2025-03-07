import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np

# Set font scale for better visibility
sns.set_context("notebook", font_scale=1.5)

# Read in the data
df_trifunctional_20 = pd.read_csv('input_data/trifunctional_20.csv', index_col=0).T
df_trifunctional_30 = pd.read_csv('input_data/trifunctional_30.csv', index_col=0).T
df_trifunctional_40 = pd.read_csv('input_data/trifunctional_40.csv', index_col=0).T
df_trifunctional_60 = pd.read_csv('input_data/trifunctional_60.csv', index_col=0).T
df_trifunctional_120 = pd.read_csv('input_data/trifunctional_120.csv', index_col=0).T

df_bi_20 = pd.read_csv('input_data/bifunctional_20.csv', index_col=0).T
df_bi_30 = pd.read_csv('input_data/bifunctional_30.csv', index_col=0).T
df_bi_40 = pd.read_csv('input_data/bifunctional_40.csv', index_col=0).T
df_bi_60 = pd.read_csv('input_data/bifunctional_60.csv', index_col=0).T
df_bi_120 = pd.read_csv('input_data/bifunctional_120.csv', index_col=0).T

# Check required columns
required_columns = ['Model Accuracy', 'Cluster Precision']
dataframes = [
    df_trifunctional_20, df_trifunctional_30, df_trifunctional_40, df_trifunctional_60, df_trifunctional_120,
    df_bi_20, df_bi_30, df_bi_40, df_bi_60, df_bi_120
]

for df in dataframes:
    for col in required_columns:
        if col not in df.columns:
            raise KeyError(f"'{col}' not found in one of the DataFrames.")

# Experimental data for Model Accuracy and Cluster Precision
experimental_data = {
    'Model Accuracy': {
        'TSTO': [15.09, 16.34, 12.30, 19.72, 15.23],
        'DSSO': [33.2]
    },
    'Cluster Precision': {
        'TSTO': [3.43, 7.8, 1.34, 3.7, 1.36],
        'DSSO': [12.0]
    }
}

# Function to create plots
def create_plots(metric):
    # Combine main plot data
    data_main = {
        'Bifunctional 20 XLs': df_bi_20[metric],
        'Trifunctional 20 XLs': df_trifunctional_20[metric],
        'Bifunctional 30 XLs': df_bi_30[metric],
        'Trifunctional 30 XLs': df_trifunctional_30[metric],
        'Bifunctional 40 XLs': df_bi_40[metric],
        'Trifunctional 40 XLs': df_trifunctional_40[metric],
        'Bifunctional 60 XLs': df_bi_60[metric],
        'Trifunctional 60 XLs': df_trifunctional_60[metric],
        'Bifunctional 120 XLs': df_bi_120[metric],
        'Trifunctional 120 XLs': df_trifunctional_120[metric],
    }

    df_main = pd.DataFrame(data_main)

    # Experimental data
    tsto_data = experimental_data[metric]['TSTO']
    dsso_data = experimental_data[metric]['DSSO']

    # Create DataFrames
    df_TSTO = pd.DataFrame({
        'Type': 'TSTO',
        metric: tsto_data,
        'Number': '83'
    })

    df_DSSO = pd.DataFrame({
        'Type': 'DSSO',
        metric: dsso_data,
        'Number': '83'
    })

    # Combine the DataFrames
    df_exp = pd.concat([df_DSSO, df_TSTO], ignore_index=True)

    # Create a figure with 3 subplots sharing the y-axis
    fig, (ax1, ax2, ax3) = plt.subplots(
        1, 3, figsize=(10, 4), sharey=True,
        gridspec_kw={'width_ratios': [3.5, 0.5, 1]}
    )

    # Define the custom color palette
    custom_palette = {
        #'Bifunctional': '#1b85b8',
        'Bifunctional': '#0000FF',
        #'Trifunctional': '#FF6961',
        'Trifunctional': '#FF0000',
        #'TSTO': '#CF9FFF',
        'TSTO': '#FF00FF',
        #'DSSO': '#2ca02c'
        'DSSO': '#2e8bc0'
    }
    
    custom_palette = {
    'Bifunctional': '#007BFF',
    'Trifunctional': '#FF4500',
    'TSTO': '#FF00FF',
    'DSSO': '#00CC99'
}

    ## Plot 1: Synthetic Data (20, 30, 40, 60 XLs)
    df_main_melted = df_main.melt(var_name='Crosslinks', value_name=metric)
    df_main_melted['Type'] = df_main_melted['Crosslinks'].apply(lambda x: x.split()[0])
    df_main_melted['Number'] = df_main_melted['Crosslinks'].apply(lambda x: x.split()[1])

    # Filter data for the first plot
    df_main_melted_1 = df_main_melted[df_main_melted['Number'].isin(['20', '30', '40', '60'])]
    df_main_melted_1.loc[:, 'Number'] = df_main_melted_1['Number'].astype(int)

    sns.stripplot(
        x='Number',
        y=metric,
        hue='Type',
        data=df_main_melted_1,
        palette=custom_palette,
        dodge=True,
        jitter=False,
        marker='o',
        edgecolor='black',
        linewidth=0.5,
        size=12,
        ax=ax1
    )
    pointprops = dict(marker='D', markerfacecolor='#1b85b8',
                      markeredgecolor='black', markersize=6, linestyle='--')
    # Overlay boxplot to show the mean
    sns.boxplot(
        x='Number',
        y=metric,
        hue='Type',
        data=df_main_melted_1,
        palette=custom_palette,
        showmeans=True,
        meanline=True,
        meanprops={'color': 'k', 'ls': '-', 'lw': 2},
        medianprops={'visible': False},
        whiskerprops={'visible': False},
        zorder=10,
        showcaps=False,
        boxprops={'visible': False},
        showfliers=False,
        dodge=True,
        ax=ax1
    )

    ax1.set_title('Synthetic Data')
    ax1.set_xlabel('No. of XL sites')
    ax1.set_ylabel(f"{metric}" r"(${\AA}$)")
    ax1.grid(True)
    # Remove the duplicated legend
    if ax1.get_legend() is not None:
        ax1.legend_.remove()

    ## Plot 2: Synthetic Data (120 XLs)
    df_main_melted_2 = df_main_melted[df_main_melted['Number'] == '120']
    df_main_melted_2.loc[:, 'Number'] = df_main_melted_2['Number'].astype(int)

    sns.stripplot(
        x='Number',
        y=metric,
        hue='Type',
        data=df_main_melted_2,
        palette=custom_palette,
        dodge=True,
        jitter=False,
        marker='o',
        edgecolor='black',
        linewidth=0.5,
        size=12,
        ax=ax2
    )

    # Overlay boxplot to show the mean
    sns.boxplot(
        x='Number',
        y=metric,
        hue='Type',
        data=df_main_melted_2,
        palette=custom_palette,
        showmeans=True,
        meanline=True,
        meanprops={'color': 'k', 'ls': '-', 'lw': 2},
        medianprops={'visible': False},
        whiskerprops={'visible': False},
        zorder=10,
        showcaps=False,
        boxprops={'visible': False},
        showfliers=False,
        dodge=True,
        ax=ax2
    )

    ax2.set_xlabel('')
    ax2.set_ylabel('')
    ax2.grid(True)
    # Remove the duplicated legend
    if ax2.get_legend() is not None:
        ax2.legend_.remove()

    # Hide the spines and ticks between the plots to create the illusion of a broken axis
    ax1.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax1.yaxis.tick_left()
    ax2.yaxis.tick_right()
    ax2.tick_params(labelleft=False)

    d = .01  # how big to make the diagonal lines in axes coordinates
    kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
    ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)
    ax1.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

    kwargs.update(transform=ax2.transAxes)
    ax2.plot((-d, +d), (-d, +d), **kwargs)
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)

    ## Plot 3: Experimental Data
    sns.stripplot(
        x='Type',
        y=metric,
        data=df_exp,
        palette=custom_palette,
        order=['DSSO', 'TSTO'],
        dodge=False,
        jitter=False,
        marker='o',
        edgecolor='black',
        linewidth=0.5,
        size=12,
        ax=ax3
    )

    # Overlay boxplot to show the mean
    sns.boxplot(
        x='Type',
        y=metric,
        data=df_exp,
        palette=custom_palette,
        order=['DSSO', 'TSTO'],
        showmeans=True,
        meanline=True,
        meanprops={'color': 'k', 'ls': '-', 'lw': 2},
        medianprops={'visible': False},
        whiskerprops={'visible': False},
        zorder=10,
        showcaps=False,
        boxprops={'visible': False},
        showfliers=False,
        ax=ax3
    )

    ax3.set_title('Experimental Data')
    ax3.set_xlabel('Crosslinker')
    ax3.set_ylabel('')
    ax3.grid(True)
    # Remove the duplicated legend, if any
    if ax3.get_legend() is not None:
        ax3.legend_.remove()

    # Adjust legends
    handles0, labels0 = ax1.get_legend_handles_labels()
    ax1.legend(
        handles0[:2], labels0[:2], title='Type',
        loc='upper right'
    )
    ax2.legend().set_visible(False)

    # Reduce space between subplots
    plt.subplots_adjust(wspace=0.04)
    plt.tight_layout()

    # Save the plot
    output_dir = 'output_files'
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(
            os.path.join(output_dir, f'stripplot_with_mean_{metric.lower().replace(" ", "_")}.svg'),
            format='svg',
            dpi=300
        )
    plt.show()

    # Calculate means and write to file and print to terminal

    # Collect all datasets and their means
    mean_results = []

    # Synthetic Data
    datasets = {
        'Bifunctional 20 XLs': df_bi_20,
        'Trifunctional 20 XLs': df_trifunctional_20,
        'Bifunctional 30 XLs': df_bi_30,
        'Trifunctional 30 XLs': df_trifunctional_30,
        'Bifunctional 40 XLs': df_bi_40,
        'Trifunctional 40 XLs': df_trifunctional_40,
        'Bifunctional 60 XLs': df_bi_60,
        'Trifunctional 60 XLs': df_trifunctional_60,
        'Bifunctional 120 XLs': df_bi_120,
        'Trifunctional 120 XLs': df_trifunctional_120
    }

    for name, df in datasets.items():
        mean_value = df[metric].mean()
        mean_results.append({
            'Dataset': name,
            'Mean': mean_value
        })

    # Experimental Data
    exp_types = ['DSSO', 'TSTO']
    for exp_type in exp_types:
        data = df_exp[df_exp['Type'] == exp_type][metric]
        mean_value = data.mean()
        mean_results.append({
            'Dataset': f'Experimental {exp_type}',
            'Mean': mean_value
        })

    # Create DataFrame of means
    df_means = pd.DataFrame(mean_results)

    # Write to file
    means_filename = os.path.join(output_dir, f'means_{metric.lower().replace(" ", "_")}.csv')
    df_means.to_csv(means_filename, index=False)

    # Print means to terminal
    print(f"\nMeans for {metric}:")
    print(df_means.to_string(index=False))

# Create plots for 'Model Accuracy' and 'Cluster Precision'
create_plots('Model Accuracy')
create_plots('Cluster Precision')