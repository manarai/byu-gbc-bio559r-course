

# A Guide to Data Visualization with Matplotlib and Seaborn

*A tutorial for undergraduates on creating and customizing plots in Python,created by T W. Terooatea.*

## Introduction

Data visualization is a critical skill in any scientific discipline. It allows us to transform raw data into insightful stories, making complex information understandable at a glance. A well-crafted plot can reveal patterns, trends, and outliers that might be missed in a table of numbers. In Python, **Matplotlib** and **Seaborn** are two of the most powerful and widely-used libraries for creating high-quality visualizations.

This tutorial is designed to provide you with the fundamental skills to start plotting with these libraries. We will not just show you how to create plots, but also how to take control of their appearance—from colors and labels to fonts and layout. The goal is to empower you to create publication-quality figures that effectively communicate your findings.

We will work through practical examples based on a real-world cell growth experiment, guiding you step-by-step so you can learn to write the code yourself.

### Learning Objectives:

- Understand the basic components of a Matplotlib plot.
- Learn how to create common statistical plots like line plots and bar plots using Seaborn.
- Gain proficiency in customizing plot aesthetics, including colors, labels, fonts, and sizes.
- Apply these skills to reproduce and enhance plots from a sample biological dataset.

### Official Documentation (References):

For more in-depth information, the official documentation is the best resource. We will refer to it throughout this tutorial.

- **Matplotlib:** [https://matplotlib.org/stable/users/index](https://matplotlib.org/stable/users/index)
- **Seaborn:** [https://seaborn.pydata.org/](https://seaborn.pydata.org/)



## Matplotlib Fundamentals: The Anatomy of a Plot

Matplotlib is the foundational plotting library in Python. While it can be complex, understanding its core components will give you the power to create almost any static visualization you can imagine. The key is to think of a plot as a collection of objects that you can manipulate individually.

### The Figure and Axes Objects

Every plot in Matplotlib resides within a **Figure** object. Think of the `Figure` as the canvas or the window that contains everything. Inside this `Figure`, you create one or more **Axes** objects. An `Axes` object is the actual plot—the area where your data is visualized with x and y axes, tick marks, and labels.

This is the most common and flexible way to use Matplotlib. By creating `Figure` and `Axes` objects explicitly, you gain fine-grained control over your plots.

Here is the standard workflow:

1.  **Create a Figure and Axes:** Use `plt.subplots()` to create a figure and one or more axes.
2.  **Plot your data:** Call plotting methods directly on the `Axes` object (e.g., `ax.plot()`, `ax.scatter()`).
3.  **Customize:** Use `Axes` methods to set titles, labels, limits, and more (e.g., `ax.set_title()`, `ax.set_xlabel()`).
4.  **Show the plot:** Use `plt.show()` to display the final figure.

### Your First Plot: A Simple Line Graph

Let's start with a basic example. Imagine you have some time-series data. The code below provides a template for creating a simple line plot. Try to understand what each line does.

```python
import matplotlib.pyplot as plt

# Sample data
x_data = [0, 1, 2, 3, 4]
y_data = [0, 2, 4, 6, 8]

# Step 1: Create a Figure and an Axes
fig, ax = plt.subplots()

# Step 2: Plot the data on the Axes
ax.plot(x_data, y_data)

# Step 3: Customize the plot
ax.set_title("A Simple Line Plot")
ax.set_xlabel("X-Axis Label")
ax.set_ylabel("Y-Axis Label")

# Step 4: Show the plot
plt.show()
```

> **Challenge:** Can you modify the code above to plot a different dataset? Try plotting `y = x^2`.

### Customizing the Basics

Matplotlib gives you control over every element of your plot. Here are a few common customizations you can apply to the `ax.plot()` method:

| Parameter | Description | Example |
| :--- | :--- | :--- |
| `color` | Sets the color of the line. | `color='red'` or `color='#FF5733'` |
| `linestyle` | Sets the style of the line. | `linestyle='--'` (dashed), `linestyle=':'` (dotted) |
| `linewidth` | Sets the thickness of the line. | `linewidth=2.5` |
| `marker` | Adds markers to your data points. | `marker='o'` (circles), `marker='x'` (x's) |

Here is how you would apply these customizations:

```python
# Inside your plotting code...
ax.plot(x_data, y_data, color='green', linestyle='--', linewidth=3, marker='o')
```

This object-oriented approach (`fig`, `ax`) is the foundation for building complex and customized visualizations. In the next section, we will see how Seaborn builds on top of Matplotlib to make statistical plotting even easier.



## Seaborn Fundamentals: Statistical Visualization Made Easy

Seaborn is a data visualization library based on Matplotlib. It provides a high-level interface for drawing attractive and informative statistical graphics. While Matplotlib is excellent for fine-grained control, Seaborn excels at making complex statistical plots simple to create, especially when working with pandas DataFrames.

### Why Use Seaborn?

-   **Integration with Pandas:** Seaborn is designed to work seamlessly with DataFrames.
-   **Statistical Functions:** It automatically handles statistical estimation (like confidence intervals) when plotting.
-   **Attractive Defaults:** Seaborn comes with several built-in themes and color palettes to make your plots look great out of the box.
-   **High-Level Functions:** Functions like `lineplot()`, `barplot()`, and `scatterplot()` make it easy to create common plot types from your data.

### The Seaborn and Matplotlib Connection

Seaborn doesn't replace Matplotlib; it complements it. Most Seaborn plotting functions return a Matplotlib `Axes` object. This means you can use Seaborn for the heavy lifting of statistical plotting and then use your Matplotlib knowledge to customize the final details.

Here is the typical workflow:

1.  **Prepare your data:** Organize your data in a pandas DataFrame.
2.  **Create a Figure and Axes:** Use `plt.subplots()` just like before.
3.  **Plot with Seaborn:** Call a Seaborn function (e.g., `sns.lineplot()`) and pass your DataFrame, column names, and the `ax` object to it.
4.  **Customize with Matplotlib:** Use the returned `ax` object to fine-tune titles, labels, and other elements.

### Your First Seaborn Plot: A Statistical Line Plot

Let's see how Seaborn simplifies plotting from a DataFrame. Notice the `hue` parameter, which automatically creates different lines for each category in a column and adds a legend.

```python
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

# Sample data in a pandas DataFrame, which is common in data analysis
data = {
    'time': [0, 1, 2, 3, 4, 0, 1, 2, 3, 4],
    'value': [5, 7, 8, 12, 15, 3, 4, 6, 5, 7],
    'treatment': ['A', 'A', 'A', 'A', 'A', 'B', 'B', 'B', 'B', 'B']
}
df = pd.DataFrame(data)

# Step 1: Create a Figure and an Axes
fig, ax = plt.subplots()

# Step 2: Plot with Seaborn
sns.lineplot(data=df, x='time', y='value', hue='treatment', ax=ax, marker='o')

# Step 3: Customize with Matplotlib
ax.set_title("Effect of Treatment Over Time")
ax.set_xlabel("Time (hours)")
ax.set_ylabel("Measured Value")

# Step 4: Show the plot
plt.show()
```

> **Challenge:** The `sns.lineplot()` function automatically calculates and displays a confidence interval around the lines (if you have replicate data). Can you find the parameter in the [Seaborn documentation](https://seaborn.pydata.org/generated/seaborn.lineplot.html) to turn this off? Hint: Look for a parameter related to 'error'.



## Creating Multi-Panel Figures: Subplots for Complex Data

One of the most powerful features of Matplotlib is the ability to create multiple plots within a single figure. This is especially useful when you want to compare different aspects of your data side-by-side or create a comprehensive dashboard-style visualization.

### Understanding Subplots

When you call `plt.subplots()`, you can specify the number of rows and columns to create a grid of plots. Each subplot is an independent `Axes` object that you can customize separately.

```python
# Create a 2x2 grid of subplots
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 10))

# Access individual subplots using array indexing
# axes[0, 0] is top-left, axes[0, 1] is top-right
# axes[1, 0] is bottom-left, axes[1, 1] is bottom-right
```

### Task: Creating a 4-Panel Dashboard

Let's recreate the multi-panel figure from the notebook that shows all four key metrics in one comprehensive view. This approach allows viewers to see all the important trends at once.

```python
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Load your data
df = pd.read_csv('hybridoma_culture_data.csv')

# Create a 2x2 subplot grid
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(15, 12))

# Plot 1: Cell Density (top-left)
sns.lineplot(data=df, x='Time_hours', y='Cell_density_per_mL', 
             hue='Treatment', ax=axes[0, 0], marker='o')
axes[0, 0].set_title('Cell Density Over Time', fontsize=14, fontweight='bold')
axes[0, 0].set_xlabel('Time (hours)', fontsize=12)
axes[0, 0].set_ylabel('Cell Density (cells/mL)', fontsize=12)
axes[0, 0].ticklabel_format(style='sci', axis='y', scilimits=(0,0))

# Plot 2: Viability (top-right)
sns.lineplot(data=df, x='Time_hours', y='Viability_percent', 
             hue='Treatment', ax=axes[0, 1], marker='o')
axes[0, 1].set_title('Cell Viability Over Time', fontsize=14, fontweight='bold')
axes[0, 1].set_xlabel('Time (hours)', fontsize=12)
axes[0, 1].set_ylabel('Viability (%)', fontsize=12)

# Plot 3: Antibody Production (bottom-left)
sns.lineplot(data=df, x='Time_hours', y='Antibody_mg_per_L', 
             hue='Treatment', ax=axes[1, 0], marker='o')
axes[1, 0].set_title('Antibody Production Over Time', fontsize=14, fontweight='bold')
axes[1, 0].set_xlabel('Time (hours)', fontsize=12)
axes[1, 0].set_ylabel('Antibody Concentration (mg/L)', fontsize=12)

# Plot 4: Glucose Consumption (bottom-right)
sns.lineplot(data=df, x='Time_hours', y='Glucose_g_per_L', 
             hue='Treatment', ax=axes[1, 1], marker='o')
axes[1, 1].set_title('Glucose Consumption Over Time', fontsize=14, fontweight='bold')
axes[1, 1].set_xlabel('Time (hours)', fontsize=12)
axes[1, 1].set_ylabel('Glucose Concentration (g/L)', fontsize=12)

# Adjust layout to prevent overlapping
plt.tight_layout()
plt.show()
```

> **Challenge:** Notice how each subplot has its own legend, which can be cluttered. Can you find a way to remove individual legends and create a single shared legend for the entire figure? Hint: Look up `fig.legend()` and `ax.get_legend().remove()`.

## Handling Replicate Data: Error Bars and Statistical Visualization

Real experimental data often includes multiple replicates (repeated measurements). Seaborn automatically handles this by calculating confidence intervals, but sometimes you want more control over how variability is displayed.

### Understanding Replicates in Your Data

In the cell growth dataset, each treatment has 3 replicates (Replicate 1, 2, and 3). When you use `sns.lineplot()` or `sns.barplot()`, Seaborn automatically:
1. Calculates the mean across replicates for each time point
2. Displays confidence intervals (the shaded area around lines or error bars on bars)
3. Shows individual data points if you specify certain parameters

### Customizing Error Representation

You can control how Seaborn handles replicates using several parameters:

| Parameter | Description | Example |
|:----------|:------------|:--------|
| `errorbar` | Type of error representation | `errorbar='sd'` (standard deviation), `errorbar='se'` (standard error) |
| `err_style` | Style of error bars | `err_style='bars'` (error bars), `err_style='band'` (shaded area) |
| `estimator` | Function to calculate central tendency | `estimator='mean'`, `estimator='median'` |

```python
# Example: Bar plot with standard deviation error bars
fig, ax = plt.subplots(figsize=(12, 8))

# Calculate maximum values for each replicate
max_antibody = df.loc[df.groupby(['Treatment', 'Replicate'])['Antibody_mg_per_L'].idxmax()]

sns.barplot(data=max_antibody, x='Treatment', y='Antibody_mg_per_L', 
            errorbar='sd', err_style='bars', ax=ax, palette='viridis')

ax.set_title('Maximum Antibody Production by Treatment\n(Error bars show standard deviation)', 
             fontsize=16, fontweight='bold')
ax.set_xlabel('Treatment', fontsize=14)
ax.set_ylabel('Maximum Antibody Concentration (mg/L)', fontsize=14)

# Rotate x-axis labels for better readability
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()
```

## Advanced Customization: Making Your Plots Your Own

Creating a plot is just the first step. To make it truly effective, you need to customize its appearance to highlight the key message and match your desired style. This section will guide you through some of the most common customizations.

### Controlling Colors with Palettes

Color is one of the most powerful tools in data visualization. Seaborn makes it easy to apply beautiful and effective color schemes.

-   **Seaborn Palettes:** Seaborn has a wide range of built-in color palettes. You can set a palette for a specific plot using the `palette` parameter.
-   **Manual Colors:** You can also pass a list of color names or hex codes to the `palette` parameter for full control.

```python
# Using a built-in Seaborn palette
sns.barplot(data=df, x='treatment', y='value', palette='viridis', ax=ax)

# Using a custom list of colors
my_colors = ['#FF5733', '#33FF57', '#3357FF']
sns.barplot(data=df, x='treatment', y='value', palette=my_colors, ax=ax)
```

> **Challenge:** Explore the [Seaborn color palette documentation](https://seaborn.pydata.org/tutorial/color_palettes.html) and try applying a "cubehelix" palette to the line plot example from the previous section.

### Customizing Fonts and Labels

Clear and readable labels are essential. Matplotlib gives you precise control over the font properties of your titles, axis labels, and tick labels.

You can set these properties using a dictionary or by passing individual arguments to the labeling functions.

| Function | Description |
| :--- | :--- |
| `ax.set_title()` | Sets the title for the Axes. |
| `ax.set_xlabel()` | Sets the label for the x-axis. |
| `ax.set_ylabel()` | Sets the label for the y-axis. |
| `ax.tick_params()` | Customizes the appearance of ticks, tick labels, and gridlines. |

Here’s how to adjust font size and weight:

```python
# ... after creating your plot with sns.lineplot() ...

# Customize title and labels with specific font properties
ax.set_title("Cell Growth Under Different Treatments", fontsize=16, fontweight='bold')
ax.set_xlabel("Time (days)", fontsize=12)
ax.set_ylabel("Cell Density (cells/mL)", fontsize=12)

# Customize the tick labels
ax.tick_params(axis='x', labelsize=10, rotation=45) # Rotate x-axis tick labels
ax.tick_params(axis='y', labelsize=10)
```

### Fine-Tuning Legends

When you use the `hue` parameter in Seaborn, a legend is automatically created. You can access and modify it through the `Axes` object.

```python
# ... after plotting ...

# Get the legend object
legend = ax.get_legend()

# Set a new title for the legend
legend.set_title("Treatment Group")

# You can also remove the legend entirely
# ax.get_legend().remove()
```

### Mastering X-Axis Control: Rotation, Ordering, and Formatting

The x-axis often contains categorical data (like treatment names) that can be difficult to read or may need specific ordering. Here's how to take full control of your x-axis appearance.

#### Rotating X-Axis Labels

Long treatment names or many categories can cause overlapping labels. Rotation is your solution:

```python
# Method 1: Using plt.xticks() (works on the current axes)
plt.xticks(rotation=45, ha='right')  # ha='right' aligns text to the right

# Method 2: Using ax.tick_params() (more explicit)
ax.tick_params(axis='x', rotation=45, labelsize=10)

# Method 3: Using ax.set_xticklabels() with rotation
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
```

#### Controlling the Order of Categories

By default, Seaborn plots categorical data in the order it appears in your DataFrame. You can control this ordering in several ways:

```python
# Method 1: Using the 'order' parameter in Seaborn functions
treatment_order = ['Control', 'Glucose_supplement', 'Serum_supplement', 
                   'IGF1_supplement', 'Glutamine_supplement', 'Ascorbic_acid', 'Sodium_butyrate']

sns.barplot(data=df, x='Treatment', y='Antibody_mg_per_L', 
            order=treatment_order, ax=ax)

# Method 2: Sort your DataFrame before plotting
df_sorted = df.sort_values('Antibody_mg_per_L', ascending=False)
sns.barplot(data=df_sorted, x='Treatment', y='Antibody_mg_per_L', ax=ax)

# Method 3: Create a custom categorical order using pandas
df['Treatment'] = pd.Categorical(df['Treatment'], categories=treatment_order, ordered=True)
```

#### Practical Example: Ordered Bar Plot with Rotated Labels

Let's create a bar plot that combines ordering by performance with readable rotated labels:

```python
# Calculate mean antibody production for ordering
mean_antibody = df.groupby('Treatment')['Antibody_mg_per_L'].mean().reset_index()
treatment_order = mean_antibody.sort_values('Antibody_mg_per_L', ascending=False)['Treatment'].tolist()

# Create the plot
fig, ax = plt.subplots(figsize=(12, 8))

sns.barplot(data=df, x='Treatment', y='Antibody_mg_per_L', 
            order=treatment_order, ax=ax, palette='plasma')

# Customize the plot
ax.set_title('Antibody Production by Treatment\n(Ordered by Performance)', 
             fontsize=16, fontweight='bold')
ax.set_xlabel('Treatment', fontsize=14)
ax.set_ylabel('Antibody Concentration (mg/L)', fontsize=14)

# Rotate and format x-axis labels
plt.xticks(rotation=45, ha='right', fontsize=11)

# Add some styling
ax.grid(axis='y', alpha=0.3)  # Add subtle horizontal gridlines
plt.tight_layout()
plt.show()
```

> **Challenge:** Try creating the same plot but with treatments ordered alphabetically instead of by performance. How would you modify the `treatment_order` list?

#### Advanced X-Axis Formatting

Sometimes you need even more control over your x-axis labels:

```python
# Replace underscores with spaces and capitalize
treatment_labels = [label.replace('_', ' ').title() for label in treatment_order]
ax.set_xticklabels(treatment_labels, rotation=45, ha='right')

# Or create a mapping dictionary for custom names
label_mapping = {
    'Control': 'Control',
    'Glucose_supplement': 'Glucose+',
    'Glutamine_supplement': 'Glutamine+',
    'Serum_supplement': 'Serum+',
    'IGF1_supplement': 'IGF-1',
    'Ascorbic_acid': 'Vitamin C',
    'Sodium_butyrate': 'Na-Butyrate'
}

# Apply the mapping
custom_labels = [label_mapping.get(label, label) for label in treatment_order]
ax.set_xticklabels(custom_labels, rotation=45, ha='right')
```

### Adjusting the Overall Style

-   **Seaborn Themes:** Seaborn has five preset themes to quickly change the overall look of your plots: `darkgrid`, `whitegrid`, `dark`, `white`, and `ticks`. You can set this at the beginning of your script with `sns.set_theme()`.

    ```python
    sns.set_theme(style="whitegrid")
    ```

-   **Fixing Layout Issues:** Sometimes, plot elements like titles and axis labels can overlap. Matplotlib's `plt.tight_layout()` function automatically adjusts plot parameters to give a tight layout. It's a good practice to call this right before `plt.show()`.

    ```python
    # ... after all your plotting and customization ...
    plt.tight_layout() # Adjusts plot to prevent labels from overlapping
    plt.show()
    ```



## Practical Application: Visualizing Cell Growth Data

Now, let's apply what we've learned to the cell growth experiment from the provided notebooks. This is where theory meets practice. We will guide you through the process of recreating the plots, focusing on the "how" and "why" of each step.

First, make sure you have the `hybridoma_culture_data.csv` file in the same directory as your script or notebook. Let's load the data and see what it looks like.

```python
import pandas as pd

# Load the dataset
df = pd.read_csv('hybridoma_culture_data.csv')

# Display the first few rows to understand its structure
print(df.head())
```

### Task 1: Creating the Line Plots

The first notebook (`01_ODE_cell_growth_lineplots.ipynb`) visualizes how different metrics (Cell Density, Viability, and Antibody Concentration) change over time for each treatment. This is a perfect use case for Seaborn's `lineplot`.

#### Guidance for Plotting Cell Density

Your goal is to create a line plot showing `Time_hours` on the x-axis and `Cell_density_per_mL` on the y-axis, with different colored lines for each `Treatment`.

1.  **Setup:** Start by importing the necessary libraries (`matplotlib.pyplot`, `seaborn`, `pandas`) and creating your `Figure` and `Axes` objects with `plt.subplots()`.
2.  **Plotting:** Use `sns.lineplot()`. Think about which columns from the DataFrame should be assigned to the `x`, `y`, and `hue` parameters.
3.  **Customization:** The plot in the notebook has a title and axis labels. Use `ax.set_title()`, `ax.set_xlabel()`, and `ax.set_ylabel()` to add these. Remember to make them descriptive!

Here is a template to get you started. Try to fill in the blanks.

```python
# Import libraries
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Load data
df = pd.read_csv('hybridoma_culture_data.csv')

# Create figure and axes
fig, ax = plt.subplots(figsize=(10, 6)) # figsize makes the plot bigger

# Use Seaborn to create the line plot
sns.lineplot(
    data=df, 
    x=..., # What goes on the x-axis?
    y=..., # What goes on the y-axis?
    hue=..., # How do we separate the lines?
    ax=ax,
    marker='o' # Adds markers to the data points
)

# Add a title and labels
ax.set_title('...', fontsize=16, fontweight='bold')
ax.set_xlabel('...', fontsize=12)
ax.set_ylabel('...', fontsize=12)

# Improve the layout and show the plot
plt.tight_layout()
plt.show()
```

> **Challenge:** The y-axis label for cell density might be better in scientific notation. The notebook uses `ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))`. Try adding this line after your customizations to see what it does.

> **Further Practice:** Repeat the process above to create separate plots for `Viability_percent` and `Antibody_mg_per_L`. How would you change the `y` parameter and the labels for each plot?

### Task 2: Creating the Bar Plots with Replicate Variation

The second notebook (`01_ODE_cell_growth_barplot.ipynb`) compares the maximum values achieved for each treatment while showing the variation across replicates. This demonstrates how to handle experimental replicates and create publication-quality bar plots.

#### Guidance for Plotting Maximum Antibody Production with Error Bars

Your goal is to create a bar plot showing the maximum `Antibody_mg_per_L` for each `Treatment`, including error bars that represent the variation across the three replicates.

1.  **Data Preparation:** First, you need to find the maximum antibody concentration for each replicate of each treatment. This preserves the replicate structure for statistical analysis.

    ```python
    # Find maximum antibody concentration for each treatment-replicate combination
    max_antibody_per_replicate = df.loc[df.groupby(['Treatment', 'Replicate'])['Antibody_mg_per_L'].idxmax()]
    
    # Take a look at the structure
    print(max_antibody_per_replicate[['Treatment', 'Replicate', 'Antibody_mg_per_L']].head())
    ```

2.  **Ordering by Performance:** Calculate the mean across replicates to determine the ranking order.

    ```python
    # Calculate mean maximum antibody production for ordering
    mean_max_antibody = max_antibody_per_replicate.groupby('Treatment')['Antibody_mg_per_L'].mean().reset_index()
    treatment_order = mean_max_antibody.sort_values('Antibody_mg_per_L', ascending=False)['Treatment'].tolist()
    ```

3.  **Creating the Bar Plot:** Use `sns.barplot()` with the replicate data. Seaborn will automatically calculate means and error bars.

Here is a comprehensive template:

```python
# ... (imports and data loading) ...

# Step 1: Find maximum values for each replicate
max_antibody_per_replicate = df.loc[df.groupby(['Treatment', 'Replicate'])['Antibody_mg_per_L'].idxmax()]

# Step 2: Calculate ordering
mean_max_antibody = max_antibody_per_replicate.groupby('Treatment')['Antibody_mg_per_L'].mean().reset_index()
treatment_order = mean_max_antibody.sort_values('Antibody_mg_per_L', ascending=False)['Treatment'].tolist()

# Step 3: Create the plot
fig, ax = plt.subplots(figsize=(14, 8))

sns.barplot(
    data=max_antibody_per_replicate,
    x='Treatment', 
    y='Antibody_mg_per_L',
    order=treatment_order,  # Use our calculated order
    errorbar='sd',  # Show standard deviation as error bars
    err_style='bars',  # Use error bars instead of bands
    ax=ax,
    palette='viridis'  # Professional color palette
)

# Step 4: Customize the appearance
ax.set_title('Maximum Antibody Production by Treatment\n(Error bars show standard deviation across replicates)', 
             fontsize=16, fontweight='bold')
ax.set_xlabel('Treatment', fontsize=14)
ax.set_ylabel('Maximum Antibody Concentration (mg/L)', fontsize=14)

# Step 5: Format x-axis labels
plt.xticks(rotation=45, ha='right', fontsize=11)

# Step 6: Add professional styling
ax.grid(axis='y', alpha=0.3)  # Subtle horizontal gridlines
ax.spines['top'].set_visible(False)  # Remove top border
ax.spines['right'].set_visible(False)  # Remove right border

plt.tight_layout()
plt.show()
```

> **Challenge:** Try modifying the code to show individual data points on top of the bars. Look up the `sns.stripplot()` function and see if you can overlay it on your bar plot using the same axes.

#### Advanced Task: Creating a Multi-Metric Comparison

Now let's create a more sophisticated visualization that compares multiple metrics side by side:

```python
# Prepare data for multiple metrics
metrics = ['Antibody_mg_per_L', 'Viability_percent', 'Cell_density_per_mL']
metric_labels = ['Max Antibody (mg/L)', 'Max Viability (%)', 'Max Cell Density (cells/mL)']

# Create subplots for each metric
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

for i, (metric, label) in enumerate(zip(metrics, metric_labels)):
    # Find max values for this metric
    max_values = df.loc[df.groupby(['Treatment', 'Replicate'])[metric].idxmax()]
    
    # Calculate order for this metric
    mean_values = max_values.groupby('Treatment')[metric].mean().reset_index()
    order = mean_values.sort_values(metric, ascending=False)['Treatment'].tolist()
    
    # Create the bar plot
    sns.barplot(data=max_values, x='Treatment', y=metric, 
                order=order, errorbar='sd', ax=axes[i], palette='Set2')
    
    # Customize each subplot
    axes[i].set_title(f'{label}', fontsize=14, fontweight='bold')
    axes[i].set_xlabel('Treatment', fontsize=12)
    axes[i].set_ylabel(label, fontsize=12)
    axes[i].tick_params(axis='x', rotation=45, labelsize=10)
    
    # Format y-axis for cell density
    if metric == 'Cell_density_per_mL':
        axes[i].ticklabel_format(style='sci', axis='y', scilimits=(0,0))

plt.tight_layout()
plt.show()
```

> **Challenge:** Can you modify the x-axis labels to use cleaner names (e.g., "Glucose+" instead of "Glucose_supplement")? Use the label mapping technique from the X-axis control section.



## Conclusion and Best Practices

Congratulations! You have now walked through the fundamentals of data visualization in Python using Matplotlib and Seaborn. You have learned how to create plots from scratch, customize their appearance, and apply these skills to a real-world dataset. The key is to see plotting not as a rigid procedure, but as a creative process of building and refining a visual story.

-   **Matplotlib** is your tool for ultimate control. Understanding the `Figure` and `Axes` objects is the key to unlocking its power.
-   **Seaborn** is your specialist for statistical plotting. It simplifies the creation of complex plots from DataFrames, while still allowing for Matplotlib-level customization.

As you continue your journey in data science, remember that effective visualization is as important as the analysis itself. The more you practice, the more intuitive this process will become.

### Best Practices for Scientific Visualization

1.  **Clarity is Key:** Always label your axes and provide a descriptive title. Your audience should be able to understand the plot without reading your paper.
2.  **Choose the Right Plot:** Use line plots for time-series data, bar plots for comparing categories, scatter plots for relationships between variables, and so on.
3.  **Use Color Purposefully:** Use color to distinguish categories, not just for decoration. Be mindful of colorblind-friendly palettes (Seaborn has many).
4.  **Keep it Clean:** Avoid clutter. Don't add unnecessary gridlines, backgrounds, or effects that distract from the data.
5.  **Ensure Readability:** Make sure your font sizes for titles, labels, and ticks are large enough to be read easily, especially if the plot will be used in a presentation.
6.  **Tell a Story:** Your plot should have a clear message. Use annotations, ordering, and color to guide the viewer's eye to the most important findings.
7.  **Cite Your Tools:** When publishing, it's good practice to mention the tools you used (e.g., "Visualizations were created using Matplotlib and Seaborn in Python.").

Happy plotting!

