

# A Practical Guide to Statistical Analysis with SciPy

*A tutorial for undergraduates on performing statistical tests in Python, created by T W.Terooatea

## Introduction

While data visualization helps us see patterns in our data, statistical analysis allows us to quantify these patterns and determine if they are significant or simply due to random chance. The **SciPy** library, a cornerstone of the scientific Python ecosystem, provides a powerful and comprehensive module called `scipy.stats` for performing a wide range of statistical tests.

This tutorial is designed to be a practical introduction to statistical analysis for undergraduate students. We will move beyond simply calculating means and standard deviations to performing formal hypothesis tests. You will learn how to ask statistical questions, choose the appropriate test, and interpret the results to make data-driven conclusions.

We will use the same cell growth dataset from our plotting tutorial to provide a familiar context. The goal is to empower you to not only create beautiful plots but also to back them up with rigorous statistical evidence.

### Learning Objectives:

- Understand the fundamentals of hypothesis testing (null hypothesis, p-value).
- Learn how to perform and interpret a t-test for comparing two groups.
- Learn how to perform and interpret an ANOVA for comparing multiple groups.
- Apply these statistical tests to a real-world biological dataset to answer scientific questions.

### Official Documentation (Reference):

SciPy is a vast library. For a complete list of functions and more detailed explanations, the official documentation is your best friend.

- **SciPy Reference Guide:** [https://docs.scipy.org/doc/scipy/reference/index.html](https://docs.scipy.org/doc/scipy/reference/index.html)



## Descriptive Statistics: The First Look at Your Data

Before diving into complex hypothesis tests, it's crucial to start with **descriptive statistics**. This involves summarizing the main features of your dataset. It helps you get a feel for the data's distribution, central tendency, and spread. This initial exploration can guide your choice of statistical tests and help you spot any potential issues, such as outliers.

### Key Descriptive Statistics

Here are some of the most common descriptive statistics you will encounter:

| Statistic | Description | How to get it in Python (with pandas) |
| :--- | :--- | :--- |
| **Mean** | The average value. | `.mean()` |
| **Median** | The middle value of a sorted dataset. Less sensitive to outliers than the mean. | `.median()` |
| **Standard Deviation (std)** | A measure of how spread out the data is from the mean. | `.std()` |
| **Variance** | The square of the standard deviation. | `.var()` |
| **Min & Max** | The minimum and maximum values in the dataset. | `.min()`, `.max()` |
| **Quartiles** | Values that divide your data into four equal parts (25th, 50th, 75th percentiles). | `.quantile(0.25)`, `.quantile(0.75)` |

### Using `.describe()` for a Quick Summary

Pandas provides a convenient method called `.describe()` that calculates several of these statistics for all numerical columns in your DataFrame at once.

```python
import pandas as pd

# Load the dataset
df = pd.read_csv('hybridoma_culture_data.csv')

# Get a statistical summary of the key metrics
summary = df[['Cell_density_per_mL', 'Viability_percent', 'Antibody_mg_per_L']].describe()

print(summary)
```

### Grouping Data for Meaningful Insights

While a general summary is useful, you'll often want to calculate statistics for different groups within your data. For our cell growth experiment, we want to see how these statistics differ between treatments. This is where `groupby()` comes in.

```python
# Calculate the mean antibody production for each treatment
mean_antibody_by_treatment = df.groupby('Treatment')['Antibody_mg_per_L'].mean()

print(mean_antibody_by_treatment)

# You can calculate multiple statistics at once
stats_by_treatment = df.groupby('Treatment')['Antibody_mg_per_L'].agg(['mean', 'std', 'median'])

print(stats_by_treatment)
```

> **Challenge:** The bar plots you created in the visualization tutorial showed the *maximum* antibody production. The code above calculates the *mean*. How do these different statistics (mean vs. max) tell different stories about the data? Which one do you think is more representative of a treatment's overall performance?

This initial analysis gives us a quantitative summary and sets the stage for the next logical step: **hypothesis testing**. For example, we can see from the descriptive statistics that the 'Sodium_butyrate' treatment appears to have the highest mean antibody production. But is this difference statistically significant, or could it have occurred by chance? To answer that, we need a formal statistical test.



## Hypothesis Testing: Asking Questions of Your Data

Descriptive statistics give us a summary, but **hypothesis testing** lets us make formal conclusions about our data. It provides a framework for asking specific questions and determining if the answers are statistically significant.

### The Core Concepts

At the heart of hypothesis testing are two competing ideas:

-   **The Null Hypothesis (H₀):** This is the default assumption, usually stating that there is **no effect** or **no difference**. For example, "There is no difference in antibody production between the Control group and the Sodium Butyrate group."
-   **The Alternative Hypothesis (H₁):** This is what you are trying to prove, stating that there **is an effect** or **a difference**. For example, "There is a difference in antibody production between the two groups."

Our goal is to see if we have enough evidence to reject the null hypothesis in favor of the alternative.

### The P-value: Your Evidence Score

The **p-value** is a probability that measures the strength of evidence against the null hypothesis. It tells you the likelihood of observing your data (or something more extreme) if the null hypothesis were actually true.

> **How to Interpret a P-value:**
> -   A **small p-value (typically ≤ 0.05)** suggests that your observed data is unlikely under the null hypothesis. This provides strong evidence to **reject the null hypothesis**.
> -   A **large p-value (> 0.05)** suggests that your data is consistent with the null hypothesis. You **fail to reject the null hypothesis** (which is not the same as proving it is true!).

The threshold of 0.05 is a common convention in many scientific fields and is known as the **significance level (alpha)**.

### Comparing Two Groups: The Independent T-test

A t-test is used to determine if there is a significant difference between the means of **two** independent groups. Let's use it to formally test if the Sodium Butyrate treatment significantly improved antibody production compared to the Control.

#### Guidance for Performing a T-test

1.  **Prepare your data:** You need to isolate the data for the two groups you want to compare.
2.  **Run the test:** Use the `scipy.stats.ttest_ind()` function.
3.  **Interpret the result:** Look at the p-value.

Here is a template to guide you. We will use the maximum antibody production for each replicate, which we calculated in the plotting tutorial.

```python
from scipy import stats
import pandas as pd

# Load data and find max antibody production per replicate
df = pd.read_csv("hybridoma_culture_data.csv")
max_antibody = df.loc[df.groupby(["Treatment", "Replicate"])["Antibody_mg_per_L"].idxmax()]

# Step 1: Isolate the data for the two groups
control_data = max_antibody[max_antibody["Treatment"] == "Control"]["Antibody_mg_per_L"]
butyrate_data = max_antibody[max_antibody["Treatment"] == "Sodium_butyrate"]["Antibody_mg_per_L"]

# Step 2: Perform the independent t-test
t_statistic, p_value = stats.ttest_ind(control_data, butyrate_data)

# Step 3: Interpret and print the results
print(f"T-test between Control and Sodium Butyrate:")
print(f"T-statistic: {t_statistic:.2f}")
print(f"P-value: {p_value:.4f}")

# Check for statistical significance
alpha = 0.05
if p_value < alpha:
    print("The difference is statistically significant (reject H₀).")
else:
    print("The difference is not statistically significant (fail to reject H₀).")
```

> **Challenge:** The `ttest_ind` function has a parameter called `equal_var`. By default, it's `True`, which assumes the two groups have the same variance. Run the test again with `equal_var=False` (this is known as Welch's t-test). Does the result change? Welch's t-test is often a safer choice when you are unsure if the variances are equal.

### Comparing Multiple Groups: Analysis of Variance (ANOVA)

What if you want to compare all seven treatments at once? Running t-tests between every possible pair would be tedious and increase the chance of a false positive. This is where **ANOVA** comes in.

ANOVA tells you if there is a statistically significant difference somewhere **among the means of three or more groups**. It tests the following hypotheses:

-   **H₀:** The means of all treatment groups are equal.
-   **H₁:** At least one treatment group mean is different from the others.

#### Guidance for Performing a One-Way ANOVA

1.  **Prepare your data:** You need to create a list where each element is an array of data for one group.
2.  **Run the test:** Use the `scipy.stats.f_oneway()` function.
3.  **Interpret the result:** A small p-value tells you that at least one group is different, but it doesn't tell you *which* one.

```python
# ... (setup code from before) ...

# Step 1: Prepare the data for all groups
treatment_groups = []
for treatment_name in max_antibody["Treatment"].unique():
    group_data = max_antibody[max_antibody["Treatment"] == treatment_name]["Antibody_mg_per_L"]
    treatment_groups.append(group_data)

# Step 2: Perform the one-way ANOVA
f_statistic, p_value = stats.f_oneway(*treatment_groups)

# Step 3: Interpret and print the results
print(f"\nOne-Way ANOVA across all treatments:")
print(f"F-statistic: {f_statistic:.2f}")
print(f"P-value: {p_value:.4f}")

if p_value < alpha:
    print("There is a significant difference among the treatment groups (reject H₀).")
else:
    print("There is no significant difference among the treatment groups (fail to reject H₀).")
```

### After ANOVA: Post-Hoc Tests

If your ANOVA test is significant (p < 0.05), you know there's a difference to be found, but where? To find out which specific pairs of groups are different from each other (e.g., is Butyrate different from Glucose? Is Glucose different from Control?), you need to perform a **post-hoc test**. A common and reliable choice is **Tukey's Honestly Significant Difference (HSD) test**, which is available in the `statsmodels` library. This is an important next step for a complete analysis.



## Correlation and Regression: Exploring Relationships

So far, we have focused on comparing groups. But what if we want to understand the relationship *between* two continuous variables? For example, as cells consume glucose, do they produce more lactate? Is there a relationship between cell density and antibody production? Correlation and regression are the tools for these questions.

### Correlation: Measuring the Strength of a Relationship

**Correlation** measures the strength and direction of a linear relationship between two variables. The result is a **correlation coefficient (r)**, which ranges from -1 to +1:

-   **r = +1:** Perfect positive linear relationship (as one variable increases, the other increases).
-   **r = -1:** Perfect negative linear relationship (as one variable increases, the other decreases).
-   **r = 0:** No linear relationship.

SciPy's `stats.pearsonr()` function calculates the Pearson correlation coefficient and also gives a p-value to test for non-correlation.

#### Guidance for Calculating Correlation

Let's test the relationship between glucose consumption and lactate production. We would expect a negative correlation (as glucose goes down, lactate goes up).

```python
from scipy import stats
import pandas as pd

# Load data
df = pd.read_csv("hybridoma_culture_data.csv")

# Isolate the two variables of interest
glucose = df["Glucose_g_per_L"]
lactate = df["Lactate_g_per_L"]

# Calculate Pearson correlation
corr_coefficient, p_value = stats.pearsonr(glucose, lactate)

# Print and interpret the results
print(f"Correlation between Glucose and Lactate:")
print(f"Pearson Correlation Coefficient (r): {corr_coefficient:.2f}")
print(f"P-value: {p_value:.4f}")

if p_value < 0.05:
    print("The correlation is statistically significant.")
else:
    print("The correlation is not statistically significant.")
```

> **Important:** Correlation does not imply causation! Just because two variables are correlated doesn't mean one causes the other. There could be a third, unmeasured variable influencing both.

### Linear Regression: Modeling the Relationship

While correlation tells you the strength of a relationship, **linear regression** goes a step further by finding the best-fit line that describes that relationship. This allows you to make predictions.

The equation for a simple linear regression is: **y = mx + c**

-   **y:** The dependent variable (what you are trying to predict).
-   **x:** The independent variable (what you are using to predict).
-   **m:** The **slope** of the line (how much y changes for a one-unit change in x).
-   **c:** The **intercept** (the value of y when x is 0).

SciPy's `stats.linregress()` function is perfect for this. It returns the slope, intercept, r-value, p-value, and standard error.

#### Guidance for Performing Linear Regression

Let's model the relationship between cell density and antibody production.

```python
# ... (setup code) ...

# Perform linear regression
# Let's see if cell density can predict antibody production
slope, intercept, r_value, p_value, std_err = stats.linregress(
    df["Cell_density_per_mL"], df["Antibody_mg_per_L"]
)

# Print and interpret the results
print(f"\nLinear Regression: Antibody Production vs. Cell Density")
print(f"Slope (m): {slope:.6f}")
print(f"Intercept (c): {intercept:.2f}")
print(f"R-squared: {r_value**2:.2f}") # R-squared is a measure of how well the model fits the data
print(f"P-value: {p_value:.4f}")

print(f"\nModel: Antibody_mg_per_L = {slope:.6f} * Cell_density_per_mL + {intercept:.2f}")
```

> **Challenge:** A good way to visualize a regression is to create a scatter plot of the two variables and then overlay the regression line. Use `seaborn.regplot()` to do this in one easy step. Try it for `Cell_density_per_mL` vs. `Antibody_mg_per_L`.



## Beyond the Basics: Non-Parametric Tests and Distribution Fitting

Many common statistical tests, like the t-test and ANOVA, are **parametric tests**. This means they make certain assumptions about the data, most notably that the data is normally distributed. But what if your data doesn't follow a normal distribution? This is common with biological data. In such cases, you should use **non-parametric tests**.

Non-parametric tests do not assume a particular distribution for the data. They often work by ranking the data instead of using the raw values, which makes them more robust to outliers and skewed distributions.

### When to Use Non-Parametric Tests

-   When your data is not normally distributed.
-   When you have a small sample size.
-   When your data contains significant outliers.

### The Mann-Whitney U Test: The Non-Parametric T-test

The **Mann-Whitney U test** is the non-parametric equivalent of the independent t-test. It is used to determine if there is a significant difference between two independent groups, but it tests for a difference in the *medians* rather than the means.

#### Guidance for Performing a Mann-Whitney U Test

The process is very similar to the t-test, but you use `stats.mannwhitneyu()`.

```python
# ... (setup code) ...

# Isolate the data for the two groups
control_data = max_antibody[max_antibody["Treatment"] == "Control"]["Antibody_mg_per_L"]
butyrate_data = max_antibody[max_antibody["Treatment"] == "Sodium_butyrate"]["Antibody_mg_per_L"]

# Perform the Mann-Whitney U test
u_statistic, p_value = stats.mannwhitneyu(control_data, butyrate_data)

# Print and interpret the results
print(f"\nMann-Whitney U test between Control and Sodium Butyrate:")
print(f"U-statistic: {u_statistic:.2f}")
print(f"P-value: {p_value:.4f}")
```

### The Kruskal-Wallis H-Test: The Non-Parametric ANOVA

The **Kruskal-Wallis H-test** is the non-parametric equivalent of the one-way ANOVA. It is used to determine if there are statistically significant differences between the medians of **three or more** independent groups.

#### Guidance for Performing a Kruskal-Wallis H-Test

Just like with ANOVA, you prepare a list of your groups and then pass them to the `stats.kruskal()` function.

```python
# ... (setup code) ...

# Prepare the data for all groups
treatment_groups = []
for treatment_name in max_antibody["Treatment"].unique():
    group_data = max_antibody[max_antibody["Treatment"] == treatment_name]["Antibody_mg_per_L"]
    treatment_groups.append(group_data)

# Perform the Kruskal-Wallis H-test
h_statistic, p_value = stats.kruskal(*treatment_groups)

# Print and interpret the results
print(f"\nKruskal-Wallis H-test across all treatments:")
print(f"H-statistic: {h_statistic:.2f}")
print(f"P-value: {p_value:.4f}")
```

> **Challenge:** How can you check if your data is normally distributed? The `scipy.stats` module has a function called `shapiro()` that performs the Shapiro-Wilk test for normality. A p-value less than 0.05 from this test suggests your data is *not* normally distributed. Try running it on the `control_data` and `butyrate_data`.

### Fitting Distributions

Another powerful feature of `scipy.stats` is the ability to fit a theoretical distribution (like a normal, log-normal, or gamma distribution) to your data. This can be useful for understanding the underlying process that generated your data or for creating models.

The `stats.norm.fit()` function, for example, will estimate the mean and standard deviation of the normal distribution that best fits your data.

```python
# Fit a normal distribution to the control group's antibody data
mu, std = stats.norm.fit(control_data)

print(f"\nFitted Normal Distribution for Control Group:")
print(f"Estimated Mean (mu): {mu:.2f}")
print(f"Estimated Standard Deviation (std): {std:.2f}")
```



## Practical Application: A Complete Statistical Workflow

Let's tie everything together by performing a complete statistical analysis on our cell growth data. Our goal is to determine which treatment is the most effective at increasing antibody production and to present our findings with statistical rigor.

### The Research Question

Does any of the supplementary treatments significantly increase the maximum antibody production compared to the control group?

### The Workflow

1.  **Explore the Data:** Start with descriptive statistics and visualization.
2.  **Global Comparison:** Use ANOVA (or its non-parametric equivalent) to see if there are *any* significant differences among the treatments.
3.  **Pairwise Comparisons:** If the global test is significant, use t-tests (or their non-parametric equivalent) to compare each treatment directly to the control.
4.  **Summarize and Report:** Present your findings in a clear and concise way.

### Step-by-Step Implementation

Here is the code to perform this full workflow. Try to follow along and understand the purpose of each step.

```python
import pandas as pd
from scipy import stats

# --- 1. Data Preparation ---
df = pd.read_csv("hybridoma_culture_data.csv")
# Find the maximum antibody production for each replicate
max_antibody = df.loc[df.groupby(["Treatment", "Replicate"])["Antibody_mg_per_L"].idxmax()]

# --- 2. Descriptive Statistics ---
print("--- Descriptive Statistics for Max Antibody Production ---")
descriptive_stats = max_antibody.groupby('Treatment')['Antibody_mg_per_L'].agg(['mean', 'std', 'median'])
print(descriptive_stats.sort_values(by='mean', ascending=False))

# --- 3. Global Comparison (ANOVA) ---
print("\n--- Global Comparison: One-Way ANOVA ---")
treatment_groups = [group["Antibody_mg_per_L"].values for name, group in max_antibody.groupby("Treatment")]
f_statistic, p_value_anova = stats.f_oneway(*treatment_groups)

print(f"ANOVA P-value: {p_value_anova:.4f}")
if p_value_anova < 0.05:
    print("Conclusion: There is a significant difference among the treatments.")
else:
    print("Conclusion: No significant difference was found among the treatments.")

# --- 4. Pairwise Comparisons to Control (T-tests) ---
if p_value_anova < 0.05:
    print("\n--- Pairwise Comparisons with Control ---")
    control_data = max_antibody[max_antibody["Treatment"] == "Control"]["Antibody_mg_per_L"]
    
    # Get a list of all other treatments
    supplement_treatments = [t for t in max_antibody["Treatment"].unique() if t != "Control"]
    
    results = []
    for treatment in supplement_treatments:
        treatment_data = max_antibody[max_antibody["Treatment"] == treatment]["Antibody_mg_per_L"]
        
        # Perform Welch's t-test (often safer)
        t_stat, p_val = stats.ttest_ind(treatment_data, control_data, equal_var=False)
        
        results.append({
            "Treatment": treatment,
            "P-value": p_val,
            "Is_Significant": "Yes" if p_val < 0.05 else "No"
        })

    # Create a summary DataFrame
    results_df = pd.DataFrame(results)
    print(results_df.sort_values(by="P-value"))

```

### Interpreting the Final Results

The output of this script gives you a complete story:

1.  The **descriptive statistics** give you a ranked list of which treatments performed best on average.
2.  The **ANOVA result** confirms that the observed differences are not just random noise; there is a real effect to be found.
3.  The **table of p-values** from the t-tests tells you exactly which treatments were significantly better than the control.

> **Challenge:** The pairwise comparison we did increases the risk of a false positive (a Type I error). To correct for this, you can use a **multiple comparisons correction**, such as the Bonferroni correction. A simple way to apply it is to divide your significance level (alpha) by the number of tests you performed. Modify the code to use a corrected alpha. Does it change which treatments you consider significant?



## Conclusion and Best Practices

This tutorial has guided you through the fundamental tools for statistical analysis in Python using the `scipy.stats` library. You have learned how to move from simple data summaries to formal hypothesis testing, enabling you to draw statistically-supported conclusions from your data. By understanding and applying these techniques, you can add a new layer of rigor to your scientific investigations.

### Best Practices for Statistical Analysis

1.  **Visualize Your Data First:** Always plot your data before running statistical tests. Visualizations can reveal the underlying distribution, potential outliers, and trends that summary statistics alone might miss. A box plot or a violin plot is an excellent companion to a t-test or ANOVA.

2.  **Know Your Test's Assumptions:** Parametric tests like t-tests and ANOVA assume that your data is normally distributed and that the groups have equal variances. Use tests like the Shapiro-Wilk test (`stats.shapiro()`) to check for normality and consider using non-parametric alternatives (like Mann-Whitney U or Kruskal-Wallis) if these assumptions are violated.

3.  **Choose the Right Test for Your Question:** Be clear about what you are asking. Are you comparing two groups? Use a t-test. Are you comparing more than two groups? Use ANOVA. Are you looking for a relationship between two continuous variables? Use correlation or regression.

4.  **Don't Just Chase P-values:** A p-value tells you whether an effect is statistically significant, but it doesn't tell you the *size* of the effect. Always report descriptive statistics (like the mean difference) alongside your p-values to provide context.

5.  **Correct for Multiple Comparisons:** When you perform many statistical tests at once (e.g., comparing multiple treatments to a control), your chance of getting a false positive increases. Be aware of this and consider using corrections like the Bonferroni correction to maintain statistical rigor.

6.  **Report Your Findings Clearly:** When reporting your results, always include the name of the test you used, the test statistic (e.g., t-statistic or F-statistic), the degrees of freedom, and the exact p-value. This transparency allows others to understand and evaluate your analysis.

By combining powerful visualization with robust statistical analysis, you can build a compelling and evidence-based narrative from your data. Happy analyzing!

