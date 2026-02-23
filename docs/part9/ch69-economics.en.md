# Chapter 69: Linear Algebra in Economics

<div class="context-flow" markdown>

**Prerequisites**: Non-negative Matrices & Perron-Frobenius (Ch17) · Matrix Equations (Ch20) · Basics of Probability and Statistics

**Chapter Outline**: Algebraic Balance of Social Production → The Leontief Input-Output Model → Consumption Matrix $A$ and the Production Equation $(I-A)x=d$ → Economic Meaning of Matrix Inversion: The Multiplier Effect → Price Evolution and Spectral Analysis → Portfolio Optimization: Covariance Matrices in the Markowitz Model → Linear Representations in Arbitrage Pricing Theory (APT) → Applications: National Industrial Planning, Financial Risk Assessment, and Inflation Prediction

**Extension**: Economics is the "resource mapping" of linear algebra; it abstracts complex transactions and supply chains into flow matrices. It proves that long-term economic equilibrium is essentially the principal eigenvector of a positive operator—the algebraic engine for understanding modern macro-regulation and financial mathematics.

</div>

In modern economics, industries are interdependent. Agriculture needs machinery, machine shops need steel, and steel mills need electricity. **Linear Algebra** provides the perfect framework for describing this complex web of supply. Through the **Leontief Input-Output Model**, we can precisely calculate the total output required from each sector to satisfy a specific set of societal demands. This chapter introduces this Nobel Prize-winning theory of economic algebra.

---

## 69.1 The Leontief Input-Output Model

!!! definition "Definition 69.1 (Consumption Matrix $A$)"
    In an economic system with $n$ sectors, $a_{ij}$ represents the value of sector $i$'s output required to produce 1 unit of sector $j$'s output.
    **Linear System**: Let $\mathbf{x}$ be the vector of total outputs and $\mathbf{d}$ be the vector of external demands. The system satisfies:
    $$(I - A)\mathbf{x} = \mathbf{d}$$

!!! note "The Multiplier Effect"
    If $(I-A)$ is invertible, then $\mathbf{x} = (I-A)^{-1}\mathbf{d}$. Using the power series $(I-A)^{-1} = I + A + A^2 + \cdots$, each term represents successive rounds of indirect demand triggered by the initial external consumption.

---

## 69.2 Portfolio Optimization (Markowitz)

!!! technique "Technique: Minimum Variance Portfolio"
    Given the covariance matrix $\Sigma \succ 0$ of $n$ assets' returns, find the weight vector $\mathbf{w}$ that minimizes risk for a target return. This is a **Quadratic Programming** problem under linear constraints, whose analytic solution involves the calculation of $\Sigma^{-1}$.

---

## 69.3 Price Equilibrium and Spectral Theory

!!! theorem "Theorem 69.1 (Price Equilibrium)"
    In long-term equilibrium, the price vector $\mathbf{p}$ satisfies $\mathbf{p}^T = \mathbf{p}^T A + \mathbf{v}^T$, where $\mathbf{v}$ is the vector of value-added. This essentially involves finding the left eigenvector associated with the production process.

---

## Exercises

**1. [Basics] Given consumption matrix $A = \begin{pmatrix} 0.5 & 0.4 \\ 0.2 & 0.5 \end{pmatrix}$ and external demand $\mathbf{d} = (1, 1)^T$, find the total output $\mathbf{x}$.**

??? success "Solution"
    **Steps:**
    1. Calculate $I - A = \begin{pmatrix} 0.5 & -0.4 \\ -0.2 & 0.5 \end{pmatrix}$.
    2. Inversion: $\det = 0.25 - 0.08 = 0.17$.
    3. $(I-A)^{-1} = \frac{1}{0.17} \begin{pmatrix} 0.5 & 0.4 \\ 0.2 & 0.5 \end{pmatrix}$.
    4. $\mathbf{x} = \frac{1}{0.17} \begin{pmatrix} 0.9 \\ 0.7 \end{pmatrix} \approx \begin{pmatrix} 5.29 \\ 4.12 \end{pmatrix}$.

**2. [Property] Why must the spectral radius $\rho(A)$ be less than 1 in the Leontief model?**

??? success "Solution"
    **Economic Reasoning:**
    1. $\rho(A) < 1$ ensures $(I-A)^{-1}$ exists and is non-negative.
    2. Physically, this means the total input (direct + indirect) required to produce 1 unit of output is less than 1 unit.
    3. If $\rho(A) \ge 1$, the system consumes resources faster than it produces them, making it impossible to satisfy any positive external demand.
    **Conclusion**: This is the algebraic criterion for a "productive" or sustainable economy.

**3. [Calculation] Two assets have correlation 0.5 and standard deviations 10% and 20%. Write the covariance matrix $\Sigma$.**

??? success "Solution"
    **Steps:**
    1. Variances: $\sigma_1^2 = 0.01$, $\sigma_2^2 = 0.04$.
    2. Covariance: $\sigma_{12} = \rho \sigma_1 \sigma_2 = 0.5 \cdot 0.1 \cdot 0.2 = 0.01$.
    **Conclusion**: $\Sigma = \begin{pmatrix} 0.01 & 0.01 \\ 0.01 & 0.04 \end{pmatrix}$.

**4. [Portfolio] Prove that the positive definiteness of the covariance matrix determines a unique risk minimum.**

??? success "Solution"
    **Reasoning:**
    Portfolio variance is $\mathbf{w}^T \Sigma \mathbf{w}$. If $\Sigma \succ 0$, this is a strictly convex hyper-paraboloid. Optimization theory states that a strictly convex function over a convex feasible set has a unique global minimum.

**5. [Application] What is "Forward Linkage" in industrial analysis?**

??? success "Solution"
    **Algebraic manifestation:**
    It corresponds to the row sums of the matrix $(I-A)^{-1}$. It describes how much a particular sector's output is utilized as input for other sectors. Linear algebra reveals the "bottleneck" sectors of a national economy through these row sum calculations.

**6. [Calculation] Is the economic system described by $A = \begin{pmatrix} 0.1 & 0.9 \\ 0.8 & 0.1 \end{pmatrix}$ feasible?**

??? success "Solution"
    **Steps:**
    Characteristic equation: $(0.1-\lambda)^2 - 0.72 = 0 \implies \lambda = 0.1 \pm \sqrt{0.72} \approx 0.1 \pm 0.85$.
    The maximum eigenvalue $\rho(A) \approx 0.95 < 1$.
    **Conclusion**: Feasible. Despite high cross-consumption, the system remains productive.

**7. [Duality] What do "Shadow Prices" in economics correspond to in Linear Programming?**

??? success "Solution"
    **Conclusion: The dual variables $\mathbf{y}$.**
    Shadow prices represent the marginal increase in the objective function (total value) for a unit increase in a constrained resource ($b_i$). This captures the marginal utility of matrix constraints.

**8. [Arbitrage] What is the linear factor model in Arbitrage Pricing Theory (APT)?**

??? success "Solution"
    It assumes returns $R$ can be modeled as $R = E + B\mathbf{f} + \epsilon$, where $B$ is the **loading matrix** (Betas) and $\mathbf{f}$ is the vector of common factors. This essentially decomposes asset fluctuations into a systemic subspace component and an idiosyncratic component.

**9. [Property] Prove: If $A \ge 0$ is a consumption matrix with $\rho(A) < 1$, then $(I-A)^{-1} \ge 0$.**

??? success "Solution"
    **Reasoning:**
    The series $(I-A)^{-1} = \sum_{k=0}^\infty A^k$ converges. Since powers of non-negative matrices are non-negative, and their sum is non-negative, the result follows. This is an application of M-matrix theory (Ch38A), ensuring that increased demand leads to increased production.

**10. [Application] Describe the role of linear algebra in analyzing "Global Value Chains" (GVC).**

??? success "Solution"
    By constructing Multi-Regional Input-Output (MRIO) tables, the matrix entries cross borders. Using Schur complements and matrix inversion on these partitioned matrices, researchers can trace how many times a product crosses borders during production and precisely quantify each nation's value-added contribution.

## Chapter Summary

Linear algebra is the "universal scale" of quantitative economics:

1.  **Production Closure**: The Leontief model proves that complex social production can be described via fixed points of linear operators, establishing the mathematical basis for sectoral coordination.
2.  **Geometrization of Risk**: Portfolio theory transforms uncertain financial fluctuations into optimization paths over convex surfaces using covariance matrices, powering modern quantitative finance.
3.  **Multiplier Effect of Value**: Inversion theory reveals how local economic shocks are amplified through industrial chains, providing the logic engine for government policy evaluation.
