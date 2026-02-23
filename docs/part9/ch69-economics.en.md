# Chapter 69: Linear Algebra in Economics

<div class="context-flow" markdown>

**Prerequisites**: Non-negative Matrices & Perron-Frobenius (Ch17) · M-matrices (Ch38A) · Optimization (Ch25) · Quadratic Forms (Ch09)

**Chapter Outline**: Leontief Input-Output Model → Consumption Matrix $A$ & Production Equation $x = Ax + d$ → Profitability Conditions & M-matrix Properties → Economic Equilibrium & Eigenvalues (Meaning of the Perron Vector) → Game Theory: Zero-Sum Games & the Minimax Theorem → Role of Linear Programming in Resource Allocation → Financial Engineering: Markowitz Portfolio Optimization & the Efficient Frontier → Risk Metrics: Spectral Decomposition of Covariance Matrices → Applications: National Industry Linkage Analysis, Financial Risk Hedging, and Pricing Models

**Extension**: Economics is linear algebra at a societal scale; it abstracts complex transactions among millions of households into large-scale systems of simultaneous equations. Through matrix inversion and spectral analysis, it reveals the fragility of production chains and the laws of capital market volatility.

</div>

How do various industries in a nation depend on one another? If the price of steel rises, how will it ripple through the entire economic system? **Economic Linear Algebra** answers these questions through matrix models. From the Nobel-winning Leontief model to modern Wall Street asset allocation algorithms, linear algebra provides the unified language for describing resource flow, risk assessment, and game equilibrium.

---

## 69.1 Leontief Input-Output Model

!!! definition "Definition 69.1 (Input-Output Equation)"
    Assume an economic system with $n$ sectors. Let $\mathbf{x}$ be the total output vector, $\mathbf{d}$ the final external demand vector, and $A$ the **consumption matrix** ($a_{ij}$ represents the units of product $i$ needed to produce one unit of product $j$). The balance equation is:
    $$\mathbf{x} = A\mathbf{x} + \mathbf{d} \quad \implies \quad (I - A)\mathbf{x} = \mathbf{d}$$

!!! theorem "Theorem 69.1 (Condition for Productivity)"
    A system can yield a positive output $\mathbf{x} > 0$ for any non-negative demand $\mathbf{d} \ge 0$ iff:
    1.  The spectral radius of the consumption matrix satisfies $\rho(A) < 1$.
    2.  $(I - A)$ is a **non-singular M-matrix** (meaning $(I-A)^{-1} \ge 0$).

---

## 69.2 Game Theory and Linear Programming

!!! technique "Zero-Sum Games and Minimax"
    In a two-person zero-sum game with payoff matrix $M$, Player 1 seeks a mixed strategy $p$ and Player 2 seeks $q$.
    **Minimax Theorem**: $\max_p \min_q p^T M q = \min_q \max_p p^T M q$.
    Solving for this equilibrium value is exactly equivalent to a pair of dual **linear programming** problems.

---

## 69.3 Financial Engineering: Markowitz Portfolio

!!! definition "Definition 69.2 (Mean-Variance Model)"
    Let the returns of $n$ assets be random variables with covariance matrix $\Sigma \succ 0$. Investment weights are $\mathbf{w}$ (with $\sum w_i = 1$).
    - **Risk** (Variance): $\sigma^2 = \mathbf{w}^T \Sigma \mathbf{w}$.
    - **Goal**: Minimize risk given a target expected return $\mu^T \mathbf{w} = E$.
    This is a **quadratic programming** problem with linear constraints, whose solutions form the "Efficient Frontier."

---

## Exercises


****
??? success "Solution"
     $I-A = \begin{pmatrix} 0.5 & -0.2 \\ -0.2 & 0.5 \end{pmatrix}$.
     $\det(I-A) = 0.25 - 0.04 = 0.21$.
     $(I-A)^{-1} = \frac{1}{0.21} \begin{pmatrix} 0.5 & 0.2 \\ 0.2 & 0.5 \end{pmatrix}$.
     $x = (I-A)^{-1} d = \frac{1}{0.21} (7, 7)^T \approx (33.3, 33.3)^T$.


****
??? success "Solution"
     If $\rho(A) \ge 1$, the internal consumption required to produce one unit of product exceeds the output itself. The system would shrink continuously and fail to satisfy any external demand.


****
??? success "Solution"
     Since $x = (I-A)^{-1} d$ and $(I-A)^{-1} \ge 0$, a change $\Delta d \ge 0$ results in $\Delta x = (I-A)^{-1} \Delta d \ge 0$.


****
??? success "Solution"
     Due to symmetry, the optimal strategy for both players is the mixed strategy $(0.5, 0.5)$.


****
??? success "Solution"
     Because portfolio risk (variance) $w^T \Sigma w$ must physically be positive unless all assets are perfectly linearly dependent (in which case it is semi-definite).


****
??? success "Solution"
     It typically represents the "Market Factor," the direction of systematic risk common to all assets.


****
??? success "Solution"
     They represent the marginal value of a scarce resource—the maximum increase in profit from adding one unit of that resource.


****
??? success "Solution"
     The value is 0 (a fair game).


****
??? success "Solution"
     In a closed model, $Ax = x$, meaning output is exactly consumed by the system. This corresponds to the eigenvector associated with the eigenvalue 1 (the Perron vector).

****
??? success "Solution"
    ## Chapter Summary

Linear algebra in economics establishes the structural logic of social systems:


****: The Leontief model proves that mutual industry consumption forms the underlying operator of the system, whose spectral radius determines if the economy can generate surplus value.

****: The Minimax theorem and linear programming reveal the dual essence of optimal decision-making in competitive environments, proving the algebraic unity of individual interest and social resource allocation.

****: The mean-variance model quantifies financial risk as the curvature of ellipsoids in high-dimensional space, proving that through matrix diagonalization (diversification), one can effectively find the optimal trade-off between return and risk.
