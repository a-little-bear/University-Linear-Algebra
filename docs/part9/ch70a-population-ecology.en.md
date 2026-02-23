# Chapter 70A: Linear Algebra in Demography and Population Ecology

<div class="context-flow" markdown>

**Prerequisites**: Non-negative Matrices & Perron-Frobenius (Ch17) · Matrix Analysis (Ch14) · Difference Equations

**Chapter Outline**: Algebraic Modeling of Population Structure → The Leslie Matrix (Age-Structured Models) → The Lefkovitch Matrix (Stage-Structured Models) → Population Evolution Equation $\mathbf{x}_{k+1} = L \mathbf{x}_k$ → Asymptotic Behavior: Growth Rates and the Perron Eigenvalue → Stable Age Distribution (Right Eigenvector) and Reproductive Value (Left Eigenvector) → Sensitivity and Elasticity Analysis → Applications: Endangered Species Conservation, Fisheries Management, and Human Demography

**Extension**: Linear algebra in ecology transforms the succession of life into matrix power iterations; it proves that the long-term survival of biological systems depends on specific spectral radii, serving as the scientific ruler for resource management and conservation strategy.

</div>

Will a species go extinct or grow infinitely? **Population Ecology Linear Algebra** answers these questions through structured matrices. By dividing a population into different age groups or life stages (e.g., juvenile, adult), we can use a non-negative matrix to describe the transitions and reproduction between these groups. This chapter demonstrates how Perron-Frobenius theory serves as the ultimate logic for predicting biological evolution.

---

## 70A.1 The Leslie Matrix Model

!!! definition "Definition 70A.1 (Leslie Matrix)"
    Suppose a population is divided into $n$ age classes. The **Leslie Matrix** $L$ describes the transition from time $k$ to $k+1$:
    $$L = \begin{pmatrix} f_1 & f_2 & \cdots & f_{n-1} & f_n \\ s_1 & 0 & \cdots & 0 & 0 \\ 0 & s_2 & \cdots & 0 & 0 \\ \vdots & \vdots & \ddots & \vdots & \vdots \\ 0 & 0 & \cdots & s_{n-1} & 0 \end{pmatrix}$$
    - $f_i \ge 0$: Fecundity (reproduction rate) of the $i$-th class.
    - $s_i \in (0, 1]$: Survival probability from the $i$-th class to the $(i+1)$-th class.

---

## 70A.2 Asymptotic Behavior and Eigenvalues

!!! theorem "Theorem 70A.1 (Fundamental Theorem of Populations)"
    For an irreducible and primitive Leslie matrix $L$:
    1.  **Growth Rate**: The long-term growth rate is determined by the Perron eigenvalue $\lambda_1 = \rho(L)$.
        - $\lambda_1 > 1$: The population grows.
        - $\lambda_1 < 1$: The population goes extinct.
    2.  **Stable Age Distribution**: The proportions of the age classes converge to the **right eigenvector** $\mathbf{v}$ associated with $\lambda_1$.
    3.  **Reproductive Value**: The contribution of each class to future population size is characterized by the **left eigenvector** $\mathbf{w}$.

---

## 70A.3 Sensitivity and Elasticity

!!! technique "Sensitivity Analysis"
    The impact of a small change in parameter $a_{ij}$ on the growth rate $\lambda_1$ is:
    $$\frac{\partial \lambda_1}{\partial a_{ij}} = \frac{w_i v_j}{\mathbf{w}^T \mathbf{v}}$$
    This informs us which life stage (e.g., improving juvenile survival vs. adult fecundity) is most effective for maintaining or increasing population size.

---

## Exercises

1.  **[Basics] Write the $2 \times 2$ Leslie matrix for a population with fecundity $(0, 2)$ and survival $0.5$.**
    ??? success "Solution"
        $L = \begin{pmatrix} 0 & 2 \\ 0.5 & 0 \end{pmatrix}$.

2.  **[Growth] Calculate the eigenvalues of the matrix above and determine population stability.**
    ??? success "Solution"
        $\det(L-\lambda I) = \lambda^2 - 1 = 0 \implies \lambda = \pm 1$.
        $\rho(L) = 1$. The population is in dynamic equilibrium (neither growing nor shrinking).

3.  **[Distribution] Find the stable age distribution for $\lambda=1$ in the previous problem.**
    ??? success "Solution"
        Solve $(L-I)v = 0 \implies \begin{pmatrix} -1 & 2 \\ 0.5 & -1 \end{pmatrix} \begin{pmatrix} v_1 \\ v_2 \end{pmatrix} = 0 \implies v_1 = 2v_2$.
        The stable ratio is $2:1$.

4.  **[Irreducibility] If all $f_i = 0$ except $f_n > 0$, is the Leslie matrix irreducible?**
    ??? success "Solution"
        Yes, because there is a cycle $1 \to 2 \to \cdots \to n \to 1$.

5.  **[Reproduction] Why is the reproductive value (left eigenvector) usually lower for the oldest age class?**
    ??? success "Solution"
        Because older classes have fewer remaining years of fecundity and cannot transition into future high-productivity classes; their "investment value" for the future is low.

6.  **[Calculation] If $\lambda_1 = 1.05$, what is the annual percentage growth rate?**
    ??? success "Solution"
        $5\%$.

7.  **[Lefkovitch] What is the main difference between a Lefkovitch matrix and a Leslie matrix?**
    ??? success "Solution"
        A Lefkovitch matrix allows for non-zero diagonal entries (representing individuals remaining in the same life stage, like perennial plants), making it more general.

8.  **[Application] How can matrix models be used to set fishing quotas?**
    ??? success "Solution"
        By determining how much of each age class can be "harvested" without dropping the Perron eigenvalue $\lambda_1$ below 1.

9.  **[Primitivity] Prove: If the first row of $L$ is strictly positive, $L$ is primitive.**
    ??? success "Solution"
        This ensures that every age class reproduces back to the first class, and because there is a 1-cycle ($f_1 > 0$), the lengths of cycles are coprime. By the Perron-Frobenius theorem, it is primitive.

10. **[Oscillation] If $L = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$, will the population converge to a stable distribution?**

   ??? success "Solution"
        No. It is a cyclic matrix (non-primitive). The population counts will oscillate indefinitely between the two classes.

## Chapter Summary

Linear algebra in ecology proves the digital laws governing life:

1.  **Algebraic Essence of Growth**: The spectral radius $\rho(L)$ is the theme of biological evolution, simplifying complex physiological processes into a single survival metric.
2.  **Inevitability of Form**: The right eigenvector proves that regardless of the initial population state, life structures eventually converge to an inherent proportional beauty as long as the environment is stable.
3.  **Precision of Intervention**: Sensitivity analysis transforms the differential tools of linear algebra into a compass for ecological protection, revealing how algebraic optimization strategies can save endangered species when resources are limited.
