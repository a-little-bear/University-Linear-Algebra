# Chapter 70A: Population Ecology

<div class="context-flow" markdown>

**Prerequisites**: Non-negative Matrices & Perron-Frobenius (Ch17) · Matrix Analysis (Ch14) · Difference Equations

**Chapter Outline**: Linear Representation of Biological Systems → Construction of the Leslie Matrix: Survival and Fecundity → Discrete Population Evolution $\mathbf{x}_{t+1} = L \mathbf{x}_t$ → Core Criterion: The Perron Eigenvalue (Dominant Growth Rate) → Stable Age Distribution → Spectral Analysis of Population Dynamics → Decay Coefficients and Extinction Criteria → Applications: Conservation Planning for Endangered Species, Setting Fishery Harvest Quotas, and Invasive Species Prediction

**Extension**: Population ecology is the "pulse of life" in linear algebra; it abstracts the biological process of reproduction and survival into iterated matrix powers. It proves that the final structure of a population is pre-determined by its physiological parameters (eigenvalues and eigenvectors)—the algebraic base for mathematical biology.

</div>

In nature, different age groups contribute differently to population growth. Juveniles do not reproduce but may have high survival, while adults possess high fecundity. **Population Ecology** utilizes the **Leslie Matrix** to integrate these physiological parameters into a unified linear model. By studying the spectral properties of this matrix, we can not only predict the future total population but also understand how the proportion of age groups converges to a "perfect equilibrium" over time. This chapter introduces this algebraic tool central to ecological resource management.

---

## 70A.1 Construction of the Leslie Matrix

!!! definition "Definition 70A.1 (Leslie Matrix $L$)"
    For a population with $n$ age groups, $L$ has the following specific structure:
    $$L = \begin{pmatrix} f_1 & f_2 & \cdots & f_{n-1} & f_n \\ s_1 & 0 & \cdots & 0 & 0 \\ 0 & s_2 & \cdots & 0 & 0 \\ \vdots & \vdots & \ddots & \vdots & \vdots \\ 0 & 0 & \cdots & s_{n-1} & 0 \end{pmatrix}$$
    - $f_i \ge 0$: The fecundity (birth rate) of the $i$-th age group (first row).
    - $s_i \in (0, 1]$: The probability of surviving from age group $i$ to $i+1$ (sub-diagonal).

---

## 70A.2 Dominant Eigenvalue and Stable Distribution

!!! theorem "Theorem 70A.1 (Application of Perron-Frobenius)"
    Since the Leslie matrix is non-negative and usually irreducible, it possesses a unique positive dominant eigenvalue $\lambda_1$.
    1.  **Growth Criterion**: If $\lambda_1 > 1$, the population grows; if $\lambda_1 < 1$, it decays.
    2.  **Stable Distribution**: The corresponding positive eigenvector $\mathbf{v}_1$ describes the **Stable Age Distribution** attained in the long term.

---

## 70A.3 Long-term Projections

!!! technique "Technique: Power Method Projection"
    For any initial population $\mathbf{x}_0$, the dominance of $\lambda_1$ implies:
    $$\mathbf{x}_t \approx c \lambda_1^t \mathbf{v}_1 \quad \text{as } t \to \infty$$
    This means that regardless of the initial state, the population will eventually grow at rate $\lambda_1$ with an age structure fixed by $\mathbf{v}_1$.

---

## Exercises

**1. [Basics] Given a $2 \times 2$ Leslie matrix $L = \begin{pmatrix} 0 & 2 \\ 0.5 & 0 \end{pmatrix}$, describe its biological meaning.**

??? success "Solution"
    **Analysis:**
    1. First row: Juveniles are sterile ($f_1=0$), while each adult produces 2 offspring ($f_2=2$) per period.
    2. Sub-diagonal: The survival rate from juvenile to adult is 0.5 ($s_1=0.5$).
    3. Adults do not survive to a next period (end of lifecycle).

**2. [Calculation] Find the growth rate $\lambda$ for the population in the previous problem.**

??? success "Solution"
    **Steps:**
    1. Characteristic equation: $\det(L - \lambda I) = \begin{vmatrix} -\lambda & 2 \\ 0.5 & -\lambda \end{vmatrix} = \lambda^2 - 1 = 0$.
    2. Solve for $\lambda$: $\lambda = \pm 1$.
    **Conclusion**: The dominant eigenvalue is 1. The population size remains **constant** in the long term.

**3. [Stable Distribution] Find the eigenvector (stable age proportion) for $\lambda=1$.**

??? success "Solution"
    **Solve $(L-I)\mathbf{v} = 0$:**
    1. $\begin{pmatrix} -1 & 2 \\ 0.5 & -1 \end{pmatrix} \begin{pmatrix} v_1 \\ v_2 \end{pmatrix} = 0 \implies v_1 = 2v_2$.
    2. Set $v_2 = 1, v_1 = 2$.
    **Conclusion**: The stable age distribution is $2:1$. In the long run, there are twice as many juveniles as adults.

**4. [Survival Analysis] If environmental degradation drops $s_1$ from 0.5 to 0.4, will the population go extinct?**

??? success "Solution"
    **Determination:**
    1. New matrix $L' = \begin{pmatrix} 0 & 2 \\ 0.4 & 0 \end{pmatrix}$.
    2. Eigenvalues: $\lambda^2 - 0.8 = 0 \implies \lambda = \sqrt{0.8} \approx 0.89$.
    3. Since $\rho(L') < 1$.
    **Conclusion**: Yes, the population will shrink by about 11% per generation and eventually vanish.

**5. [Properties] Why is the first row of a Leslie matrix not necessarily all positive?**

??? success "Solution"
    **Reasoning:**
    Most organisms have a juvenile stage where they are not yet capable of reproduction. Thus, the first few entries $f_1, f_2, \ldots$ are typically zero. As long as at least one age group eventually reproduces, the system can still have a positive Perron eigenvalue.

**6. [Calculation] Show that for a Leslie matrix, the characteristic equation can be written as $\sum_{i=1}^n f_i \left( \prod_{j=1}^{i-1} s_j \right) \lambda^{-i} = 1$.**

??? success "Solution"
    **Proof Sketch:**
    Use Laplace expansion or observe the state transitions: each term represents the probability of an individual surviving to age $i$ and producing offspring, discounted by the growth rate. This is the matrix version of the **Euler-Lotka equation**.

**7. [Application] Briefly state how "Fishery Harvest Quotas" are set using matrix eigenvalues.**

??? success "Solution"
    **Strategy:**
    1. Assume the natural growth rate is $\lambda > 1$.
    2. Harvesting acts as a reduction factor on the survival rates: $H = \operatorname{diag}(1-h_1, 1-h_2, \ldots)$.
    3. The goal is to find a harvest vector $\mathbf{h}$ such that the principal eigenvalue of the adjusted matrix $L_{new}$ is exactly 1.
    This yields the maximum sustainable yield while keeping the population at equilibrium.

**8. [Primitivity] If a Leslie matrix is primitive, will the population oscillate?**

??? success "Solution"
    **Conclusion: No.**
    **Reasoning**: According to primitivity theory (Ch17), the dominant eigenvalue $\lambda_1$ is strictly greater than the magnitude of any other eigenvalue. Any initial fluctuations are exponentially dampened, and the population smoothly converges to the stable age distribution.

**9. [Calculation] If $\mathbf{x}_0 = (100, 0)^T$, use the matrix from Problem 1 to find $\mathbf{x}_1$ and $\mathbf{x}_2$.**

??? success "Solution"
    **Iteration:**
    1. $\mathbf{x}_1 = L \mathbf{x}_0 = \begin{pmatrix} 0 & 2 \\ 0.5 & 0 \end{pmatrix} \begin{pmatrix} 100 \\ 0 \end{pmatrix} = \begin{pmatrix} 0 \\ 50 \end{pmatrix}$.
    2. $\mathbf{x}_2 = L \mathbf{x}_1 = \begin{pmatrix} 0 & 2 \\ 0.5 & 0 \end{pmatrix} \begin{pmatrix} 0 \\ 50 \end{pmatrix} = \begin{pmatrix} 100 \\ 0 \end{pmatrix}$.
    **Note**: In this specific case, $\lambda_2 = -1$, causing the population to oscillate between generations (it is not primitive).

**10. [Application] What is a "Lefkovitch Matrix" and how does it differ from a Leslie Matrix?**

??? success "Solution"
    **Comparison:**
    - **Leslie Matrix**: Based strictly on **age**. Individuals must move to the next group each year.
    - **Lefkovitch Matrix**: Based on **stage** (e.g., larva, pupa, adult). Individuals can stay in the same stage or skip stages.
    Algebraically, a Lefkovitch matrix can have non-zero diagonals, providing a more flexible linear model for populations with variable development rates.

## Chapter Summary

Linear algebra is the "digital DNA" parsing the logic of biological reproduction:

1.  **Spectral Analysis of Life**: Via the dominant eigenvalue $\lambda_1$, it distills complex life-and-death data into a single growth metric, establishing the mathematical criterion for species survival.
2.  **Inevitability of Equilibrium**: Eigenvector theory proves that a stable age structure is the necessary destination under environmental and physiological constraints, showcasing the natural beauty of systems converging to attractors.
3.  **Science of Intervention**: From fishery quotas to forest management, matrix models provide the only analytic path for quantitative natural resource management—the rational guide for harmony between humanity and ecosystems.
