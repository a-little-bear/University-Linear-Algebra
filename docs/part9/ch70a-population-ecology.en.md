# Chapter 70A: Matrix Models in Population Ecology

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues (Ch6) · Non-negative Matrices (Ch17) · Difference Equations

**Chapter Outline**: Leslie Matrix Model → Structured Populations → Age Distribution → Survival and Fertility Rates → Long-term Growth Rate (Perron Root) → Steady-state Age Distribution → Reproductive Value → Lefkovitch Matrix → Sensitivity and Elasticity Analysis

**Extension**: Matrix models are the quantitative decision-making foundation for modern wildlife conservation, fisheries management, and pest control.

</div>

Biological systems, though complex, have basic generational evolution laws that can be accurately characterized using linear algebra. The Leslie matrix model partitions a population by age and reveals the "long-term rhythm" of life through its spectral structure.

---

## 70A.1 Leslie Matrix Model and Evolutionary Asymptotics

!!! definition "Definition 70A.1 (Leslie Transition Operator)"
    Suppose a population is divided into $n$ discrete age groups. The Leslie matrix $L$ encodes the fertility rate $f_i$ and survival rate $s_i$ of each group. The discrete dynamical equation for the population vector $x_k$ is $x_{k+1} = L x_k$. The special structure of this matrix ensures the existence of its principal eigenvalue.

!!! theorem "Theorem 70A.1 (Principal Eigenvalue and Long-term Growth)"
    If the Leslie matrix $L$ is primitive, then according to the Perron-Frobenius theorem, there exists a unique principal eigenvalue $\lambda_1 > 0$. The asymptotic growth rate of the total population $N_k$ is $\lambda_1$. If $\lambda_1 > 1$, the population grows geometrically; if $\lambda_1 < 1$, the population faces asymptotic extinction.

---

## Exercises

1. **[Survival Rate Constraints] Explain why the sub-diagonal elements $s_i$ of a Leslie matrix must satisfy $0 \le s_i \le 1$.**
   ??? success "Solution"
       $s_i$ represents the probability of an individual surviving from group $i$ to group $i+1$. By the axioms of probability, a probability cannot be negative or exceed 1. This ensures the operator norm of $L$ is consistent with physiological limits.

2. **[Zero Fertility] Analyze the final fate of a population if the first row of its Leslie matrix is all zeros.**
   ??? success "Solution"
       In this case, $L$ is a strictly lower triangular matrix, making it nilpotent ($L^n = 0$). All eigenvalues are 0. Physically, this means that with no new offspring being produced, the existing population will age and die out within at most $n$ time periods.

3. **[Growth Rate] Given $L = \begin{pmatrix} 0 & 2 \\ 0.5 & 0 \end{pmatrix}$. Calculate the eigenvalues and determine if the population is stable.**
   ??? success "Solution"
       Characteristic equation: $\det(\lambda I - L) = \lambda^2 - 1 = 0 \implies \lambda = \pm 1$. The principal eigenvalue is $\lambda_1 = 1$. The population size remains constant in the long term (replacement level).

4. **[Steady-state Structure] Prove: The steady-state age distribution corresponds to the right eigenvector $\mathbf{w}$ of $L$ associated with $\lambda_1$.**
   ??? success "Solution"
       As $k \to \infty$, the state vector $x_k = L^k x_0$ is dominated by the term $c_1 \lambda_1^k \mathbf{w}$. The relative proportions of the age groups $x_{i,k} / \sum x_{j,k}$ converge to $w_i / \sum w_j$, which is independent of the initial population distribution.

5. **[Reproductive Value] Interpret the components $v_i$ of the left eigenvector $\mathbf{v}$ associated with $\lambda_1$.**
   ??? success "Solution"
       The left eigenvector $\mathbf{v}^T L = \lambda_1 \mathbf{v}^T$ is known as the reproductive value vector. Component $v_i$ quantifies the relative contribution of a single individual in age group $i$ to the future total population size compared to other groups.

6. **[Stage-structured] Contrast the Lefkovitch matrix with the Leslie matrix in terms of development.**
   ??? success "Solution"
       The Leslie matrix assumes all survivors move to the next age class. The Lefkovitch matrix allows individuals to remain in the same developmental stage (e.g., adult), corresponding to non-zero diagonal entries $l_{ii} > 0$. This provides flexibility for modeling species like trees or insects.

7. **[Sensitivity] Provide the formula for the sensitivity of $\lambda_1$ to changes in $l_{ij}$.**
   ??? success "Solution"
       The sensitivity is given by $s_{ij} = \frac{\partial \lambda_1}{\partial l_{ij}} = \frac{v_i w_j}{\mathbf{v}^T \mathbf{w}}$. This derivative identifies which vital rate (like juvenile survival) has the most impact on the population's overall growth rate.

8. **[Elasticity] Show that the sum of elasticities $e_{ij} = \frac{l_{ij}}{\lambda_1} s_{ij}$ for a Leslie matrix equals 1.**
   ??? success "Solution"
       Since $\lambda_1$ is a homogeneous function of degree 1 in the matrix entries, Euler's homogeneous function theorem implies $\lambda_1 = \sum l_{ij} \frac{\partial \lambda_1}{\partial l_{ij}}$. Dividing by $\lambda_1$ yields $\sum e_{ij} = 1$.

9. **[Primitivity] State the graph-theoretic condition for a Leslie matrix to be primitive.**
   ??? success "Solution"
       The directed graph associated with the matrix must be strongly connected and its cycle lengths must be coprime. This ensures that the population structure eventually settles into a unique stationary distribution.

10. **[Spectral Radius] Relate $\rho(L)$ to the dissipativity of the ecological system.**
    ??? success "Solution"
        $\rho(L)$ measures the per-period expansion factor of "life energy" in the system. As a non-negative operator, its principal eigenvalue determines the asymptotic growth rate in the $\ell_1$ norm (total population count).

## Chapter Summary

This chapter explores the linear modeling of population evolution:

1. **Transition Operators**: Established the Leslie matrix as the standard linear model for structured population growth.
2. **Asymptotic Analysis**: Utilized eigenvalues and eigenvectors to define growth rates and steady-state structures.
3. **Weighting Metrics**: Introduced reproductive value to quantify individual contribution to group dynamics.
4. **Perturbation Analysis**: Developed sensitivity and elasticity tools to support resource allocation in wildlife management.
