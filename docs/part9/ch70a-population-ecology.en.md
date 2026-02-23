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

1. **[Survival Rate Constraints] Explain the biological-physical basis for the Leslie matrix sub-diagonal elements $s_i$ satisfying $0 \le s_i \le 1$ and its significance for the matrix norm.**
   ??? success "Solution"
       $s_i$ represents the probability of surviving from group $i$ to $i+1$. By the axioms of probability, its value must be in $[0, 1]$. This ensures the induced norm of $L$ is physiologically limited, reflecting the dissipative nature of material flow across generations.

2. **[Zero Fertility Analysis] If the first row of the Leslie matrix (fertility row) is all zeros, analyze the spectral structure of the operator and the final fate of the population.**
   ??? success "Solution"
       In this case, $L$ is a strictly lower triangular matrix, and all its eigenvalues are 0. The system is nilpotent, $L^n = 0$. Physically, this means that since no new individuals are produced, the existing population will be exhausted within $n$ periods.

3. **[Calculation] Given the Leslie matrix $L = \begin{pmatrix} 0 & 2 \ 0.5 & 0 \end{pmatrix}$. Calculate its eigenvalues and determine whether the population is growing, declining, or in dynamic equilibrium.**
   ??? success "Solution"
       The characteristic equation is $\det(\lambda I - L) = \lambda^2 - 1 = 0 \implies \lambda = \pm 1$. The principal eigenvalue is $\lambda_1 = 1$. The population size remains constant in the long term, i.e., it is at a stable replacement level.

4. **[Steady-state Age Structure] Prove: The steady-state age distribution of a population corresponds to the right eigenvector $\mathbf{w}$ of $L$ associated with $\lambda_1$.**
   ??? success "Solution"
       Since $L\mathbf{w} = \lambda_1 \mathbf{w}$, as $k 	o \infty$, $x_k \propto \lambda_1^k \mathbf{w}$. The ratios of age groups $x_{i,k} / \sum x_{j,k}$ converge to $w_i / \sum w_j$, meaning the age structure no longer changes over time.

5. **[Reproductive Value] Explain the algebraic meaning of the components $v_i$ of the left eigenvector $\mathbf{v}$ in population contribution analysis.**
   ??? success "Solution"
       $\mathbf{v}$ is the reproductive value vector. Since $\mathbf{v}^T L = \lambda_1 \mathbf{v}^T$, the component $v_i$ quantifies the expected contribution weight of an individual in group $i$ to the future total population size.

6. **[Stage-structured Models] Analyze the algebraic differences between the Lefkovitch matrix (non-zero diagonal) and the Leslie matrix in describing developmental arrest.**
   ??? success "Solution"
       The Lefkovitch matrix allows individuals to remain in the same developmental stage, corresponding to diagonal entries $l_{ii} > 0$. This makes the cyclic structure of the matrix more complex and may alter the damping characteristics of convergence to the principal eigenvalue.

7. **[Parameter Sensitivity] Derive the formula for the partial derivative of the eigenvalue $\lambda_1$ with respect to the matrix element $l_{ij}$, and explain its use in ecological decision-making.**
   ??? success "Solution"
       $\frac{\partial \lambda_1}{\partial l_{ij}} = \frac{v_i w_j}{\mathbf{v}^T \mathbf{w}}$. it identifies which life-history parameter (e.g., juvenile survival or adult fertility) contributes most significantly to the total growth rate, thereby guiding the allocation of conservation resources.

8. **[Elasticity Analysis] Prove: For a Leslie matrix, the sum of the elasticities $e_{ij} = \frac{l_{ij}}{\lambda_1} \frac{\partial \lambda_1}{\partial l_{ij}}$ of all elements equals 1.**
   ??? success "Solution"
       This is determined by the homogeneity of eigenvalues (Euler's theorem): $\lambda_1(k L) = k \lambda_1(L)$. A total elasticity of 1 means that changes in the growth rate can be fully decomposed into the percentage contributions of each parameter.

9. **[Primitivity Determination] Analyze the graph-theoretic criteria for a Leslie matrix to satisfy the requirement of being a primitive matrix.**
   ??? success "Solution"
       Its directed graph must be strongly connected and have a period of 1. Usually, as long as there are two coprime age groups with non-zero fertility, this criterion is met, ensuring unique convergence to the steady-state structure.

10. **[Ecological Stability] Explain the relationship between the matrix spectral radius $ho(L)$ and system dissipation.**
    ??? success "Solution"
        $ho(L)$ measures the expansion factor of energy (living individuals) in the system. As a non-negative operator, its principal eigenvalue characterizes the asymptotic growth rate of the discrete dynamical system in the sense of the $\ell_1$ norm.

## Chapter Summary

This chapter systematically discusses the application of linear algebra in structured population dynamics:

1. **Transition Operator Theory**: Established the Leslie matrix as the core linear model for describing population evolution.
2. **Asymptotic Structural Analysis**: Utilized eigenvalues and eigenvectors to formulate the algebraic representation of population growth rates and steady-state age distributions.
3. **Value and Contribution**: Introduced reproductive value (left eigenvector) to quantify the dynamic weight of individuals within the group.
4. **Perturbation Analysis**: Demonstrated the role of sensitivity and elasticity theory in supporting decisions by identifying key physiological parameters of the system.
