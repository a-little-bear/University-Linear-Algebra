# Chapter 17: Non-negative Matrices and Perron-Frobenius Theory

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues (Ch06) · Matrix Analysis (Ch14) · Graph Theory Basics (Ch27)

**Chapter Outline**: Definition of Non-negative Matrices → Positive vs. Non-negative → Irreducibility (Graph Criteria) → Perron's Theorem (Positive Matrices) → Frobenius Extension (Irreducible Matrices) → Uniqueness of Perron Eigenvalues and Vectors → Stochastic Matrices and Markov Chains → Collatz-Wielandt Formula → Application: Google's PageRank Algorithm

**Extension**: Perron-Frobenius theory is the ultimate tool for analyzing steady states in discrete-time systems (e.g., population models, economic input-output, web ranking); it reveals how the physical constraint of "non-negativity" leads to exceptionally strong algebraic laws.

</div>

In the real world, most physical quantities (e.g., probability, population, money, pixels) are non-negative. Non-negative matrices provide the mathematical models for the evolution of these quantities. Perron-Frobenius theory not only guarantees the existence of a dominant eigenvalue but also establishes a "physically meaningful" positive eigenvector associated with it.

---

## 17.1 Positive Matrices and Perron's Theorem

!!! definition "Definition 17.1 (Positive and Non-negative Matrices)"
    1.  **Non-negative Matrix $A \ge 0$**: All entries $a_{ij} \ge 0$.
    2.  **Positive Matrix $A > 0$**: All entries $a_{ij} > 0$.

!!! theorem "Theorem 17.1 (Perron's Theorem)"
    If $A > 0$, then:
    1.  **Dominant Eigenvalue**: The spectral radius $r = \rho(A)$ is an eigenvalue of $A$ and $r > 0$.
    2.  **Simplicity**: $r$ is a simple root of the characteristic equation (algebraic multiplicity 1).
    3.  **Positive Vector**: The eigenvector $\mathbf{v}$ associated with $r$ can be chosen strictly positive ($v_i > 0$).
    4.  **Dominance**: For any other eigenvalue $\lambda$, $|\lambda| < r$.

---

## 17.2 Irreducibility and the Frobenius Extension

!!! definition "Definition 17.2 (Irreducible Matrix)"
    A matrix $A$ is **irreducible** if it cannot be transformed into a block upper triangular form via permutation. In graph theory, this is equivalent to its associated directed graph being strongly connected.

!!! theorem "Theorem 17.2 (Frobenius Theorem)"
    For an irreducible non-negative matrix $A \ge 0$, most of Perron's conclusions hold, but "dominance" is weakened to $|\lambda| \le r$ (there may be multiple eigenvalues with modulus $r$, as seen in cyclic matrices).

---

## 17.3 Stochastic Matrices and Markov Chains

!!! definition "Definition 17.3 (Stochastic Matrix)"
    A non-negative matrix $A$ is **row-stochastic** if the sum of elements in each row is 1.
    **Properties**: A stochastic matrix has $\rho(A) = 1$, and the all-ones vector $\mathbf{1}$ is the eigenvector associated with the eigenvalue 1.

---

## Exercises

1. **[Criteria] Determine if $A = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$ is a positive matrix.**
   ??? success "Solution"
       No. It is a non-negative matrix but contains zero entries.

2. **[Irreducibility] Is the matrix $A$ from the previous exercise irreducible?**
   ??? success "Solution"
       Yes. Its associated graph $1 \leftrightarrow 2$ is strongly connected.

3. **[Eigenvalue] Find the Perron eigenvalue of $\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$.**
   ??? success "Solution"
       The eigenvalues are 2 and 0. Thus $r = \rho(A) = 2$.

4. **[Eigenvector] Find the positive eigenvector for $r=2$ from the matrix above.**
   ??? success "Solution"
       $\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix} \begin{pmatrix} 1 \\ 1 \end{pmatrix} = 2 \begin{pmatrix} 1 \\ 1 \end{pmatrix}$. Thus $\mathbf{v} = (1, 1)^T$.

5. **[Stochastic] Prove a stochastic matrix always has an eigenvalue of 1.**
   ??? success "Solution"
       The $i$-th component of $A \mathbf{1}$ is the sum of the $i$-th row. By definition, this sum is 1, so $A \mathbf{1} = 1 \cdot \mathbf{1}$.

6. **[Graph] Why is $A = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$ reducible?**
   ??? success "Solution"
       It is already upper triangular (or observe the graph: vertex 2 cannot reach vertex 1).

7. **[Frobenius] Give an example of an irreducible non-negative matrix with $|\lambda| = r$ but $\lambda \neq r$.**
   ??? success "Solution"
       $A = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$. Eigenvalues are 1 and -1. $|-1| = 1$.

8. **[PageRank] Why is a small positive perturbation added to the link matrix in PageRank?**
   ??? success "Solution"
       To make the matrix strictly positive, ensuring the uniqueness of the Perron eigenvalue 1 and the convergence of the power method.

9. **[Monotonicity] If $A \ge 0$, prove $\rho(A)$ is a non-decreasing function of the entries of $A$.**
   ??? success "Solution"
       From the Collatz-Wielandt formula $\rho(A) = \sup_{x>0} \min_i \frac{(Ax)_i}{x_i}$, increasing an entry of $A$ clearly increases the supremum.

10. **[Limit] If $A$ is primitive, what is $\lim_{k \to \infty} (A/r)^k$?**
    ??? success "Solution"
        It equals $\mathbf{v}\mathbf{w}^T$, where $\mathbf{v}$ and $\mathbf{w}$ are the normalized left and right Perron eigenvectors.

## Chapter Summary

Non-negative matrix theory establishes a harmonic resonance between physical laws and algebraic structures:

1.  **Preservation of Positivity**: Perron-Frobenius theory establishes the positivity of the principal eigenvector, providing unique and valid steady-state solutions for problems in probability and resource allocation.
2.  **Structural Connectivity**: Irreducibility links matrix algebra to graph topology, proving that "global circulation" is a prerequisite for a system to converge to a unique stable state.
3.  **Computational Convergence**: Stochastic matrix theory lays the foundation for all equilibrium evolution processes (e.g., Markov chains), revealing the final destination of information flow in complex networks.
