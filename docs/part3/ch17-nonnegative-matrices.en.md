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

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

## Chapter Summary

Non-negative matrix theory establishes a harmonic resonance between physical laws and algebraic structures:

1.  **Preservation of Positivity**: Perron-Frobenius theory establishes the positivity of the principal eigenvector, providing unique and valid steady-state solutions for problems in probability and resource allocation.
2.  **Structural Connectivity**: Irreducibility links matrix algebra to graph topology, proving that "global circulation" is a prerequisite for a system to converge to a unique stable state.
3.  **Computational Convergence**: Stochastic matrix theory lays the foundation for all equilibrium evolution processes (e.g., Markov chains), revealing the final destination of information flow in complex networks.
