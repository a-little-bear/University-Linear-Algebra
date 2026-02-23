# Chapter 14: Matrix Analysis

<div class="context-flow" markdown>

**Prerequisites**: Matrix Functions (Ch13) · Calculus Series Theory · Norm Basics (Ch15)

**Chapter Outline**: From Scalar to Matrix Analysis → Matrix Sequences and Convergence → Convergence of Matrix Series → The Central Role of Spectral Radius $\rho(A)$ → Gelfand's Formula → Criteria for Convergence of Matrix Powers $A^k$ → Gershgorin Circle Theorem → Applications: Iterative Methods for Linear Systems (Jacobi, Gauss-Seidel)

**Extension**: Matrix analysis is the intersection of numerical analysis and continuous dynamical systems; it introduces continuous tools from calculus into the discrete world of matrices and is key to analyzing the stability of large-scale systems.

</div>

Matrix analysis studies the analytical properties of matrices when treated as variables. The core problems lie in defining what it means for matrices to be "close" and how to handle the long-term evolution of matrix sequences. Unlike elementary algebra, matrix analysis focuses on limits, rates of convergence, and the distribution of the spectrum.

---

## 14.1 Matrix Sequences and Convergence

!!! definition "Definition 14.1 (Matrix Sequence Convergence)"
    A matrix sequence $\{A_k\}$ converges to $A$, denoted $\lim_{k \to \infty} A_k = A$, if for every entry of the matrices, $\lim_{k \to \infty} (A_k)_{ij} = a_{ij}$.
    **Equivalence**: This is equivalent to saying that for any matrix norm $\|\cdot\|$, $\lim_{k \to \infty} \|A_k - A\| = 0$.

---

## 14.2 Spectral Radius and Convergence Criteria

!!! definition "Definition 14.2 (Spectral Radius $\rho(A)$)"
    The **spectral radius** of a square matrix $A$ is the maximum of the absolute values of its eigenvalues:
    $$\rho(A) = \max \{|\lambda| : \lambda \in \sigma(A)\}$$

!!! theorem "Theorem 14.1 (Convergence of Powers)"
    The matrix sequence of powers $\lim_{k \to \infty} A^k = O$ if and only if $\rho(A) < 1$.

!!! theorem "Theorem 14.2 (Gelfand's Formula)"
    For any matrix norm $\|\cdot\|$, we have:
    $$\rho(A) = \lim_{k \to \infty} \|A^k\|^{1/k}$$
    This establishes a profound link between the purely algebraic property of an operator (spectral radius) and its analytical property (norm).

---

## 14.3 Gershgorin Circle Theorem

!!! theorem "Theorem 14.3 (Gershgorin Circle Theorem)"
    All eigenvalues of a matrix $A$ lie within the union of the following $n$ disks in the complex plane:
    $$R_i = \left\{ z \in \mathbb{C} : |z - a_{ii}| \le \sum_{j \neq i} |a_{ij}| \right\}$$
    **Application**: This provides a way to quickly estimate the general range of eigenvalues without computing them explicitly.

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

Matrix analysis endows static algebra with the dimension of time:

1.  **Convergence Criteria**: The spectral radius is the ultimate judge of evolution, determining whether discrete dynamical systems stabilize or explode.
2.  **Analytical Tools**: Matrix calculus provides the symbolic engine for optimization algorithms and variational problems, allowing for dynamic analysis of complex matrix models.
3.  **Global Perspective**: Gelfand's formula and the Circle Theorem unify local entry-wise information with global spectral properties, providing efficient paths for numerical stability estimation.
