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

1. **[Convergence] Let $A = \begin{pmatrix} 0.5 & 1 \\ 0 & 0.5 \end{pmatrix}$. Does the sequence $A^k$ converge to the zero matrix?**
   ??? success "Solution"
       Yes. The eigenvalues are 0.5, so the spectral radius $\rho(A) = 0.5 < 1$.

2. **[Series] Given $A = \begin{pmatrix} 0.1 & 0 \\ 0 & 0.2 \end{pmatrix}$, compute the series $\sum_{k=0}^\infty A^k$.**
   ??? success "Solution"
       Since $\rho(A) < 1$, the series converges to $(I-A)^{-1} = \begin{pmatrix} 1/0.9 & 0 \\ 0 & 1/0.8 \end{pmatrix} = \operatorname{diag}(1.11, 1.25)$.

3. **[Spectral Radius] Provide an example of a matrix $A$ such that $\|A\| > 1$ for some norm but $\rho(A) < 1$.**
   ??? success "Solution"
       For instance, $A = \begin{pmatrix} 0 & 10 \\ 0 & 0 \end{pmatrix}$. Its eigenvalues are 0, so $\rho(A)=0 < 1$. However, its spectral norm (2-norm) is 10.

4. **[Circle Theorem] Estimate the range of eigenvalues for $A = \begin{pmatrix} 4 & 1 \\ 1 & 3 \end{pmatrix}$ using Gershgorin's theorem.**
   ??? success "Solution"
       The eigenvalues lie in the union of a disk centered at 4 with radius 1 and a disk centered at 3 with radius 1, which corresponds to the real interval $[2, 5]$.

5. **[Derivative] Compute the derivative $A'(t)$ for $A(t) = \begin{pmatrix} t & t^2 \\ 1 & e^t \end{pmatrix}$.**
   ??? success "Solution"
       $A'(t) = \begin{pmatrix} 1 & 2t \\ 0 & e^t \end{pmatrix}$.

6. **[Integral] Evaluate $\int_0^1 \begin{pmatrix} x & 1 \\ 0 & 3x^2 \end{pmatrix} dx$.**
   ??? success "Solution"
       $\begin{pmatrix} 1/2 & 1 \\ 0 & 1 \end{pmatrix}$.

7. **[Gelfand] Prove that for any induced norm, $\rho(A) \le \|A\|$.**
   ??? success "Solution"
       Let $Ax = \lambda x$. Then $\|\lambda x\| = | \lambda | \|x\| = \|Ax\| \le \|A\|\|x\|$. Since $x \neq 0$, we have $|\lambda| \le \|A\|$.

8. **[Inverse Derivative] Using $A(t)A^{-1}(t) = I$, prove $(A^{-1})' = -A^{-1} A' A^{-1}$.**
   ??? success "Solution"
       Differentiate both sides: $A' A^{-1} + A (A^{-1})' = 0 \implies A (A^{-1})' = -A' A^{-1}$. Multiplying by $A^{-1}$ on the left gives the result.

9. **[Trace Derivative] Prove $\frac{d}{dt} \operatorname{tr}(A(t)) = \operatorname{tr}(A'(t))$.**
   ??? success "Solution"
       The trace is a linear combination of entries; differentiation and summation commute.

10. **[Iterative] Why is the spectral radius of $B$ crucial when studying the iteration $x_{k+1} = Bx_k + f$?**
    ??? success "Solution"
        Error evolution follows $e_{k+1} = B e_k$. The error vanishes (convergence) if and only if $B^k \to 0$, which is equivalent to $\rho(B) < 1$.

## Chapter Summary

Matrix analysis endows static algebra with the dimension of time:

1.  **Convergence Criteria**: The spectral radius is the ultimate judge of evolution, determining whether discrete dynamical systems stabilize or explode.
2.  **Analytical Tools**: Matrix calculus provides the symbolic engine for optimization algorithms and variational problems, allowing for dynamic analysis of complex matrix models.
3.  **Global Perspective**: Gelfand's formula and the Circle Theorem unify local entry-wise information with global spectral properties, providing efficient paths for numerical stability estimation.
