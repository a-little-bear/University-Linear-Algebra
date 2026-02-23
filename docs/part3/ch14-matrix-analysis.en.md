# Chapter 14: Matrix Analysis

<div class="context-flow" markdown>

**Prerequisites**: Matrix Functions (Ch13) · Matrix Norms (Ch15) · Series Theory in Calculus

**Chapter Outline**: From Scalar to Matrix Analysis → Matrix Sequences and Convergence (Entry-wise and Norm-wise) → Convergence Criteria for Matrix Series → The Central Role of Spectral Radius $\rho(A)$ → Gelfand’s Formula (The Bridge between Spectrum and Norm) → Convergence of Matrix Powers $A^k$ → Neumann Series and Inverse Matrix Approximation → Introduction to Matrix Calculus → Applications: Iterative Solvers for Linear Systems, Stability Criteria for Dynamical Systems

**Extension**: Matrix analysis is the intersection of numerical analysis and continuous dynamical systems; it introduces the continuous tools of calculus into the discrete world of matrices, serving as the ultimate key to analyzing the long-term stability of massive systems.

</div>

Matrix analysis studies the analytical properties of matrices when treated as variables. Unlike elementary algebra, which focuses on exact equality, matrix analysis emphasizes "closeness," "limits," and "trends." This chapter establishes the standards for determining whether matrix evolutions converge and reveals how the spectral radius acts as the "vital sign" of a system, dictating the long-term behavior of operators.

---

## 14.1 Matrix Sequences and Series

!!! definition "Definition 14.1 (Convergence of Matrix Sequences)"
    A matrix sequence $\{A_k\}$ converges to $A$ if for every entry, $\lim_{k \to \infty} (A_k)_{ij} = a_{ij}$.
    **Equivalence**: This is equivalent to saying that for any matrix norm $\|\cdot\|$, $\lim_{k \to \infty} \|A_k - A\| = 0$.

!!! theorem "Theorem 14.1 (Neumann Series)"
    The series $\sum_{k=0}^\infty A^k$ converges iff the spectral radius $\rho(A) < 1$.
    If it converges, the sum is $(I - A)^{-1}$.

---

## 14.2 Spectral Radius $\rho(A)$ and Gelfand's Formula

!!! definition "Definition 14.2 (Spectral Radius)"
    The **spectral radius** of a square matrix $A$ is the maximum absolute value of its eigenvalues:
    $$\rho(A) = \max \{ |\lambda| : \lambda \in \sigma(A) \}$$

!!! theorem "Theorem 14.2 (Gelfand’s Formula)"
    For any matrix norm $\|\cdot\|$, we have:
    $$\rho(A) = \lim_{k \to \infty} \|A^k\|^{1/k}$$
    This profound conclusion indicates that while a single norm value might be large, the long-term average growth rate of matrix powers is entirely determined by its eigenvalues.

---

## 14.3 Convergence of Matrix Powers $A^k$

!!! theorem "Theorem 14.3 (Convergence Criteria for Powers)"
    1.  $\lim_{k \to \infty} A^k = O \iff \rho(A) < 1$.
    2.  $\{A^k\}$ is bounded $\iff \rho(A) \le 1$ and every eigenvalue with $|\lambda|=1$ has a Jordan block of size 1 (i.e., its algebraic multiplicity equals its geometric multiplicity).

---

## Exercises

**1. [Convergence] Let $A = \begin{pmatrix} 0.5 & 1 \\ 0 & 0.5 \end{pmatrix}$. Does the sequence $A^k$ converge to the zero matrix?**

??? success "Solution"
    **Analysis Steps:**
    1. Extract eigenvalues: Since it is upper triangular, the eigenvalues are $\lambda_1 = 0.5, \lambda_2 = 0.5$.
    2. Calculate spectral radius: $\rho(A) = \max(|0.5|) = 0.5$.
    3. Apply criterion: Since $\rho(A) < 1$, the power sequence must converge.
    **Conclusion**: Yes, $\lim_{k \to \infty} A^k = O$.

**2. [Series] Given $A = \begin{pmatrix} 0.1 & 0 \\ 0 & 0.2 \end{pmatrix}$, compute the series $\sum_{k=0}^\infty A^k$.**

??? success "Solution"
    **Calculation:**
    1. Check convergence: $\rho(A) = 0.2 < 1$. It converges.
    2. Use Neumann series formula: $S = (I - A)^{-1}$.
    3. $I - A = \begin{pmatrix} 0.9 & 0 \\ 0 & 0.8 \end{pmatrix}$.
    4. Inverse: $S = \begin{pmatrix} 1/0.9 & 0 \\ 0 & 1/0.8 \end{pmatrix} = \begin{pmatrix} 1.111 & 0 \\ 0 & 1.25 \end{pmatrix}$.

**3. [Spectral Radius] Give an example of a matrix $A$ such that the spectral norm $\|A\|_2 > 1$ but $\rho(A) < 1$.**

??? success "Solution"
    **Typical Counter-example:**
    Let $A = \begin{pmatrix} 0 & 10 \\ 0 & 0 \end{pmatrix}$.
    1. Eigenvalues are $0, 0$, so $\rho(A) = 0$.
    2. 2-norm (maximum singular value): $A^* A = \begin{pmatrix} 0 & 0 \\ 0 & 100 \end{pmatrix}$, so $\|A\|_2 = 10$.
    **Insight**: Spectral radius only reflects long-term behavior; the norm reflects single-step behavior. This matrix expands in the first step but vanishes from the second step onward.

**4. [Gelfand] If $\|A\| = 0.9$ for some norm, prove that $\|A^k\|^{1/k} \le 0.9$ for all $k \ge 1$.**

??? success "Solution"
    **Proof:**
    1. Use sub-multiplicativity: $\|A^k\| \le \|A\|^k$.
    2. Take the $k$-th root: $\|A^k\|^{1/k} \le (\|A\|^k)^{1/k} = \|A\|$.
    3. Substitute known value: $\|A^k\|^{1/k} \le 0.9$.
    By Gelfand’s Formula, the spectral radius $\rho(A)$ is the limit of this sequence and must satisfy $\rho(A) \le \|A\|$.

**5. [Calculus] Calculate the derivative $A'(t)$ for $A(t) = \begin{pmatrix} e^t & \sin t \\ 0 & 1 \end{pmatrix}$.**

??? success "Solution"
    **Application:**
    The matrix derivative is the matrix of component derivatives.
    $A'(t) = \begin{pmatrix} \frac{d}{dt} e^t & \frac{d}{dt} \sin t \\ \frac{d}{dt} 0 & \frac{d}{dt} 1 \end{pmatrix} = \begin{pmatrix} e^t & \cos t \\ 0 & 0 \end{pmatrix}$.

**6. [Integration] Evaluate $\int_0^1 \begin{pmatrix} x & 1 \\ 0 & 3x^2 \end{pmatrix} dx$.**

??? success "Solution"
    **Calculation:**
    Integrate entry-wise.
    $\begin{pmatrix} \int_0^1 x dx & \int_0^1 1 dx \\ 0 & \int_0^1 3x^2 dx \end{pmatrix} = \begin{pmatrix} 1/2 & 1 \\ 0 & 1 \end{pmatrix}$.

**7. [Inverse Derivative] Using $A(t)A^{-1}(t) = I$, prove $(A^{-1})' = -A^{-1} A' A^{-1}$.**

??? success "Solution"
    **Proof:**
    1. Differentiate $A(t)A^{-1}(t) = I$ using the product rule.
    2. $A'(t)A^{-1}(t) + A(t)(A^{-1}(t))' = O$.
    3. Rearrange: $A(t)(A^{-1}(t))' = -A'(t)A^{-1}(t)$.
    4. Left-multiply by $A^{-1}(t)$: $(A^{-1})' = -A^{-1} A' A^{-1}$.
    This shows how calculus laws balance elegantly in non-commutative algebra.

**8. [Boundedness] If $\rho(A) = 1$, why is the sequence $\{A^k\}$ not necessarily bounded? Give an example.**

??? success "Solution"
    **Reasoning**: Spectral radius only bounds the magnitude of eigenvalues. If there is a non-trivial Jordan block for an eigenvalue with magnitude 1, the power operation will produce polynomial growth.
    **Example**: $A = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$. The eigenvalue is 1, so $\rho(A)=1$.
    However, $A^k = \begin{pmatrix} 1 & k \\ 0 & 1 \end{pmatrix}$. As $k \to \infty$, the entry $k$ grows without bound.

**9. [Trace Derivative] Prove $\frac{d}{dt} \operatorname{tr}(A(t)) = \operatorname{tr}(A'(t))$.**

??? success "Solution"
    **Proof:**
    Trace is a linear summation of components. Differentiation is also linear.
    $\frac{d}{dt} \sum A_{ii}(t) = \sum \frac{d}{dt} A_{ii}(t) = \operatorname{tr}(A'(t))$.
    Since summation and differentiation commute for finite terms, the result holds.

**10. [Application] In the iterative method $x_{k+1} = Bx_k + f$, why is $\rho(B) < 1$ required?**

??? success "Solution"
    **Stability Analysis:**
    1. Define the error $e_k = x_k - x^*$ (where $x^*$ is the exact solution).
    2. Substitute into the iteration: $e_{k+1} = B e_k = B^2 e_{k-1} = \cdots = B^k e_0$.
    3. To guarantee that the error vanishes for any $e_0$ (i.e., $\lim e_k = 0$), we must have $\lim B^k = O$.
    **Conclusion**: According to matrix analysis, this is equivalent to $\rho(B) < 1$. This establishes the spectral radius as the "life-or-death" boundary for iterative numerical solvers.

## Chapter Summary

Matrix analysis endows static algebra with the dimension of time:

1.  **Convergence Criteria**: The spectral radius is the ultimate judge of evolution, determining whether discrete dynamical systems stabilize or explode.
2.  **Analytical Tools**: Matrix calculus and series provide the symbolic engine for non-linear approximation and sensitivity analysis, allowing for dynamic analysis of matrix models.
3.  **Global Perspective**: Gelfand’s Formula reveals the asymptotic consistency between local norms and global spectral structures, proving that operators are essentially defined by their deepest eigen-structures.
