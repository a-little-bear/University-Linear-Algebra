# Chapter 61: Inverse Eigenvalue Problems (IEP)

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues (Ch06) · Matrix Analysis (Ch14) · Positive Definite Matrices (Ch16)

**Chapter Outline**: From Forward to Inverse Problems → General Definition of Inverse Eigenvalue Problems (IEP) → Importance of Structural Constraints (Symmetric, Jacobi, Toeplitz) → IEP for Jacobi Matrices → Matrix Construction from Spectral Data (Lanczos Methods) → Criteria for Existence and Uniqueness → Numerical Algorithms: Alternating Projections and Newton's Method → Applications: Structural Design (Bridge Vibration Control), Geophysical Exploration (Determining Stratum Stiffness via Seismic Waves), and Pole Placement in Control Systems

**Extension**: The Inverse Eigenvalue Problem is the "design model" of linear algebra; it elevates us from "analyzing operator properties" to "constructing operators that satisfy specific physical characteristics," serving as the core link between mathematical analysis and engineering synthesis.

</div>

In standard linear algebra textbooks, we learn how to compute the eigenvalues of a given matrix. In engineering practice, however, the problem is often the reverse: we pre-specify a set of ideal eigenvalues (such as the vibration frequencies of a bridge) and must construct a matrix with a specific structure (such as a stiffness matrix) that possesses exactly those eigenvalues. This is known as the **Inverse Eigenvalue Problem** (IEP). This chapter discusses the algebraic challenges and constructive algorithms for these "custom-made" matrices.

---

## 61.1 Definition and Classification

!!! definition "Definition 61.1 (Inverse Eigenvalue Problem)"
    Given a set of scalars $\Lambda = \{\lambda_1, \lambda_2, \ldots, \lambda_n\}$, find a matrix $A$ belonging to a specific class $\mathcal{S}$ such that the spectrum of $A$ is exactly $\Lambda$.
    - **Structural Constraint $\mathcal{S}$**: Typically requires the matrix to be symmetric, tridiagonal (a Jacobi matrix), or have a specific zero pattern.

---

## 61.2 IEP for Jacobi Matrices

!!! theorem "Theorem 61.1 (Uniqueness of Jacobi IEP)"
    Given $n$ distinct real numbers $\{\lambda_i\}$ as eigenvalues for the full matrix, and $n-1$ distinct real numbers $\{\mu_j\}$ as eigenvalues for its $(n-1) \times (n-1)$ leading principal submatrix. If the interlacing property holds:
    $$\lambda_1 < \mu_1 < \lambda_2 < \mu_2 < \cdots < \mu_{n-1} < \lambda_n$$
    then there exists a unique Jacobi matrix (tridiagonal with positive off-diagonals) satisfying these conditions.

---

## 61.3 Numerical Methods

!!! technique "Technique: Alternating Projections"
    Search for a point in the intersection of two sets:
    1.  **Set 1**: Matrices with the target eigenvalues (isospectral orbit).
    2.  **Set 2**: Matrices with the required structure (a linear subspace).
    By repeatedly projecting between these sets, the algorithm can converge to an optimal or exact construction.

---

## Exercises

**1. [Basics] Construct a $2 \times 2$ real symmetric matrix with eigenvalues $\{1, 3\}$.**

??? success "Solution"
    **Construction:**
    The simplest solution is the diagonal matrix:
    $A = \begin{pmatrix} 1 & 0 \\ 0 & 3 \end{pmatrix}$.
    Note: Any matrix of the form $Q \operatorname{diag}(1, 3) Q^T$ for orthogonal $Q$ is also a solution.

**2. [Constraint] Using the eigenvalues from problem 1, find a matrix where the off-diagonal entries must be 1.**

??? success "Solution"
    **Steps:**
    1. Let $A = \begin{pmatrix} a & 1 \\ 1 & c \end{pmatrix}$.
    2. Trace property: $a + c = 1 + 3 = 4$.
    3. Determinant property: $ac - 1 = 1 \cdot 3 = 3 \implies ac = 4$.
    4. Solve: $c = 4-a \implies a(4-a) = 4 \implies a^2 - 4a + 4 = 0$.
    5. $a = 2, c = 2$.
    **Conclusion**: $A = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$.

**3. [Interlacing] Can a Jacobi matrix be constructed with eigenvalues $\{1, 10\}$ and a $1 \times 1$ submatrix eigenvalue of $\{12\}$?**

??? success "Solution"
    **Conclusion: No.**
    **Reasoning**: By the Cauchy Interlacing Theorem, the submatrix eigenvalue $\mu$ must satisfy $\lambda_1 \le \mu \le \lambda_2$, i.e., $1 \le \mu \le 10$. Since $12 > 10$, the problem is unsolvable under the Jacobi constraint.

**4. [Uniqueness] Why is the solution to an IEP non-unique without structural constraints?**

??? success "Solution"
    **Algebraic Intuition:**
    Eigenvalues only determine the "scale" of the matrix, not its "direction." For any orthogonal matrix $Q$, the transform $Q \Lambda Q^T$ preserves the spectrum. This rotational freedom leads to an infinite number of solutions. Only strict structural constraints (like tridiagonalization) consume these degrees of freedom to yield uniqueness.

**5. [Application] What task does the IEP correspond to in bridge seismic design?**

??? success "Solution"
    **Explanation:**
    1. A bridge has inherent resonant frequencies (eigenvalues).
    2. If external noise (wind, traffic) matches these frequencies, the bridge collapses.
    3. Engineers modify the mass and stiffness (matrix entries) so that the eigenvalues of the dynamic system matrix avoid dangerous frequency bands. This is a constrained IEP.

**6. [Calculation] Construct a $2 \times 2$ matrix with trace 0 and determinant -4.**

??? success "Solution"
    **Construction:**
    $A = \begin{pmatrix} 0 & 2 \\ 2 & 0 \end{pmatrix}$.
    Verification: $\operatorname{tr}=0, \det=-4$. Characteristic equation: $\lambda^2 - 4 = 0 \implies \lambda = \pm 2$.

**7. [Control] Is "Pole Placement" in control theory an IEP?**

??? success "Solution"
    **Yes.**
    By state feedback $u = -Kx$, the closed-loop eigenvalues become $\sigma(A-BK)$. The goal is to choose the gain matrix $K$ (the structural constraint) such that the eigenvalues are at specified locations in the left half-plane.

**8. [Property] Prove that if $A$ is a real symmetric Jacobi matrix with non-zero off-diagonals, its eigenvalues must be distinct.**

??? success "Solution"
    **Proof Sketch:**
    1. The matrix $A-\lambda I$ has a principal submatrix of order $n-1$ that is non-singular if off-diagonals are non-zero.
    2. This implies $\operatorname{rank}(A-\lambda I) \ge n-1$ for any $\lambda$.
    3. The geometric multiplicity $\gamma = n - \operatorname{rank} \le 1$.
    4. Since $A$ is diagonalizable, algebraic multiplicity equals geometric multiplicity.
    **Conclusion**: All eigenvalues must be simple (multiplicity 1).

**9. [Numerical] Briefly describe the use of Newton's method for solving IEPs.**

??? success "Solution"
    Treat the eigenvalues as a non-linear function of matrix entries: $f(\mathbf{a}) = \Lambda$. Using the eigenvalue derivative formulas (Ch42), one computes the Jacobian matrix and iteratively updates the entries to minimize the discrepancy from the target spectrum.

**10. [Complexity] How does the difficulty of an IEP with a specific zero-pattern change as $n$ increases?**

??? success "Solution"
    **Conclusion: It increases exponentially.**
    Determining if an IEP with a fixed sparsity pattern has a solution is a highly non-linear problem in algebraic geometry. For large systems, we often solve **least-squares IEPs** to find the structured matrix closest to the target spectrum.

## Chapter Summary

The Inverse Eigenvalue Problem is the "Reverse Engineering" of linear algebra:

1.  **From Observation to Creation**: It breaks the passive analysis of existing systems, establishing algebraic criteria for designing physical structures based on target performance.
2.  **Constraint of Constraints**: IEP proves the deep antagonism between structural sparsity and spectral freedom, revealing the mathematical boundaries of physical feasibility.
3.  **Hub of Computation**: Linking system identification, pole placement, and structural optimization, IEP algorithms are the core logic engine of modern high-precision engineering design.
