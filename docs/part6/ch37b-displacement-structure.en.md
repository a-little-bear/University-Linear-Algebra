# Chapter 37B: Displacement Structure and Fast Algorithms

<div class="context-flow" markdown>

**Prerequisites**: Toeplitz Matrices (Ch37A) · Matrix Decompositions (Ch10) · Matrix Analysis (Ch14)

**Chapter Outline**: From Exact to Approximate Structures → Definition of Displacement Operators → The Concept of Displacement Rank → Toeplitz-like and Hankel-like Matrices → Core Technique: The Schur Algorithm (Layered Elimination and Generators) → Displacement Structure of Inverse Matrices → Generalized Levinson Algorithms → Applications: Ultra-fast Inversion in Radar Signal Processing, Variance Updates in Large-scale Kalman Filtering, and Coding Theory

**Extension**: Displacement structure is the unified language for handling "near-structured" matrices; it proves that even if a matrix is not strictly Toeplitz, as long as its displacement rank is low, it can still be decomposed and inverted with near-linear complexity.

</div>

In Ch37A, we studied matrices with exact shift-invariance. However, in reality, many matrices (such as the inverse of a Toeplitz matrix or a circulant matrix with small perturbations) do not possess this perfect diagonal equality but still retain some "structural signature." **Displacement Structure** quantifies the deviation of a matrix from an ideal template by defining specific difference operators. This chapter introduces how to utilize this deep algebraic trait to achieve ultra-fast algorithms orders of magnitude faster than the traditional $O(n^3)$.

---

## 37B.1 Displacement Operators and Rank

!!! definition "Definition 37B.1 (Displacement Operator)"
    Let $Z$ be the lower-shift matrix. For an $n \times n$ matrix $A$, define its **Stein Displacement** as:
    $$\nabla_{Z, Z^T}(A) = A - Z A Z^T$$
    If the rank of this difference matrix is $r$, we say the **displacement rank** of $A$ is $r$.

!!! note "Defining Structured Matrices"
    - A strictly Toeplitz matrix typically has a displacement rank of at most 2.
    - Any matrix where $r \ll n$ is termed a **Toeplitz-like matrix**.

---

## 37B.2 Schur Algorithm and Generators

!!! technique "Generator Representation"
    If $\operatorname{rank}(\nabla A) = r$, then $A$ can be represented using two sets of vectors $\{g_i, h_i\}$ (known as generators). This means the entire matrix is encoded using only $O(rn)$ elements.

!!! algorithm "Algorithm 37B.1 (The Schur Elimination Algorithm)"
    The Schur algorithm allows elimination to be performed directly on the generators without explicitly forming the full $n \times n$ matrix. Each elimination step corresponds to a fractional linear transformation, resulting in a total complexity of $O(rn^2)$.

---

## 37B.3 Preservation of Structure in Inversion

!!! theorem "Theorem 37B.1 (Invariance of Displacement Rank)"
    If a matrix $A$ has displacement rank $r$, then its inverse $A^{-1}$ also has displacement rank $r$ with respect to the corresponding adjoint operator.
    **Significance**: This result explains why the inverse of a Toeplitz matrix, though not Toeplitz itself, can still be computed extremely efficiently.

---

## Exercises

**1. [Basics] Verify the displacement rank of the identity matrix $I$ with respect to $\nabla(A) = A - ZAZ^T$.**

??? success "Solution"
    **Calculation:**
    1. $ZIZ^T = ZZ^T$.
    2. For an $n \times n$ matrix, $ZZ^T = \operatorname{diag}(0, 1, 1, \ldots, 1)$ (since the last row shifts out and the first row is zeroed).
    3. $\nabla(I) = I - ZZ^T = \operatorname{diag}(1, 0, 0, \ldots, 0)$.
    **Conclusion**: The rank is 1. The identity matrix is one of the simplest matrices in terms of displacement.

**2. [Generator] If $\nabla A = \mathbf{g}\mathbf{h}^T$, write the explicit recovery formula for $A$.**

??? success "Solution"
    **Derivation:**
    Using the geometric series expansion: $A = \sum_{k=0}^{n-1} Z^k (\mathbf{g}\mathbf{h}^T) (Z^T)^k$.
    This shows that $A$ can be written as a sum of low-rank outer products formed by shifted vectors.

**3. [Rank Check] Prove that the displacement rank of an $n \times n$ Toeplitz matrix is $\le 2$.**

??? success "Solution"
    **Proof:**
    Observe the entries of $A - ZAZ^T$. Except for the first row and column, the internal entries $(A)_{i,j} - (A)_{i-1,j-1}$ are zero due to the diagonal constancy. The difference matrix is only non-zero on the boundary, so its rank is at most 2.

**4. [Complexity] For a Toeplitz-like matrix with $r=2$ and $n=1000$, approximately what is the inversion cost?**

??? success "Solution"
    **Estimation:**
    Using the generalized Schur or Levinson algorithm, the complexity is $O(rn^2)$.
    Calculation: $2 \cdot (1000)^2 = 2 \times 10^6$.
    Contrast with $O(n^3)$: $10^9$. The efficiency is increased by approximately 500 times.

**5. [Reciprocity] What is the "reciprocity" of displacement operators?**

??? success "Solution"
    It refers to the relationship that if $\nabla_{F, G}(A) = 0$, then $\nabla_{G^*, F^*}(A^{-1}) = 0$. It guarantees the duality of structure under the inversion operation.

**6. [Reduction] Why is periodic re-orthogonalization of generators needed in the Schur algorithm?**

??? success "Solution"
    Due to rounding errors in floating-point arithmetic, the generator vectors can lose their algebraic rank accuracy. Re-compressing generators via QR or SVD maintains the "compactness" of the displacement structure and improves algorithm stability.

**7. [Application] Briefly describe the role of displacement structure in fast Kalman filtering.**

??? success "Solution"
    In steady-state stochastic processes, the state transition and covariance matrices often exhibit Toeplitz-like structures. Utilizing displacement structure allows the update of the Kalman gain to be reduced from $O(n^3)$ to $O(n^2)$, enabling real-time high-dimensional filtering.

**8. [Sum] What is the maximum displacement rank of the sum of two matrices with displacement ranks $r_1$ and $r_2$?**

??? success "Solution"
    **Conclusion**: $r_1 + r_2$. Rank satisfies sub-additivity under matrix addition.

**9. [Limits] Is there an $O(n \log^2 n)$ algorithm for inverting Toeplitz-like matrices as $n \to \infty$?**

??? success "Solution"
    **Yes.**
    By combining displacement structure with the Fast Fourier Transform (FFT) via a "super-Schur" algorithm or divide-and-conquer methods, the complexity can be reduced to near-linear levels.

**10. [Identity] What is the Gohberg-Semencul formula?**

??? success "Solution"
    It is a formula expressing the inverse of a Toeplitz matrix explicitly as the product of two triangular Toeplitz matrices. It is the most famous closed-form result in displacement theory, revealing how inversion transforms diagonal structures into convolutional ones.

## Chapter Summary

Displacement structure is a modern powerhouse for handling "quasi-structured" data:

1.  **From Exact to Approximate**: It breaks the limit of strict equality, defining the algebraic distance between a matrix and an ideal template via "displacement rank," greatly expanding the reach of fast algorithms.
2.  **Generator Magic**: By compressing massive matrices into a small set of generator vectors, displacement theory achieves exponential leaps in storage efficiency and computational throughput.
3.  **Structural Heredity**: The stability of displacement rank under inversion and Schur complement operations provides a solid theoretical foundation for a unified system of structured numerical linear algebra.
