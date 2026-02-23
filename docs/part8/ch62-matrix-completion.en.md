# Chapter 62: Matrix Completion

<div class="context-flow" markdown>

**Prerequisites**: SVD (Ch11) · Optimization Foundations (Ch25) · Matrix Norms (Ch15) · Probability Theory

**Chapter Outline**: Motivation for Matrix Completion (The Netflix Problem) → The Central Role of the Low-rank Assumption → NP-hardness of Rank Minimization → Nuclear Norm Convex Relaxation → Completeness Condition: Incoherence → Core Algorithms: Singular Value Thresholding (SVT) and Alternating Least Squares (ALS) → Candès-Recht Theorem (Sample Complexity) → Applications: Recommender Systems, Image Inpainting, and Sensor Network Localization

**Extension**: Matrix completion is the algebraic magic of inferring global truth from local observations; it exploits the "low-rank structure" in high-dimensional space to break the limits of sampling theorems, proving that information is often more compact than the data itself.

</div>

Imagine you have a massive matrix, but 99% of its entries are missing. Is there any way to correctly fill in the rest? This is the task of **Matrix Completion**. Under the powerful constraint of being Low-rank, local information is sufficient to "spill over" into the entire space through the correlations inherent in the operator structure. This chapter demonstrates this miracle of modern data science.

---

## 62.1 Problem Definition and the Low-Rank Assumption

!!! definition "Definition 62.1 (Matrix Completion Problem)"
    Given a matrix $M$ with observed entries $M_{ij}$ on a sample set $(i,j) \in \Omega$, find the matrix $X$ with the minimum rank that satisfies the observed entries:
    $$\min \operatorname{rank}(X) \quad \text{subject to } P_\Omega(X) = P_\Omega(M)$$
    **Challenge**: Directly minimizing rank is a combinatorial optimization problem, proven to be NP-hard.

---

## 62.2 Convex Relaxation and the Nuclear Norm

!!! technique "Nuclear Norm Relaxation"
    To make the problem tractable, we replace $\operatorname{rank}(X)$ with the **Nuclear Norm** $\|X\|_*$ (the sum of its singular values).
    $$\min \|X\|_* \quad \text{subject to } P_\Omega(X) = P_\Omega(M)$$
    This is a convex optimization problem that can be solved efficiently via Semidefinite Programming (SDP) or specialized iterative algorithms.

---

## 62.3 Conditions for Exact Recovery

!!! theorem "Theorem 62.1 (Candès-Recht Theorem)"
    If an $n \times n$ matrix satisfies the **Incoherence Condition** (i.e., its singular vectors are not aligned with the coordinate axes), and the number of random samples $m$ satisfies:
    $$m \ge C \cdot n r \log^2 n$$
    then nuclear norm minimization will recover the original matrix exactly with very high probability. Here $r$ is the rank of the matrix.

---

## 62.4 Core Algorithms

!!! algorithm "Algorithm 62.1 (Singular Value Thresholding - SVT)"
    1.  Initialize $Y_0 = 0$.
    2.  Compute $X_k = \mathcal{D}_\tau(Y_k)$, where $\mathcal{D}_\tau$ is the singular value shrinkage operator (keeping and reducing singular values $> \tau$).
    3.  Update $Y_{k+1} = Y_k + \delta P_\Omega(M - X_k)$.
    4.  Repeat until convergence.

---

## Exercises

1.  **[Basics] Why can't a rank-1 matrix with only one non-zero row be completed?**
    ??? success "Solution"
        Because this matrix is highly "coherent." If an element in that non-zero row is not sampled, no other row provides information to recover it.

2.  **[Nuclear Norm] Calculate the nuclear norm of $\operatorname{diag}(3, 4, 0)$.**
    ??? success "Solution"
        The singular values are 3, 4, and 0. The nuclear norm is $3 + 4 + 0 = 7$.

3.  **[Sampling] For a $1000 \times 1000$ matrix of rank 10, approximately how many samples are theoretically needed?**
    ??? success "Solution"
        Based on $n r \log n$, roughly $1000 \times 10 \times \ln(1000) \approx 7 \times 10^4$ samples (about 7% of the data).

4.  **[Comparison] What is the link between NMF and matrix completion?**
    ??? success "Solution"
        NMF can also be used for completion; it enforces low-rank through the $WH$ form and adds non-negativity, which is useful in recommendation systems where scores are non-negative.

5.  **[Property] Prove that the nuclear norm is the dual of the operator norm.**
    ??? success "Solution"
        Derived from $\|A\|_* = \sup \{ \operatorname{tr}(A^T B) : \|B\|_2 \le 1 \}$. This generalizes the $L_1/L_\infty$ vector norm duality to matrices.

6.  **[SVT] How does the shrinkage operator $\mathcal{D}_\tau$ handle a singular value of 2 if $\tau=1$?**
    ??? success "Solution"
        It becomes $2 - 1 = 1$. It acts like the soft-thresholding operator in $L_1$ optimization.

7.  **[Application] In the Netflix problem, what do missing values represent?**
    ??? success "Solution"
        Movies that a user has not yet watched or rated. The completion predicts user interest.

8.  **[Uniqueness] Is the completion result unique if samples are too sparse?**
    ??? success "Solution"
        No. Infinitely many high-rank matrices satisfy the observations; the low-rank constraint is the physical assumption used to pick the unique "correct" solution.

9.  **[Noise] How should the constraint $P_\Omega(X) = P_\Omega(M)$ be modified for noisy data?**
    ??? success "Solution"
        Relax it to an inequality: $\|P_\Omega(X - M)\|_F \le \delta$.

10. **[Limits] Why can't gradient descent be applied directly to the rank minimization problem?**

   ??? success "Solution"
        Because the rank function is step-like; its gradient is zero almost everywhere and discontinuous at transitions, providing no directional information for descent.

## Chapter Summary

Matrix completion is the pinnacle application of modern sparsity theory:

1.  **Victory of Correlation**: Proved that in a low-rank context, data are not isolated entries but coupled parts of a whole; local observations suffice to infer the global state via operator consistency.
2.  **Convex Bridge**: The nuclear norm, as the "best convex approximation" of the rank function, transforms unsolvable combinatorial puzzles into tractable convex optimization tasks.
3.  **Structural Information**: The introduction of incoherence provides the algebraic criterion for "high-quality data," establishing the central value of randomized sampling in information recovery.
