# Chapter 62: Matrix Completion Problems

<div class="context-flow" markdown>

**Prerequisites**: SVD (Ch11) · Matrix Norms (Ch15) · Convex Optimization (Ch25)

**Chapter Outline**: From Missing Data to Low-rank Recovery → Mathematical Definition of Matrix Completion → Importance of the Low-rank Assumption → The Core Challenge: Non-convexity of Rank Minimization → Convex Relaxation: Nuclear Norm Minimization → Sampling Rates and Recoverability (Phase Transitions) → Algorithmic Implementation: Singular Value Thresholding (SVT) → Applications: Netflix Recommendation Engine, Image Inpainting in Computer Vision, and Localization in Sensor Networks

**Extension**: Matrix completion reveals the immense power of "informational redundancy" in high-dimensional data; it proves that if the underlying structure is low-rank, we can perfectly reconstruct the whole from a minimal set of random observations—an elegant extension of Compressed Sensing into the two-dimensional domain.

</div>

In modern data science, we often face incomplete observations. For example, in a recommendation system, only a tiny fraction of users rate a small subset of movies. How can we guess the entire rating matrix from these sparse numbers? This is the **Matrix Completion Problem**. Utilizing the **Low-rank** nature of data, we can mathematically prove that, under certain conditions, the missing information can be precisely "calculated." This chapter introduces the theory at the heart of big data intelligence.

---

## 62.1 Mathematical Definition

!!! definition "Definition 62.1 (Matrix Completion)"
    Let $M$ be an unknown $m \times n$ matrix. We observe a subset of entries $M_{ij}$ for $(i,j) \in \Omega$. The goal is to find the matrix $X$ with minimum rank that satisfies the observations:
    $$\min \operatorname{rank}(X) \quad \text{s.t. } X_{ij} = M_{ij}, \forall (i,j) \in \Omega$$

---

## 62.2 Convex Relaxation and Nuclear Norm

!!! note "Computational Barrier"
    Directly minimizing $\operatorname{rank}(X)$ is a combinatorial problem and is NP-hard.

!!! technique "Nuclear Norm Minimization"
    Since the sum of singular values (the **Nuclear Norm** $\|X\|_*$) is the best convex envelope of the rank function, we relax the problem to:
    $$\min \|X\|_* \quad \text{s.t. } X_{ij} = M_{ij}, \forall (i,j) \in \Omega$$
    This is a **convex optimization** problem that can be solved efficiently using Semidefinite Programming (SDP).

---

## 62.3 Recovery Bounds and Algorithms

!!! theorem "Theorem 62.1 (Exact Recovery Theorem)"
    If $M$ satisfies "incoherence" (energy is not concentrated in a few entries) and has rank $r$, then for a sample size $|\Omega| \ge C n r \log^2 n$, $M$ can be reconstructed exactly with high probability via nuclear norm minimization.

---

## Exercises

**1. [Basics] Why is the "low-rank" assumption necessary for matrix completion?**

??? success "Solution"
    **Reasoning:**
    1. A general $m \times n$ matrix has $mn$ degrees of freedom.
    2. Without structure, every missing entry could be any value, making prediction impossible.
    3. **Low-rankness** implies the matrix has only $r(m+n-r)$ independent variables.
    4. When this is much smaller than $mn$, redundancy allows the whole to be determined by a part.

**2. [Calculation] Let $M = \mathbf{u}\mathbf{v}^T$ be a rank-1 matrix. If $M_{11}=1, M_{12}=2, M_{21}=3$, find $M_{22}$.**

??? success "Solution"
    **Derivation:**
    1. A rank-1 matrix satisfies the determinantal condition $M_{11}M_{22} - M_{12}M_{21} = 0$.
    2. Substitute: $1 \cdot M_{22} = 2 \cdot 3$.
    3. $M_{22} = 6$.
    **Conclusion**: Under the low-rank assumption, the missing item is fixed by a unique algebraic constraint.

**3. [Norms] Calculate the nuclear norm $\|A\|_*$ and the operator norm $\|A\|_2$ for $A = \operatorname{diag}(3, 4, 0)$.**

??? success "Solution"
    **Calculation:**
    1. Nuclear norm is the sum of singular values: $\|A\|_* = 3 + 4 + 0 = 7$.
    2. Operator norm is the maximum singular value: $\|A\|_2 = 4$.
    **Significance**: The nuclear norm acts like the $L_1$ norm for singular values, promoting sparsity (low rank).

**4. [Recoverability] Give an example of a low-rank matrix that cannot be completed.**

??? success "Solution"
    **Example: The Sparse Impulse** $M = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$.
    **Reason**: Although rank-1, its energy is perfectly concentrated at $(1,1)$. If $\Omega$ does not include $(1,1)$, we can never know it is non-zero. This is called high **Coherence**.

**5. [Application] Briefly describe matrix completion in the Netflix problem.**

??? success "Solution"
    The matrix has users as rows and movies as columns. Entries are ratings. Since users only see a few movies, 99% of entries are missing. Matrix completion uncovers the "taste dimensions" (low-rank space), allowing the algorithm to predict ratings for un-watched movies and provide personalized recommendations.

**6. [Relaxation] Why use the nuclear norm instead of the rank function?**

??? success "Solution"
    **Reasoning:**
    1. The rank function is non-continuous and non-convex with zero gradients everywhere, making it impossible to optimize directly.
    2. The nuclear norm is the convex envelope of the rank on the unit ball.
    3. It tends to push singular values to zero, enabling the use of mature convex optimization algorithms to solve what was originally a combinatorial nightmare.

**7. [Numerical] Describe one step of the Singular Value Thresholding (SVT) algorithm.**

??? success "Solution"
    **Steps:**
    1. **Update**: $X_{k} = X_{k-1} + \delta \mathcal{P}_\Omega(M - X_{k-1})$ (correcting errors at observed positions).
    2. **Shrink**: $X_{k+1} = \mathcal{D}_\tau(X_k)$ (perform SVD on $X_k$, subtract threshold $\tau$ from singular values, and truncate at 0).
    This iteratively fits data and compresses rank.

**8. [Sampling] How does the required observation ratio change as rank $r$ increases?**

??? success "Solution"
    **Conclusion: It increases linearly.**
    Theoretical bounds show $|\Omega|$ must scale with $nr$. If the system becomes more complex (higher rank), more data must be collected to maintain reconstruction accuracy.

**9. [Properties] Does the nuclear norm satisfy the triangle inequality?**

??? success "Solution"
    **Yes.**
    The nuclear norm is a valid norm (it is the dual of the spectral norm). Thus $\|A+B\|_* \le \|A\|_* + \|B\|_*$, which ensures global convexity of the objective function.

**10. [Application] How does matrix completion solve sensor network localization?**

??? success "Solution"
    In sensor networks, often only distances between nearby nodes are known. The squared distance matrix $D$ has a very low rank (rank 5 in 3D). By completing the missing long-distance entries based on local ones, the global coordinates of all nodes can be determined.

## Chapter Summary

Matrix completion theory is an algebraic miracle in the era of sparse data:

1.  **Low-dimensional Projection**: It reveals the underlying simplicity of high-dimensional phenomena, proving that low-rankness is the bridge across the gap of missing data.
2.  **Wisdom of Relaxation**: The introduction of the nuclear norm transforms an intractable non-convex chasm into a solvable convex path, establishing the algorithmic paradigm for modern statistical learning.
3.  **Philosophy of Reconstruction**: From local fragments to global wholes, matrix completion is not just a tool but a profound law describing the redundancy of information structures and the possibility of recovery.
