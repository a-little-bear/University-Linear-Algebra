# Chapter 42: Invariant Subspaces and Perturbations

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues (Ch06) · Matrix Norms and Perturbations (Ch15) · Matrix Equations (Ch20)

**Chapter Outline**: Definition of Invariant Subspaces → Block Triangularization and Subspace Representation → Distance between Subspaces (Gap Metric) → Stability of Invariant Subspaces → Key Formula: Sylvester Equations and Perturbation Bounds → The Role of Riccati Equations in Estimating Subspace Shifts → Eigenvalue Sensitivity and the Angle between Eigenvectors → Applications: Model Order Reduction (MOR), Decoupling in Control Theory, and Convergence of Numerical Eigen-algorithms (QR Algorithm)

**Extension**: If Chapter 15 focuses on the perturbation of individual eigenvalues, this chapter focuses on the stability of an entire "set of directions" (subspaces); it reveals the ability of a system to maintain its critical dynamical features when the operator is disturbed, forming the basis for modern robust control.

</div>

In advanced applications of linear algebra, we are often more concerned with **Invariant Subspaces** spanned by a set of eigenvectors rather than just individual eigenvalues. For instance, in model reduction, we aim to preserve the most critical modes of a system. However, when a matrix is perturbed, these subspaces "drift." This chapter establishes mathematical measures to quantify this drift and explores why some subspaces are inherently more stable than others.

---

## 42.1 Invariant Subspaces and Block Structures

!!! definition "Definition 42.1 (Invariant Subspace)"
    A subspace $\mathcal{X} \subseteq V$ is **invariant** under $A$ if for every $x \in \mathcal{X}$, we have $Ax \in \mathcal{X}$.
    **Matrix Representation**: In an appropriate basis, $A$ takes a block upper triangular form:
    $$A = \begin{pmatrix} A_{11} & A_{12} \\ 0 & A_{22} \end{pmatrix}$$
    where $A_{11}$ characterizes the dynamics within the subspace.

---

## 42.2 Subspace Distance and Stability

!!! definition "Definition 42.2 (Gap Metric)"
    The distance between two subspaces $\mathcal{X}$ and $\mathcal{Y}$ is typically measured by the distance between their corresponding orthogonal projectors $P_{\mathcal{X}}$ and $P_{\mathcal{Y}}$:
    $$\operatorname{dist}(\mathcal{X}, \mathcal{Y}) = \|P_{\mathcal{X}} - P_{\mathcal{Y}}\|_2$$

!!! theorem "Theorem 42.1 (Stewart’s Perturbation Bound)"
    Let $\mathcal{X}$ be an invariant subspace of $A$ corresponding to the block $A_{11}$. If $A$ is perturbed by $E$, the shift of the resulting subspace $\mathcal{X}'$ depends on the **separation** between $A_{11}$ and $A_{22}$:
    $$\operatorname{sep}(A_{11}, A_{22}) = \inf_{\|X\|=1} \|A_{11}X - X A_{22}\|$$
    **Conclusion**: The smaller the spectral separation, the more sensitive the subspace is to perturbations.

---

## 42.3 Sensitivity Analysis

!!! technique "Eigenvalue Sensitivity"
    For a simple eigenvalue $\lambda$, the perturbation $\Delta \lambda$ is approximately $w^* E v / (w^* v)$, where $v$ and $w$ are the right and left eigenvectors.
    The factor $1/|w^* v|$ is called the **condition number** of the eigenvalue.

---

## Exercises

**1. [Basics] Determine if the subspace spanned by $(1, 0)^T$ in $\mathbb{R}^2$ is invariant under $A = \begin{pmatrix} 2 & 1 \\ 0 & 3 \end{pmatrix}$.**

??? success "Solution"
    **Steps:**
    1. Take a general vector in the subspace: $x = (k, 0)^T$.
    2. Multiply by $A$: $Ax = \begin{pmatrix} 2 & 1 \\ 0 & 3 \end{pmatrix} \begin{pmatrix} k \\ 0 \end{pmatrix} = \begin{pmatrix} 2k \\ 0 \end{pmatrix}$.
    3. The result is still only non-zero in the first component, so it remains in the original subspace.
    **Conclusion**: Yes, it is an invariant subspace.

**2. [Property] Prove: If $\lambda$ is an eigenvalue of $A$, then the eigenspace $E_\lambda$ is an invariant subspace of $A$.**

??? success "Solution"
    **Proof:**
    1. Let $v \in E_\lambda$, so $Av = \lambda v$.
    2. Since $\lambda v$ is just a scaling of $v$, it remains in $E_\lambda$.
    3. Thus $A(E_\lambda) \subseteq E_\lambda$.

**3. [Calculation] For $A = \begin{pmatrix} 1 & 0 \\ 0 & 1.01 \end{pmatrix}$, calculate the separation $\operatorname{sep}$ between $A_{11}=1$ and $A_{22}=1.01$.**

??? success "Solution"
    **Steps:**
    1. In the scalar case, $\operatorname{sep}(a, b) = |a - b|$.
    2. $\operatorname{sep} = |1 - 1.01| = 0.01$.
    **Significance**: Because the separation is tiny, even a small perturbation $E = \begin{pmatrix} 0 & \epsilon \\ \epsilon & 0 \end{pmatrix}$ can cause a large rotation of the eigenvectors.

**4. [Eigenvectors] Why are subspaces of normal matrices usually more stable than those of non-normal matrices?**

??? success "Solution"
    **Analysis:**
    1. For a normal matrix, left and right eigenvectors are the same, so $|w^* v| = \|v\|^2 = 1$.
    2. The eigenvalue condition number is always 1.
    3. At the subspace level, invariant subspaces of normal matrices are orthogonal, maximizing the separation between blocks (in the unitary equivalence sense).

**5. [Gap] If two subspaces coincide, what is their gap distance? What if they are orthogonal?**

??? success "Solution"
    **Conclusion:**
    1. Coincident: Distance is 0.
    2. Orthogonal: Distance is 1 (for the 2-norm).

**6. [Riccati] The shift of an invariant subspace is related to which type of matrix equation?**

??? success "Solution"
    **Conclusion:**
    It is related to the **Non-linear Algebraic Riccati Equation**.
    In deriving subspace perturbation bounds, the offset $P$ of the basis vectors is expressed as the root of a quadratic matrix equation. For small perturbations, this is often linearized to a Sylvester equation.

**7. [Application] Briefly state the significance of "Eigenvalue Clusters" in perturbation theory.**

??? success "Solution"
    If a set of eigenvalues is very close to each other but far from the rest of the spectrum, the **total subspace** spanned by this cluster is very stable, even if the individual eigenvectors within the cluster are highly unstable. This is key to preserving "slow modes" in model reduction.

**8. [Stability] If $A$ has eigenvalues 10 and 0.1, and due to rounding errors it becomes $A+E$ with eigenvalues 10.001 and 0.099, assess the stability of the subspace.**

??? success "Solution"
    **Assessment:**
    Since $\operatorname{sep}(10, 0.1) \approx 9.9$ is large, the subspace shift will be very small (proportional to roughly $\|E\|/9.9$). This is a very stable subspace configuration.

**9. [Intersection] Prove that the intersection of two invariant subspaces is also an invariant subspace.**

??? success "Solution"
    **Proof:**
    1. Let $x \in \mathcal{X} \cap \mathcal{Y}$.
    2. Then $x \in \mathcal{X}$ and $x \in \mathcal{Y}$.
    3. By invariance, $Ax \in \mathcal{X}$ and $Ax \in \mathcal{Y}$.
    4. Thus $Ax \in \mathcal{X} \cap \mathcal{Y}$.

**10. [Numerical] Why is computing the Schur form (quasi-triangular) more robust than computing the eigenvector basis directly?**

??? success "Solution"
    **Reasoning:**
    Eigenvector bases (especially for non-normal matrices) can be near-linearly dependent, leading to a huge condition number for the transformation matrix $V$. The Schur form uses unitary transformations ($Q$ has condition number 1), maintaining numerical stability and revealing the hierarchical structure of invariant subspaces through diagonal blocks.

## Chapter Summary

The perturbation theory of invariant subspaces is the "stress test" for linear system structures:

1.  **Directional Inertia**: It proves that the core evolutionary directions of a linear system are not only defined by the operator but are also closely guarded by the gaps in the spectral distribution (Separation).
2.  **Robustness of Decoupling**: Subspace stability is the criterion for determining whether a complex system can be safely decomposed into multiple low-dimensional subsystems, providing precision guarantees for industrial-grade model reduction.
3.  **Geometric Metrics**: Via the gap metric and operator equations, this chapter transforms the intuitive "drift of direction" into rigorous algebraic inequalities, establishing the reliability boundaries for numerical eigen-analysis.
