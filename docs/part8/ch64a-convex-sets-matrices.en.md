# Chapter 64A: Convex Sets in Matrix Spaces

<div class="context-flow" markdown>

**Prerequisites**: Positive Definite Matrices (Ch16) · Convex Optimization Basics (Ch25) · Matrix Norms (Ch15)

**Chapter Outline**: From Euclidean to Matrix Convex Sets → Definition of Matrix Convex Sets → Core Object: The Positive Semidefinite (PSD) Cone $\mathcal{S}_n^+$ → Matrix Convex Combinations and Completely Positive Maps → Extreme Points and Extreme Rays → Self-dual Properties of Cones → Intersections and Direct Sums of Convex Sets → Applications: Feasible Regions of Linear Matrix Inequalities (LMI), Density Matrix Sets in Quantum Mechanics, and Structural Optimization

**Extension**: The theory of convex sets in matrix spaces is the geometric cornerstone of modern Semidefinite Programming (SDP); it proves that even in high-dimensional operator spaces, simple "linear mixing" preserves connectivity and optimality—key to understanding the stability of complex constrained systems.

</div>

In linear algebra, we typically deal with flat linear spaces. In optimization and quantum mechanics, however, we often study sets of matrices with constraints, such as the set of all positive definite matrices. **Matrix Convex Sets** require not only closure under ordinary linear combinations but also stability under "matrix-weighted" combinations. This deep geometric structure provides the convexity guarantees needed for handling complex non-linear constraints. This chapter introduces the most important convex structure in matrix space: the PSD cone.

---

## 64A.1 Fundamental Definitions

!!! definition "Definition 64A.1 (Matrix Convex Set)"
    A subset $\mathcal{K}$ of the matrix space $M_n$ is **convex** if for any $A, B \in \mathcal{K}$ and $\lambda \in [0, 1]$:
    $$\lambda A + (1-\lambda)B \in \mathcal{K}$$

!!! definition "Definition 64A.2 (Positive Semidefinite Cone $\mathcal{S}_n^+$)"
    The set of all $n \times n$ symmetric positive semidefinite matrices. It is a closed convex cone that plays a role in modern optimization similar to that of the non-negative orthant in traditional programming.

---

## 64A.2 Extreme Points and Geometry

!!! theorem "Theorem 64A.1 (Extreme Rays of the PSD Cone)"
    The **extreme rays** of the PSD cone $\mathcal{S}_n^+$ are exactly the set of all rank-1 matrices $\mathbf{vv}^T$.
    **Physical Meaning**: This means any complex positive operator can be decomposed into a sum of "pure state" projection operators.

---

## 64A.3 Self-Duality

!!! theorem "Theorem 64A.2 (Self-dual Property)"
    The PSD cone is **self-dual**, meaning:
    $(\mathcal{S}_n^+)^* = \{ Y : \operatorname{tr}(YX) \ge 0, \forall X \in \mathcal{S}_n^+ \} = \mathcal{S}_n^+$
    This beautiful symmetry is the underlying reason why SDP duality theory (Ch25) holds.

---

## Exercises

**1. [Basics] Determine if the unit ball $\mathcal{B} = \{ X : \|X\|_2 \le 1 \}$ is a convex set in matrix space.**

??? success "Solution"
    **Proof:**
    1. Take $A, B \in \mathcal{B}$, so $\|A\|_2 \le 1$ and $\|B\|_2 \le 1$.
    2. Consider the convex combination $C = \lambda A + (1-\lambda)B$ for $\lambda \in [0, 1]$.
    3. Use the triangle inequality for norms: $\|C\|_2 \le \lambda \|A\|_2 + (1-\lambda) \|B\|_2$.
    4. Substitute the bounds: $\|C\|_2 \le \lambda(1) + (1-\lambda)(1) = 1$.
    **Conclusion**: It satisfies the definition of convexity, so the unit ball is a convex set.

**2. [Traces] Prove that the set of all PSD matrices with trace 1 is convex.**

??? success "Solution"
    **Proof:**
    1. Let $\operatorname{tr}(A)=1$ and $\operatorname{tr}(B)=1$.
    2. By linearity of the trace: $\operatorname{tr}(\lambda A + (1-\lambda)B) = \lambda \operatorname{tr}(A) + (1-\lambda) \operatorname{tr}(B)$.
    3. $= \lambda(1) + (1-\lambda)(1) = 1$.
    **Conclusion**: The subset is convex. In quantum mechanics, this is exactly the set of **density matrices**.

**3. [Calculation] Determine if the $(x, y)$ region defined by $\begin{pmatrix} x & 1 \\ 1 & y \end{pmatrix} \in \mathcal{S}_2^+$ is convex.**

??? success "Solution"
    **Transformation:**
    1. The matrix belongs to the PSD cone iff $x \ge 0, y \ge 0$, and the determinant $xy - 1 \ge 0$.
    2. This implies $y \ge 1/x$, which is the region above a hyperbola in the first quadrant.
    **Conclusion**: Yes. This region is defined by a Linear Matrix Inequality (LMI), and all LMI feasible sets are convex.

**4. [Extreme Points] Is the identity matrix $I$ an extreme point of the set of trace-1 PSD matrices?**

??? success "Solution"
    **Conclusion: No (for $n > 1$).**
    **Reasoning**: Extreme points must be rank-1 matrices (pure states). The identity matrix has rank $n$ and can be written as a convex combination $\sum \frac{1}{n} e_i e_i^T$. Only rank-1 projectors are extreme points of this set.

**5. [Inclusion] Prove the intersection of two convex sets is convex.**

??? success "Solution"
    **Proof Sketch:**
    If $A, B$ are in the intersection $\mathcal{K}_1 \cap \mathcal{K}_2$, they are in both $\mathcal{K}_1$ and $\mathcal{K}_2$. Since both are convex individually, their combination must be in both, and thus in the intersection.

**6. [Self-duality] Why is self-duality so important in optimization?**

??? success "Solution"
    It allows us to place the primal constraints (e.g., $X \succeq 0$) and dual variables (e.g., $Y \succeq 0$) in the same space and within the same cone. This greatly simplifies the analysis of the Duality Gap and the design of interior-point algorithms.

**7. [Combinations] What is a "Matrix Convex Combination"?**

??? success "Solution"
    **Definition**: $\sum V_i^* A_i V_i$, where $\sum V_i^* V_i = I$.
    **Difference**: It allows coefficients to be matrices (operators) rather than just scalars. This is an advanced convexity structure unique to operator algebras.

**8. [Application] Briefly state the geometric meaning of a Linear Matrix Inequality (LMI).**

??? success "Solution"
    An LMI defines a "slice" or "cross-section" of a high-dimensional space cut out by a PSD cone. Since the cone is convex, the resulting feasible set is always convex. This ensures that finding stabilizing gains in control systems is a well-posed convex optimization problem.

**9. [Property] Prove that the interior of the PSD cone consists of all strictly positive definite matrices $A \succ 0$.**

??? success "Solution"
    **Proof Summary:**
    A strictly PD matrix has all eigenvalues $> 0$. For a sufficiently small perturbation $E$, the shift in eigenvalues (by Weyl's Inequality) is not enough to make them zero or negative. Thus, there exists a neighborhood around $A \succ 0$ entirely contained within the cone.

**10. [Quantum] How do convex sets describe noisy channels in quantum information?**

??? success "Solution"
    A quantum channel can be viewed as a completely positive map from one set of density matrices (a convex set) to another. Because the map preserves convexity, mixed states (uncertainty) can only be maintained or increased during transmission, aligning with the laws of quantum entropy.

## Chapter Summary

Convex sets in matrix spaces are the geometric soul of modern analysis and optimization:

1.  **Order of Constraints**: They prove that complex matrix inequality constraints possess elegant convex shapes, providing a "theoretical pass" for finding global optima.
2.  **Extremality and Synthesis**: Extreme ray theory (rank-1 decomposition) reveals how complex operators are synthesized from fundamental physical entities, establishing atomic analysis methods for systems.
3.  **Symmetry of Duality**: The self-dual property constructs a perfect mirror relationship between primal and dual spaces, supporting the entire architecture of contemporary Semidefinite Programming.
