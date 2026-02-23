# Chapter 64A: Convex Sets of Matrices

<div class="context-flow" markdown>

**Prerequisites**: Positive Definite Matrices (Ch16) · Optimization (Ch25) · Vector Spaces (Ch04)

**Chapter Outline**: Definition of Convexity in Matrix Spaces → The Space of Symmetric Matrices $S_n$ and the Positive Definite Cone $S_n^+$ → Properties of Convex Cones (Self-duality, Closure) → Definition of Spectrahedra → Linear Matrix Inequalities (LMI) & Convex Representations → Structure of Extreme Points and Faces → Matrix version of Carathéodory’s Theorem → Applications: Geometric Foundations of Semidefinite Programming (SDP) and the Set of Quantum States (Density Matrix Space)

**Extension**: Convex sets of matrices are the geometric foundations of modern optimization theory; they generalize classical linear programming to smooth, non-linear domains with spectral constraints, providing the ultimate perspective for studying linear operator constraints.

</div>

In classical geometry, we study convex sets in $\mathbb{R}^n$. In modern control theory and quantum mechanics, however, the variables are often matrices. **Convex Sets of Matrices** investigates sets within matrix spaces that satisfy convexity axioms. The core object is the **Positive Definite Cone**, which is not only the algebraic generalization of positive numbers but the generator of all convex matrix inequalities.

---

## 64A.1 The Positive Definite Cone $S_n^+$

!!! definition "Definition 64A.1 (The PSD Cone)"
    The set of all $n \times n$ real symmetric positive semi-definite matrices $S_n^+$ is a **closed convex cone** in $S_n$:
    1.  **Convexity**: If $A, B \in S_n^+$, then for any $\lambda \in [0, 1]$, $\lambda A + (1-\lambda) B \in S_n^+$.
    2.  **Conicity**: If $A \in S_n^+$ and $c \ge 0$, then $cA \in S_n^+$.

!!! theorem "Theorem 64A.1 (Self-Duality)"
    The PSD cone is **self-dual** under the Frobenius inner product: $S_n^+ = \{ X \in S_n : \langle X, A \rangle \ge 0, \forall A \in S_n^+ \}$. This property is the root of duality theory in Semidefinite Programming.

---

## 64A.2 Spectrahedra

!!! definition "Definition 64A.2 (Spectrahedron)"
    A **spectrahedron** is the intersection of the positive definite cone with an affine subspace. It can be expressed as the set of points satisfying a **Linear Matrix Inequality (LMI)**:
    $$\{ x \in \mathbb{R}^m : A_0 + x_1 A_1 + \cdots + x_m A_m \succeq 0 \}$$
    **Status**: Spectrahedra are the feasible regions of SDPs, just as polytopes are the feasible regions of LPs.

---

## 64A.3 Extreme Points and Facial Structure

!!! technique "Extreme Points"
    - For the set of positive definite matrices with unit trace (the space of quantum states), the extreme points are precisely all **rank-1 projection matrices** (pure states).
    - The facial structure of a spectrahedron is determined by the **kernels** of its matrices. Increasing the rank of a matrix corresponds to moving into the interior of higher-dimensional faces.

---

## Exercises

1.  **[Basics] Prove that the convex combination of two PSD matrices is PSD.**
    ??? success "Solution"
        For any vector $v$, $v^T(\lambda A + (1-\lambda)B)v = \lambda v^T A v + (1-\lambda) v^T B v$. Since $v^T A v \ge 0$ and $v^T B v \ge 0$ and the weights are non-negative, the sum is non-negative.

2.  **[Dual] Verify: If $A, B \in S_n^+$, then $\operatorname{tr}(AB) \ge 0$.**
    ??? success "Solution"
        Let $B = \sum \lambda_i q_i q_i^T$ be its spectral decomposition. Then $\operatorname{tr}(AB) = \sum \lambda_i q_i^T A q_i$. Since each $\lambda_i \ge 0$ and $q_i^T A q_i \ge 0$, the sum is non-negative.

3.  **[Spectrahedron] Is a disk $\{ (x, y) : x^2 + y^2 \le 1 \}$ a spectrahedron?**
    ??? success "Solution"
        Yes. Using the Schur complement, it can be written as the LMI: $\begin{pmatrix} 1+x & y \\ y & 1-x \end{pmatrix} \succeq 0$.

4.  **[Dimension] What is the dimension of the space of symmetric matrices $S_n$?**
    ??? success "Solution"
        $n(n+1)/2$. This is the dimension of the linear space containing the PSD cone.

5.  **[Extreme] Why is $I/n$ not an extreme point of the set of unit-trace PSD matrices?**
    ??? success "Solution"
        Because $I/n = \frac{1}{n} \sum e_i e_i^T$, it can be written as a convex combination of other points (pure states). Extreme points must have rank 1.

6.  **[Interior] How do you identify an interior point of the PSD cone?**
    ??? success "Solution"
        A matrix is an interior point iff it is **strictly positive definite** ($\operatorname{rank}(A) = n$).

7.  **[Compactness] Is the PSD cone $S_n^+$ compact?**
    ??? success "Solution"
        No, it is unbounded. However, its intersection with the hyperplane $\operatorname{tr}(A)=1$ is compact (the set of density matrices).

8.  **[Separation] Given $A \notin S_n^+$, does there exist a symmetric matrix $H$ such that $\operatorname{tr}(AH) < 0$ but $\operatorname{tr}(BH) \ge 0$ for all $B \in S_n^+$?**
    ??? success "Solution"
        Yes. By the hyperplane separation theorem and the self-duality of the PSD cone, such an $H$ must exist.

9.  **[Control] Why are LMIs used to describe stability regions in control theory?**
    ??? success "Solution"
        Because the Lyapunov stability condition ($A^T P + PA \prec 0$) is itself a linear matrix inequality in the variable $P$, defining a convex region of stable parameters.

****

??? success "Solution"
    

## Chapter Summary

Convex sets of matrices define the physical boundaries of operator space:

1.  **Dominance of the Cone**: The PSD cone, as the natural order structure in matrix space, elevates scalar inequality comparisons to inclusion relations between operators, establishing the geometric framework for stability analysis.
2.  **Expressive Power of Spectrahedra**: Through LMI representations, spectrahedra unify seemingly non-linear geometric bodies (spheres, ellipsoids, semi-algebraic sets) under the language of linear algebra, greatly expanding the scope of computable problems.
3.  **Purity of Extrema**: The correspondence between extreme points and rank proves that "simple structures" (low-rank matrices) are the atoms used to build complex convex forms—an insight supporting modern matrix completion and quantum information theory.
