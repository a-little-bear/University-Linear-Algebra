# Chapter 64A: Convex Sets in Matrix Spaces

<div class="context-flow" markdown>

**Prerequisites**: Positive Definite Matrices (Ch16) · Matrix Inequalities (Ch18) · Optimization (Ch25) · Copositive Matrices (Ch45b)

**Chapter Outline**: Dimensions and Trace Inner Product → PSD Cone → Birkhoff Polytope → Correlation Matrices (Elliptope) → CP and COP Cones → Spectrahedra → Separation Theorems → S-Lemma → Minimax Theorem → SDP Duality

**Extension**: Convex set theory in matrix spaces is the mathematical foundation for Semidefinite Programming (SDP) and quantum information theory (entanglement witness, state space structure).

</div>

Convexity is the central pillar of optimization theory. When generalized from $\mathbb{R}^n$ to matrix spaces, convexity reveals rich geometric structures. The semi-positive definite cone $S_n^+$ is the most fundamental matrix convex cone, and its self-duality enables the powerful framework of Semidefinite Programming.

---

## 64A.1 Fundamental Matrix Convex Sets

!!! definition "Definition 64A.3 (The PSD Cone)"
    $S_n^+ = \{A \in S_n(\mathbb{R}) : A \succeq 0\}$ is a closed convex cone in the space of real symmetric matrices. Its self-duality implies $\langle A, B \rangle = \operatorname{tr}(AB) \ge 0$ for all $A, B \in S_n^+$.

!!! theorem "Theorem 64A.6 (Birkhoff's Theorem)"
    The set of $n \times n$ doubly stochastic matrices, $\mathcal{B}_n$, is a convex polytope whose vertices are exactly the $n!$ permutation matrices.

!!! definition "Definition 64A.10 (Spectrahedron)"
    A spectrahedron is the feasible set of a Linear Matrix Inequality (LMI): $\{x \in \mathbb{R}^d : A_0 + \sum x_i A_i \succeq 0\}$. It is the fundamental object of study in SDP.

---

## Exercises

1. **[Fundamentals] Prove that the PSD cone $S_n^+$ is convex.**
   ??? success "Solution"
       Let $A, B \in S_n^+$ and $t \in [0, 1]$. For any vector $x$, $x^T(tA + (1-t)B)x = t(x^T Ax) + (1-t)(x^T Bx)$. Since $x^T Ax \ge 0$ and $x^T Bx \ge 0$, and the weights $t, 1-t \ge 0$, the sum remains non-negative. Thus the convex combination is in $S_n^+$.

2. **[Trace Inner Product] Compute the inner product $\langle A, B \rangle$ for $A = \begin{pmatrix} 1 & 2 \\ 2 & 1 \end{pmatrix}$ and $B = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$.**
   ??? success "Solution"
       $\langle A, B \rangle = \operatorname{tr}(A^T B) = \operatorname{tr}(AB) = 1(0) + 2(1) + 2(1) + 1(0) = 4$.

3. **[Birkhoff] Describe the general form of all $2 \times 2$ doubly stochastic matrices and verify its dimension.**
   ??? success "Solution"
       Constraints: $a_{11}+a_{12}=1, a_{21}+a_{22}=1, a_{11}+a_{21}=1, a_{12}+a_{22}=1$. Letting $a_{11}=t \in [0,1]$, we have $A = \begin{pmatrix} t & 1-t \\ 1-t & t \end{pmatrix}$. The dimension is $(2-1)^2 = 1$.

4. **[Separation] If a symmetric matrix $A$ has a negative eigenvalue $-1$, find $H \succeq 0$ such that $\operatorname{tr}(AH) < 0$.**
   ??? success "Solution"
       Let $v$ be the unit eigenvector associated with $-1$. Take $H = vv^T \succeq 0$. Then $\operatorname{tr}(AH) = v^T Av = -1 < 0$. Geometrically, $H$ represents a hyperplane separating $A$ from the PSD cone.

5. **[S-Lemma] State the significance of the S-procedure in control theory.**
   ??? success "Solution"
       The S-lemma allows for converting a condition where one quadratic constraint implies another into a single LMI. This enables the verification of stability for non-linear systems using efficient convex optimization solvers.

6. **[Spectrahedra] Identify the geometric shape of $\{ (x, y) : \begin{pmatrix} 1+x & y \\ y & 1-x \end{pmatrix} \succeq 0 \}$.**
   ??? success "Solution"
       The PSD condition requires $1+x \ge 0, 1-x \ge 0$ and the determinant $(1+x)(1-x) - y^2 \ge 0 \implies 1 - x^2 - y^2 \ge 0$. This is a unit disk in $\mathbb{R}^2$.

7. **[Duality] What condition is typically required for strong duality in SDP?**
   ??? success "Solution"
       Strict feasibility (Slater's condition) is required: there must exist a feasible point in the interior of the cone (i.e., $X \succ 0$ or the dual slack $Z \succ 0$).

8. **[Self-duality] Prove that $S_n^+$ is self-dual under the trace inner product.**
   ??? success "Solution"
       The dual cone $K^* = \{Y : \operatorname{tr}(XY) \ge 0, \forall X \succeq 0\}$. Taking $X=vv^T$ shows $v^T Y v \ge 0$, so $Y \succeq 0$. Conversely, the trace of the product of two PSD matrices is always non-negative, so $S_n^+ \subseteq K^*$.

9. **[Minimax] What does the Von Neumann Minimax Theorem guarantee for matrix games?**
   ??? success "Solution"
       It guarantees that for any two-person zero-sum game, there exists a value $V$ and mixed strategies such that neither player can improve their expected outcome by changing only their own strategy.

10. **[Complexity] Contrast the complexity of checking membership in the PSD cone versus the Copositive cone.**
    ??? success "Solution"
        Checking if $A \succeq 0$ is a polynomial-time problem (via eigenvalues). Checking if $A$ is copositive ($x^T Ax \ge 0, \forall x \ge 0$) is NP-hard, showcasing the computational difficulty introduced by adding non-negativity constraints to vectors.

## Chapter Summary

This chapter examines the geometry of matrix spaces through the lens of convexity:

1. **Cone Dominance**: Established the PSD cone as the primary regular cone in matrix analysis, detailing its self-duality and extreme rays.
2. **Combinatorial Convexity**: Linked doubly stochastic matrices to permutations via Birkhoff's theorem, bridging discrete and continuous symmetry.
3. **Spectrahedral Geometry**: Defined the feasible regions of LMIs, showcasing spectrahedra as the generalized "polyhedra" of modern optimization.
4. **Foundations of SDP**: Formulated the duality theory for semidefinite programming, providing the analytical framework for solving matrix-valued constraints.
