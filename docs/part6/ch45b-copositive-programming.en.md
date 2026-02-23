# Chapter 45B: Copositive Programming

<div class="context-flow" markdown>

**Prerequisites**: Positive Definite Matrices (Ch16) · Convex Optimization (Ch25) · Completely Positive Matrices (Ch45A)

**Chapter Outline**: Definition of Copositive (COP) Matrices → Decompositions of the COP Cone → Characterization of Copositivity → Duality between CP and COP → COP Extreme Rays → Copositive Formulation of NP-hard Problems (Max Clique, Independent Set, MAX-CUT) → Computational Complexity → SOS Hierarchies (Parrilo)

**Extension**: Copositive programming provides a unified framework for combinatorial optimization; every NP-hard problem can be expressed exactly as a copositive program.

</div>

A symmetric matrix $A$ is **copositive** if the quadratic form $x^T A x$ is non-negative for all non-negative vectors $x \ge 0$. While this seems like a minor tweak to positive semi-definiteness, it captures immense computational complexity. The **copositive cone** $\mathcal{COP}_n$ is the dual of the completely positive cone $\mathcal{CP}_n$.

---

## 45B.1 Core Concepts

!!! definition "Definition 45B.1 (Copositive Matrix)"
    A matrix $A \in \mathcal{S}_n$ is **copositive** if $x^T A x \ge 0$ for all $x \in \mathbb{R}^n_{\ge 0}$.

!!! theorem "Theorem 45B.1 (Cone Inclusion)"
    For all $n$, $\mathcal{S}_n^+ + \mathcal{N}_n \subseteq \mathcal{COP}_n$. For $n \le 4$, equality holds. For $n \ge 5$, the inclusion is strict (e.g., the Horn matrix).

---

## Exercises

1. **[Fundamentals] Determine if $A = \begin{pmatrix} 1 & -1 \ -1 & 2 \end{pmatrix}$ is copositive.**

   ??? success "Solution"
       Since $A$ is positive definite ($\lambda \approx 0.38, 2.62$), $x^T A x \ge 0$ for all $x$, including $x \ge 0$. Thus $A$ is copositive despite having negative entries.

2. **[Duality] Prove that the dual of the completely positive cone $\mathcal{CP}_n$ is the copositive cone $\mathcal{COP}_n$.**

   ??? success "Solution"
       The generators of $\mathcal{CP}_n$ are $xx^T$ for $x \ge 0$. By definition, $C \in \mathcal{CP}_n^*$ iff $\langle xx^T, C angle = x^T C x \ge 0$ for all $x \ge 0$. This is exactly the definition of copositivity.

3. **[Complexity] What is the computational complexity of checking if a matrix is copositive?**

   ??? success "Solution"
       It is **co-NP-complete**. This reflects the inherent difficulty of the combinatorial optimization problems that copositive programming can encode.

4. **[Independent Set] Express the independence number $\alpha(G)$ as a copositive program.**

   ??? success "Solution"
       $\frac{1}{\alpha(G)} = \min \{ t : tI + A_G \in \mathcal{COP}_n \}$, where $A_G$ is the adjacency matrix. This formula provides an exact value, unlike the SDP theta function which provides an upper bound.

5. **[Small Dimensions] Why does $\mathcal{COP}_n = \mathcal{S}_n^+ + \mathcal{N}_n$ for $n \le 4$?**

   ??? success "Solution"
       For small dimensions, the geometry of the non-negative orthant is simple enough that any copositive quadratic form can be decomposed into a sum of a PSD form and a non-negative form. This breaks down at $n=5$ due to the appearance of pentagonal structures (odd cycles).

6. **[Horn Matrix] Provide an example of a copositive matrix that is not in $\mathcal{S}_5^+ + \mathcal{N}_5$.**

   ??? success "Solution"
       The Horn matrix $H$ (with 1s on diagonal, -1s on super-diagonals, and 1s elsewhere) is the classic counterexample. It is copositive but cannot be written as a sum of a PSD matrix and a non-negative matrix.

7. **[Boundary] Characterize the boundary of the copositive cone.**

   ??? success "Solution"
       A matrix $A$ is on the boundary of $\mathcal{COP}_n$ if and only if there exists a non-zero $x \ge 0$ such that $x^T A x = 0$.

8. **[MAX-CUT] Can MAX-CUT be formulated exactly as a copositive program?**

   ??? success "Solution"
       Yes. All standard NP-hard problems, including MAX-CUT and graph coloring, have exact formulations as linear programs over the CP or COP cones.

9. **[SOS Hierarchy] Briefly describe Parrilo's SDP hierarchy for copositivity.**

   ??? success "Solution"
       It approximates the COP cone by a sequence of nested cones $\mathcal{K}_n^{(r)}$ defined by sum-of-squares (SOS) conditions. As $r 	o \infty$, the sequence converges to the copositive cone.

10. **[LCP] Relate copositive matrices to the Linear Complementarity Problem (LCP).**

   ??? success "Solution"
        If $M$ is copositive, the LCP$(M, q)$ is guaranteed to have a solution whenever it is feasible.

## Chapter Summary

Copositive programming serves as a bridge between continuous and discrete optimization:

1. **Restricted Positivity**: Defined the COP cone as the relaxation of PSD requirements to the non-negative orthant.
2. **Exact Encoding**: Established that COP/CP programs provide exact (non-relaxed) formulations for the hardest problems in computer science.
3. **Hierarchical Computation**: Introduced SOS methods to provide computable, monotonic approximations to the intractable COP constraints.
4. **Duality Anchor**: Leveraged the CP-COP duality to unify the study of non-negative matrix factorization and global optimization.
