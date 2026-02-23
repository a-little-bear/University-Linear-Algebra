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

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

## Chapter Summary

Copositive programming serves as a bridge between continuous and discrete optimization:

1. **Restricted Positivity**: Defined the COP cone as the relaxation of PSD requirements to the non-negative orthant.
2. **Exact Encoding**: Established that COP/CP programs provide exact (non-relaxed) formulations for the hardest problems in computer science.
3. **Hierarchical Computation**: Introduced SOS methods to provide computable, monotonic approximations to the intractable COP constraints.
4. **Duality Anchor**: Leveraged the CP-COP duality to unify the study of non-negative matrix factorization and global optimization.
