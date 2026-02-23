# Chapter 45B: Copositive Matrices and Copositive Programming

<div class="context-flow" markdown>

**Prerequisites**: Positive Definite Matrices (Ch16) · Completely Positive Matrices (Ch45A) · Convex Optimization (Ch25)

**Chapter Outline**: From Global Positivity to Local Sector Positivity → Definition of Copositive Matrices → The Dual Relationship between the Copositive Cone $\mathcal{C}_n$ and the Completely Positive Cone $\mathcal{CP}_n$ → Computational Hardness of Testing Copositivity → Mathematical Form of Copositive Programming → Key Theorem: Exact Relaxations of Combinatorial Optimization Problems (Stability Number, Hamiltonian Cycles) → Applications: Conformation Search in Computational Biology, Non-convex Quadratic Programming in Operations Research

**Extension**: Copositive matrices are the "dual soul" of completely positive matrices; by restricting the positive definiteness criterion to the non-negative orthant, they capture the continuous extrema features of discrete combinatorial problems, serving as a cutting-edge "continuation" tool for solving traditional NP-hard problems.

</div>

While completely positive matrices require an "internal non-negative construction," **Copositive Matrices** require an "external non-negative action." They only require the quadratic form to be non-negative on non-negative vectors, not necessarily on the entire space. This seemingly relaxed condition introduces immense non-convexity and combinatorial complexity due to the "orthant restriction." This chapter introduces this key structure, existing as the dual to completely positive matrices, and its central role in cracking combinatorial optimization challenges.

---

## 45B.1 Definition and Duality

!!! definition "Definition 45B.1 (Copositive Matrix)"
    A symmetric matrix $A \in S_n$ is **Copositive** if for all non-negative vectors $\mathbf{x} \ge 0$:
    $$\mathbf{x}^T A \mathbf{x} \ge 0$$
    The set of all $n \times n$ copositive matrices is denoted by $\mathcal{C}_n$.

!!! theorem "Theorem 45B.1 (Dual Relationship)"
    The copositive cone $\mathcal{C}_n$ and the completely positive cone $\mathcal{CP}_n$ are **dual cones** to each other:
    $$\mathcal{C}_n = (\mathcal{CP}_n)^*$$
    This means that determining if a matrix is completely positive is equivalent to determining if a linear functional is non-negative over the copositive cone.

---

## 45B.2 Hardness and Structure

!!! note "The Challenge of Determination"
    Unlike positive definite matrices, there is no simple polynomial-time algorithm (like leading principal minors) to determine copositivity. The problem has been proven to be **co-NP-complete**.

!!! technique "Standard Copositive Matrices"
    1.  **Non-negative Matrices**: If $A \ge 0$, it is trivially copositive.
    2.  **Positive Definite Matrices**: If $A \succeq 0$, it is trivially copositive.
    3.  **Sum Matrices**: $A = P + N$ where $P \succeq 0$ and $N \ge 0$. Such matrices are called "standard copositive matrices."

---

## 45B.3 Copositive Programming

!!! definition "Definition 45B.2 (Copositive Program)"
    Copositive programming is the optimization of a matrix variable over the copositive cone (or its dual) under a set of linear constraints.
    **Value**: It can represent many NP-hard discrete optimization problems (like the Independent Set or Graph Coloring) as **convex optimization** problems, albeit over a computationally intractable cone.

---

## Exercises

**1. [Basics] Determine if $\begin{pmatrix} 1 & -1 \\ -1 & 1 \end{pmatrix}$ is a copositive matrix.**

??? success "Solution"
    **Verification Steps:**
    1. Calculate the quadratic form: $f(x, y) = x^2 - 2xy + y^2 = (x-y)^2$.
    2. For any real numbers $(x, y)$, $(x-y)^2 \ge 0$.
    3. Especially for non-negative $x, y$, the condition holds.
    **Conclusion**: The matrix is not only copositive but also positive definite.

**2. [Comparison] Determine if $A = \begin{pmatrix} 1 & 2 \\ 2 & 1 \end{pmatrix}$ is copositive.**

??? success "Solution"
    **Calculation:**
    1. $f(x, y) = x^2 + 4xy + y^2$.
    2. Since $x, y \ge 0$, all three terms are non-negative.
    3. For non-zero non-negative vectors, $f(x, y) > 0$.
    **Conclusion**: It is copositive. Note: Its determinant is $-3 < 0$, so it is not positive definite. This illustrates that copositive matrices are a broader class than PD matrices.

**3. [Counter-example] Give an example of a symmetric matrix that is not copositive.**

??? success "Solution"
    **Example:** $A = \begin{pmatrix} -1 & 0 \\ 0 & 1 \end{pmatrix}$.
    **Proof**: Take $\mathbf{x} = (1, 0)^T \ge 0$. Then $\mathbf{x}^T A \mathbf{x} = -1 < 0$. This violates the definition for non-negative vectors.

**4. [Properties] Prove: The diagonal entries of a copositive matrix must be non-negative.**

??? success "Solution"
    **Proof:**
    1. Take the standard basis vector $\mathbf{e}_i \ge 0$.
    2. $\mathbf{e}_i^T A \mathbf{e}_i = a_{ii}$.
    3. By definition of copositivity, this value must be $\ge 0$.

**5. [Closure] Is the sum of two copositive matrices always copositive?**

??? success "Solution"
    **Yes.**
    **Reasoning**: $\mathbf{x}^T (A+B) \mathbf{x} = \mathbf{x}^T A \mathbf{x} + \mathbf{x}^T B \mathbf{x}$. If each term is non-negative in the non-negative orthant, the sum is too. This proves the set of copositive matrices forms a **convex cone**.

**6. [Duality] Why is the dual of the Completely Positive (CP) cone the Copositive cone?**

??? success "Solution"
    **Algebraic Logic:**
    1. A CP matrix is a sum of $BB^T$, each term being $\mathbf{b}\mathbf{b}^T$ with $\mathbf{b} \ge 0$.
    2. By dual cone definition, $A \in \mathcal{CP}^*$ iff $\langle A, X \rangle \ge 0$ for all $X \in \mathcal{CP}$.
    3. $\langle A, \mathbf{b}\mathbf{b}^T \rangle = \mathbf{b}^T A \mathbf{b}$.
    4. This requires $\mathbf{b}^T A \mathbf{b} \ge 0$ for all $\mathbf{b} \ge 0$, which is the exact definition of copositivity.

**7. [Graph Theory] What is the "copositive representation" of a graph $G$?**

??? success "Solution"
    The stability number (maximum independent set) $\alpha(G)$ of any graph $G$ can be expressed exactly as a linear program over copositive matrices. This demonstrates the power of copositivity theory in capturing local non-connectivity in graphs.

**8. [Calculation] How do you determine if a $2 \times 2$ symmetric matrix is copositive?**

??? success "Solution"
    **Criterion:**
    For $\begin{pmatrix} a & b \\ b & c \end{pmatrix}$, it is copositive iff $a \ge 0, c \ge 0$ and ($b \ge 0$ or $b \ge -\sqrt{ac}$). This means the diagonals must be positive, and the off-diagonal must not be "too negative."

**9. [Stability] How are copositive matrices applied in population dynamics?**

??? success "Solution"
    In studying Lotka-Volterra competition models, the global stability of equilibrium points often depends on the behavior of the interaction matrix in the positive orthant. If the interaction matrix is copositive, one can guarantee that population levels (which are always non-negative) will not explode or collapse without bound.

**10. [Complexity] Why is copositive programming called "Hard Convex"?**

??? success "Solution"
    **Deep Contradiction:**
    1. The constraint set is convex and the objective is linear, so it is theoretically convex optimization.
    2. However, the "boundary" of the copositive cone is extremely complex (co-NP-complete) and cannot be characterized by simple polynomial inequalities.
    3. This shows that convexity does not always imply computational ease; certain algebraic structures hide combinatorial properties under a convex shell.

## Chapter Summary

Copositive matrices are the extension of positivity theory into constrained spaces:

1.  **Orthant Intelligence**: They prove that by narrowing the focus from the entire space to the non-negative orthant, we can capture local stability and combinatorial traits that positive definite matrices cannot describe.
2.  **Harmonious Duality**: The dual relationship between the copositive and completely positive cones forms one of the most symmetric architectures in convex analysis, allowing the same problem to be attacked from different angles.
3.  **New Paradigm in Optimization**: Copositive programming blurs the line between continuous optimization and discrete combinatorics. Despite its hardness, it provides a unified and exact algebraic path for solving NP-hard problems.
