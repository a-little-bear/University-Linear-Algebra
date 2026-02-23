# Chapter 25: Linear Algebra in Optimization

<div class="context-flow" markdown>

**Prerequisites**: Linear Equations (Ch01) · Positive Definite Matrices (Ch16) · Matrix Manifolds (Ch24) · Convex Sets (Ch64A)

**Chapter Outline**: The Linear Algebraic Essence of Optimization → Standard Form and Matrix Representation of Linear Programming (LP) → Algebraic Steps of the Simplex Method → Strong Duality and Farkas’s Lemma → Quadratic Programming (QP) and KKT Conditions → Core Advancement: Semidefinite Programming (SDP) → Algebraic Logic of Interior Point Methods → Applications: Portfolio Optimization, Support Vector Machines (SVM), Compressed Sensing, and Sparse Coding

**Extension**: Optimization is the "action guide" of linear algebra; it elevates systems of equations from "finding equality" to "seeking optimality." It proves that matrix inequality constraints define the convex geometric boundaries of high-dimensional space, serving as the mathematical engine for parameter updates in all AI algorithms.

</div>

In previous chapters, we focused on solving $Ax = b$. In the real world, however, resources are finite, and our goal is often to maximize profit or minimize cost subject to a set of linear constraints. This is the domain of **Optimization Theory**. Linear algebra provides the core language for optimization: the basis rotations in the simplex method are essentially base changes, while semidefinite programming is the ultimate application of eigenvalue theory in constrained spaces. This chapter demonstrates how linear algebra serves as the mathematical tool for extracting maximum value from limited resources.

---

## 25.1 Linear Programming and Simplex Algebra

!!! definition "Definition 25.1 (Standard LP Form)"
    $$\min \mathbf{c}^T \mathbf{x} \quad \text{s.t. } A\mathbf{x} = \mathbf{b}, \mathbf{x} \ge 0$$
    **Simplex Logic**: Start from a basic feasible solution (a vertex) of the augmented matrix and move to an adjacent vertex that decreases the objective function via elementary row operations (pivoting).

---

## 25.2 Duality Theory and Farkas's Lemma

!!! theorem "Theorem 25.1 (Strong Duality)"
    If the primal problem has an optimal solution, then its dual problem $\max \mathbf{b}^T \mathbf{y}$ also has an optimal solution, and their optimal values are equal.
    **Algebraic Essence**: This is a direct consequence of Farkas’s Lemma, revealing the deep symmetry between "feasibility" and "hyperplane separation."

---

## 25.3 Semidefinite Programming (SDP)

!!! definition "Definition 25.2 (Semidefinite Programming)"
    $$\min \operatorname{tr}(CX) \quad \text{s.t. } \operatorname{tr}(A_i X) = b_i, X \succeq 0$$
    **Significance**: SDP is the pinnacle of convex optimization, extending the variables of linear programming from vectors to matrices, enabling precise handling of complex matrix eigenvalue constraints.

---

## Exercises

**1. [Basics] Transform $\min x_1 - x_2$ subject to $x_1+x_2 \le 4, x_1, x_2 \ge 0$ into standard form.**

??? success "Solution"
    **Conversion Steps:**
    1. Introduce a slack variable $s_1 \ge 0$.
    2. Convert the inequality to an equality: $x_1 + x_2 + s_1 = 4$.
    3. The objective vector is $\mathbf{c} = (1, -1, 0)^T$.
    **Standard Form**: $\min (1, -1, 0) \begin{pmatrix} x_1 \\ x_2 \\ s_1 \end{pmatrix}$ subject to $(1, 1, 1) \mathbf{x} = 4$ and $\mathbf{x} \ge 0$.

**2. [Simplex] What is the algebraic basis for choosing an entering variable in the Simplex method?**

??? success "Solution"
    **Reasoning:**
    Calculate the **Reduced Cost** $\bar{c}_j = c_j - \mathbf{c}_B^T B^{-1} A_j$. If the reduced cost of a non-basic variable is negative, increasing that variable (letting it enter the basis) will decrease the objective value. This is essentially finding the steepest descent direction along the edges of the feasible region.

**3. [Duality] Write the dual of $\min \mathbf{c}^T \mathbf{x}, A\mathbf{x} \ge \mathbf{b}, \mathbf{x} \ge 0$.**

??? success "Solution"
    **Conclusion:**
    $\max \mathbf{b}^T \mathbf{y}, A^T \mathbf{y} \le \mathbf{c}, \mathbf{y} \ge 0$.
    **Physical Meaning**: The primal problem seeks the lowest cost under constraints, while the dual seeks the highest revenue under resource pricing constraints.

**4. [SDP] Is Linear Programming (LP) a special case of Semidefinite Programming (SDP)?**

??? success "Solution"
    **Yes.**
    **Reasoning**: If we restrict the matrix variable $X$ in SDP to be diagonal, the constraint $X \succeq 0$ reduces to $x_i \ge 0$ for its diagonal entries. The trace operations become standard dot products. Thus, LP is simply SDP performed on the subspace of diagonal matrices.

**5. [KKT] In quadratic programming $\min \frac{1}{2}\mathbf{x}^T Q \mathbf{x} + \mathbf{c}^T \mathbf{x}$, what property of $Q$ guarantees a global optimum?**

??? success "Solution"
    **Conclusion: $Q$ must be positive semi-definite ($Q \succeq 0$).**
    **Reasoning**: Positive semi-definiteness ensures the objective function is convex. In convex optimization, any local minimum is a global minimum. If $Q$ had negative eigenvalues, the problem would be non-convex and finding a global solution would be extremely difficult.

**6. [Calculation] Using duality: if a primal feasible solution has value 10 and a dual feasible solution has value 12, is this possible for a minimization primal?**

??? success "Solution"
    **Conclusion: No.**
    **Reasoning**: By the **Weak Duality Theorem**, the value of any dual feasible solution (maximization) is always less than or equal to the value of any primal feasible solution (minimization). Since $12 > 10$, this contradicts the duality principle.

**7. [Application] How does the "Kernel Trick" in SVM demonstrate the power of linear algebra?**

??? success "Solution"
    The primal SVM problem seeks a hyperplane in a high-dimensional space. Through duality, the problem is transformed into one involving only inner products $x_i^T x_j$ between samples. Using **Positive Definiteness**, we can replace these inner products with a kernel function $K(x_i, x_j)$, enabling non-linear classification without ever explicitly computing high-dimensional mappings.

**8. [Interior Point] What is the core advantage of Interior Point methods over the Simplex method?**

??? success "Solution"
    **Comparison:**
    - **Simplex**: Moves along the boundary (vertices); worst-case complexity is exponential (Klee-Minty examples).
    - **Interior Point**: Moves through the interior of the feasible set using Newton's method along a "central path."
    **Conclusion**: Interior Point methods have **polynomial-time complexity**, making them superior for very large-scale problems or those with matrix constraints like SDP.

**9. [Farkas] What does Farkas’s Lemma describe in terms of matrix equations?**

??? success "Solution"
    It describes an "either-or" alternative: Either there exists $x \ge 0$ such that $Ax=b$ (the target point is in the cone), or there exists $y$ such that $A^T y \ge 0$ and $b^T y < 0$ (a hyperplane separates the target from the cone). This is the algebraic template for existence proofs in constrained optimization.

**10. [Application] What is $\ell_1$-minimization in "Compressed Sensing"?**

??? success "Solution"
    **Explanation:**
    We seek the sparsest solution ($\min \|\mathbf{x}\|_0$) subject to $A\mathbf{x} = \mathbf{b}$. Since $\ell_0$ optimization is NP-hard, linear algebra proves that under certain conditions (the RIP property), it can be exactly replaced by **linear programming** ($\min \|\mathbf{x}\|_1$). This reduction to convex optimization is a milestone in modern signal processing.

## Chapter Summary

Linear algebra is the engine and the measure of optimization theory:

1.  **Geometric Perspective**: It transforms dry numerical constraints into polyhedra and cones in high-dimensional space, establishing geometric paths for finding optima.
2.  **Duality Harmony**: Strong duality reveals the deep algebraic balance within optimization, proving that the lower bound of cost and the upper bound of value must coincide at the optimum.
3.  **From Vectors to Operators**: The leap from LP to SDP marks the transition of linear algebra from handling simple discrete allocations to complex system stability and structural matrix optimization—the necessary path to modern AI algorithms.
