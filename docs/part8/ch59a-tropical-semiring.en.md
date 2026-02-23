# Chapter 59A: Tropical Semiring and Matrix Theory

<div class="context-flow" markdown>

**Prerequisites**: Matrix Operations (Ch2) · Eigenvalues (Ch6) · Graph Theory (Ch27) · Perron-Frobenius Theorem (Ch17)

**Chapter Outline**: Tropical Semiring $(\mathbb{R}\cup\{-\infty\}, \max, +)$ → Min-plus Variant → Tropical Matrix Multiplication and Weighted Directed Graphs → Tropical Determinant = Optimal Assignment → Tropical Singularity → Tropical Eigenvalues → Tropical Perron-Frobenius Theorem → Karp's Algorithm → Kleene Star and Convergence → Path Problems

**Extension**: Tropical linear algebra provides a unified framework for shortest/longest path problems, scheduling, and discrete event systems.

</div>

Tropical linear algebra is the study of linear algebra over the **tropical semiring**, where addition is replaced by the maximum operator and multiplication is replaced by standard addition. This transformation linearizes certain classes of optimization problems, mapping matrix operations directly to path-finding algorithms on graphs.

---

## 59A.1 The Tropical Semiring

!!! definition "Definition 59A.1 (Tropical Semiring)"
    The **tropical semiring** $\mathbb{T} = \mathbb{R} \cup \{-\infty\}$ is equipped with:
    - **Tropical Addition**: $a \oplus b = \max(a, b)$
    - **Tropical Multiplication**: $a \odot b = a + b$
    The additive identity is $\mathbb{0} = -\infty$, and the multiplicative identity is $\mathbb{1} = 0$.

!!! theorem "Theorem 59A.5 (Tropical Perron-Frobenius Theorem)"
    If $A \in \mathbb{T}^{n 	imes n}$ is irreducible (strongly connected), it has a unique tropical eigenvalue $\lambda$ equal to the **maximum cycle mean** of its associated weighted directed graph.

---

## Exercises

1. **[Fundamentals] Compute $(2 \oplus 5) \odot (3 \oplus 4)$ in the tropical semiring.**
   ??? success "Solution"
       $(2 \oplus 5) = \max(2, 5) = 5$. $(3 \oplus 4) = \max(3, 4) = 4$. $5 \odot 4 = 5 + 4 = 9$.

2. **[Matrix Multiplication] Let $A = \begin{pmatrix} 1 & 2 \ 3 & 4 \end{pmatrix}, B = \begin{pmatrix} 0 \ 1 \end{pmatrix}$. Compute the tropical product $A \odot B$.**
   ??? success "Solution"
       $(A \odot B)_1 = \max(1+0, 2+1) = 3$. $(A \odot B)_2 = \max(3+0, 4+1) = 5$. Thus $A \odot B = \begin{pmatrix} 3 \ 5 \end{pmatrix}$.

3. **[Determinant] Compute the tropical determinant of $A = \begin{pmatrix} 3 & 1 \ 2 & 4 \end{pmatrix}$. Is it tropically singular?**
   ??? success "Solution"
       $\operatorname{tdet}(A) = \max(3+4, 1+2) = 7$. Since the maximum is attained by a unique permutation (the identity), the matrix is tropically non-singular.

4. **[Eigenvalues] Define the tropical eigenvalue and its relation to graph cycles.**
   ??? success "Solution"
       A tropical eigenvalue $\lambda$ satisfies $A \odot \mathbf{x} = \lambda \odot \mathbf{x}$. For a strongly connected graph, $\lambda$ is exactly the maximum average weight among all directed cycles in the graph.

5. **[Calculation] Find the tropical eigenvalue of $A = \begin{pmatrix} 2 & 5 \ 3 & 1 \end{pmatrix}$.**
   ??? success "Solution"
       Cycles: $(1,1)$ weight 2; $(2,2)$ weight 1; $(1,2,1)$ weight $(5+3)/2 = 4$. The maximum cycle mean is 4, so $\lambda = 4$.

6. **[Kleene Star] Interpret the Kleene star operator $A^* = I \oplus A \oplus A^{\odot 2} \oplus \dots$ in graph theory.**
   ??? success "Solution"
       $A^*$ represents the all-pairs longest path matrix (in max-plus) or shortest path matrix (in min-plus). $(A^*)_{ij}$ is the maximum weight of any path from vertex $i$ to $j$.

7. **[Linear Equations] Contrast tropical linear equations $A \odot \mathbf{x} = \mathbf{b}$ with the classical case.**
   ??? success "Solution"
       Tropical equations rarely have exact solutions. We typically seek the "maximal sub-solution" $\mathbf{x}^* = A \backslash \mathbf{b}$ (satisfying $A \odot \mathbf{x} \le \mathbf{b}$). An exact solution exists if and only if $A \odot \mathbf{x}^* = \mathbf{b}$.

8. **[Scheduling] What does the tropical eigenvalue represent in a train timetable model?**
   ??? success "Solution"
       It represents the **minimum stable cycle time** of the system. It determines the maximum frequency at which the system can operate without accumulating delays.

9. **[Rank] Why are there multiple definitions of "rank" for tropical matrices?**
   ??? success "Solution"
       Because the equivalence between determinantal rank, row/column rank, and factor rank fails in the tropical setting. Tropical rank (largest non-singular submatrix) is generally less than or equal to factor rank (Barvinok rank).

10. **[Discrete Events] How does the tropical linear recurrence $x(k+1) = A \odot x(k)$ describe system dynamics?**
    ??? success "Solution"
        It models systems where state updates depend on the completion of preceding tasks (e.g., assembly lines). The tropical eigenvalue gives the steady-state growth rate, while the tropical eigenvector gives the stable timing offset between tasks.

## Chapter Summary

This chapter transplants linear algebra to the idempotent structure of the tropical semiring:

1. **Tropical Arithmetic**: Established the $(\max, +)$ system to convert non-linear optimization into linear-style algebraic problems.
2. **Combinatorial Determinants**: Proved that tropical determinants correspond to optimal assignment problems, linking algebra to matching theory.
3. **Spectral Graph Theory**: Revealed tropical eigenvalues as maximum cycle means, providing a bottleneck analysis for scheduling systems.
4. **Path Closures**: Utilized the Kleene star to unify path optimization algorithms like Floyd-Warshall within an algebraic framework.
