# Chapter 45A: Completely Positive Matrices

<div class="context-flow" markdown>

**Prerequisites**: Positive Definite Matrices (Ch16) · Non-negative Matrices (Ch17) · Matrix Factorization (Ch10)

**Chapter Outline**: From Positive Definiteness to Non-negative Factorization → Definition of Completely Positive (CP) Matrices → Distinguishing CP from General Non-negative PSD Matrices → Core Criterion: Existence of Factorization $A = BB^T$ with $B \ge 0$ → Definition and Bounds of CP-rank → Relationship with Non-negative Matrix Factorization (NMF) → Applications: Clique Number Estimation in Graph Theory, Semidefinite Programming (SDP) Relaxations in Combinatorial Optimization, and Mixture Models of Probability Distributions

**Extension**: Completely positive matrices are a subset of positive definite matrices with a "non-negative soul"; they require the operator to be not only energetic but also constructible from purely non-negative basis vectors—a perfect bridge between continuous convex optimization and discrete graph theory.

</div>

In Ch16, we learned that if $A = BB^T$, then $A$ is positive definite. If we go a step further and require that every element of $B$ must be non-negative (i.e., $B \ge 0$), we obtain a special class known as **Completely Positive** (CP) matrices. Although this seems like a minor restriction on the factorization, its implications are vastly complex and link deeply to hard problems in graph theory, such as finding the maximum clique. This chapter explores the ultimate structural constraints of positive operators.

---

## 45A.1 Definition and Core Criteria

!!! definition "Definition 45A.1 (Completely Positive Matrix)"
    A square matrix $A$ is **Completely Positive** (CP) if there exists an $n \times k$ **non-negative matrix** $B$ such that:
    $$A = B B^T, \quad B_{ij} \ge 0$$

!!! note "Nuance: CP vs. PSD + Non-negative"
    - **Non-negative PSD Matrix**: A matrix that is both positive semi-definite and entry-wise non-negative ($A \succeq 0$ and $A \ge 0$).
    - **Completely Positive Matrix**: Not only $A \succeq 0$ and $A \ge 0$, but its factors must also be non-negative.
    **Conclusion**: Every CP matrix is PSD and non-negative, but the converse is not true for $n \ge 5$.

---

## 45A.2 CP-rank

!!! definition "Definition 45A.2 (CP-rank)"
    The minimum number of columns $k$ required for a non-negative matrix $B$ such that $A = BB^T$ is called the **CP-rank** of $A$.
    **Bounds**: Generally, $r \le \text{cp-rank}(A) \le n(n+1)/2$, where $r$ is the standard rank. Computing the exact CP-rank is an NP-hard problem.

---

## 45A.3 Applications in Graph Theory

!!! technique "Application: Graph Total Positivity"
    Every graph $G$ is associated with a class of matrices. If all non-negative PSD matrices associated with a graph are completely positive, the graph is called a **CP-graph**.
    **Key Result**: A graph is a CP-graph iff it does not contain odd cycles of length greater than 4 (i.e., it is bipartite or specific types of chordal graphs).

---

## Exercises

**1. [Basics] Determine if $\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$ is a completely positive matrix.**

??? success "Solution"
    **Construct Factorization:**
    1. Observe $A = \begin{pmatrix} 1 \\ 1 \end{pmatrix} \begin{pmatrix} 1 & 1 \end{pmatrix}$.
    2. The factor $B = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$.
    3. All elements of $B$ are non-negative.
    **Conclusion**: Yes, it is CP. Its CP-rank is 1.

**2. [Counter-example] Give an example of a matrix that is PSD and non-negative but not CP (for $n=5$).**

??? success "Solution"
    **Explanation:**
    This is one of the most famous discoveries in CP theory. For $n \le 4$, $A \succeq 0$ and $A \ge 0$ implies $A$ is CP. However, for $n=5$, there are variants like the **Horn matrix**. While every entry is positive and eigenvalues are non-negative, its associated graph contains a 5-cycle (pentagon), making a non-negative $BB^T$ factorization impossible.

**3. [Properties] Prove: The diagonal entries of a CP matrix must be non-negative.**

??? success "Solution"
    **Proof:**
    1. By definition, $a_{ii} = \sum_{j=1}^k B_{ij}^2$.
    2. Since squares of real numbers are non-negative, $a_{ii} \ge 0$.
    Indeed, since CP matrices are PSD, this property is naturally satisfied.

**4. [CP-rank] Estimate the CP-rank of $A = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$.**

??? success "Solution"
    **Analysis:**
    1. Trying a $k=2$ factorization: $A = \begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix} \begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}^T$ contains negative numbers.
    2. However, for $2 \times 2$ matrices, $A \ge 0$ and $A \succeq 0$ is sufficient for CP-ness.
    3. A valid CP factor is $B = \begin{pmatrix} \sqrt{1.5} & \sqrt{0.5} \\ \sqrt{0.5} & \sqrt{1.5} \end{pmatrix}$.
    **Conclusion**: The CP-rank equals the standard rank 2.

**5. [Determination] Prove: If $A$ is CP, then $A$ must be positive semi-definite.**

??? success "Solution"
    **Proof:**
    1. For any vector $x$, $x^T A x = x^T (BB^T) x = (x^T B)(B^T x) = \|B^T x\|^2$.
    2. The square of a norm is always $\ge 0$.
    **Conclusion**: CP matrices are a proper subset of PSD matrices.

**6. [Hadamard] Is the Hadamard product of two CP matrices always a CP matrix?**

??? success "Solution"
    **Yes.**
    **Reasoning**: If $A = BB^T$ and $C = DD^T$, then $A \circ C$ can be factored using the Kronecker product of rows of $B$ and $D$. This ensures the CP cone is closed under Hadamard operations.

**7. [Application] Why are CP matrices related to the "Maximum Clique Problem"?**

??? success "Solution"
    **Explanation:**
    The clique number $\omega(G)$ of a graph $G$ can be exactly obtained by solving an optimization problem over CP matrices:
    $\omega(G) = \max \{ \mathbf{1}^T X \mathbf{1} : \operatorname{tr}(X)=1, X \in CP, X_{ij}=0 \text{ if } (i,j) \notin E \}$.
    This equivalence transforms a combinatorial explosion problem into a matrix determination problem in continuous space.

**8. [NMF] What is the link between CP factorization and Non-negative Matrix Factorization (NMF)?**

??? success "Solution"
    **Connection:**
    CP factorization is a special symmetric form of NMF. While NMF decomposes $V \approx WH$, CP requires the symmetry $V = BB^T$. CP theory provides the deep algebraic foundation for existence, uniqueness, and rank estimation in NMF.

**9. [Convexity] Prove that the set of all $n \times n$ CP matrices is a convex cone.**

??? success "Solution"
    **Proof:**
    1. If $A, C$ are CP, then $A = BB^T$ and $C = DD^T$.
    2. For $\alpha, \beta > 0$, $\alpha A + \beta C = \begin{pmatrix} \sqrt{\alpha}B & \sqrt{\beta}D \end{pmatrix} \begin{pmatrix} \sqrt{\alpha}B & \sqrt{\beta}D \end{pmatrix}^T$.
    3. The new factor is clearly non-negative.
    **Conclusion**: The set is closed under addition and positive scalar multiplication.

**10. [Complexity] How does the computational complexity of determining CP-ness change as $n$ increases?**

??? success "Solution"
    **Conclusion**: It becomes **NP-hard**.
    While simple analytic criteria exist for small dimensions, for general dimensions, there is no known polynomial-time algorithm to accurately determine if a non-negative PSD matrix possesses a non-negative factorization.

## Chapter Summary

Completely positive matrices are among the most structurally constrained operators in linear algebra:

1.  **Factorization Purity**: They require the positivity of the operator to be traceable back to non-negative atomic components, a "purity of origin" that yields a special status in probability mixing and combinatorial counting.
2.  **Algebraic Graph Theory**: CP theory proves that graph topology (such as the presence of odd cycles) can transcend domains, directly dictating the algebraic factorization properties of associated matrices.
3.  **Ultimate Challenge in Optimization**: As one of the most attractive yet difficult cones to handle in semidefinite programming, the study of CP matrices marks the transition from "calculating matrices" to "understanding the laws of matrix composition."
