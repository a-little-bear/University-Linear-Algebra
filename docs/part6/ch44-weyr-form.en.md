# Chapter 44: Weyr Canonical Form

<div class="context-flow" markdown>

**Prerequisites**: Jordan Canonical Form (Ch12) · Nilpotent Matrices · Invariant Subspaces (Ch42)

**Chapter Outline**: Another Choice for Canonical Forms → Definition of Weyr Blocks → Structure of the Weyr Canonical Form (WCF) as a Large Block Diagonal Matrix → Core Concept: The Weyr Characteristic → Duality with Jordan Characteristic (Conjugate Partitions) → Existence and Uniqueness of WCF → Advantages in Commutative Algebra of Operators → Applications: Commutativity Testing, Dimension of Matrix Algebras, and Simultaneous Canonical Forms

**Extension**: The Weyr Canonical Form is the "dual" of the Jordan Canonical Form; while the Jordan form is more intuitive for differential equations, the Weyr form provides a clearer hierarchical structure for matrix commutativity and advanced problems in algebraic geometry—a powerful tool for investigating the commutants of operators.

</div>

In discussions of similarity canonical forms, the Jordan form almost entirely dominates the scene. However, there exists another structure, equally important and superior in certain algebraic aspects: the **Weyr Canonical Form** (WCF). Unlike the Jordan form, which decomposes space into independent chains, the Weyr form arranges the successive kernels of the generalized eigenspaces "horizontally." This structure makes the relationship between a matrix and its commutant immediately apparent.

---

## 44.1 Definition of Weyr Blocks and WCF

!!! definition "Definition 44.1 (Weyr Block)"
    A **Weyr block** corresponding to eigenvalue $\lambda$ is a block upper triangular matrix:
    $$W = \begin{pmatrix} \lambda I_{n_1} & E_1 & 0 & \cdots \\ 0 & \lambda I_{n_2} & E_2 & \cdots \\ \vdots & \vdots & \ddots & \ddots \\ 0 & 0 & \cdots & \lambda I_{n_k} \end{pmatrix}$$
    where $E_i = \begin{pmatrix} I_{n_{i+1}} \\ 0 \end{pmatrix}$ are identity-like blocks of appropriate dimensions. The sequence $(n_1, n_2, \ldots, n_k)$ is called the **Weyr Characteristic**.

---

## 44.2 Weyr vs. Jordan Characteristics

!!! theorem "Theorem 44.1 (Duality)"
    1.  **Jordan Characteristic**: The partition formed by the sizes of the Jordan blocks.
    2.  **Weyr Characteristic**: The **conjugate partition** of the Jordan characteristic.
    Example: If the Jordan blocks have sizes $(2, 1)$, the Weyr characteristic is $(2, 1)$. If the Jordan block is $(3)$, the Weyr characteristic is $(1, 1, 1)$.

---

## 44.3 Applications in Commutativity

!!! technique "Commutant Criterion"
    A matrix $X$ commutes with a matrix $A$ in Weyr Canonical Form iff $X$ is also block upper triangular and satisfies specific compatibility conditions between its blocks.
    **Significance**: WCF makes finding common invariant subspaces and computing the dimension of the centralizer algebra extremely straightforward.

---

## Exercises

**1. [Basics] Find the Weyr Canonical Form of $J_3(0)$.**

??? success "Solution"
    **Steps:**
    1. $J_3(0)$ has a single block of size 3. Its partition is $(3)$.
    2. Find the conjugate partition: Rotate a row of 3 squares into a column of 3 squares.
    3. The result is $(1, 1, 1)$.
    **Construct Weyr Block**:
    $W = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{pmatrix}$.
    **Conclusion**: For a single Jordan block, the Weyr form is identical to the Jordan form.

**2. [Calculation] If $A = \operatorname{diag}(0, 0)$, what is its Weyr characteristic?**

??? success "Solution"
    **Analysis:**
    1. The Jordan characteristic is $(1, 1)$.
    2. Conjugate partition: Two columns of height 1 combine into one column of height 2.
    3. The result is $(2)$.
    **Conclusion**: The Weyr characteristic is $(2)$. WCF groups the eigenspace into a single diagonal block.

**3. [Structure] Describe the relationship between the dimensions $n_i$ of the identity blocks $I_{n_i}$ in a Weyr block.**

??? success "Solution"
    **Property:**
    Following the properties of Weyr characteristics, $n_1 \ge n_2 \ge \cdots \ge n_k$. This reflects that the increments in the kernel dimensions $\ker(A-\lambda I)^i$ are monotonically non-increasing.

**4. [Duality] Given Jordan block sizes $(2, 2, 1)$, calculate the Weyr characteristic.**

??? success "Solution"
    **Calculation:**
    1. Draw the Young Tableau:
       Row 1: 2 squares
       Row 2: 2 squares
       Row 3: 1 square
    2. Count squares in each column:
       Column 1: 3 squares
       Column 2: 2 squares
    **Conclusion**: The Weyr characteristic is $(3, 2)$.

**5. [Uniqueness] Why is the WCF unique?**

??? success "Solution"
    **Reasoning:**
    1. The Weyr characteristic is uniquely determined by the sequence of ranks $\operatorname{rank}(A-\lambda I)^k$.
    2. Since rank is an invariant under similarity, the Weyr characteristic is an invariant.
    3. Once the order of eigenvalues is fixed, the WCF is perfectly unique.

**6. [Commutativity] If $A$ is in WCF and has distinct eigenvalues, what is the structure of its commutant $X$?**

??? success "Solution"
    **Conclusion:**
    $X$ must be a **diagonal matrix**.
    When eigenvalues are distinct, every Weyr block is $1 \times 1$, so the whole matrix is diagonal. Any matrix commuting with a diagonal matrix (with distinct entries) is itself diagonal.

**7. [Comparison] What is the main morphological difference between Jordan and Weyr forms?**

??? success "Solution"
    **Core Difference:**
    - **Jordan form**: Emphasizes "vertical" depth (the length of Jordan chains).
    - **Weyr form**: Emphasizes "horizontal" width (the dimension of each generalized eigenspace).
    WCF groups generalized eigenvectors of the same rank together, making the hierarchical projection structure of the operator clearer.

**8. [Nilpotent] Prove: If $A$ is nilpotent, the diagonal entries of its WCF are all 0.**

??? success "Solution"
    **Proof:**
    1. All eigenvalues of a nilpotent matrix are 0.
    2. The diagonal blocks of WCF are $\lambda I_{n_i}$.
    3. Substituting $\lambda = 0$ results in diagonal blocks of zero matrices.

**9. [Dimension] If the Weyr characteristic of $A$ is $(n_1, \ldots, n_k)$, what is the dimension of its commutant algebra?**

??? success "Solution"
    **Formula:**
    $\dim \mathcal{C}(A) = \sum_{i=1}^k (2i-1) n_i$.
    This demonstrates the convenience of the Weyr characteristic for precisely calculating the dimension of an operator's centralizer.

**10. [Application] Why is WCF preferred over Jordan form in matrix algebra studies?**

??? success "Solution"
    **Reasoning:**
    The WCF is **block upper triangular** with **scalar** diagonal blocks. This ensures that interactions between blocks follow very simple algebraic rules when dealing with matrix polynomials and equations. It provides a clearer basis for studying the algebraic structures generated by a set of operators.

## Chapter Summary

The Weyr Canonical Form is the perfect representation of operator hierarchies:

1.  **Conjugate Features**: By mapping to the Jordan characteristic via conjugation, WCF proves that matrix similarity classes possess two complementary geometric descriptions, broadening our understanding of operator degeneracy.
2.  **Map of Commutativity**: The block triangular structure of WCF provides intuitive criteria for matrix commutativity, making it the tool of choice for problems involving families of commuting operators.
3.  **Algebraic Clarity**: By grouping generalized eigenvectors of the same rank, WCF simplifies the calculation of dimensions for commutants and sub-algebras, occupying an irreplaceable position in advanced representation theory and matrix analysis.
