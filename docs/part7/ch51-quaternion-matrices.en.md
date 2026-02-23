# Chapter 51: Quaternion Matrices

<div class="context-flow" markdown>

**Prerequisites**: Matrix Algebra (Ch2) · Quaternions (Ch50) · Eigenvalues (Ch6) · SVD (Ch11)

**Chapter Outline**: Quaternion Algebra $\mathbb{H}$ → Matrices over $\mathbb{H}$ → Complex Representation of Quaternion Matrices → Left and Right Eigenvalues → Canonical Forms → Quaternion Singular Value Decomposition (QSVD) → Applications in Signal Processing

**Extension**: Quaternion matrices are used in color image processing (representing RGB as a pure quaternion) and attitude control systems.

</div>

Matrices with quaternion entries generalize complex matrices. Due to the non-commutativity of quaternions ($ij = -ji$), the theory of quaternion matrices differs significantly from standard linear algebra. Specifically, one must distinguish between **left** and **right** eigenvalues, and standard determinantal definitions must be modified (e.g., Dieudonné determinant).

---

## 51.1 Quaternion Algebra and Matrix Representation

!!! definition "Definition 51.1 (Complex Representation)"
    A quaternion matrix $A = A_1 + A_2 j$ (where $A_1, A_2 \in M_n(\mathbb{C})$) can be represented by a $2n 	imes 2n$ complex matrix:
    $$\chi(A) = \begin{pmatrix} A_1 & A_2 \ -\bar{A}_2 & \bar{A}_1 \end{pmatrix}$$
    This mapping preserves matrix addition and multiplication.

!!! theorem "Theorem 51.1 (Right Eigenvalues)"
    Every square quaternion matrix $A$ has at least one **right eigenvalue** $\lambda \in \mathbb{H}$ satisfying $Aq = q\lambda$ for some $q 
eq 0$.

---

## Exercises

1. **[Non-commutativity] Show that if $\lambda$ is a right eigenvalue of $A$, then any $s^{-1}\lambda s$ is also a right eigenvalue.**
   ??? success "Solution"
       Let $Aq = q\lambda$. For any $s 
eq 0$, $A(qs) = (Aq)s = (q\lambda)s = (qs)(s^{-1}\lambda s)$. Thus $s^{-1}\lambda s$ is an eigenvalue with eigenvector $qs$. This implies right eigenvalues exist in "equivalence classes" (spheres in $\mathbb{H}$).

2. **[Complex Mapping] Map the quaternion $q = 1 + i + j + k$ to its $2 	imes 2$ complex representation.**
   ??? success "Solution"
       $q = (1+i) + (1+i)j$. $\chi(q) = \begin{pmatrix} 1+i & 1+i \ -1+i & 1-i \end{pmatrix}$.

3. **[Standard Eigenvalues] Why do we typically focus on the "complex part" of the right eigenvalues?**
   ??? success "Solution"
       Each equivalence class of right eigenvalues contains exactly one pair of complex conjugate numbers $\alpha \pm \beta i$ (with $\beta \ge 0$). These are called the **standard eigenvalues** and uniquely represent the class.

4. **[Determinant] Why is the standard definition of a determinant problematic for quaternion matrices?**
   ??? success "Solution"
       The Leibniz formula $\sum \operatorname{sgn}(\sigma) \prod a_{i,\sigma(i)}$ depends on the order of multiplication. Swapping entries changes the value, making it non-invariant. The **Dieudonné determinant** or the determinant of the complex representation $\det(\chi(A))$ are used instead.

5. **[QSVD] State the form of the Quaternion Singular Value Decomposition.**
   ??? success "Solution"
       For any $A \in M_{m 	imes n}(\mathbb{H})$, there exist unitary quaternion matrices $U, V$ such that $U^* A V = \Sigma$, where $\Sigma$ is a real diagonal matrix of singular values.

6. **[Hermitian] Define a Hermitian quaternion matrix and its spectral properties.**
   ??? success "Solution"
       $A^* = A$, where $A^*$ is the conjugate transpose (using quaternion conjugation). Hermitian quaternion matrices have real standard eigenvalues and are unitarily diagonalizable.

7. **[Trace] Is the trace $\operatorname{tr}(AB)$ equal to $\operatorname{tr}(BA)$ for quaternion matrices?**
   ??? success "Solution"
       Generally no, due to $q_1 q_2 
eq q_2 q_1$. However, the real part of the trace $\operatorname{Re}(\operatorname{tr}(AB)) = \operatorname{Re}(\operatorname{tr}(BA))$ remains valid.

8. **[Inversion] How can you compute the inverse of a quaternion matrix $A$ using its complex representation $\chi(A)$?**
   ??? success "Solution"
       Compute the complex inverse $B = (\chi(A))^{-1}$. The quaternion inverse $A^{-1}$ is then recovered from the $n 	imes n$ blocks of $B$ using the inverse of the complex representation mapping.

9. **[Rank] Define the rank of a quaternion matrix.**
   ??? success "Solution"
       The rank is the maximum number of left (or right) linearly independent rows (or columns). For quaternion matrices, row rank always equals column rank.

10. **[Signal Processing] Explain how a color pixel (RGB) can be represented as a pure quaternion.**
    ??? success "Solution"
        A pixel is represented as $q = R i + G j + B k$. This allows for processing all color channels simultaneously using quaternion operators, preserving the correlation between channels during filtering or rotation.

## Chapter Summary

This chapter extends linear algebra to the non-commutative domain of quaternions:

1. **Algebraic Complexity**: Highlighted the distinctions between left and right linear structures necessitated by non-commutativity.
2. **Representation Theory**: Utilized complex matrix embeddings to solve quaternion problems using standard numerical tools.
3. **Spectral Classes**: Defined standard eigenvalues as representatives of similarity orbits in the quaternion spectrum.
4. **Unitary Geometry**: Established the existence of QSVD and the spectral theorem for Hermitian matrices over $\mathbb{H}$.
