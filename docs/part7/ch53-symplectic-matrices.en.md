# Chapter 53: Symplectic Matrices and Symplectic Geometry

<div class="context-flow" markdown>

**Prerequisites**: Matrix Algebra (Ch2) · Orthogonality (Ch7) · Characteristic Polynomial (Ch6) · Quadratic Forms (Ch9)

**Chapter Outline**: Symplectic Form $J$ → Symplectic Group $Sp(2n, \mathbb{R})$ → Symplectic Matrices → Eigenvalue Properties → Williamson's Theorem → Hamiltonian Matrices → Symplectic Basis → Application in Classical Mechanics

**Extension**: Symplectic matrices are the mathematical foundation of Hamiltonian mechanics (preserving phase space volume) and linear optics.

</div>

Symplectic matrices preserve a non-degenerate skew-symmetric bilinear form, often represented by the matrix $J = \begin{pmatrix} 0 & I \ -I & 0 \end{pmatrix}$. While orthogonal matrices preserve the Euclidean dot product (distances), symplectic matrices preserve the "symplectic area" in phase space. This preservation is fundamental to the laws of motion in classical mechanics.

---

## 53.1 The Symplectic Group and Symplectic Forms

!!! definition "Definition 53.1 (Symplectic Matrix)"
    A $2n 	imes 2n$ real matrix $M$ is **symplectic** if it satisfies:
    $$M^T J M = J, \quad 	ext{where } J = \begin{pmatrix} 0 & I_n \ -I_n & 0 \end{pmatrix}$$
    The set of all such matrices forms the **symplectic group** $Sp(2n, \mathbb{R})$.

!!! theorem "Theorem 53.1 (Eigenvalue Symmetry)"
    If $\lambda$ is an eigenvalue of a symplectic matrix $M$, then $1/\lambda$, $\bar{\lambda}$, and $1/\bar{\lambda}$ are also eigenvalues with the same multiplicity.

---

## Exercises

1. **[Fundamentals] Show that every symplectic matrix $M$ has $\det M = 1$.**
   ??? success "Solution"
       From $M^T J M = J$, taking the determinant yields $(\det M)^2 \det J = \det J$. Since $J$ is non-singular, $(\det M)^2 = 1$. A more subtle argument (e.g., using the Pfaffian) shows that $\det M$ must be $+1$, never $-1$.

2. **[Dimension] Calculate the dimension of the Lie algebra $\mathfrak{sp}(2n, \mathbb{R})$.**
   ??? success "Solution"
       The Lie algebra consists of Hamiltonian matrices $A$ satisfying $AJ + J A^T = 0$. This condition implies $JA$ is symmetric. The dimension of symmetric $2n 	imes 2n$ matrices is $n(2n+1)$.

3. **[Hamiltonian Matrices] Define a Hamiltonian matrix and its relation to symplectic matrices.**
   ??? success "Solution"
       A matrix $A$ is Hamiltonian if $AJ$ is symmetric. Hamiltonian matrices are the generators of the symplectic group: if $A$ is Hamiltonian, then $e^{At}$ is symplectic for all $t$.

4. **[Spectral Structure] Analyze the spectrum of $M$ if $\lambda = 2+i$ is an eigenvalue.**
   ??? success "Solution"
       The spectrum must contain $\{2+i, 2-i, \frac{2-i}{5}, \frac{2+i}{5}\}$. Eigenvalues of symplectic matrices always appear in quadruplets (or pairs if they lie on the unit circle or real axis).

5. **[Williamson's Theorem] State the significance of Williamson's Theorem for positive definite matrices.**
   ??? success "Solution"
       Every $2n 	imes 2n$ positive definite matrix $V$ can be "symplectically diagonalized": $M^T V M = \operatorname{diag}(\lambda_1, \dots, \lambda_n, \lambda_1, \dots, \lambda_n)$ for some $M \in Sp(2n, \mathbb{R})$. The $\lambda_i$ are called the **symplectic eigenvalues**.

6. **[Composition] Prove that the product of two symplectic matrices is symplectic.**
   ??? success "Solution"
       Let $M_1, M_2$ be symplectic. $(M_1 M_2)^T J (M_1 M_2) = M_2^T (M_1^T J M_1) M_2 = M_2^T J M_2 = J$. Thus the symplectic property is preserved under composition.

7. **[Orthogonal Symplectic] Show that a matrix is both orthogonal and symplectic if and only if it has the form $\begin{pmatrix} A & B \ -B & A \end{pmatrix}$ where $A+iB$ is a complex unitary matrix.**
   ??? success "Solution"
       This reflects the isomorphism $Sp(2n, \mathbb{R}) \cap O(2n) \cong U(n)$. It shows that unitary transformations are exactly those that preserve both the Euclidean and symplectic structures.

8. **[Inversion] Prove that if $M$ is symplectic, then its inverse $M^{-1}$ is also symplectic.**
   ??? success "Solution"
       $M^T J M = J \implies J = (M^T)^{-1} J M^{-1} = (M^{-1})^T J M^{-1}$. Thus $M^{-1}$ satisfies the defining identity.

9. **[Phase Space] In classical mechanics, why is the Jacobian of a Hamiltonian flow always a symplectic matrix?**
   ??? success "Solution"
       Hamilton's equations $\dot{z} = J 
abla H(z)$ describe a flow that preserves the symplectic form $dz_1 \wedge dz_2 + \dots$. By Liouville's theorem, the linearization of this flow (the Jacobian) must preserve $J$, hence it is a symplectic matrix.

10. **[Optics] How are symplectic matrices used in Gaussian beam optics?**
    ??? success "Solution"
        The transfer matrix (ABCD matrix) of an optical system operating on position and momentum (angle) coordinates must be symplectic to satisfy the conservation of etendue (phase space volume).

## Chapter Summary

This chapter explores matrices that preserve the non-Euclidean geometry of phase space:

1. **Defining Symmetries**: Formulated symplectic matrices through the preservation of the standard skew-symmetric form $J$.
2. **Spectral Quadruplets**: Analyzed the unique reciprocal and conjugate symmetry of eigenvalues in $Sp(2n)$.
3. **Hamiltonian Connection**: Linked symplectic groups to their generators (Hamiltonian matrices) used in physics.
4. **Symplectic Normalization**: Utilized Williamson's theorem to extend spectral analysis to positive definite matrices under symplectic transformations.
