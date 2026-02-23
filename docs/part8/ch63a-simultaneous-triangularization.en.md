# Chapter 63A: Simultaneous Triangularization and Diagonalization

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues (Ch6) · Schur Decomposition (Ch14) · Lie Algebra Foundations (Ch58)

**Chapter Outline**: Definition → McCoy's Theorem → Triangularizing Commuting Matrices → Lie's Theorem → Criteria for Simultaneous Diagonalization → Simultaneous Unitary Diagonalization of Normal Matrices → Common Invariant Subspaces

**Extension**: Simultaneous triangularization is the foundation for representation theory and joint spectral analysis (Ch63B); common invariant subspaces link operator algebras to invariant subspace theory.

</div>

When moving from a single matrix to a **set of matrices**, the central question is: Under what conditions can a collection of matrices be brought to a triangular or diagonal form by the *same* similarity transformation? A single matrix over an algebraically closed field can always be triangularized (Schur decomposition), but for a set of matrices, profound algebraic constraints are required.

---

## 63A.1 Criteria for Simultaneous Form

!!! definition "Definition 63A.1 (Simultaneous Triangularization)"
    A set of matrices $\{A_1, \dots, A_m\}$ is simultaneously triangularizable if there exists an invertible matrix $P$ such that $P^{-1} A_i P$ is upper triangular for all $i$. This is equivalent to the existence of a common **complete flag** of invariant subspaces.

!!! theorem "Theorem 63A.1 (McCoy's Theorem, 1936)"
    Over an algebraically closed field, a set of matrices is simultaneously triangularizable if and only if for every non-commutative polynomial $p$, the commutator $[A_i, p(A_1, \dots, A_m)]$ is nilpotent for all $i$.

!!! theorem "Theorem 63A.2 (Commuting Case)"
    Any collection of pairwise commuting matrices ($A_i A_j = A_j A_i$) over an algebraically closed field is simultaneously triangularizable.

---

## Exercises

1. **[Fundamentals] Is the product $AB$ of two upper triangular matrices necessarily upper triangular? What about their commutator $[A, B]$?**
   ??? success "Solution"
       Yes, upper triangular matrices form a subalgebra. For the commutator $[A, B] = AB - BA$, since both $AB$ and $BA$ are upper triangular with diagonal entries $a_{ii}b_{ii}$, the diagonal of the commutator is zero. Thus, $[A, B]$ is **strictly upper triangular** and hence nilpotent.

2. **[Diagonalization] Prove: If a set of matrices is simultaneously diagonalizable, they must pairwise commute.**
   ??? success "Solution"
       Let $P^{-1}A_i P = D_i$. Since diagonal matrices commute ($D_i D_j = D_j D_i$), we have $A_i A_j = (P D_i P^{-1})(P D_j P^{-1}) = P D_i D_j P^{-1} = P D_j D_i P^{-1} = A_j A_i$. Commutativity is a necessary condition.

3. **[McCoy's Theorem] State the core criterion of McCoy's Theorem.**
   ??? success "Solution"
       A set of matrices is simultaneously triangularizable iff all commutators of the form $[A_i, p(A_1, \dots, A_m)]$ are nilpotent. This essentially requires that the non-commutativity of the family is "weak enough" to be compressed into a single triangular framework.

4. **[Lie's Theorem] What is a "Solvable Lie Algebra" and how does it relate to triangularization?**
   ??? success "Solution"
       A solvable Lie algebra is one whose derived series terminates at zero. Lie's Theorem states that all matrices in a solvable Lie algebra over $\mathbb{C}$ can be simultaneously triangularized, generalizing the result for commuting matrices (Abelian algebras).

5. **[Quantum Mechanics] Explain why "commuting observables can be measured simultaneously" using matrix language.**
   ??? success "Solution"
       Commuting Hermitian matrices (observables) can be simultaneously diagonalized. This means there exists a common basis of eigenvectors (states) in which all physical quantities have definite, non-interfering values (eigenvalues).

6. **[Calculation] Determine if $A = \begin{pmatrix} 1 & 1 \ 0 & 1 \end{pmatrix}$ and $B = \begin{pmatrix} 2 & 3 \ 0 & 2 \end{pmatrix}$ are simultaneously diagonalizable.**
   ??? success "Solution"
       They are both already upper triangular, so they are simultaneously triangularizable. They also commute ($AB = BA$). However, since they both have non-trivial Jordan blocks, they are not diagonalizable individually, let alone simultaneously.

7. **[Normal Matrices] Prove: A family of pairwise commuting normal matrices can be simultaneously unitarily diagonalized.**
   ??? success "Solution"
       This is the Spectral Theorem for commuting families. Since normal matrices are diagonalizable and commute, they share common eigenvectors. One can recursively decompose common eigenspaces using the unitary invariance of the normal property.

8. **[Invariant Subspaces] Let $AB=BA$. Prove that the eigenspaces of $A$ are invariant under $B$.**
   ??? success "Solution"
       Let $Av = \lambda v$. Then $A(Bv) = B(Av) = B(\lambda v) = \lambda (Bv)$, so $Bv$ is also an eigenvector of $A$ for eigenvalue $\lambda$.

9. **[Burnside] What is an "Irreducible" set of matrices? Can such a set be simultaneously triangularized?**
   ??? success "Solution"
       An irreducible set has no common non-trivial invariant subspaces. Since simultaneous triangularization requires a full chain of common invariant subspaces, an irreducible set (for $n>1$) can never be simultaneously triangularized.

10. **[Pole Placement] How does simultaneous triangularization relate to the stability of switched linear systems?**
    ??? success "Solution"
        If a set of system matrices $\{A_i\}$ is simultaneously triangularizable, then the stability of the switched system depends only on the eigenvalues of the matrices, significantly simplifying the analysis of the joint spectral radius (Ch63B).

## Chapter Summary

This chapter examines the consistency of structure across sets of matrices:

1. **McCoy's Criterion**: Established the nilpotency of commutators as the ultimate test for joint triangular form.
2. **Commutativity vs. Diagonalization**: Clarified that commutativity is the algebraic equivalent of sharing a common eigenbasis for diagonalizable matrices.
3. **Lie Algebraic Perspective**: Extended triangularization theory to solvable algebras, linking infinitesimal generators to global triangular forms.
4. **Subspace Synchronicity**: Defined simultaneous forms through the existence of common invariant subspace chains (flags).
