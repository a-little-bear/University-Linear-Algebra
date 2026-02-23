# Chapter 60: Linear Preserver Problems

<div class="context-flow" markdown>

**Prerequisites**: Linear Transformations (Ch5) · Determinants (Ch3) · Eigenvalues (Ch6) · Matrix Operations (Ch2)

**Chapter Outline**: Framework of LPP → Preserving Determinants (Frobenius Theorem) → Preserving Rank → Preserving Spectrum → Positive Maps and Completely Positive Maps → Preserving Invertibility → Preserving Commutativity → Modern Directions

**Extension**: Linear preserver problems (LPP) characterize the symmetries of matrix algebras; they have profound implications in quantum information (CPTP maps) and operator theory.

</div>

Linear Preserver Problems (LPP) involve characterizing linear maps $\phi: M_n 	o M_n$ that leave certain properties, subsets, or relations invariant. This field, initiated by Frobenius in 1897, reveals the "rigidity" of matrix structures: maps that preserve core properties like the determinant or rank must typically be of a very restricted form, such as similarity or transposition.

---

## 60.1 The Frobenius and Marcus-Moyls Theorems

!!! theorem "Theorem 60.1 (Frobenius Determinant Preserver Theorem, 1897)"
    Let $\phi: M_n(\mathbb{C}) 	o M_n(\mathbb{C})$ be a linear map such that $\det(\phi(A)) = \det(A)$ for all $A$. Then there exist invertible matrices $M, N$ with $\det(MN) = 1$ such that:
    $$\phi(A) = MAN \quad 	ext{or} \quad \phi(A) = MA^T N$$

!!! theorem "Theorem 60.3 (Marcus-Moyls Rank-1 Preserver Theorem)"
    A linear map that preserves rank-1 matrices must be of the form $\phi(A) = MAN$ or $\phi(A) = MA^T N$. Since rank-1 matrices are the "building blocks" of the matrix space, many LPP results reduce to this theorem.

---

## Exercises

1. **[Concept] Define the "Standard Forms" in linear preserver problems.**
   ??? success "Solution"
       Standard forms refer to maps of the type $\phi(A) = MAN$ (Type I) and $\phi(A) = MA^T N$ (Type II). Most LPP results state that maps preserving a specific property must be of one of these two forms.

2. **[Trace Preserver] Show that the transposition map $\phi(A) = A^T$ preserves the trace and the determinant.**
   ??? success "Solution"
       $\operatorname{tr}(A^T) = \sum (A^T)_{ii} = \sum a_{ii} = \operatorname{tr}(A)$. For the determinant, $\det(A^T) = \det(A)$ is a fundamental property of the determinant expansion.

3. **[Spectrum] Prove that a similarity transformation $\phi(A) = PAP^{-1}$ is a spectrum preserver.**
   ??? success "Solution"
       The characteristic polynomial $\det(\lambda I - PAP^{-1}) = \det(P(\lambda I - A)P^{-1}) = \det(P) \det(\lambda I - A) \det(P)^{-1} = \det(\lambda I - A)$. Since the characteristic polynomials are identical, the set of eigenvalues (spectrum) is preserved.

4. **[Rank-1] Why is the preservation of rank-1 matrices considered a crucial step in solving LPPs?**
   ??? success "Solution"
       Rank-1 matrices are the extreme rays of the cone of PSD matrices and the simplest elements of the matrix manifold. Many properties (rank, determinantal degree, etc.) are determined by the behavior of the map on rank-1 elements.

5. **[Invertibility] State the Dieudonné Theorem regarding invertibility preservers.**
   ??? success "Solution"
       Any linear bijection that maps the set of singular matrices to itself (or preserves invertibility) must be of the standard form $\phi(A) = MAN$ or $\phi(A) = MA^T N$.

6. **[Positive Maps] Define a "Positive Map" and provide an example that is not a similarity.**
   ??? success "Solution"
       A map $\phi$ is positive if $A \succeq 0 \implies \phi(A) \succeq 0$. The transposition map $\phi(A) = A^T$ is a positive map but is not a similarity transformation.

7. **[Complete Positivity] Distinguish between positive maps and completely positive (CP) maps.**
   ??? success "Solution"
       A map is completely positive if its extension $\phi \otimes \operatorname{id}_k$ is positive for all $k$. Transposition is positive but not completely positive, while $\phi(A) = VAV^*$ is completely positive.

8. **[Choi's Theorem] What is the Kraus representation of a CP map?**
   ??? success "Solution"
       $\phi(A) = \sum V_i A V_i^*$. This representation guarantees positivity regardless of the dimension of the environment to which the system is coupled.

9. **[Commutativity] Describe the form of a linear map that preserves commutativity ($AB=BA \implies \phi(A)\phi(B)=\phi(B)\phi(A)$).**
   ??? success "Solution"
       For $n \ge 3$, such maps must be of the form $\phi(A) = \alpha P A P^{-1} + f(A)I$ or $\phi(A) = \alpha P A^T P^{-1} + f(A)I$, where $f$ is a linear functional.

10. **[Quantum Mechanics] How are LPP results applied in quantum information theory?**
    ??? success "Solution"
        Quantum channels are represented as CPTP (Completely Positive Trace Preserving) maps. LPP results characterize the allowed physical transformations of quantum states (density matrices) and define the boundaries of entanglement witnesses.

## Chapter Summary

This chapter examines the rigidity of matrix algebra through Linear Preserver Problems:

1. **Foundational Theorems**: Frobenius and Marcus-Moyls theorems established the standard forms for maps preserving the determinant and rank.
2. **Spectral Rigidity**: Demonstrated that spectrum-preserving maps are essentially restricted to similarity and transposition.
3. **Positivity Classes**: Explored the hierarchy between positive and completely positive maps, providing the mathematical basis for quantum operations.
4. **Structural Symmetry**: Revealed that the preservation of essential algebraic relations (like commutativity or invertibility) forces the map into a highly structured geometric form.
