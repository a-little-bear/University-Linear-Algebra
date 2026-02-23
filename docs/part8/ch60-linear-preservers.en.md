# Chapter 60: Linear Preserver Problems

<div class="context-flow" markdown>

**Prerequisites**: Linear Transformations (Ch5) · Determinants (Ch3) · Eigenvalues (Ch6) · Matrix Operations (Ch2)

**Chapter Outline**: Framework of LPP → Preserving Determinants (Frobenius Theorem) → Preserving Rank → Preserving Spectrum → Positive Maps and Completely Positive Maps → Preserving Invertibility → Preserving Commutativity → Modern Directions

**Extension**: Linear preserver problems (LPP) characterize the symmetries of matrix algebras; they have profound implications in quantum information (CPTP maps) and operator theory.

</div>

Linear Preserver Problems (LPP) involve characterizing linear maps $\phi: M_n \to M_n$ that leave certain properties, subsets, or relations invariant. This field, initiated by Frobenius in 1897, reveals the "rigidity" of matrix structures: maps that preserve core properties like the determinant or rank must typically be of a very restricted form, such as similarity or transposition.

---

## 60.1 The Frobenius and Marcus-Moyls Theorems

!!! theorem "Theorem 60.1 (Frobenius Determinant Preserver Theorem, 1897)"
    Let $\phi: M_n(\mathbb{C}) \to M_n(\mathbb{C})$ be a linear map such that $\det(\phi(A)) = \det(A)$ for all $A$. Then there exist invertible matrices $M, N$ with $\det(MN) = 1$ such that:
    $$\phi(A) = MAN \quad \text{or } \quad \phi(A) = MA^T N$$

!!! theorem "Theorem 60.3 (Marcus-Moyls Rank-1 Preserver Theorem)"
    A linear map that preserves the set of rank-1 matrices must be of the form $\phi(A) = MAN$ or $\phi(A) = MA^T N$. Since rank-1 matrices are the "building blocks" of the matrix space, many LPP results reduce to this theorem.

---

## Exercises

1. **[Concept] Define the "Standard Forms" in linear preserver problems.**
   ??? success "Solution"
       Standard forms refer to maps of the type $\phi(A) = MAN$ (similarity-like) and $\phi(A) = MA^T N$ (transposition-like). Most LPP results prove that maps preserving specific algebraic invariants must fall into these two categories.

2. **[Trace Preserver] Show that the transposition map $\phi(A) = A^T$ preserves the trace and the determinant.**
   ??? success "Solution"
       $\operatorname{tr}(A^T) = \sum (A^T)_{ii} = \sum a_{ii} = \operatorname{tr}(A)$. For the determinant, $\det(A^T) = \det(A)$ is an identity from determinant theory. Transposition is the prototype for non-trivial LPPs.

3. **[Spectrum] Prove that a similarity transformation $\phi(A) = PAP^{-1}$ is a spectrum preserver.**
   ??? success "Solution"
       The characteristic polynomial satisfies $\det(\lambda I - PAP^{-1}) = \det(P(\lambda I - A)P^{-1}) = \det(P) \det(\lambda I - A) \det(P)^{-1} = \det(\lambda I - A)$. Identical characteristic polynomials imply identical sets of eigenvalues.

4. **[Rank-1] Why is the preservation of rank-1 matrices considered a crucial step in solving many LPPs?**
   ??? success "Solution"
       Rank-1 matrices are the extreme rays of the PSD cone and are the simplest non-zero elements of the matrix manifold. Many properties (like rank, determinantal degree, and image subspaces) are determined by the behavior of the map on these fundamental elements.

5. **[Invertibility] State the Dieudonné Theorem regarding the preservation of singular matrices.**
   ??? success "Solution"
       A linear bijection that maps the set of singular matrices onto itself must be of the form $\phi(A) = MAN$ or $\phi(A) = MA^T N$. This shows that invertibility is a "rigid" property that restricts the symmetry group of the space.

6. **[Positive Maps] Define a "Positive Map" and provide an example that is not a similarity transformation.**
   ??? success "Solution"
       A map $\phi$ is positive if $A \succeq 0 \implies \phi(A) \succeq 0$. The transposition map $\phi(A) = A^T$ is a positive map but is not a similarity transformation.

7. **[Complete Positivity] Distinguish between positive maps and completely positive (CP) maps.**
   ??? success "Solution"
       A map is completely positive if its extension $\phi \otimes \operatorname{id}_k$ is positive for all $k \ge 1$. Transposition is positive but not CP (it creates negative eigenvalues in entangled spaces), while $\phi(A) = VAV^*$ is CP.

8. **[Choi's Theorem] What is the Kraus representation of a CP map?**
   ??? success "Solution"
       $\phi(A) = \sum V_i A V_i^*$. This sum-of-operators form ensures that the mapping remains positive even when the system is coupled to an arbitrary environment.

9. **[Commutativity] Describe the form of a linear map that preserves commutativity for $n \ge 3$.**
   ??? success "Solution"
       Such maps must satisfy $\phi(A) = \alpha P A P^{-1} + f(A)I$ or $\phi(A) = \alpha P A^T P^{-1} + f(A)I$, where $f$ is a linear functional. The term $f(A)I$ reflects that the identity commutes with all matrices and thus acts as a center.

10. **[Quantum Mechanics] Explain the relevance of LPP results to quantum channels.**
    ??? success "Solution"
        Quantum channels are mathematically modeled as CPTP (Completely Positive Trace Preserving) maps. LPP theory characterizes the allowed physical transformations of density matrices and defines the structure of entanglement-preserving operators.

## Chapter Summary

This chapter examines the rigidity of matrix algebra through Linear Preserver Problems:

1. **Foundational Theorems**: Frobenius and Marcus-Moyls theorems established the standard forms for maps preserving the determinant and rank.
2. **Spectral Rigidity**: Demonstrated that spectrum-preserving maps are essentially restricted to similarity and transposition.
3. **Positivity Classes**: Explored the hierarchy between positive and completely positive maps, providing the mathematical basis for quantum operations.
4. **Structural Symmetry**: Revealed that the preservation of essential algebraic relations (like commutativity or invertibility) forces the map into a highly structured geometric form.
