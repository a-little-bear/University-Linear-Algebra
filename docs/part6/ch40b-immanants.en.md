# Chapter 40B: Immanants

<div class="context-flow" markdown>

**Prerequisites**: Determinants (Ch03) · Permanents (Ch40A) · Group Representation Theory (Ch55)

**Chapter Outline**: Group Theoretic Motivation for Immanants (as generalizations of Determinants and Permanents) → Definition via Characters of the Symmetric Group $S_n$ → Types of Immanants corresponding to different characters → Spectral perspectives on Determinants and Permanents → Schur’s Inequality for Immanants → Key Theorem: The Schur Power Series Identity → Applications: Description of State Spaces for Identical Particles in Quantum Mechanics, and S-functions in Algebraic Combinatorics

**Extension**: Immanants are the projection of the representation theory of the symmetric group onto matrix functions; they prove that determinants and permanents are not isolated constructs, but correspond to the two most extreme one-dimensional representations (antisymmetric and trivial), with a vast family of functions characterized by Young Tableaux lying in between.

</div>

In previous sections, we saw the determinant with alternating signs and the permanent with all positive signs. The **Immanant** is their ultimate unification. By utilizing the **characters** of the symmetric group $S_n$, we can define a series of matrix functions spanning the spectrum between the determinant and the permanent. Immanants are not only natural products of group representation theory but also universal algebraic tools for describing quantum many-body systems and combinatorial symmetries.

---

## 40B.1 Definition and Group Background

!!! definition "Definition 40B.1 (Immanant)"
    Let $\chi$ be a complex character of the symmetric group $S_n$. The **immanant** of a square matrix $A$ with respect to $\chi$ is defined as:
    $$d_\chi(A) = \sum_{\sigma \in S_n} \chi(\sigma) a_{1,\sigma(1)} a_{2,\sigma(2)} \cdots a_{n,\sigma(n)}$$

!!! note "Special Cases"
    1.  If $\chi$ is the **alternating character** $\epsilon(\sigma) = \operatorname{sgn}(\sigma)$, then $d_\epsilon(A) = \det(A)$.
    2.  If $\chi$ is the **trivial character** $1(\sigma) = 1$, then $d_1(A) = \operatorname{perm}(A)$.
    3.  For other irreducible characters, the resulting functions lie between these two extremes.

---

## 40B.2 Schur’s Inequality

!!! theorem "Theorem 40B.1 (Schur’s Inequality)"
    For any positive semi-definite matrix $A$ and any irreducible character $\chi$ of $S_n$:
    $$\det(A) \le \frac{d_\chi(A)}{\chi(id)} \le \operatorname{perm}(A)$$
    where $\chi(id)$ is the dimension of the representation.
    **Significance**: This chain of inequalities proves that the determinant is the lower bound for all immanants, while the permanent is the upper bound.

---

## 40B.3 Application: Quantum Mechanics

!!! technique "Identical Particle States"
    In quantum mechanics, the wavefunctions of fermions are described by determinants (satisfying the Pauli Exclusion Principle), while the wavefunctions of bosons are described by permanents. Immanants appear in describing quasiparticle systems with more complex exchange statistics (such as fractional statistics or anyons).

---

## Exercises

**1. [Basics] For $2 \times 2$ matrices, how many irreducible characters does $S_2$ have? What are the corresponding immanants?**

??? success "Solution"
    **Analysis:**
    $S_2$ has only two elements: the identity $(1,2)$ and the swap $(2,1)$. It has two irreducible characters:
    1. Trivial character $\chi_1$: $\chi_1(id)=1, \chi_1(\text{swap})=1$. This corresponds to the **permanent**.
    2. Alternating character $\chi_2$: $\chi_2(id)=1, \chi_2(\text{swap})=-1$. This corresponds to the **determinant**.
    **Conclusion**: For order 2, the immanants are restricted to the determinant and the permanent.

**2. [Calculation] For $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$, calculate all immanants of $S_2$.**

??? success "Solution"
    **Calculation:**
    1. $\operatorname{perm}(A) = 1\cdot 4 + 2\cdot 3 = 10$.
    2. $\det(A) = 1\cdot 4 - 2\cdot 3 = -2$.
    **Conclusion**: These are the two immanant values for the matrix.

**3. [Dimension] If $\chi$ is an irreducible character of $S_n$ with representation dimension $f$, what is $d_\chi(I)$?**

??? success "Solution"
    **Derivation:**
    1. For the identity matrix $I$, the product term is non-zero only for the identity permutation $\sigma=id$, where the value is 1.
    2. The definition simplifies to $d_\chi(I) = \chi(id) \cdot 1$.
    3. Character theory states $\chi(id)$ equals the dimension $f$.
    **Conclusion**: $d_\chi(I) = f$.

**4. [Schur Inequality] Verify for the PD matrix $A = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$ that the determinant $\le$ the permanent.**

??? success "Solution"
    **Calculation:**
    1. $\det(A) = 4 - 1 = 3$.
    2. $\operatorname{perm}(A) = 4 + 1 = 5$.
    **Conclusion**: $3 \le 5$, verification successful.

**5. [Transpose] Prove: Immanants are invariant under matrix transposition.**

??? success "Solution"
    **Proof:**
    1. Both $\operatorname{perm}$ and $\det$ satisfy $\operatorname{perm}(A^T) = \operatorname{perm}(A)$.
    2. For a general immanant, the sum corresponds to the same product terms. Since characters are class functions, they take the same value on conjugate elements (and the permutation for $A^T$ is the inverse of the permutation for $A$, which belongs to the same class).
    **Conclusion**: The property holds for all immanants.

**6. [Normalization] Why divide by $\chi(id)$ in Schur’s Inequality?**

??? success "Solution"
    **Reason:**
    To eliminate the scaling effect of the representation dimension. After dividing by $\chi(id)$, the value of the function at the identity matrix $I$ is normalized to 1, allowing for a direct comparison between functions generated by characters of different dimensions.

**7. [Multilinearity] Are immanants linear with respect to each row of the matrix?**

??? success "Solution"
    **Yes.**
    **Reason**: From the definition, each product term contains exactly one element from row $i$. Thus, immanants inherit the multilinear property of determinants and permanents.

**8. [Young Tableaux] What is the immanant generated by the character corresponding to the Young diagram $(n-1, 1)$?**

??? success "Solution"
    It is often referred to as the immanant of the "standard representation" of the symmetric group. It is the simplest branch beyond the trivial and alternating representations.

**9. [Hardness] How difficult is it to compute an irreducible immanant in general?**

??? success "Solution"
    **Conclusion**: Typically **#P-complete**.
    **Reason**: Since the permanent is a special case of an immanant, computing a general immanant is at least as hard as computing a permanent, except for the rare case of the determinant.

**10. [Combinatorics] How do immanants relate to S-functions (Schur functions)?**

??? success "Solution"
    **Connection:**
    If the entries of matrix $A$ are treated as variables of symmetric polynomials, immanants serve as the algebraic basis for constructing S-functions. S-functions are the ultimate family of functions for describing partitions and representation dimensions, and immanants provide their matrix characterizations (generalizing the Jacobi-Trudi identity).

## Chapter Summary

The immanant is the master of algebraic symmetry:

1.  **Unification of the Spectrum**: It proves that the determinant and the permanent are not isolated phenomena but endpoints of the representation spectrum of the symmetric group, establishing a deep link between matrix functions and character theory.
2.  **Ladder of Energy**: Schur’s inequality reveals how the degree of sign alternation monotonically affects the global "measure" of an operator, providing a unified algebraic perspective on quantum interference and statistical exclusion.
3.  **Projection of Groups**: By mapping the irreducible representations of the symmetric group into matrix space, immanants become the essential language for parsing combinatorial structures and particle statistical properties.
