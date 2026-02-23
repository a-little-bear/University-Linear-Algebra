# Chapter 55B: Lie Algebras and Infinitesimal Generators

<div class="context-flow" markdown>

**Prerequisites**: Matrix Groups (Ch55A) · Matrix Exponential (Ch13) · Matrix Calculus (Ch47A)

**Chapter Outline**: From Continuous Groups to Infinitesimal Transformations → Abstract Definition of a Lie Algebra → The Lie Bracket and the Commutator $[A, B]$ → The Jacobi Identity → The Exponential Map for Matrix Groups → Lie Algebras of Typical Lie Groups ($\mathfrak{gl}, \mathfrak{sl}, \mathfrak{so}, \mathfrak{su}, \mathfrak{sp}$) → Generators and Structure Constants → Applications: Angular Momentum Operators in Quantum Mechanics, Reachability Analysis in Control Theory, and Velocity Mappings in Robot Kinematics

**Extension**: A Lie algebra is the "linearization" of a Lie group at the identity; it simplifies complex group multiplication into linear addition and bracket operations in a vector space. It proves that all local information of a continuous symmetry is contained within its tangent space, serving as the core link between symmetry and dynamics.

</div>

An extremely powerful method for studying continuous transformation groups (Lie groups) is to examine their "infinitesimal" behavior. A **Lie Algebra** is the mathematical manifestation of this microscopic local structure. Through the matrix exponential map, we can amplify linear components in the Lie algebra into global transformations within the Lie group. This chapter introduces how to use the commutator as an algebraic tool to move the study of symmetry from curved manifolds to flat tangent spaces.

---

## 55B.1 Lie Algebras and the Lie Bracket

!!! definition "Definition 55B.1 (Lie Algebra)"
    A vector space $\mathfrak{g}$ equipped with a binary operation $[\cdot, \cdot]$ (the Lie bracket) is a **Lie Algebra** if it satisfies:
    1.  **Bilinearity**.
    2.  **Antisymmetry**: $[X, Y] = -[Y, X]$.
    3.  **Jacobi Identity**: $[X, [Y, Z]] + [Y, [Z, X]] + [Z, [X, Y]] = 0$.
    For matrix Lie algebras, the Lie bracket is the **Commutator**: $[A, B] = AB - BA$.

---

## 55B.2 The Exponential Map and Tangent Space

!!! theorem "Theorem 55B.1 (Lie Algebra as Tangent Space)"
    The tangent space of a Lie group $G$ at the identity $I$ is its associated Lie algebra $\mathfrak{g}$. For any $X \in \mathfrak{g}$, the curve $\gamma(t) = e^{tX}$ lies entirely within $G$.
    **Intuition**: $X$ is the "velocity" or "generator" that produces the group's evolution.

---

## 55B.3 Typical Matrix Lie Algebras

!!! note "Correspondences"
    1.  **$\mathfrak{sl}(n)$**: Trace-zero matrices (corresponding to $SL(n)$ with det=1).
    2.  **$\mathfrak{so}(n)$**: Skew-symmetric matrices $X^T = -X$ (corresponding to rotation groups).
    3.  **$\mathfrak{su}(n)$**: Skew-Hermitian trace-zero matrices $X^* = -X, \operatorname{tr}(X)=0$.

---

## Exercises

**1. [Basics] Verify if the matrix commutator $[A, B] = AB - BA$ satisfies antisymmetry.**

??? success "Solution"
    **Verification:**
    $[B, A] = BA - AB = -(AB - BA) = -[A, B]$.
    **Conclusion**: Yes. This reflects that at the infinitesimal scale, reversing the order of operations results in a sign flip.

**2. [Jacobi] Prove that for any matrices $A, B, C$, the identity $[A, [B, C]] + [B, [C, A]] + [C, [A, B]] = O$ holds.**

??? success "Solution"
    **Proof:**
    Expand each term:
    1. $[A, BC-CB] = ABC - ACB - BCA + CBA$.
    2. Similarly write out the other two nested brackets.
    3. Summing all 12 triple products reveals that every term (e.g., $ABC$) cancels with its negative counterpart from another part of the sum.
    **Conclusion**: The Jacobi identity always holds for matrix commutators.

**3. [Generator] What is the basis for the Lie algebra $\mathfrak{so}(2)$?**

??? success "Solution"
    **Analysis:**
    1. $\mathfrak{so}(2)$ consists of $2 \times 2$ skew-symmetric matrices.
    2. They take the form $\begin{pmatrix} 0 & -\theta \\ \theta & 0 \end{pmatrix}$.
    **Basis**: $X = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}$.
    **Verification**: $e^{\theta X} = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}$, which is an element of the rotation group $SO(2)$.

**4. [Trace] Prove: If $e^X \in SL(n)$, then $\operatorname{tr}(X) = 0$.**

??? success "Solution"
    **Proof:**
    1. Use Jacobi's identity for determinants: $\det(e^X) = e^{\operatorname{tr}(X)}$.
    2. $SL(n)$ requires the determinant to be 1.
    3. $e^{\operatorname{tr}(X)} = 1 \implies \operatorname{tr}(X) = 0$.
    **Conclusion**: The Lie algebra of the special linear group consists of trace-zero matrices.

**5. [Calculation] Compute $[ \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}, \begin{pmatrix} 0 & 0 \\ 1 & 0 \end{pmatrix} ]$.**

??? success "Solution"
    **Steps:**
    1. $AB = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix} \begin{pmatrix} 0 & 0 \\ 1 & 0 \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$.
    2. $BA = \begin{pmatrix} 0 & 0 \\ 1 & 0 \end{pmatrix} \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix} = \begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}$.
    3. $[A, B] = AB - BA = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}$.
    **Conclusion**: These two nilpotent elements generate a diagonal (weight) operator through the bracket, a standard feature of the $SL(2)$ Lie algebra structure.

**6. [Constants] What are the "Structure Constants" of a Lie algebra?**

??? success "Solution"
    **Definition:**
    If $\{T_i\}$ is a basis, the bracket of any two basis elements is a linear combination of the basis: $[T_i, T_j] = \sum_k f_{ij}^k T_k$. The coefficients $f_{ij}^k$ are the **structure constants**. They contain all the local geometric information of the Lie group.

**7. [Unitary] Why do physicists often insert an '$i$' before Lie algebra generators?**

??? success "Solution"
    In quantum mechanics, physical observables must be **Hermitian** (real eigenvalues). However, the Lie algebra $\mathfrak{su}(n)$ consists of **skew-Hermitian** matrices. By defining $X = iH$, the skew-Hermitian operator $X$ is mapped to a Hermitian operator $H$, aligning the generators with physical concepts like energy or momentum.

**8. [Properties] Prove: If $X \in \mathfrak{so}(n)$, its eigenvalues are purely imaginary or zero.**

??? success "Solution"
    **Proof:**
    1. Skew-symmetry implies $X^T = -X$.
    2. This means $iX$ is a Hermitian matrix.
    3. The eigenvalues $\mu$ of a Hermitian matrix are real.
    4. The eigenvalues $\lambda$ of $X$ are $\mu/i = -i\mu$.
    **Conclusion**: Eigenvalues lie on the imaginary axis, corresponding to oscillatory behavior without decay or growth.

**9. [Application] Briefly state the meaning of the Lie bracket in robot reachability.**

??? success "Solution"
    If a robot has only two control directions $X$ and $Y$, by rapidly switching between them, it can effectively move in the $[X, Y]$ direction. If the Lie brackets of all available controls span the entire space dimension (Chow's Theorem), the robot is **completely controllable**.

**10. [Limit] Write the Taylor approximation of $e^{tX}$ as $t \to 0$.**

??? success "Solution"
    **Conclusion: $e^{tX} \approx I + tX$.**
    This proves that the Lie algebra element $X$ is indeed the "first derivative" or "tangent vector" of the group evolution.

## Chapter Summary

Lie algebras are the infinitesimal microscope for studying the evolution of symmetry:

1.  **Power of Linearization**: They compress complex non-linear group structures into vector operations in flat spaces, greatly simplifying classification and spectral analysis.
2.  **Essence of Dynamics**: Via the exponential map, Lie algebras reveal that all continuous transformations evolve from a set of "generators" through integration, establishing the mapping between conservation laws and symmetry.
3.  **Code of Structure**: Lie brackets and structure constants form the algebraic cipher for describing physical interactions. From angular momentum addition to the Standard Model gauge fields, Lie algebras are the singular pedestal for modern Grand Unified Theories.
