# Chapter 51: Quaternion Matrices

<div class="context-flow" markdown>

**Prerequisites**: Clifford Algebra (Ch50) · Matrix Algebra (Ch02) · Eigenvalues (Ch06)

**Chapter Outline**: From Complex Numbers to Quaternions → Non-commutativity of the Quaternion Algebra $\mathbb{H}$ → Definitions and Operational Laws of Quaternion Matrices → The Divergence between Left and Right Eigenvalues → The Core Challenge: Non-commutative Determinants (Dieudonné Determinant) → Complex Representation of Quaternion Matrices → Spectral Structure Theorem → Unitary Quaternion Matrices ($Sp(n)$) → Applications: Robotic Pose Control, Quaternionic Quantum Mechanics, and Color Image Representation in Signal Processing

**Extension**: Quaternion matrices are an attempt to place linear algebra on a "non-commutative field"; by introducing three imaginary units $i, j, k$, they perfectly encapsulate 3D rotation information—the ultimate mathematical framework for processing data with strong rotational correlation.

</div>

When we replace the entries of a matrix with **Quaternions**, we enter a non-commutative algebraic world. Since $ij \neq ji$, traditional determinant definitions and eigenvalue theories undergo drastic changes. **Quaternion Matrices** are not only theoretically challenging in pure mathematics but also demonstrate incredible elegance in describing robotic joint movements, spacecraft attitude, and color image processing. This chapter introduces how to re-establish order in this non-commutative domain.

---

## 51.1 Foundations of Quaternion Algebra $\mathbb{H}$

!!! definition "Definition 51.1 (Quaternion)"
    A quaternion $q \in \mathbb{H}$ is of the form $q = a + bi + cj + dk$, satisfying Hamilton's rules:
    $$i^2 = j^2 = k^2 = ijk = -1$$
    **Key Conflict**: $ij = k$, but $ji = -k$. Multiplication is not commutative.

---

## 51.2 Eigenvalues: Left vs. Right

!!! definition "Definition 51.2 (Left and Right Eigenvalues)"
    For a square quaternion matrix $A$:
    1.  **Left Eigenvalue**: $A\mathbf{x} = \lambda \mathbf{x}$.
    2.  **Right Eigenvalue**: $A\mathbf{x} = \mathbf{x} \lambda$.
    **Difference**: Due to non-commutativity, right eigenvalues are far more research-worthy and exhibit "spectral" properties similar to complex eigenvalues.

---

## 51.3 Complex Representation

!!! technique "Technique: Mapping to Complex Matrices"
    Every $n \times n$ quaternion matrix $A$ can be mapped to a $2n \times 2n$ complex matrix $\mathcal{X}(A)$:
    $$q = \alpha + j\beta \to \begin{pmatrix} \alpha & \beta \\ -\bar{\beta} & \bar{\alpha} \end{pmatrix}$$
    This allows us to leverage classical complex linear algebra tools to solve quaternionic problems.

---

## Exercises

**1. [Basics] Calculate the quaternion product $(1+i)j$.**

??? success "Solution"
    **Steps:**
    1. Apply distributivity: $1 \cdot j + i \cdot j$.
    2. Apply Hamilton's rules: $ij = k$.
    **Conclusion**: The result is $j + k$. Note: $j(1+i) = j + ji = j - k$, which is different.

**2. [Complex Rep] Represent $q = i$ as a $2 \times 2$ complex matrix.**

??? success "Solution"
    **Construction:**
    1. Write as $q = \alpha + j\beta$: here $\alpha = i, \beta = 0$.
    2. Apply the formula: $\begin{pmatrix} i & 0 \\ 0 & -i \end{pmatrix}$.
    **Verification**: The square of this matrix is $-I$, consistent with $i^2 = -1$.

**3. [Determinant] Why can't the Leibniz formula be used directly for quaternion determinants?**

??? success "Solution"
    **Reasoning:**
    The Leibniz formula involves products of elements (e.g., $a_{1,\sigma(1)} a_{2,\sigma(2)} \cdots$). Since quaternion multiplication is non-commutative, changing the order of terms changes the value of the product. Thus, the traditional determinant loses its uniqueness. Mathematicians introduced the **Dieudonné determinant** to solve this.

**4. [Right Eigenvalues] Prove: If $\lambda$ is a right eigenvalue of $A$, then $q^{-1}\lambda q$ is also an eigenvalue for any $q \in \mathbb{H}, q \neq 0$.**

??? success "Solution"
    **Proof:**
    1. Given $A\mathbf{x} = \mathbf{x}\lambda$.
    2. Consider vector $\mathbf{y} = \mathbf{x}q$.
    3. $A\mathbf{y} = A(\mathbf{x}q) = (A\mathbf{x})q = (\mathbf{x}\lambda)q = \mathbf{x}q(q^{-1}\lambda q) = \mathbf{y}(q^{-1}\lambda q)$.
    **Conclusion**: Eigenvalues form "conjugacy classes." This shows that every non-real eigenvalue of a quaternion matrix actually corresponds to a sphere in the 3D imaginary space.

**5. [Transpose] Define the conjugate transpose $A^*$ for quaternion matrices.**

??? success "Solution"
    **Definition:**
    $(A^*)_{ij} = \overline{a_{ji}}$.
    Where the quaternion conjugate $\bar{q} = a - bi - cj - dk$. Note: $\overline{pq} = \bar{q}\bar{p}$.

**6. [Unitary] What is a unitary quaternion matrix (Symplectic matrix)?**

??? success "Solution"
    **Definition:**
    A quaternion matrix satisfying $A^* A = I$. These matrices form the **Symplectic group** $Sp(n)$, describing symmetries that preserve a specific skew-symmetric structure in physics.

**7. [Calculation] Find the conjugate transpose of $\begin{pmatrix} i & j \\ k & 1 \end{pmatrix}$.**

??? success "Solution"
    **Steps:**
    1. Transpose and conjugate: $\begin{pmatrix} \bar{i} & \bar{k} \\ \bar{j} & \bar{1} \end{pmatrix}$.
    2. Substitute conjugates: $\begin{pmatrix} -i & -k \\ -j & 1 \end{pmatrix}$.

**8. [Real Spectra] Prove: Right eigenvalues of a Hermitian quaternion matrix ($A^*=A$) must be real.**

??? success "Solution"
    **Proof Sketch:**
    Use the complex representation $\mathcal{X}(A)$. If $A$ is Hermitian, $\mathcal{X}(A)$ is a complex Hermitian matrix. Since complex Hermitian matrices have real eigenvalues, and quaternion right eigenvalues are contained within the spectrum of $\mathcal{X}(A)$, the result follows.

**9. [Inversion] What is the condition for a quaternion matrix to be invertible?**

??? success "Solution"
    **Conclusion:**
    Its complex representation $\mathcal{X}(A)$ must be non-singular. In the quaternion domain, this corresponds to its Dieudonné determinant being non-zero.

**10. [Application] Briefly state the advantage of quaternion matrices in color image processing.**

??? success "Solution"
    Color images have R, G, B channels. Using pure quaternions $q = ri + gj + bk$ to represent a pixel treats the image as a quaternion matrix.
    **Advantage**: When performing rotations or filtering, quaternion algebra automatically preserves the correlation between R, G, and B components, avoiding color shifts that occur when processing channels independently.

## Chapter Summary

Quaternion matrices are the pinnacle of non-commutative algebra applications:

1.  **Fusion of Dimensions**: Quaternions are not just numbers but geometric entities with internal rotation; quaternion matrices weave 3D rotation perfectly into the framework of linear operators.
2.  **Spherical Spectra**: The phenomenon of right eigenvalue conjugacy classes reveals that the spectrum of a non-commutative operator is no longer a set of isolated points but symmetric orbits, greatly enriching the geometry of operator theory.
3.  **Representation Switching**: The complex representation technique bridges non-commutative algebra and classical complex analysis, proving that even in environments without commutativity, we can regain control through dimension doubling.
