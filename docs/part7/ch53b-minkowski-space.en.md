# Chapter 53B: Minkowski Space and Lorentz Groups

<div class="context-flow" markdown>

**Prerequisites**: Inner Product Spaces (Ch08) · Symplectic Matrices (Ch53A) · Basics of Special Relativity

**Chapter Outline**: From Euclidean to Minkowski Space → Indefinite Inner Product and the Metric Matrix $\eta$ → Matrix Definition of Lorentz Transformations → The Lorentz Group $O(1, 3)$ and its Branches → Proper Lorentz Transformations, Rotations, and Boosts → Classification of Vectors: Time-like, Space-like, and Null (Light-like) → Causal Structure and the Light Cone → Polar Decomposition of Lorentz Matrices → Applications: Spacetime Coordinate Transforms, Doppler Effect, and Momentum Tensors in Particle Physics

**Extension**: Minkowski space is the algebraic vehicle for the "unification of spacetime" in physics; it weaves time and space into a four-dimensional manifold with a negative-signature metric, proving that the invariance of physical laws is essentially the manifestation of the Lorentz group at the operator level. It is the only mathematical foundation for understanding relativity.

</div>

In everyday experience, distance is always positive. However, in special relativity, the interval between time and space requires a specific type of combination. **Minkowski Space** introduces a metric with a negative sign, allowing "distances" to be zero or even negative. The matrices that preserve this interval form the **Lorentz Group**. This chapter introduces the linear algebraic language describing the four-dimensional skeleton of the universe.

---

## 53B.1 Minkowski Metric and Indefinite Inner Product

!!! definition "Definition 53B.1 (Minkowski Metric $\eta$)"
    In four-dimensional spacetime, the metric matrix is typically taken as:
    $$\eta = \operatorname{diag}(-1, 1, 1, 1)$$
    The **Minkowski Inner Product** of two four-vectors $u, v$ is defined as:
    $$\langle u, v \rangle_\eta = u^T \eta v = -u_0 v_0 + u_1 v_1 + u_2 v_2 + u_3 v_3$$

!!! note "Classification of Vectors"
    1.  **Time-like**: $\langle v, v \rangle_\eta < 0$.
    2.  **Space-like**: $\langle v, v \rangle_\eta > 0$.
    3.  **Null / Light-like**: $\langle v, v \rangle_\eta = 0$ (the trajectory of a light signal).

---

## 53B.2 Matrix Representation of Lorentz Transformations

!!! definition "Definition 53B.2 (Lorentz Matrix)"
    A square matrix $\Lambda \in M_4$ is a **Lorentz Transformation** if it preserves the Minkowski inner product:
    $$\Lambda^T \eta \Lambda = \eta$$
    The group of all such matrices is denoted by $O(1, 3)$.

---

## 53B.3 Boosts and Rotations

!!! technique "Lorentz Boost"
    A boost along the $x$-axis with velocity $v$ is represented by the matrix:
    $$\Lambda = \begin{pmatrix} \gamma & -\beta\gamma & 0 & 0 \\ -\beta\gamma & \gamma & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}$$
    where $\beta = v/c$ and $\gamma = 1/\sqrt{1-\beta^2}$. This corresponds to a "rotation" (hyperbolic rotation) between time and space coordinates.

---

## Exercises

**1. [Basics] Calculate the norm squared of vector $v = (1, 1, 0, 0)^T$ under the metric $\eta = \operatorname{diag}(-1, 1, 1, 1)$.**

??? success "Solution"
    **Steps:**
    1. $\langle v, v \rangle_\eta = v^T \eta v = -v_0^2 + v_1^2 + v_2^2 + v_3^2$.
    2. Substitute components: $-1^2 + 1^2 + 0^2 + 0^2 = -1 + 1 = 0$.
    **Conclusion**: The norm squared is 0. This is a **Null vector** (light-like), representing a path traveling at the speed of light.

**2. [Lorentz Group] Prove that the determinant of a Lorentz matrix must be $\pm 1$.**

??? success "Solution"
    **Proof:**
    1. Given $\Lambda^T \eta \Lambda = \eta$.
    2. Take the determinant of both sides: $\det(\Lambda^T) \det(\eta) \det(\Lambda) = \det(\eta)$.
    3. Since $\det(\eta) = -1 \neq 0$, we divide to get $(\det \Lambda)^2 = 1$.
    **Conclusion**: $\det \Lambda = \pm 1$. Matrices with $\det = 1$ are called proper Lorentz transformations.

**3. [Calculation] Verify that the $2 \times 2$ boost matrix $\begin{pmatrix} \cosh \phi & -\sinh \phi \\ -\sinh \phi & \cosh \phi \end{pmatrix}$ satisfies the Lorentz definition.**

??? success "Solution"
    **Steps:**
    1. Let $\eta = \operatorname{diag}(-1, 1)$.
    2. Multiply: $\Lambda^T \eta \Lambda = \begin{pmatrix} \cosh \phi & -\sinh \phi \\ -\sinh \phi & \cosh \phi \end{pmatrix} \begin{pmatrix} -1 & 0 \\ 0 & 1 \end{pmatrix} \begin{pmatrix} \cosh \phi & -\sinh \phi \\ -\sinh \phi & \cosh \phi \end{pmatrix}$.
    3. The top-left entry is $-\cosh^2 \phi + \sinh^2 \phi = -1$.
    4. The bottom-right entry is $-\sinh^2 \phi + \cosh^2 \phi = 1$.
    5. Cross terms: $\cosh\sinh - \sinh\cosh = 0$.
    **Conclusion**: The result is exactly $\eta$. The parameter $\phi$ is known as **rapidity**.

**4. [Causality] Why is the angle between time-like vectors undefined, but they still have "future" and "past" designations?**

??? success "Solution"
    **Algebraic Reason:**
    For a time-like vector $v$ ($v_0^2 > v_1^2+v_2^2+v_3^2$), the sign of the $v_0$ component is invariant under continuous Lorentz transformations. This splits spacetime into two disconnected regions: $v_0 > 0$ (the future light cone) and $v_0 < 0$ (the past light cone). This forms the mathematical basis for causality in physics.

**5. [Invariants] Does the rest mass of an object (the norm of the momentum four-vector) change under a Lorentz transformation?**

??? success "Solution"
    **Conclusion: No.**
    **Reasoning**: The squared norm of the momentum four-vector $P^T \eta P = -E^2 + p^2$ corresponds to $-m^2 c^2$. Since Lorentz transformations preserve the Minkowski inner product, the rest mass is a **Lorentz Invariant** agreed upon by all observers.

**6. [Property] Prove that the inverse of a Lorentz transformation is $\Lambda^{-1} = \eta \Lambda^T \eta$.**

??? success "Solution"
    **Proof:**
    1. From $\Lambda^T \eta \Lambda = \eta$, right-multiply by $\Lambda^{-1}$.
    2. $\Lambda^T \eta = \eta \Lambda^{-1}$.
    3. Since $\eta^{-1} = \eta$, left-multiply by $\eta$.
    4. $\eta \Lambda^T \eta = \Lambda^{-1}$.

**7. [Classification] If a vector has only spatial components $(0, 1, 0, 0)$ in a given frame, what type of vector is it?**

??? success "Solution"
    **Determination:**
    $\langle v, v \rangle = -0^2 + 1^2 + 0^2 + 0^2 = 1 > 0$.
    **Conclusion**: It is a **Space-like vector**. Such vectors represent points that cannot have a causal connection in any reference frame.

**8. [Application] What is a "Hyperbolic Rotation" in spacetime?**

??? success "Solution"
    Standard rotations preserve $x^2+y^2$ (circles). Lorentz boosts preserve $-t^2+x^2$ (hyperbolas). Thus, a boost is mathematically equivalent to a **hyperbolic rotation** in the spacetime plane, with the rotation "angle" being the rapidity $\phi$.

**9. [Spectral] What are the characteristics of the eigenvalues of a Lorentz matrix?**

??? success "Solution"
    **Conclusion:**
    1. If $\lambda$ is an eigenvalue, $1/\lambda$ must also be an eigenvalue (similar to symplectic matrices).
    2. For a boost, the eigenvalues are $e^\phi$ and $e^{-\phi}$. These correspond to the Doppler shift factors.

**10. [Application] Briefly state the role of the Lorentz group in particle collider data analysis.**

??? success "Solution"
    In collision experiments, detectors are in the lab frame, while physics is simplest in the center-of-mass frame. Using Lorentz matrices to transform measured four-momenta between frames is the only algebraic path to calculating invariant masses and confirming the discovery of new particles (like the Higgs boson).

## Chapter Summary

Minkowski space and Lorentz groups reconstruct our understanding of reality:

1.  **Absoluteness of Interval**: By introducing an indefinite inner product, the spacetime interval replaces Euclidean distance as the only physical entity independent of the observer's state.
2.  **Unity of Symmetry**: The Lorentz group unifies spatial rotations and velocity boosts into a single algebraic framework, proving that time is just a special direction with a negative metric in the spacetime manifold.
3.  **Boundaries of Causality**: The null vectors and light-cone structure establish the ultimate limit for information propagation, proving that the topological causality of the universe is a direct consequence of the Minkowski metric signature.
