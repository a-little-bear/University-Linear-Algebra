# Chapter 53B: Minkowski Space and Lorentz Transformations

<div class="context-flow" markdown>

**Prerequisites**: Inner Product Spaces (Ch08) · Orthogonality (Ch07) · Matrix Groups (Ch55)

**Chapter Outline**: The Minkowski Metric $\eta = \operatorname{diag}(-1, 1, 1, 1)$ → Definition of Indefinite Inner Product Spaces → Lorentz Transformations $L^T \eta L = \eta$ → Classification of Vectors: Time-like, Space-like, and Light-like (Null) → Structure of the Light Cone → The Lorentz Group $O(1, 3)$ and its Components → Boosts and Rotations → Applications: Time Dilation, Length Contraction, and Invariant Spacetime Intervals in Special Relativity → Connection to Hyperbolic Geometry

**Extension**: Minkowski space provides the local geometric background for General Relativity and Quantum Field Theory; it demonstrates how linear algebra, by modifying the definition of "length," can perfectly describe the four-dimensional structure of the universe, serving as the essential link between geometry and modern physics.

</div>

In Euclidean geometry, the square of the distance is the sum of the squares of the components. However, in Einstein's Special Relativity, time and space are intertwined into a four-dimensional space with an "indefinite metric"—**Minkowski Space**. Transformations in this space preserve the square of the spacetime interval rather than simple spatial distance. This chapter utilizes the language of linear algebra to reinterpret the core geometric structures of relativity.

---

## 53B.1 The Minkowski Metric and Indefinite Inner Product

!!! definition "Definition 53B.1 (Minkowski Metric)"
    In the four-dimensional real vector space $\mathbb{R}^4$, define the **Minkowski metric tensor** $\eta$:
    $$\eta = \begin{pmatrix} -1 & 0 & 0 & 0 \ 0 & 1 & 0 & 0 \ 0 & 0 & 1 & 0 \ 0 & 0 & 0 & 1 \end{pmatrix}$$
    The **Minkowski inner product** of two four-vectors $u, v$ is defined as:
    $$\langle u, v angle_\eta = u^T \eta v = -u_0 v_0 + u_1 v_1 + u_2 v_2 + u_3 v_3$$

---

## 53B.2 Vector Classification and the Light Cone

!!! definition "Definition 53B.2 (Causal Classification of Vectors)"
    Based on the sign of the squared length $s^2 = \langle v, v angle_\eta$, vectors are classified as:
    1.  **Time-like**: $s^2 < 0$. Corresponds to possible worldlines of massive particles.
    2.  **Space-like**: $s^2 > 0$. Corresponds to pairs of events that cannot be causally linked.
    3.  **Light-like (Null)**: $s^2 = 0$. Corresponds to the trajectories of photons.

!!! technique "Geometry: The Light Cone"
    The set of all light-like vectors ($s^2 = 0$) forms a double cone in spacetime known as the **Light Cone**. it defines the geometric boundary of causality.

---

## 53B.3 Lorentz Transformations and the Lorentz Group

!!! definition "Definition 53B.3 (Lorentz Transformation)"
    A linear transformation $L$ is a **Lorentz Transformation** if it preserves the Minkowski inner product:
    $$L^T \eta L = \eta$$
    The group of all such transformations is the **Lorentz Group** $O(1, 3)$.

!!! theorem "Theorem 53B.1 (Components of Lorentz Transformations)"
    Every proper Lorentz transformation can be decomposed into:
    1.  **Spatial Rotations**: Rotations among the three spatial dimensions.
    2.  **Boosts**: "Hyperbolic rotations" between time and spatial dimensions, corresponding to relative motion between inertial frames.

---

## Exercises

1.  **[Basics] Calculate the Minkowski squared length of the vector $v = (1, 1, 0, 0)$. What type of vector is it?**
    ??? success "Solution"
        $s^2 = -1^2 + 1^2 + 0^2 + 0^2 = 0$. It is a **light-like (null) vector**.

2.  **[Determinant] Prove: If $L$ is a Lorentz transformation, then $\det(L) = \pm 1$.**
    ??? success "Solution"
        From $L^T \eta L = \eta$, taking the determinant yields $\det(L)^2 \det(\eta) = \det(\eta)$. Since $\det(\eta) = -1 
eq 0$, it follows that $\det(L)^2 = 1 \implies \det(L) = \pm 1$.

3.  **[Boost] Write the Lorentz boost matrix for a velocity $v$ along the $x$-axis (set $c=1$).**
    ??? success "Solution"
        Let $\gamma = 1/\sqrt{1-v^2}$. The matrix is $\begin{pmatrix} \gamma & -\gamma v & 0 & 0 \ -\gamma v & \gamma & 0 & 0 \ 0 & 0 & 1 & 0 \ 0 & 0 & 0 & 1 \end{pmatrix}$.

4.  **[Hyperbolic] Prove the determinant of the boost matrix is 1.**
    ??? success "Solution"
        $\gamma^2 - (-\gamma v)^2 = \gamma^2(1-v^2) = 1$.

5.  **[Classification] Is the vector $(2, 1, 1, 1)$ time-like or space-like?**
    ??? success "Solution"
        $s^2 = -2^2 + 1^2 + 1^2 + 1^2 = -4 + 3 = -1$. It is a **time-like vector**.

6.  **[Invariance] Prove that Lorentz transformations preserve the "speed of light."**
    ??? success "Solution"
        If $v$ is a light-like vector (representing light speed), then $\langle v, v angle_\eta = 0$. After transformation, $\langle Lv, Lv angle_\eta = \langle v, v angle_\eta = 0$. Thus $Lv$ remains light-like.

7.  **[Topology] How many connected components does the Lorentz group $O(1, 3)$ have?**
    ??? success "Solution"
        There are 4 components, categorized by the sign of $\det(L)$ and the sign of $L_{00}$.

8.  **[Trace] How is the trace of a boost matrix related to velocity $v$?**
    ??? success "Solution"
        $\operatorname{tr}(L) = 2\gamma + 2 = 2(1 + 1/\sqrt{1-v^2})$.

9.  **[Interval] What is the Spacetime Interval?**
    ??? success "Solution"
        It is the Minkowski inner product $\Delta s^2 = -\Delta t^2 + \Delta x^2 + \Delta y^2 + \Delta z^2$. It is an objective quantity upon which all observers agree.

10. **[Geometry] What surface is formed by the set of all time-like unit vectors in Minkowski space?**

   ??? success "Solution"
        The set of points satisfying $-t^2 + x^2 + y^2 + z^2 = -1$. This is a **hyperboloid of two sheets**, explaining the intrinsic link between relativity and hyperbolic geometry.

## Chapter Summary

Minkowski space reshapes the causality of the physical world by modifying the algebraic structure of distance:

1.  **Fusion of Dimensions**: The Minkowski metric treats time as a fourth dimension with a negative sign, proving the essential unity of spacetime and establishing the principle of covariance for physical laws.
2.  **Boundaries of Causality**: The introduction of the light cone structure, via the sign of vector lengths, sets absolute physical limits (the speed of light) for information transmission in the universe.
3.  **Essence of Transformation**: Lorentz transformations prove that "time dilation" and "length contraction" are merely hyperbolic rotations within Minkowski space—manifestations of objective intervals projected onto different frames.
