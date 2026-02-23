# Chapter 59B: Tropical Geometry

<div class="context-flow" markdown>

**Prerequisites**: Tropical Semiring (Ch59A) · Linear Transformations (Ch05) · Convex Sets (Ch25) · Polynomial Algebra (Ch00)

**Chapter Outline**: From Algebraic Varieties to Tropical Varieties → Geometric Representation of Tropical Polynomials → Tropical Hypersurfaces (Sets of Non-smoothness) → Amoebas and Tropicalization → Fundamental Theorem: Tropicalization as the Limit of Log Maps → Combinatorial Structure of Tropical Lines and Conics → Kapranov’s Theorem → Applications: Enumerative Geometry (Gromov-Witten Invariants), Phylogenetics (Space of Trees), and Geometry of Optimization

**Extension**: Tropical geometry is the "skeletonization" of algebraic geometry; it transforms complex algebraic curves into piecewise-linear graphs composed of line segments and rays, enabling the solution of profound geometric counting problems via combinatorial methods.

</div>

In the tropical semiring, addition is taking the minimum and multiplication is standard addition. When we define polynomials and study their zeros over this peculiar field, **Tropical Geometry** is born. Here, smooth algebraic curves turn into sharp, "spiky" linear graphs. This chapter demonstrates how tropical geometry provides a powerful computational toolkit for modern geometry by translating algebraic problems into problems of graph theory and convex geometry.

---

## 59B.1 Tropical Polynomials and Hypersurfaces

!!! definition "Definition 59B.1 (Tropical Polynomial)"
    A tropical polynomial in variables $\mathbf{x} = (x_1, \ldots, x_n)$ is an expression of the form $F(\mathbf{x}) = \bigoplus_i (c_i \otimes \mathbf{x}^{\mathbf{a}_i})$. In terms of classical arithmetic, it is a piecewise-linear convex function:
    $$f(x_1, \ldots, x_n) = \min_i \{ c_i + a_{i1}x_1 + \cdots + a_{in}x_n \}$$

!!! definition "Definition 59B.2 (Tropical Hypersurface)"
    The **Tropical Hypersurface** $V(F)$ defined by a tropical polynomial $F$ is the set of points where the function $f(\mathbf{x})$ is **not differentiable** (i.e., the points where the minimum is attained by at least two terms).

---

## 59B.2 Amoebas and Tropicalization

!!! technique "The Amoeba"
    For a complex algebraic variety $V \subset (\mathbb{C}^*)^n$, its **Amoeba** is its image under the map $\operatorname{Log}(z_1, \ldots, z_n) = (\ln|z_1|, \ldots, \ln|z_n|)$.
    **Tropicalization**: As the base of the logarithm tends to infinity, the amoeba shrinks into a skeleton made of line segments and rays. This skeleton is the **Tropical Variety**.

---

## 59B.3 Fundamental Theorem: Kapranov's Theorem

!!! theorem "Theorem 59B.1 (Kapranov’s Theorem)"
    For a hypersurface defined by a Laurent polynomial $f = \sum c_w z^w$, its tropicalization is exactly the tropical hypersurface defined by the tropical polynomial $F = \bigoplus \operatorname{val}(c_w) \otimes x^w$. This establishes a precise correspondence between classical and tropical algebraic geometry.

---

## 59B.4 Tropical Lines

!!! example "The Tropical Line"
    In $\mathbb{R}^2$, the tropical line $L(x, y) = 0 \oplus x \oplus y$ consists of three rays starting from the origin in directions $(-1, -1)$, $(1, 0)$, and $(0, 1)$. These rays satisfy the **Balancing Condition**: the weighted sum of their direction vectors is zero.

---

## Exercises

1.  **[Basics] Draw the tropical hypersurface (zero set) for $f(x) = 0 \oplus x$.**
    ??? success "Solution"
        $f(x) = \min(0, x)$. The non-differentiable point is where $0 = x$. Thus, the hypersurface is the single point $\{0\}$ on the real line.

2.  **[2D Line] List the linear regions for the tropical line $x \oplus y \oplus 5$.**
    ??? success "Solution"
        The regions are defined by: $x \le y, x \le 5$; $y \le x, y \le 5$; and $5 \le x, 5 \le y$. The tropical line is the set of boundaries where these regions meet.

3.  **[Dimension] Prove that a tropical hypersurface in $n$-dimensional space is an $(n-1)$-dimensional balanced polyhedral fan.**
    ??? success "Solution"
        This follows from the structure of piecewise-linear convex functions. Each piece corresponds to a region where one monomial dominates; the boundaries are where at least two monomials meet, resulting in a codimension-1 complex.

4.  **[Balancing] Explain the "Balancing Condition" for tropical varieties.**
    ??? success "Solution"
        At every ridge (face of codimension 2), the sum of the normal vectors to the facets meeting at that ridge, weighted by their multiplicities, must be zero. This is analogous to force balance in physics.

5.  **[Tropicalization] Tropicalize the classical line $z_1 + z_2 + 1 = 0$.**
    ??? success "Solution"
        Assuming the valuation of the coefficients is 0, the tropical polynomial is $x \oplus y \oplus 0$. Its graph is the standard three-way "tripod" centered at the origin.

6.  **[Degree] How many "infinite" directions does a tropical curve of degree $d$ have?**
    ??? success "Solution"
        Typically $3d$ directions, counting multiplicities.

7.  **[Application] Why can tropical geometry be used to count intersections of algebraic curves?**
    ??? success "Solution"
        Because tropicalization preserves topological invariants like intersection numbers (the tropical version of Bézout's Theorem). We can count intersections of tropical lines to infer the intersection count of complex curves.

8.  **[Valuations] What are Puiseux series and their role in tropical geometry?**
    ??? success "Solution"
        Puiseux series serve as the standard scalar field for tropical geometry. The valuation map (taking the lowest exponent) maps the series to real numbers, providing the mathematical mechanism for the tropicalization process.

9.  **[Geometry] What does a tropical conic (circle) look like?**
    ??? success "Solution"
        It typically looks like a hexagonal skeleton with six rays extending outwards.

****

??? success "Solution"
    

## Chapter Summary

Tropical geometry is the "low-energy" projection of algebraic geometry:

1.  **Simplification of Form**: It collapses exponentially complex algebraic equations into linear skeletons while preserving topological cores (like intersection counts and genus), greatly simplifying computations.
2.  **Link between Continuous and Discrete**: Through the limit of amoebas, tropical geometry proves that discrete combinatorial structures (graphs, polyhedral fans) are the ultimate essence of continuous algebraic forms.
3.  **New Paradigm of Computation**: As a tool for solving enumerative geometry problems, tropical geometry demonstrates how to solve "equation problems" by "counting lines," providing new intuition for modern mathematical physics.
