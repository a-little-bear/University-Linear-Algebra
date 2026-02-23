# Chapter 59B: Tropical Geometry and Applications

<div class="context-flow" markdown>

**Prerequisites**: Tropical Semiring (Ch59A) · Linear Equations (Ch1) · Convex Geometry (Ch64A) · Graph Theory (Ch27)

**Chapter Outline**: Tropical Hyperplanes → Tropical Linear Spaces → Amoebas and Tropical Varieties → Log-Glass Transform → Tropical Eigenvalues and Newton Polygons → Tropical Matrix Factorization → Phylogeny and Tropical Geometry → Application in Auction Theory

**Extension**: Tropical geometry is the "combinatorial skeleton" of algebraic geometry; it converts polynomial equations into piecewise-linear polyhedral complexes.

</div>

Tropical geometry is the study of algebraic varieties over the tropical semiring. It serves as a bridge between the continuous world of complex algebraic geometry and the discrete world of polyhedral combinatorics. By applying the "tropicalization" map, complex curves are transformed into graphs, and high-dimensional spaces are mapped to piecewise-linear complexes.

---

## 59B.1 Tropical Hyperplanes and Linear Spaces

!!! definition "Definition 59B.1 (Tropical Hyperplane)"
    A tropical hyperplane defined by $L(x) = \bigoplus_{i=1}^n a_i \odot x_i$ is the set of points where the maximum in $\max(a_1+x_1, \dots, a_n+x_n)$ is achieved by **at least two** terms. This creates a polyhedral structure in $\mathbb{R}^n$.

!!! theorem "Theorem 59B.3 (Fundamental Theorem of Tropical Geometry)"
    The tropicalization of an algebraic variety $V \subseteq (\mathbb{C}^*)^n$ is exactly the set of component-wise logarithms of points in $V$ in the limit of a degenerate valuation.

---

## Exercises

1. **[Tropical Line] Draw the tropical line $L(x, y) = 0 \oplus x \oplus y$ in $\mathbb{R}^2$.**
   ??? success "Solution"
       The set where $\max(0, x, y)$ is attained at least twice consists of three rays emanating from the origin: one going left ($x=y \ge 0$), one going down ($x=0 \ge y$), and one going diagonal-left-down ($y=0 \ge x$). This is the "tropical Y-shape."

2. **[Matrix Relation] How does the tropical determinant relate to the Newton polygon of a polynomial?**
   ??? success "Solution"
       The tropical eigenvalues of a matrix are the slopes of the lower boundary of the Newton polygon of the characteristic polynomial. This links the discrete spectrum to the facial geometry of the polynomial's convex hull.

3. **[Amoebas] Define the amoeba of a complex variety and its relation to tropical geometry.**
   ??? success "Solution"
       The amoeba is the image of the variety under the log-absolute-value map $\mathbb{C}^* 	o \mathbb{R}$. As the base of the logarithm goes to infinity, the amoeba shrinks and converges to the tropical variety (the "skeleton").

4. **[Phylogeny] Why are phylogenetic trees modeled as points in a tropical Grassmannian?**
   ??? success "Solution"
       The four-point condition for tree metrics (Ch70B) is equivalent to the tropical Plücker relations. Thus, the space of all possible trees on $n$ taxa forms a tropical linear space inside the tropical Grassmannian $\operatorname{TrGr}(2, n)$.

5. **[Auction Theory] Describe how tropical matrix multiplication models a first-price auction.**
   ??? success "Solution"
       Let $V_{ij}$ be the value of item $j$ to bidder $i$. The optimal assignment (maximizing total value) is given by the tropical determinant. The stable prices correspond to the tropical dual variables (shadow prices) in the dual linear programming problem.

6. **[Calculation] Compute the tropical hyperplane $H$ for $2 \odot x \oplus 3 \odot y \oplus 0 = \max(x+2, y+3, 0)$. Where is the vertex?**
   ??? success "Solution"
       The vertex is where all three terms are equal: $x+2 = y+3 = 0 \implies x=-2, y=-3$. The rays extend from $(-2, -3)$ in the directions $(-1, 0), (0, -1),$ and $(1, 1)$.

7. **[Tropical Rank] Contrast Kapranov rank with tropical determinantal rank.**
   ??? success "Solution"
       Kapranov rank is the minimum rank of a classical matrix that tropicalizes to the given matrix. Tropical determinantal rank is based on the vanishing of tropical minors. Unlike the classical case, these two definitions do not always coincide.

8. **[Linear Spaces] What is a "Tropical Linear Space" geometrically?**
   ??? success "Solution"
       It is a balanced polyhedral complex that satisfies the tropical exchange axiom. Geometrically, it looks like a collection of polyhedral cells glued together such that the "flow" of orientations is conserved at each ridge (balancing condition).

9. **[Viterbi Algorithm] How is the Viterbi algorithm in HMMs related to tropical matrix-vector multiplication?**
   ??? success "Solution"
       The Viterbi algorithm finds the most likely path in a trellis. By taking the negative log-probabilities, the product of probabilities becomes a sum, and the max-probability becomes a min-sum. This is exactly tropical multiplication in the min-plus semiring.

10. **[Discrete Convexity] Explain the link between tropical geometry and M-convex/L-convex functions in discrete optimization.**
    ??? success "Solution"
        Tropical linear spaces are the dual objects to matroids. M-convex and L-convex functions are the discrete analogs of convex functions that behave "linearly" in the tropical sense, providing the foundation for greedy algorithms and discrete duality.

## Chapter Summary

This chapter explores the combinatorial skeleton of algebraic structures:

1. **Skeletonization**: Defined tropical hyperplanes and varieties as the piecewise-linear limits of algebraic curves and surfaces.
2. **Spectral Geometry**: Linked tropical eigenvalues to Newton polygons, bridging algebraic degree and polyhedral slope.
3. **Data Metric Spaces**: Formulated the space of phylogenetic trees as a tropical Grassmannian, unifying evolutionary biology and discrete geometry.
4. **Optimization Duality**: Demonstrated the application of tropical algebra in auctions and path-finding algorithms via the log-glass transform.
