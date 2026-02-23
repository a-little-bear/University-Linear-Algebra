# Chapter 70B: Linear Algebra in Molecular Biology

<div class="context-flow" markdown>

**Prerequisites**: Graph Theory (Ch27) · Matrix Decompositions (Ch10) · Non-negative Matrices (Ch17) · Probability Theory

**Chapter Outline**: Scoring Matrices for Sequence Alignment (BLOSUM, PAM) → Markov Models of Evolution (Jukes-Cantor, Kimura) → Phylogenetic Trees & Distance Matrices → Gene Expression Profiling (PCA and NMF) → Metabolic Network Analysis: The Stoichiometric Matrix $S$ → Flux Balance Analysis (FBA) and the Null Space of $S$ → Protein Folding & Graph Laplacians → Applications: Drug Discovery, Cancer Subtyping, and Gene Regulatory Networks

**Extension**: Molecular biology is the digital encoding of life; it uses linear algebra to link microscopic chemical reactions to macroscopic biological phenotypes, proving that metabolic processes are essentially linear programs under stoichiometric constraints.

</div>

How can we find evidence of species evolution within massive DNA sequences? How can gene chip data predict the occurrence of cancer? **Molecular Biology Linear Algebra** transforms the interactions of biological macromolecules into matrix operations. From scoring matrices describing amino acid substitution probabilities to null-space analysis of metabolic flows, linear algebra provides precise mathematical tools for decoding the blueprint of life.

---

## 70B.1 Sequence Alignment and Scoring Matrices

!!! definition "Definition 70B.1 (Substitution Matrix)"
    In sequence alignment, a **substitution matrix** (e.g., BLOSUM62 or PAM) is a $20 \times 20$ symmetric matrix whose entries $s_{ij}$ represent the log-likelihood score of amino acid $i$ being replaced by amino acid $j$ through evolution.
    **Function**: it transforms biological "similarity" into computable algebraic sums.

---

## 70B.2 Metabolic Networks and Stoichiometry

!!! definition "Definition 70B.2 (Stoichiometric Matrix $S$)"
    For a system with $m$ metabolites and $n$ reactions, the **Stoichiometric Matrix** $S$ is an $m \times n$ matrix where $s_{ij}$ is the participating coefficient of metabolite $i$ in reaction $j$ (positive for products, negative for reactants).

!!! theorem "Theorem 70B.1 (Steady-State Flux)"
    At metabolic steady state, the rate of change of metabolite concentrations is zero:
    $$S \mathbf{v} = 0$$
    where $\mathbf{v}$ is the vector of reaction rates (fluxes). This means all possible steady-state fluxes lie in the **null space of $S$**, $N(S)$.

---

## 70B.3 Flux Balance Analysis (FBA)

!!! technique "FBA Optimization"
    Since the dimension of $N(S)$ is typically large (the system is redundant), organisms choose fluxes to optimize specific objectives, such as maximizing biomass yield.
    $$\max \mathbf{c}^T \mathbf{v} \quad \text{subject to } S\mathbf{v} = 0, \quad \mathbf{v}_{\min} \le \mathbf{v} \le \mathbf{v}_{\max}$$
    This is a classic **linear programming** problem.

---

## Exercises

1.  **[Basics] Write the stoichiometric column vector for the reaction $A + B \to C$.**
    ??? success "Solution"
        Assuming the metabolite order is $(A, B, C)^T$, the vector is $(-1, -1, 1)^T$.

2.  **[Null Space] If $S$ has 10 metabolites and 15 reactions, and $\operatorname{rank}(S)=8$, what is the dimension of the steady-state flux space?**
    ??? success "Solution"
        $\dim N(S) = 15 - 8 = 7$.

3.  **[Alignment] In a log-odds scoring matrix, what do positive and negative values represent?**
    ??? success "Solution"
        Positive values mean the substitution occurs more frequently in evolution than expected by chance (conserved); negative values mean the substitution is rare and likely disrupts protein function.

4.  **[PCA] Why is PCA commonly used in analyzing microarray gene chip data?**
    ??? success "Solution"
        Gene data is extremely high-dimensional (tens of thousands of genes). PCA extracts the principal components that define cell states (e.g., healthy vs. diseased), significantly reducing noise.

5.  **[Evolution] Briefly describe the characteristics of the transition matrix in the Jukes-Cantor model.**
    ??? success "Solution"
        It is a symmetric stochastic matrix where the eigenvalues reflect the rate at which nucleotide distributions tend toward uniformity over time.

6.  **[Network] What is an "essential reaction" in a metabolic network?**
    ??? success "Solution"
        A reaction such that setting its rate to zero in the $S\mathbf{v}=0$ constraint forces the objective function (e.g., growth rate) to zero.

7.  **[Graph Theory] What does the Fiedler eigenvalue $\lambda_2$ of a protein contact map represents?**
    ??? success "Solution"
        It represents the structural compactness or "folding rate" of the protein; a smaller $\lambda_2$ corresponds to more loosely connected, flexible regions.

8.  **[Calculation] Find a basis vector for $S = \begin{pmatrix} -1 & 1 \\ 1 & -1 \end{pmatrix}$.**
    ??? success "Solution"
        $\mathbf{v} = (1, 1)^T$. Physical meaning: This is a reversible reaction ($A \leftrightarrow B$) where rates are equal at steady state.

9.  **[NMF] State the advantage of NMF in gene clustering over PCA.**
    ??? success "Solution"
        NMF results are non-negative, allowing the decomposition of gene expression into additive "modules," which is more consistent with biological intuitions of gene co-regulation.

10. **[Limits] As $t \to \infty$, what effect does a Markov evolution matrix have on a DNA sequence?**
    ??? success "Solution"
        The sequence loses all its original information, and the nucleotide distribution converges to the Perron vector (steady-state distribution) of the transition matrix.

## Chapter Summary

Linear algebra in molecular biology establishes the rigorous format of life processes:

1.  **Algebraic Metric of Evolution**: Through scoring matrices and Markov chains, this theory quantifies the passage of time and the accumulation of mutations as spectral properties, providing the mathematical skeleton for reconstructing the Tree of Life.
2.  **Structural Constraints of Metabolism**: Stoichiometric matrices prove that life is not just a biochemical coincidence but a necessity under mass conservation and flux balance, establishing the computational cornerstone of systems biology.
3.  **Automated Pattern Extraction**: Matrix decomposition techniques peel away the ordered regulatory modules from chaotic genomic data, proving that high-dimensional life signals follow surprisingly simple linear superposition logic at their core.
