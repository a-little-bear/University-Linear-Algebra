window.MathJax = {
  loader: {load: ['[tex]/boldsymbol']},
  tex: {
    packages: {'[+]': ['boldsymbol']},
    inlineMath: [["\\(", "\\)"]],
    displayMath: [["\\[", "\\]"]],
    processEscapes: true,
    processEnvironments: true,
    tags: 'ams',
    macros: {
      R: "\\mathbb{R}",
      C: "\\mathbb{C}",
      F: "\\mathbb{F}",
      N: "\\mathbb{N}",
      Z: "\\mathbb{Z}",
      rank: "\\operatorname{rank}",
      tr: "\\operatorname{tr}",
      diag: "\\operatorname{diag}",
      dim: "\\operatorname{dim}",
      ker: "\\operatorname{ker}",
      im: "\\operatorname{im}",
      span: "\\operatorname{span}",
      proj: "\\operatorname{proj}",
      adj: "\\operatorname{adj}",
      sgn: "\\operatorname{sgn}",
      det: "\\operatorname{det}",
      vec: "\\operatorname{vec}",
      svd: "\\operatorname{SVD}",
      spec: "\\operatorname{spec}",
      Re: "\\operatorname{Re}",
      Im: "\\operatorname{Im}"
    }
  },
  options: {
    ignoreHtmlClass: ".*|",
    processHtmlClass: "arithmatex"
  }
};

document$.subscribe(() => {
  MathJax.typesetPromise()
})
