# Piecewise linear manifolds
_PL_ is a _GAP_ library for work with piecewise linear manifolds. We define a
manifold as _CW_-complex with extra axiom:

> Boundary of every n-dimentional cell (n > 1) is a (n-1)-dimensional sphere.

# Usage
Please see the file _doc/manual.pdf_ for manual of commands. We provide the
major commands for cell (re)bilding and some algebraic and combinatorial ideas
connected with cell-composition.

Load package by next command:

```
    SetPackagePath("PL", "path/to/package/gapPL");
    LoadPackage("PL");
```
