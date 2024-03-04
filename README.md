# Skoll, devourer of the Ricci

*Skoll* is a program to compute Ricci-flat metrics on nilpotent Lie algebras, based on the methods of

[D. Conti, The Ricci-flatness that lurks in weight, arXiv:2403.00697](https://arxiv.org/abs/2403.00697)

# Requireements

*Skoll* requires [cmake](https://cmake.org/) and [Wedge](https://github.com/diego-conti/wedge)

To compile, run 

    mkdir build
    cd build
    cmake ..
    cmake --build .
You may need to set the `WEDGE_PATH` environment variable to point to your installation of Wedge, e.g. 

    export WEDGE_PATH=/home/user/wedge

# Usage

*Skoll* can operate on a single Lie algebra passed on the command line or iterate through Lie algebras appearing in one of the classifications.

## Single Lie-algebra

Run *Skoll* as

    skoll sigma-diagonal|graded|filtered|any|derivations --lie-algebra <struture_constants> --all|--nice|--non-nice --columns cols
to study a single Lie algebra. If the flag `--nice` is indicated, *Skoll* assumes that the Lie algebra is nice, and only considers the grading induced by the split torus of diagonal derivations.

## Classifications 

Run *Skoll* as 

    skoll sigma-diagonal|graded|filtered|any --dimension n --all|--nice|--non-nice --columns cols
to print a table of all classified nilpotent Lie algebras in dimension 9. Use the flag `--all` to consider generic nilpotent Lie algebras; `--nice` to restrict to nice nilpotent Lie algebras; `--non-nice` to restrict to non-nice nilpotent Lie algebras (only allowed in dimension 7).

Notice that nilpotent Lie algebras of dimension 8 and 9 are not classified, so if n=8,9 only nice nilpotent Lie algebras are considered regardless of whether `--nice` is indicated.

## Modes of use

The first argument controls behaviour as follows:

- `sigma-diagonal`: Try to compute a sigma-diagonal Ricci-flat metric g_i e^i\otimes e^{\sigma_i}. Preference is given to metrics which are Ricci-flat regardless of the parameters g_i. Failing that, the polynomial equations that the parameters g_i must satisfy are printed.

- `graded`: Try to determine a weight sequence that satisfies (G1)-(G5) associated to a fixed grading. If the flag `--nice` is indicated, Skoll only considers the grading induced by the split torus of diagonal derivations. Otherwise, Skoll tries to compute a maximal split torus acting diagonally; an error is issued if this fails.

- `filtered`: Try to determine a weight sequence that satisfies (F1)-(F5) relative to the fixed basis.

- `any`: Try to determine a Ricci-flat metric using gradings, or failing that, filtrations, or failing that, a sigma-diagonal ansatz.

- `derivations`: Print out the space of derivations, its nilradical, and try to decompose a complement of the nilradical as the sum of a compact and a split torus.

## Output

The output is meant to be included in a LaTeX file (see e.g. the ancillary file in [arXiv:2403.00697](https://arxiv.org/abs/2403.00697)). The parameter `--columns` controls how many columns should be occupied by the Lie algebra in the output. Set it to 1 for low dimensions, and 2 or 3 for higher dimensions, which will result in the structure constants taking a separate line in the resulting table.



