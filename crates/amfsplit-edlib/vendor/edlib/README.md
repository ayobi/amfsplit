# Vendored edlib

`edlib.cpp` and `edlib.h` are copied unchanged from the PyPI source distribution of the
Python `edlib` package, version 1.3.9.post1 (`edlib-1.3.9.post1.tar.gz`, SHA-256
`b0fb6e85882cab02208ccd6daa46f80cb9ff1d05764e91bf22920a01d7a6fbfa`), paths
`edlib/src/edlib.cpp` and `edlib/include/edlib.h`. `SHA256SUMS` holds the hashes of the two
files as copied.

Why this copy and not the `edlib_rs` crate: the Stage 1 prototype (`amfsplit.py` v0.3.1)
calls edlib through that Python package, and the parity tests require identical alignment
results. The `edlib_rs` crate bundles edlib 1.2.6, whose `edlib.cpp` differs from this one
(traceback bookkeeping and an alphabet-loop type). Compiling the very same source removes
that question.

edlib is by Martin Šošić, MIT licence (`LICENSE`, copied from the edlib repository).
