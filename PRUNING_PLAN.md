# pyProximation — Pruning Plan

**Goal:** Land on a lean Hilbert-space approximation library that depends only on
`sympy` (plus `numpy`/`scipy` for numerics), with no graphics handling and no
collocation/subregion machinery.

---

## Dependency Analysis

| Module | Imports core? | Imported by core? | Decision |
|---|---|---|---|
| `graphics.py` | No | `__init__.py` only | **DELETE** |
| `collocation.py` | Yes (`measure`, `orthsys`) | `__init__.py` + `subregion` only | **DELETE** |
| `subregion.py` | Yes (all core + collocation) | Not in public API | **DELETE** |
| `rational.py` | Yes (`measure`, `orthsys`) | Not exported | **KEEP & CLEAN** |

Surviving public API (`__init__.py`): `Measure`, `OrthSystem`, `Interpolation`.

---

## Changes Made

### Phase 1 — Deleted files

| File | Reason |
|---|---|
| `pyProximation/graphics.py` | Removed graphics handling |
| `pyProximation/collocation.py` | No reverse-deps in core; out of scope |
| `pyProximation/subregion.py` | Depends entirely on collocation |
| `doc/graphics.rst` | Docs for removed module |
| `doc/collocation.rst` | Docs for removed module |
| `doc/subregion.rst` | Docs for removed module |
| `doc/symbolic.rst` | Documents sage/symengine integration (removed) |
| `examples/Exm01.py` | Uses Collocation + Graphics |
| `examples/Exm02.py` | Uses Collocation + Graphics |
| `examples/Exm03.py` | Uses Collocation + Graphics |

### Phase 2 — Source code changes

**`pyProximation/__init__.py`**
- Removed `from .collocation import Collocation`
- Removed `from .graphics import Graphics`

**`pyProximation/base.py`**
- `DetSymEnv()`: now returns `['sympy']` when sympy is loaded; sage/symengine
  detection removed.
- `CommonSymFuncs()`: collapsed to unconditional `from sympy import …`; sage and
  symengine branches removed.

**`pyProximation/orthsys.py`**
- `__init__()`: removed sage/symengine env validation; raises on anything other
  than `'sympy'`.
- `inner()`: removed `sage.fast_callable` and `symengine.Lambdify` branches;
  unconditionally uses `sympy.lambdify`.

**`pyProximation/interpolation.py`**
- `__init__()`: raises `ValueError` for `env != 'sympy'`; updated docstring.
- `Delta()`: removed `symengine.Matrix`/`Subs` branch; unconditionally uses
  `sympy.Matrix`/`Subs`.

**`pyProximation/rational.py`**
- Removed gratuitous `from sympy import *` (no sympy names used directly in the
  class body).

**`setup.py`**
- `install_requires`: removed `matplotlib` and `itertools` (stdlib); now
  `['sympy', 'numpy', 'scipy']`.

### Phase 3 — Documentation changes

**`doc/index.rst`** — removed `symbolic`, `collocation`, `subregion`, `graphics`
from the TOC tree.

**`doc/code.rst`** — removed `automodule` entries for `collocation`, `subregion`,
`graphics`.

**`doc/introduction.rst`**
- Rewritten opening paragraph to reflect current Hilbert-space approximation focus.
- Dependencies section: removed sage, symengine, matplotlib, mayavi.

**`doc/interpolation.rst`**
- Removed Graphics-dependent 3D visualisation block and its screenshot.
- Replaced `print p` (Python 2) with `print(p)`.
- Removed plotting tail of the L²-approximation example; replaced with
  `print(intr)`.

### Phase 4 — Example update

**`examples/Exm04.py`**
- Replaced `from sympy import *` / `from pyProximation import *` with explicit
  imports.
- Removed `Graphics` usage; result is shown via `sympy.pprint`.

---

## Deferred Decisions

- `rational.py` implements useful rational Least-Squares approximation and retains
  access to the core, but is **not exported** via `__init__.py`. The maintainer
  may decide to add it to the public API or document it separately.
- `measure.py` had no sage/symengine conditional code; no changes were required.

---

## Verification Checklist

```bash
# 1. Core imports work
python -c "from pyProximation import Measure, OrthSystem, Interpolation; print('OK')"

# 2. Public API contains only the 3 core classes
python -c "import pyProximation; print([x for x in dir(pyProximation) if not x.startswith('_')])"

# 3. Surviving example runs end-to-end
python examples/Exm04.py

# 4. Documentation builds (requires sphinx)
cd doc && make html
```
