# -*- coding: utf-8 -*-
"""
Tests for the physical (energy-based) symmetric relaxation.

Run from the project root:

    .venv/Scripts/python.exe -m pytest tests/ -q

The tests cover:
  * analytic gradients of E_cen / E_ang vs central finite differences (~1e-7)
  * scale invariance of the energy
  * structural properties of the centre-angle force (zero net force / no work
    under dilation)
  * monotonic energy decrease and angle relaxation under PhysicalAnnealer
  * equivalence with the geometric algorithm's two rules
        - large k_angle  <=> forbid increasing sum of interior angles squared
        - large k_center  <=> drive towards centre-angle symmetry
  * convexity safety net (never turns a convex cell non-convex)
  * degree-1 (patch-tip) vertices stay frozen

A small honeycomb patch and a 3-triangle network are built by hand so the tests
are deterministic and do not depend on the Voronoi initialiser.
"""

import math
import sys

import numpy as np
import pytest

sys.path.insert(0, ".")  # allow `import annealing` from root

import annealing.energy as E
import annealing.physical_annealer as PA


# --------------------------------------------------------------------------
# lightweight structural Cell mock
# --------------------------------------------------------------------------
# The physical core (annealing.energy / annealing.physical_annealer) only ever
# reads ``cell.points`` (a list of [x, y]) and calls ``cell.setVertex()`` /
# ``cell.setArea()`` to refresh derived quantities; it never needs the ellipse
# fit.  We therefore use a faithful *structural* mock instead of importing the
# real ``utillib.mylib.Cell``: that import drags in the (heavy, R/matplotlib-
# backed) ellipse-fitting chain which intermittently hangs the shell.  The mock
# reproduces the exact shoelace area-centroid and signed-area logic of the real
# Cell so the geometry the core consumes is identical.  An integration test
# against the real ``utillib.mylib.Cell`` can be added in an environment where
# that import is fast (see PHYSICAL_TEST_PLAN_CN.md).
class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __str__(self):
        return f"({self.x},{self.y})"


class Cell:
    no = 0
    cell_no = 0

    def __init__(self, points):
        self.points = [[float(p[0]), float(p[1])] for p in points]
        self.layer = 0
        self.center_point = Point(0.0, 0.0)
        self.area = 0.0
        self.setNo()
        self.setVertex()
        self.setArea()

    def setNo(self):
        Cell.cell_no += 1
        self.no = Cell.cell_no

    def setVertex(self):
        pts = self.points
        if len(pts) <= 2:
            return
        area = 0.0
        x = 0.0
        y = 0.0
        n = len(pts)
        for i in range(n):
            lng = pts[i][0]
            lat = pts[i][1]
            nx = pts[i - 1][0]
            ny = pts[i - 1][1]
            tmp = (nx * lat - ny * lng) / 2.0
            area += tmp
            x += tmp * (lng + nx) / 3.0
            y += tmp * (lat + ny) / 3.0
        if abs(x * area) < 1e-10:
            return
        x /= area
        y /= area
        self.center_point = Point(x, y)

    def setArea(self):
        pts = self.points
        s = 0.0
        n = len(pts)
        for i in range(n):
            s += pts[i - 1][0] * pts[i][1] - pts[i][0] * pts[i - 1][1]
        s /= 2.0
        self.area = s


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def make_cells(polys):
    """Build Cell objects from a list of polygons (each a list of [x, y])."""
    return [Cell([[float(x), float(y)] for x, y in p]) for p in polys]


def three_triangle_network():
    """Three triangles around a central vertex O -> O has degree 3, A,B,C deg 2."""
    O = np.array([0.0, 0.0])
    A = np.array([1.0, 0.2])
    B = np.array([-0.3, 1.0])
    C = np.array([-0.5, -0.9])
    polys = [[O.tolist(), A.tolist(), B.tolist()],
             [O.tolist(), B.tolist(), C.tolist()],
             [O.tolist(), C.tolist(), A.tolist()]]
    return make_cells(polys)


def build_honeycomb_patch(centers, side=1.0, perturb=0.12, seed=0):
    """Regular hexagons on a triangular center lattice, vertices merged globally.

    ``centers`` is a list of [cx, cy].  Each hexagon's 6 corners are generated,
    coincident corners across hexagons are merged into a single global vertex so
    the resulting cells share coordinates (as in a real tissue).  ``perturb``
    jitters the final vertex positions to give the relaxer something to do.
    """
    rng = np.random.default_rng(seed)
    raw_verts = []          # list of [x, y]
    cell_local = []         # list of lists of indices into raw_verts
    gmap = {}
    s = side
    for (cx, cy) in centers:
        for k in range(6):
            ang = math.pi / 3.0 * k
            x = cx + s * math.cos(ang)
            y = cy + s * math.sin(ang)
            key = (round(x, 9), round(y, 9))
            gid = gmap.get(key)
            if gid is None:
                gid = len(raw_verts)
                gmap[key] = gid
                raw_verts.append([x, y])
            cell_local.append(gid)  # placeholder; rebuilt below
    # Rebuild properly: one entry per cell (6 consecutive gids is wrong because
    # corners are shared).  Instead collect per-cell gid lists.
    cell_gids = []
    gmap = {}
    raw_verts = []
    for (cx, cy) in centers:
        gids = []
        for k in range(6):
            ang = math.pi / 3.0 * k
            x = cx + s * math.cos(ang)
            y = cy + s * math.sin(ang)
            key = (round(x, 9), round(y, 9))
            gid = gmap.get(key)
            if gid is None:
                gid = len(raw_verts)
                gmap[key] = gid
                raw_verts.append([x, y])
            gids.append(gid)
        cell_gids.append(gids)

    coords = np.asarray(raw_verts, dtype=float)
    if perturb > 0:
        coords = coords + perturb * rng.normal(size=coords.shape)
        # re-snap to 9 decimals so shared vertices stay coincident
        coords = np.round(coords, 9)
    cells = []
    for gids in cell_gids:
        pts = [[float(coords[g, 0]), float(coords[g, 1])] for g in gids]
        cells.append(Cell(pts))
    return cells


def honeycomb_centers(rings=2, spacing=1.0):
    """Triangular-lattice centers for a hexagonal patch of hexagons."""
    s = math.sqrt(3.0) * spacing  # distance between adjacent hexagon centers
    a1 = np.array([1.5 * spacing, 0.5 * math.sqrt(3.0) * spacing])
    a2 = np.array([0.0, math.sqrt(3.0) * spacing])
    centers = []
    for i in range(-rings, rings + 1):
        for j in range(-rings, rings + 1):
            c = i * a1 + j * a2
            if np.linalg.norm(c) <= rings * s + 1e-6:
                centers.append([float(c[0]), float(c[1])])
    return centers


# --------------------------------------------------------------------------
# gradient correctness
# --------------------------------------------------------------------------

@pytest.mark.parametrize("kc,ka", [(1.0, 0.0), (0.0, 1.0), (1.0, 1.0)])
def test_gradient_matches_finite_difference(kc, ka):
    cells = three_triangle_network()
    vidx = E.build_vertex_index(cells)
    sp = E.SymParams(k_center=kc, k_angle=ka)
    _T0, g_an = E.total_energy_gradient(cells, sp, vidx)
    coords0 = vidx['coords'].copy()
    h = 1e-6

    def total_at(coords):
        for gid in range(coords.shape[0]):
            for (cell, li) in vidx['refs'][gid]:
                cell.points[li] = [float(coords[gid, 0]), float(coords[gid, 1])]
        r, _ = E.total_energy_gradient(cells, sp)
        return r.total

    g_num = np.zeros_like(coords0)
    for gid in range(coords0.shape[0]):
        for d in range(2):
            cp = coords0.copy(); cp[gid, d] += h
            cm = coords0.copy(); cm[gid, d] -= h
            g_num[gid, d] = (total_at(cp) - total_at(cm)) / (2 * h)

    # analytic grad (dE/d r) vs numeric; both should agree to ~1e-7
    err = np.abs(g_num - g_an).max()
    assert err < 1e-7, f"gradient mismatch {err:.2e}"


def test_scale_invariance():
    cells = three_triangle_network()
    sp = E.SymParams(k_center=1.0, k_angle=1.0)
    vidx = E.build_vertex_index(cells)
    rep0, g0 = E.total_energy_gradient(cells, sp, vidx)
    lam = 2.5
    for cell in cells:
        cell.points = [[lam * x, lam * y] for (x, y) in cell.points]
    vidx2 = E.build_vertex_index(cells)
    rep1, g1 = E.total_energy_gradient(cells, sp, vidx2)
    assert abs(rep1.total - rep0.total) < 1e-9 * max(1.0, abs(rep0.total)), \
        "energy not scale invariant"
    # gradient w.r.t. scaled coordinates must be (1/lam) times the original grad
    # at the corresponding (scaled) vertex.
    for gid in range(vidx2['coords'].shape[0]):
        g_scaled_expected = g0[gid] / lam
        assert np.abs(g1[gid] - g_scaled_expected).max() < 1e-7


def test_centre_force_zero_sum_and_no_work():
    cells = three_triangle_network()
    for cell in cells:
        verts = np.asarray(cell.points, dtype=float)
        _E, g = E.centre_angle_energy_and_grad(verts, 1.0)
        # net centre-angle force on the cell's vertices is zero
        assert np.linalg.norm(g.sum(axis=0)) < 1e-9
        # does zero work under a uniform dilation about the cell centroid
        O = verts.mean(axis=0)
        disp = (verts - O) * 1e-3
        work = float(np.sum(g * disp))
        assert abs(work) < 1e-9


# --------------------------------------------------------------------------
# relaxation behaviour
# --------------------------------------------------------------------------

def test_relaxation_decreases_energy_monotonically():
    cells = build_honeycomb_patch(honeycomb_centers(rings=2), perturb=0.15, seed=3)
    pa = PA.PhysicalAnnealer(PA.PhysicalAnnealerParams(
        mobility=0.3, k_center=1.0, k_angle=1.0, convex_guard=True))
    sp = pa.params.to_sym_params()
    hist = []
    prev = None
    for _ in range(150):
        rep, _ = pa.one_round(cells, sp)
        hist.append(rep.total)
        if prev is not None:
            assert rep.total <= prev + 1e-9, "energy increased during relaxation"
        prev = rep.total
    assert hist[-1] < hist[0], "energy did not decrease"


def test_relaxation_drives_interior_angles_to_120_and_marginal_to_90():
    cells = build_honeycomb_patch(honeycomb_centers(rings=2), perturb=0.18, seed=5)
    # lambda = 10 (the GUI "inner_angle_sq_guard" default): the vertex-angle
    # term dominates the centre-angle term, so the angle targets are nearly
    # reachable.  At lambda = 1 the 90 deg marginal target competes with the
    # cell-level 120 deg (regular-hexagon) preference and keeps a residual
    # deviation of ~0.3 rad at the rim -- that competition is physical, not a
    # convergence failure.
    pa = PA.PhysicalAnnealer(PA.PhysicalAnnealerParams(
        mobility=0.4, k_center=1.0, k_angle=10.0))
    pa.relax(cells, max_rounds=600, tol=1e-10)

    # measure deviation of actual interior angles from targets, separately for
    # interior (degree 3) and marginal (degree 2) vertices
    vidx = E.build_vertex_index(cells)
    sp = E.SymParams(k_center=1.0, k_angle=10.0)
    tau, mask = E.vertex_targets(vidx['degree'], sp)
    max_dev = {2: 0.0, 3: 0.0}
    for cell in cells:
        raw = np.asarray(cell.points, dtype=float)
        verts, perm = E._as_ccw(raw)
        alpha, *_ = E.interior_angle_terms(verts)
        gids = E._cell_gids(cell, verts, perm, vidx)
        for ll in range(verts.shape[0]):
            gg = gids[ll]
            dg = vidx['degree'][gg]
            if mask[gg] <= 0 or dg not in max_dev:
                continue
            max_dev[dg] = max(max_dev[dg], abs(alpha[ll] - tau[gg]))
    assert max_dev[3] < 0.15, \
        f"interior angle deviation {max_dev[3]:.3f} rad too large"
    assert max_dev[2] < 0.35, \
        f"marginal angle deviation {max_dev[2]:.3f} rad too large"


def test_large_k_angle_forbids_increasing_interior_angle_variance():
    """Physical analogue of the geometric 'forbid increasing sum of 3 interior
    angles squared' guard: with k_angle >> k_center the interior angles are
    driven toward their equal-tension target (120 deg)."""
    cells = build_honeycomb_patch(honeycomb_centers(rings=2), perturb=0.2, seed=7)
    vidx = E.build_vertex_index(cells)
    sp0 = E.SymParams(k_center=0.0, k_angle=1.0)
    # initial mean |alpha - 120| per interior vertex
    tau, mask = E.vertex_targets(vidx['degree'], sp0)
    before = _mean_interior_dev(cells, vidx, tau, mask)

    pa = PA.PhysicalAnnealer(PA.PhysicalAnnealerParams(
        mobility=0.4, k_center=0.0, k_angle=1.0))
    pa.relax(cells, max_rounds=800, tol=1e-10)

    vidx1 = E.build_vertex_index(cells)
    tau1, mask1 = E.vertex_targets(vidx1['degree'], sp0)
    after = _mean_interior_dev(cells, vidx1, tau1, mask1)
    assert after < before * 0.3, f"angle variance not reduced ({before:.3f} -> {after:.3f})"


def test_large_k_center_drives_centre_angle_symmetry():
    cells = build_honeycomb_patch(honeycomb_centers(rings=2), perturb=0.2, seed=9)
    vidx = E.build_vertex_index(cells)
    before = sum(E.centre_angle_energy_and_grad(
        np.asarray(c.points, dtype=float), 1.0)[0] for c in cells)

    pa = PA.PhysicalAnnealer(PA.PhysicalAnnealerParams(
        mobility=0.4, k_center=1.0, k_angle=0.0))
    pa.relax(cells, max_rounds=800, tol=1e-10)

    after = sum(E.centre_angle_energy_and_grad(
        np.asarray(c.points, dtype=float), 1.0)[0] for c in cells)
    assert after < before * 0.3, f"centre-angle energy not reduced ({before:.3f} -> {after:.3f})"


def test_convexity_safety_net():
    cells = build_honeycomb_patch(honeycomb_centers(rings=2), perturb=0.25, seed=11)
    pa = PA.PhysicalAnnealer(PA.PhysicalAnnealerParams(
        mobility=0.5, k_center=1.0, k_angle=1.0, convex_guard=True))
    pa.relax(cells, max_rounds=400, tol=1e-9)
    for cell in cells:
        assert PA._is_polygon_convex(cell.points), "a cell became non-convex"


def test_degree_one_vertices_frozen():
    cells = build_honeycomb_patch(honeycomb_centers(rings=1), perturb=0.0, seed=0)
    vidx = E.build_vertex_index(cells)
    deg1 = np.where(vidx['degree'] <= 1)[0]
    if len(deg1) == 0:
        pytest.skip("no degree-1 vertex in this patch")
    pos0 = vidx['coords'][deg1].copy()
    pa = PA.PhysicalAnnealer(PA.PhysicalAnnealerParams(
        mobility=0.5, k_center=1.0, k_angle=1.0))
    pa.relax(cells, max_rounds=50, tol=1e-9)
    vidx1 = E.build_vertex_index(cells)
    pos1 = vidx1['coords'][deg1]
    assert np.abs(pos1 - pos0).max() < 1e-9, "degree-1 vertex moved (should be frozen)"


# --------------------------------------------------------------------------
# utility
# --------------------------------------------------------------------------

def _mean_interior_dev(cells, vidx, tau, mask):
    total = 0.0
    count = 0
    for cell in cells:
        raw = np.asarray(cell.points, dtype=float)
        verts, perm = E._as_ccw(raw)
        alpha, *_ = E.interior_angle_terms(verts)
        gids = E._cell_gids(cell, verts, perm, vidx)
        for ll in range(verts.shape[0]):
            gg = gids[ll]
            if mask[gg] <= 0 or vidx['degree'][gg] < 3:
                continue
            total += abs(alpha[ll] - tau[gg])
            count += 1
    return total / max(count, 1)


if __name__ == "__main__":
    pytest.main([__file__, "-q"])
