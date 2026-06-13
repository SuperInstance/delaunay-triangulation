# Delaunay Triangulation

**A computational geometry library that computes the Delaunay triangulation of 2D point sets** using the Bowyer-Watson incremental insertion algorithm — the standard method for generating optimal triangular meshes from scattered points.

## Why It Matters

Delaunay triangulation is one of the most important algorithms in computational geometry. Given a set of points, it produces a triangulation where no point falls inside the circumcircle of any triangle — the **Delaunay property**. This produces the "nicest" possible triangles (avoiding skinny, degenerate ones), which is critical for:

- **Finite element analysis (FEA)** — Mesh quality directly affects numerical stability
- **Terrain modeling** — GIS systems triangulate elevation points into surfaces
- **Computer graphics** — Procedural mesh generation, level of detail rendering
- **Nearest-neighbor interpolation** — Delaunay triangulation defines natural neighbor regions
- **Voronoi diagrams** — The dual graph of a Delaunay triangulation is the Voronoi diagram

**The Delaunay property:** A triangulation is Delaunay if for every triangle, no other input point lies inside its circumcircle (the circle passing through all three vertices). This maximizes the minimum angle across all triangles — no algorithm can produce a triangulation with larger minimum angles.

**The Bowyer-Watson algorithm** runs in O(n²) worst case but O(n log n) expected with random insertion order. For n > ~1000 points, divide-and-conquer or sweep-line algorithms are theoretically faster but much harder to implement correctly.

## How It Works

The algorithm proceeds in four phases:

**1. Super-triangle initialization:** Compute the bounding box of all input points and create a "super-triangle" large enough to contain all of them (20× the bounding box diameter). This simplifies the algorithm — every inserted point is guaranteed to be inside some existing triangle.

**2. Incremental insertion:** For each input point P:
- Find all triangles whose circumcircle contains P ("bad triangles")
- The union of bad triangles forms a polygonal hole
- Find the boundary edges of this hole (edges shared by exactly one bad triangle)
- Remove the bad triangles
- Create new triangles connecting P to each boundary edge

**3. Circumcircle test:** The `circumcircle_contains` function computes the determinant of a 3×3 matrix formed by the vectors from P to each triangle vertex. If the determinant is positive (for counter-clockwise triangles), P is inside the circumcircle. This is equivalent to the in-circle predicate from computational geometry, and the sign depends on triangle orientation.

**4. Super-triangle removal:** After all points are inserted, remove any triangle that shares a vertex with the super-triangle. These are artificial triangles that were only needed during construction.

**Correctness:** The key invariant is that after inserting point i, the triangulation of {p₁, ..., pᵢ} plus the super-triangle vertices is Delaunay. By induction, after inserting all points and removing super-triangle triangles, the result is the Delaunay triangulation of the input.

**Numerical robustness:** The circumcircle test uses `f64` arithmetic. For production use with adversarial inputs, exact arithmetic predicates (like Shewchuk's adaptive precision) would be needed to avoid robustness issues. The `det.abs() < 1e-12` check handles degenerate (collinear) cases.

## Quick Start

```rust
use delaunay_triangulation::{Point, triangulate};

// Define points scattered in 2D space
let points = vec![
    Point::new(0.0, 0.0),
    Point::new(1.0, 0.0),
    Point::new(0.5, 1.0),
    Point::new(0.5, 0.5), // interior point
];

// Compute Delaunay triangulation
let triangles = triangulate(&points);

for tri in &triangles {
    println!("Triangle: ({}, {}, {}) → ({:.1},{:.1}) ({:.1},{:.1}) ({:.1},{:.1})",
        tri.0, tri.1, tri.2,
        points[tri.0].x, points[tri.0].y,
        points[tri.1].x, points[tri.1].y,
        points[tri.2].x, points[tri.2].y,
    );
}
```

## API

### `Point`
- `new(x: f64, y: f64) -> Self` — Create a 2D point
- `x: f64, y: f64` — Cartesian coordinates

### `Triangle(usize, usize, usize)`
- A triangle as three indices into a points slice
- `contains_vertex(&self, v: usize) -> bool` — Check if vertex index is part of this triangle

### `triangulate(points: &[Point]) -> Vec<Triangle>`
- Compute Delaunay triangulation using Bowyer-Watson algorithm
- Returns triangles as vertex-index triples
- O(n²) worst case, O(n log n) expected with random insertion order
- Returns empty vec for fewer than 3 points

## Architecture Notes

This geometry library is used in SuperInstance's spatial analysis tools for mesh generation, spatial clustering, and visualization. The triangulation output feeds into mesh-based simulations and terrain rendering pipelines.

See the full architecture: [ARCHITECTURE.md](https://github.com/SuperInstance/SuperInstance/blob/main/ARCHITECTURE.md)

## License

MIT
