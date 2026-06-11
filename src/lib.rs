//! Delaunay triangulation using the Bowyer-Watson algorithm.

use std::collections::HashSet;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point {
    pub x: f64,
    pub y: f64,
}

impl Point {
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Triangle(pub usize, pub usize, pub usize);

impl Triangle {
    pub fn contains_vertex(&self, v: usize) -> bool {
        self.0 == v || self.1 == v || self.2 == v
    }
}

fn circumcircle_contains(a: &Point, b: &Point, c: &Point, p: &Point) -> bool {
    let ax = a.x - p.x;
    let ay = a.y - p.y;
    let bx = b.x - p.x;
    let by = b.y - p.y;
    let cx = c.x - p.x;
    let cy = c.y - p.y;

    let det = ax * (by - cy) + bx * (cy - ay) + cx * (ay - by);

    // If det <= 0, the triangle is degenerate or points are collinear
    // We consider the point inside in that case
    if det.abs() < 1e-12 {
        return true;
    }

    let asq = ax * ax + ay * ay;
    let bsq = bx * bx + by * by;
    let csq = cx * cx + cy * cy;

    let circ_x = (asq * (by - cy) + bsq * (cy - ay) + csq * (ay - by)) / (2.0 * det);
    let circ_y = (asq * (cx - bx) + bsq * (ax - cx) + csq * (bx - ax)) / (2.0 * det);

    let radius_sq = (ax - circ_x) * (ax - circ_x) + (ay - circ_y) * (ay - circ_y);
    let dist_sq = circ_x * circ_x + circ_y * circ_y;

    dist_sq <= radius_sq
}

/// Compute Delaunay triangulation for a set of 2D points.
/// Returns a list of triangles as indices into the points slice.
pub fn triangulate(points: &[Point]) -> Vec<Triangle> {
    if points.len() < 3 {
        return Vec::new();
    }

    // Find bounding box and create super-triangle
    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;

    for p in points {
        min_x = min_x.min(p.x);
        min_y = min_y.min(p.y);
        max_x = max_x.max(p.x);
        max_y = max_y.max(p.y);
    }

    let dx = max_x - min_x;
    let dy = max_y - min_y;
    let dmax = dx.max(dy);
    let mid_x = (min_x + max_x) / 2.0;
    let mid_y = (min_y + max_y) / 2.0;

    // Super-triangle vertices appended at end
    let mut all_points = points.to_vec();
    let st0 = all_points.len();
    all_points.push(Point::new(mid_x - 20.0 * dmax, mid_y - dmax));
    let st1 = all_points.len();
    all_points.push(Point::new(mid_x, mid_y + 20.0 * dmax));
    let st2 = all_points.len();
    all_points.push(Point::new(mid_x + 20.0 * dmax, mid_y - dmax));

    let mut triangles = vec![Triangle(st0, st1, st2)];

    for i in 0..points.len() {
        let mut bad = Vec::new();

        for (idx, tri) in triangles.iter().enumerate() {
            let a = &all_points[tri.0];
            let b = &all_points[tri.1];
            let c = &all_points[tri.2];
            if circumcircle_contains(a, b, c, &all_points[i]) {
                bad.push(idx);
            }
        }

        // Find boundary of the polygonal hole
        let mut polygon: Vec<(usize, usize)> = Vec::new();
        for (j, &idx) in bad.iter().enumerate() {
            let tri = triangles[idx];
            let edges = [(tri.0, tri.1), (tri.1, tri.2), (tri.2, tri.0)];
            for edge in edges {
                let mut shared = false;
                for (k, &other_idx) in bad.iter().enumerate() {
                    if j != k {
                        let other = triangles[other_idx];
                        let other_edges = [(other.0, other.1), (other.1, other.2), (other.2, other.0)];
                        for oe in &other_edges {
                            if (edge.0 == oe.0 && edge.1 == oe.1) || (edge.0 == oe.1 && edge.1 == oe.0) {
                                shared = true;
                                break;
                            }
                        }
                    }
                    if shared { break; }
                }
                if !shared {
                    polygon.push(edge);
                }
            }
        }

        // Remove bad triangles (in reverse order to keep indices valid)
        let mut bad_sorted = bad;
        bad_sorted.sort_unstable_by(|a, b| b.cmp(a));
        for idx in bad_sorted {
            triangles.remove(idx);
        }

        // Re-triangulate the polygonal hole
        for edge in polygon {
            triangles.push(Triangle(edge.0, edge.1, i));
        }
    }

    // Remove triangles that share vertices with super-triangle
    triangles.retain(|t| !t.contains_vertex(st0) && !t.contains_vertex(st1) && !t.contains_vertex(st2));

    triangles
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_three_points() {
        let pts = vec![Point::new(0.0, 0.0), Point::new(1.0, 0.0), Point::new(0.5, 1.0)];
        let tris = triangulate(&pts);
        assert_eq!(tris.len(), 1);
    }

    #[test]
    fn test_four_points() {
        let pts = vec![
            Point::new(0.0, 0.0),
            Point::new(1.0, 0.0),
            Point::new(1.0, 1.0),
            Point::new(0.0, 1.0),
        ];
        let tris = triangulate(&pts);
        assert_eq!(tris.len(), 2);
    }
}
