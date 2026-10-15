# Layout conventions and geometry (claw-plot-grid)

Read this when the user's field book uses different words for the layout, or when a grid lands in the wrong place.

## Corners, ranges and columns

```
            planting direction (bearing) ->
   c2 +--------+--------+--------+  ...         ranges increase from c1 toward c2
      | R6 C1  | R6 C2  | R6 C3  |              columns increase from c1 toward c4
      +--------+--------+--------+              (here to the RIGHT of the planting
        alley                                    direction when you stand at c1
      +--------+--------+--------+               facing c2: --across right)
      | R2 C1  | R2 C2  | R2 C3  |
      +--------+--------+--------+
        alley
      +--------+--------+--------+
      | R1 C1  | R1 C2  | R1 C3  |
   c1 +--------+--------+--------+ c4
```

The sketch is drawn with the planting direction pointing up the page. Each plot is `--plot-length` long (along the planting direction, alley excluded) and `--rows-per-plot x --row-spacing` wide.

Field books use other names. Map them with `--join-by`:

| Field book says | Means here |
|---|---|
| range, tier, or row (in field books that call each tier a row) | `range` |
| pass, column, col, or row (in field books that call each planter pass a row) | `column` |
| plot, PLOT, plot_number | `plot` |

When unsure, join on `plot` only and use `--fb-coords` with the field book's two position columns: the run stops if they do not map one-to-one onto the grid, which catches a wrong serpentine direction or start corner.

## Numbering schemes (4 ranges x 3 columns)

```
serpentine, --sweep along, start c1     row-wise, --sweep along, start c1
 range 4 |  4   5  12                     range 4 |  4   8  12
 range 3 |  3   6  11                     range 3 |  3   7  11
 range 2 |  2   7  10                     range 2 |  2   6  10
 range 1 |  1   8   9                     range 1 |  1   5   9
           C1  C2  C3                               C1  C2  C3

serpentine, --sweep across, start c1    serpentine, --sweep along, start c3
 range 4 | 12  11  10                     range 4 |  9   8   1
 range 3 |  7   8   9                     range 3 | 10   7   2
 range 2 |  6   5   4                     range 2 | 11   6   3
 range 1 |  1   2   3                     range 1 | 12   5   4
           C1  C2  C3                               C1  C2  C3
```

`--sweep along` follows a planter pass through all ranges of one column; `--sweep across` numbers across one range at a time. Programs that number plots as range x 100 + column (101, 102, ...) should join on `ROW=range,COL=column` instead of `plot`.

## Geometry

- Local layout coordinates: `s` along the planting direction from c1, `t` across. Plot (range i, column j) spans `s` from `(i - 1)(L + A)` to `+ L` and `t` from `(j - 1)(W + G)` to `+ W`.
- Map coordinates: `x = O + k (s a + t c)`, where `O` is c1, `a` the unit vector of the grid bearing, `c` the unit vector 90 deg to the right (or left), and `k` the point scale factor of the output CRS at the field.
- `k` and the grid bearing of true north are computed numerically from PROJ: two short steps east-west and north-south around the field are projected, and their map lengths are divided by their ground lengths on the WGS84 ellipsoid (radii of curvature M and N). For UTM this reproduces Snyder (1987, eq. 8-11) to about 1e-9.
- Rigid corner fit: the layout's four corners are fitted to the given corners by least squares with rotation and translation only (2-D Procrustes). The handedness (`--across`) comes from the corner order.
- Bilinear corner fit (`--fit-corners`): a point at fractional position (u, v) of the block maps to `(1-u)(1-v) c1 + u(1-v) c2 + uv c3 + (1-u)v c4`. Lines of constant u or v stay straight, so plots stay quadrilaterals.

## Typical mistakes

| Symptom | Likely cause |
|---|---|
| Block about 3.3 times too long or short | feet read as meters or the reverse |
| Grid shifted by one plot length | plot length includes the alley, or the corner is a stake in the alley |
| Grid rotated a few degrees | magnetic bearing used as true; or corners digitized in a different CRS |
| Numbers mirror the planting map | wrong `--start-corner` or `--across` |
| Alternate columns numbered backward | `--numbering rowwise` versus `serpentine` |

## References

- Snyder JP (1987) Map projections: a working manual. USGS Professional Paper 1395. https://doi.org/10.3133/pp1395
