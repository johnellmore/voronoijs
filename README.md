# voronoijs

A Javascript implementation of Steven J. Fortune's algorithm to efficiently compute Voronoi diagrams. There are no production dependencies--only dev dependencies for development usage. It also contains no rendering code: that is left to the user of the library.

This implementation is forked from [Raymond Hill's excellent JS library](https://github.com/gorhill/Javascript-Voronoi). The goal with this fork is to support the modified version of the Fortune algorithm which can construct an additively weighted Voronoi diagram. Also, this refactors much of the implementation to use modern JavaScript syntax (ES modules, classes, etc.), and adds unit tests
around the library.

VERY MUCH A WORK IN PROGRESS. THIS PROBABLY ISN'T WHAT YOU WANT RIGHT NOW.

## Usage

```javascript
const voronoi = new Voronoi();
const bbox = { xl: 0, xr: 800, yt: 0, yb: 600 };
// xl is x-left, xr is x-right, yt is y-top, and yb is y-bottom
const sites = [
  { x: 200, y: 200 },
  { x: 50, y: 250 } /* , ... */
];

// A vertex is an object exhibiting `x` and `y` properties. The Voronoi object
// will add a unique `voronoiId` property to all sites. The `voronoiId` can be
// used as a key to lookup the associated cell in `diagram.cells`.

const diagram = voronoi.compute(sites, bbox);
```

The returned 'diagram' variable is a Javascript object with the following
properties:

#### `diagram.vertices`

An array of unordered, unique `Vertex` objects making up the Voronoi diagram.
Each `Vertex` object in the list is shared by many `Edge` objects.

#### `diagram.edges`

An array of unordered, unique `Edge` objects making up the Voronoi diagram.
`Edges` are defined by two vertices, `va` and `vb`, which vertices are shared
by connected edges. This means that if you change one vertex belonging to an
edge, other connected edges will also be changed.

#### `diagram.cells`

An array of `Cell` objects making up the Voronoi diagram. A `Cell` object might
have an empty array of `halfedges`, meaning no Voronoi cell could be computed
for a particular cell.

#### `diagram.execTime`

The time it took to compute the Voronoi diagram, in milliseconds.

---

Added on October 12, 2013: In order to help improve performance,
`Voronoi.recycle()` has been added to allow the recycling of a returned Voronoi
diagram. Usage:

```javascript
var diagram;

// some kind of loop starting here (whether outright or through a timer)

voronoi.recycle(diagram);

// diagram.vertices, diagram.edges and diagram.cells can no longer be used!

diagram = voronoi.compute(sites, bbox);

// do stuff with content of `diagram`
```

This new method helps performance significantly when re-computing a Voronoi
diagram, as it saves on memory allocation, and associated garbage collection.

## License

This software is released under the MIT License (like the original implementation).

## Further reading

- [Voronoi Diagram](https://en.wikipedia.org/wiki/Voronoi_diagram)
- [Fortune's Algorithm](https://en.wikipedia.org/wiki/Fortune's_algorithm)
- [Weighted Voronoi Algorithm](https://en.wikipedia.org/wiki/Weighted_Voronoi_diagram)
