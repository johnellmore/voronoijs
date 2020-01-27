/**
 * An object which describes the computed Voronoi diagram.
 *
 * @param {Vertex[]} vertices - An array of unordered, unique `Vertex` objects
 * making up the Voronoi diagram. Each `Vertex` object in the list is shared by
 * many `Edge` objects.
 * @param {Edge[]} edges - An array of unordered, unique `Edge` objects making
 * up the Voronoi diagram. `Edge`s are defined by two vertices, `va` and `vb`,
 * which vertices are shared by connected edges. This mean that if you change
 * one vertex belonging to an edge, other connected edges will also be changed.
 * @param {Cell[]} cells - An array of `Cell` objects making up the Voronoi
 * diagram. A `Cell` object might have an empty array of halfedges, meaning no
 * Voronoi cell could be computed for a particular cell.
 * @param {Number} execTime - The time it took to compute the Voronoi diagram,
 * in milliseconds.
 */
export default class Diagram {
  constructor() {
    this.cells = [];
    this.edges = [];
    this.vertices = [];
    this.execTime = 0;
  }
}
