/**
 * Represents the edge between two sites (or one site and the edge of the box).
 *
 * @param {object|null} lSite - The Voronoi site object at the left of this
 * edge. The site object is just a reference to a site in the array of sites
 * supplied by the user when Voronoi.compute() was called. Can be null when this
 * is a border edge.
 * @param {object|null} rSite - The Voronoi site object at the right of this
 * edge. The site object is just a reference to a site in the array of sites
 * supplied by the user when Voronoi.compute() was called. Can be null when this
 * is a border edge.
 * @param {object} va - A Vertex object with an x and a y property defining the
 * start point (relative to the Voronoi site on the left) of this edge object.
 * @param {object} vb - A Vertex object with an x and a y property defining the
 * end point (relative to the Voronoi site on the left) of this edge object.
 */
export default class Edge {
  constructor(lSite, rSite) {
    this.lSite = lSite;
    this.rSite = rSite;
    this.va = null;
    this.vb = null;
  }

  setEdgeStartpoint(lSite, rSite, vertex) {
    if (!this.va && !this.vb) {
      this.va = vertex;
      this.lSite = lSite;
      this.rSite = rSite;
    } else if (this.lSite === rSite) {
      this.vb = vertex;
    } else {
      this.va = vertex;
    }
  }
}
