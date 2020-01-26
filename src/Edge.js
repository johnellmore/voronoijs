// ---------------------------------------------------------------------------
// Edge methods
//
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
