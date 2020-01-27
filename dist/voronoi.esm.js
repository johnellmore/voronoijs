// rhill 2011-06-07: For some reasons, performance suffers significantly
// when instanciating a literal object instead of an empty ctor
class Beachsection {}

// ---------------------------------------------------------------------------
// Cell methods

class Cell {
  constructor(site) {
    this.halfedges = [];
    this.closeMe = false;
    this.site = site;
  }

  prepareHalfedges() {
    const { halfedges } = this;
    let iHalfedge = halfedges.length;
    let edge;
    // get rid of unused halfedges
    // rhill 2011-05-27: Keep it simple, no point here in trying
    // to be fancy: dangling edges are a typically a minority.
    while (iHalfedge--) {
      edge = halfedges[iHalfedge].edge;
      if (!edge.vb || !edge.va) {
        halfedges.splice(iHalfedge, 1);
      }
    }

    // rhill 2011-05-26: I tried to use a binary search at insertion
    // time to keep the array sorted on-the-fly (in Cell.addHalfedge()).
    // There was no real benefits in doing so, performance on
    // Firefox 3.6 was improved marginally, while performance on
    // Opera 11 was penalized marginally.
    halfedges.sort((a, b) => {
      return b.angle - a.angle;
    });
    return halfedges.length;
  }

  getNeighborIds() {
    const neighbors = [];
    let iHalfedge = this.halfedges.length;
    let edge;
    while (iHalfedge--) {
      edge = this.halfedges[iHalfedge].edge;
      if (edge.lSite !== null && edge.lSite.voronoiId !== this.site.voronoiId) {
        neighbors.push(edge.lSite.voronoiId);
      } else if (
        edge.rSite !== null &&
        edge.rSite.voronoiId !== this.site.voronoiId
      ) {
        neighbors.push(edge.rSite.voronoiId);
      }
    }
    return neighbors;
  }

  /**
   * Compute the bounding box.
   */
  getBbox() {
    const { halfedges } = this;
    let iHalfedge = halfedges.length;
    let xmin = Infinity;
    let ymin = Infinity;
    let xmax = -Infinity;
    let ymax = -Infinity;
    let v;
    let vx;
    let vy;
    while (iHalfedge--) {
      v = halfedges[iHalfedge].getStartpoint();
      vx = v.x;
      vy = v.y;
      if (vx < xmin) {
        xmin = vx;
      }
      if (vy < ymin) {
        ymin = vy;
      }
      if (vx > xmax) {
        xmax = vx;
      }
      if (vy > ymax) {
        ymax = vy;
      }
      // we dont need to take into account end point,
      // since each end point matches a start point
    }
    return {
      x: xmin,
      y: ymin,
      width: xmax - xmin,
      height: ymax - ymin
    };
  }

  /**
   * Return whether a point is inside, on, or outside the cell:
   * -1: point is outside the perimeter of the cell
   *  0: point is on the perimeter of the cell
   *  1: point is inside the perimeter of the cell
   */
  pointIntersection(x, y) {
    // Check if point in polygon. Since all polygons of a Voronoi
    // diagram are convex, then:
    // http://paulbourke.net/geometry/polygonmesh/
    // Solution 3 (2D):
    //   "If the polygon is convex then one can consider the polygon
    //   "as a 'path' from the first vertex. A point is on the interior
    //   "of this polygons if it is always on the same side of all the
    //   "line segments making up the path. ...
    //   "(y - y0) (x1 - x0) - (x - x0) (y1 - y0)
    //   "if it is less than 0 then P is to the right of the line segment,
    //   "if greater than 0 it is to the left, if equal to 0 then it lies
    //   "on the line segment"
    const { halfedges } = this;
    let iHalfedge = halfedges.length;
    let halfedge;
    let p0;
    let p1;
    let r;
    while (iHalfedge--) {
      halfedge = halfedges[iHalfedge];
      p0 = halfedge.getStartpoint();
      p1 = halfedge.getEndpoint();
      r = (y - p0.y) * (p1.x - p0.x) - (x - p0.x) * (p1.y - p0.y);
      if (!r) {
        return 0;
      }
      if (r > 0) {
        return -1;
      }
    }
    return 1;
  }
}

// ---------------------------------------------------------------------------
// Circle event methods

// rhill 2011-06-07: For some reasons, performance suffers significantly
// when instanciating a literal object instead of an empty ctor
class CircleEvent {
  constructor() {
    // rhill 2013-10-12: it helps to state exactly what we are at ctor time.
    this.arc = null;
    this.rbLeft = null;
    this.rbNext = null;
    this.rbParent = null;
    this.rbPrevious = null;
    this.rbRed = false;
    this.rbRight = null;
    this.site = null;
    this.x = 0;
    this.y = 0;
    this.ycenter = 0;
  }
}

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
class Diagram {
  constructor() {
    this.cells = [];
    this.edges = [];
    this.vertices = [];
    this.execTime = 0;
  }
}

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
class Edge {
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

class Halfedge {
  constructor(edge, lSite, rSite) {
    this.site = lSite;
    this.edge = edge;
    // 'angle' is a value to be used for properly sorting the
    // halfsegments counterclockwise. By convention, we will
    // use the angle of the line defined by the 'site to the left'
    // to the 'site to the right'.
    // However, border edges have no 'site to the right': thus we
    // use the angle of line perpendicular to the halfsegment (the
    // edge should have both end points defined in such case.)
    if (rSite) {
      this.angle = Math.atan2(rSite.y - lSite.y, rSite.x - lSite.x);
    } else {
      const { va, vb } = edge;
      // rhill 2011-05-31: used to call getStartpoint()/getEndpoint(),
      // but for performance purpose, these are expanded in place here.
      if (edge.lSite === lSite) {
        this.angle = Math.atan2(vb.x - va.x, va.y - vb.y);
      } else {
        this.angle = Math.atan2(va.x - vb.x, vb.y - va.y);
      }
    }
  }

  getStartpoint() {
    return this.edge.lSite === this.site ? this.edge.va : this.edge.vb;
  }

  getEndpoint() {
    return this.edge.lSite === this.site ? this.edge.vb : this.edge.va;
  }
}

/* eslint-disable no-param-reassign */

/**
 * Red-Black tree code (based on C version of "rbtree" by Franck Bui-Huu
 * https://github.com/fbuihuu/libtree/blob/master/rb.c
 */
class RBTree {
  constructor() {
    this.root = null;
  }

  rbInsertSuccessor(node, successor) {
    let parent;
    if (node) {
      // >>> rhill 2011-05-27: Performance: cache previous/next nodes
      successor.rbPrevious = node;
      successor.rbNext = node.rbNext;
      if (node.rbNext) {
        node.rbNext.rbPrevious = successor;
      }
      node.rbNext = successor;
      // <<<
      if (node.rbRight) {
        // in-place expansion of node.rbRight.getFirst();
        node = node.rbRight;
        while (node.rbLeft) {
          node = node.rbLeft;
        }
        node.rbLeft = successor;
      } else {
        node.rbRight = successor;
      }
      parent = node;
    }
    // rhill 2011-06-07: if node is null, successor must be inserted
    // to the left-most part of the tree
    else if (this.root) {
      node = this.getFirst();
      // >>> Performance: cache previous/next nodes
      successor.rbPrevious = null;
      successor.rbNext = node;
      node.rbPrevious = successor;
      // <<<
      node.rbLeft = successor;
      parent = node;
    } else {
      // >>> Performance: cache previous/next nodes
      successor.rbPrevious = null;
      successor.rbNext = null;
      // <<<
      this.root = successor;
      parent = null;
    }
    successor.rbLeft = null;
    successor.rbRight = null;
    successor.rbParent = parent;
    successor.rbRed = true;
    // Fixup the modified tree by recoloring nodes and performing
    // rotations (2 at most) hence the red-black tree properties are
    // preserved.
    let grandpa;
    let uncle;
    node = successor;
    while (parent && parent.rbRed) {
      grandpa = parent.rbParent;
      if (parent === grandpa.rbLeft) {
        uncle = grandpa.rbRight;
        if (uncle && uncle.rbRed) {
          parent.rbRed = false;
          uncle.rbRed = false;
          grandpa.rbRed = true;
          node = grandpa;
        } else {
          if (node === parent.rbRight) {
            this.rbRotateLeft(parent);
            node = parent;
            parent = node.rbParent;
          }
          parent.rbRed = false;
          grandpa.rbRed = true;
          this.rbRotateRight(grandpa);
        }
      } else {
        uncle = grandpa.rbLeft;
        if (uncle && uncle.rbRed) {
          parent.rbRed = false;
          uncle.rbRed = false;
          grandpa.rbRed = true;
          node = grandpa;
        } else {
          if (node === parent.rbLeft) {
            this.rbRotateRight(parent);
            node = parent;
            parent = node.rbParent;
          }
          parent.rbRed = false;
          grandpa.rbRed = true;
          this.rbRotateLeft(grandpa);
        }
      }
      parent = node.rbParent;
    }
    this.root.rbRed = false;
  }

  rbRemoveNode(node) {
    // >>> rhill 2011-05-27: Performance: cache previous/next nodes
    if (node.rbNext) {
      node.rbNext.rbPrevious = node.rbPrevious;
    }
    if (node.rbPrevious) {
      node.rbPrevious.rbNext = node.rbNext;
    }
    node.rbNext = null;
    node.rbPrevious = null;
    // <<<
    let parent = node.rbParent;
    const left = node.rbLeft;
    const right = node.rbRight;
    let next;
    if (!left) {
      next = right;
    } else if (!right) {
      next = left;
    } else {
      next = RBTree.getFirstOfTree(right);
    }
    if (parent) {
      if (parent.rbLeft === node) {
        parent.rbLeft = next;
      } else {
        parent.rbRight = next;
      }
    } else {
      this.root = next;
    }
    // enforce red-black rules
    let isRed;
    if (left && right) {
      isRed = next.rbRed;
      next.rbRed = node.rbRed;
      next.rbLeft = left;
      left.rbParent = next;
      if (next !== right) {
        parent = next.rbParent;
        next.rbParent = node.rbParent;
        node = next.rbRight;
        parent.rbLeft = node;
        next.rbRight = right;
        right.rbParent = next;
      } else {
        next.rbParent = parent;
        parent = next;
        node = next.rbRight;
      }
    } else {
      isRed = node.rbRed;
      node = next;
    }
    // 'node' is now the sole successor's child and 'parent' its
    // new parent (since the successor can have been moved)
    if (node) {
      node.rbParent = parent;
    }
    // the 'easy' cases
    if (isRed) {
      return;
    }
    if (node && node.rbRed) {
      node.rbRed = false;
      return;
    }
    // the other cases
    let sibling;
    do {
      if (node === this.root) {
        break;
      }
      if (node === parent.rbLeft) {
        sibling = parent.rbRight;
        if (sibling.rbRed) {
          sibling.rbRed = false;
          parent.rbRed = true;
          this.rbRotateLeft(parent);
          sibling = parent.rbRight;
        }
        if (
          (sibling.rbLeft && sibling.rbLeft.rbRed) ||
          (sibling.rbRight && sibling.rbRight.rbRed)
        ) {
          if (!sibling.rbRight || !sibling.rbRight.rbRed) {
            sibling.rbLeft.rbRed = false;
            sibling.rbRed = true;
            this.rbRotateRight(sibling);
            sibling = parent.rbRight;
          }
          sibling.rbRed = parent.rbRed;
          parent.rbRed = false;
          sibling.rbRight.rbRed = false;
          this.rbRotateLeft(parent);
          node = this.root;
          break;
        }
      } else {
        sibling = parent.rbLeft;
        if (sibling.rbRed) {
          sibling.rbRed = false;
          parent.rbRed = true;
          this.rbRotateRight(parent);
          sibling = parent.rbLeft;
        }
        if (
          (sibling.rbLeft && sibling.rbLeft.rbRed) ||
          (sibling.rbRight && sibling.rbRight.rbRed)
        ) {
          if (!sibling.rbLeft || !sibling.rbLeft.rbRed) {
            sibling.rbRight.rbRed = false;
            sibling.rbRed = true;
            this.rbRotateLeft(sibling);
            sibling = parent.rbLeft;
          }
          sibling.rbRed = parent.rbRed;
          parent.rbRed = false;
          sibling.rbLeft.rbRed = false;
          this.rbRotateRight(parent);
          node = this.root;
          break;
        }
      }
      sibling.rbRed = true;
      node = parent;
      parent = parent.rbParent;
    } while (!node.rbRed);
    if (node) {
      node.rbRed = false;
    }
  }

  rbRotateLeft(node) {
    const p = node;
    const q = node.rbRight; // can't be null
    const parent = p.rbParent;
    if (parent) {
      if (parent.rbLeft === p) {
        parent.rbLeft = q;
      } else {
        parent.rbRight = q;
      }
    } else {
      this.root = q;
    }
    q.rbParent = parent;
    p.rbParent = q;
    p.rbRight = q.rbLeft;
    if (p.rbRight) {
      p.rbRight.rbParent = p;
    }
    q.rbLeft = p;
  }

  rbRotateRight(node) {
    const p = node;
    const q = node.rbLeft; // can't be null
    const parent = p.rbParent;
    if (parent) {
      if (parent.rbLeft === p) {
        parent.rbLeft = q;
      } else {
        parent.rbRight = q;
      }
    } else {
      this.root = q;
    }
    q.rbParent = parent;
    p.rbParent = q;
    p.rbLeft = q.rbRight;
    if (p.rbLeft) {
      p.rbLeft.rbParent = p;
    }
    q.rbRight = p;
  }

  getFirst() {
    return RBTree.getFirstOfTree(this.root);
  }

  static getFirstOfTree(rootNode) {
    let maybeFirst = rootNode;
    while (maybeFirst.rbLeft) {
      maybeFirst = maybeFirst.rbLeft;
    }
    return maybeFirst;
  }
}

class Vertex {
  constructor(x, y) {
    this.x = x;
    this.y = y;
  }
}

/* eslint-disable no-plusplus, no-param-reassign */

const ε = 1e-9;

function equalWithEpsilon(a, b) {
  return Math.abs(a - b) < ε;
}

function greaterThanWithEpsilon(a, b) {
  return a - b > ε;
}

function lessThanWithEpsilon(a, b) {
  return b - a > ε;
}

class Voronoi {
  constructor() {
    this.vertices = null;
    this.edges = null;
    this.cells = null;
    this.toRecycle = null;
    this.beachsectionJunkyard = [];
    this.circleEventJunkyard = [];
    this.vertexJunkyard = [];
    this.edgeJunkyard = [];
  }

  reset() {
    if (!this.beachline) {
      this.beachline = new RBTree();
    }
    // Move leftover beachsections to the beachsection junkyard.
    if (this.beachline.root) {
      let beachsection = this.beachline.getFirst();
      while (beachsection) {
        this.beachsectionJunkyard.push(beachsection); // mark for reuse
        beachsection = beachsection.rbNext;
      }
    }
    this.beachline.root = null;
    if (!this.circleEvents) {
      this.circleEvents = new RBTree();
    }
    this.circleEvents.root = null;
    this.firstCircleEvent = null;
    this.vertices = [];
    this.edges = [];
    this.cells = [];
  }

  /**
   * This creates and adds a vertex to the internal collection.
   */
  createVertex(x, y) {
    let v = this.vertexJunkyard.pop();
    if (!v) {
      v = new Vertex(x, y);
    } else {
      v.x = x;
      v.y = y;
    }
    this.vertices.push(v);
    return v;
  }

  /**
   * This creates and adds an edge to internal collection, and also creates two
   * halfedges which are added to each site's counterclockwise array of
   * halfedges.
   */
  createEdge(lSite, rSite, va, vb) {
    let edge = this.edgeJunkyard.pop();
    if (!edge) {
      edge = new Edge(lSite, rSite);
    } else {
      edge.lSite = lSite;
      edge.rSite = rSite;
      edge.va = null;
      edge.vb = null;
    }

    this.edges.push(edge);
    if (va) {
      edge.setEdgeStartpoint(lSite, rSite, va);
    }
    if (vb) {
      edge.setEdgeStartpoint(rSite, lSite, vb);
    }
    this.cells[lSite.voronoiId].halfedges.push(
      new Halfedge(edge, lSite, rSite)
    );
    this.cells[rSite.voronoiId].halfedges.push(
      new Halfedge(edge, rSite, lSite)
    );
    return edge;
  }

  createBorderEdge(lSite, va, vb) {
    let edge = this.edgeJunkyard.pop();
    if (!edge) {
      edge = new Edge(lSite, null);
    } else {
      edge.lSite = lSite;
      edge.rSite = null;
    }
    edge.va = va;
    edge.vb = vb;
    this.edges.push(edge);
    return edge;
  }

  // rhill 2011-06-02: A lot of Beachsection instanciations
  // occur during the computation of the Voronoi diagram,
  // somewhere between the number of sites and twice the
  // number of sites, while the number of Beachsections on the
  // beachline at any given time is comparatively low. For this
  // reason, we reuse already created Beachsections, in order
  // to avoid new memory allocation. This resulted in a measurable
  // performance gain.
  createBeachsection(site) {
    let beachsection = this.beachsectionJunkyard.pop();
    if (!beachsection) {
      beachsection = new Beachsection();
    }
    beachsection.site = site;
    return beachsection;
  }

  // calculate the left break point of a particular beach section,
  // given a particular sweep line
  static leftBreakPoint(arc, directrix) {
    // http://en.wikipedia.org/wiki/Parabola
    // http://en.wikipedia.org/wiki/Quadratic_equation
    // h1 = x1,
    // k1 = (y1+directrix)/2,
    // h2 = x2,
    // k2 = (y2+directrix)/2,
    // p1 = k1-directrix,
    // a1 = 1/(4*p1),
    // b1 = -h1/(2*p1),
    // c1 = h1*h1/(4*p1)+k1,
    // p2 = k2-directrix,
    // a2 = 1/(4*p2),
    // b2 = -h2/(2*p2),
    // c2 = h2*h2/(4*p2)+k2,
    // x = (-(b2-b1) + Math.sqrt((b2-b1)*(b2-b1) - 4*(a2-a1)*(c2-c1))) / (2*(a2-a1))
    // When x1 become the x-origin:
    // h1 = 0,
    // k1 = (y1+directrix)/2,
    // h2 = x2-x1,
    // k2 = (y2+directrix)/2,
    // p1 = k1-directrix,
    // a1 = 1/(4*p1),
    // b1 = 0,
    // c1 = k1,
    // p2 = k2-directrix,
    // a2 = 1/(4*p2),
    // b2 = -h2/(2*p2),
    // c2 = h2*h2/(4*p2)+k2,
    // x = (-b2 + Math.sqrt(b2*b2 - 4*(a2-a1)*(c2-k1))) / (2*(a2-a1)) + x1

    // change code below at your own risk: care has been taken to
    // reduce errors due to computers' finite arithmetic precision.
    // Maybe can still be improved, will see if any more of this
    // kind of errors pop up again.
    let { site } = arc;
    const rfocx = site.x;
    const rfocy = site.y;
    const pby2 = rfocy - directrix;
    // parabola in degenerate case where focus is on directrix
    if (!pby2) {
      return rfocx;
    }
    const lArc = arc.rbPrevious;
    if (!lArc) {
      return -Infinity;
    }
    site = lArc.site;
    const lfocx = site.x;
    const lfocy = site.y;
    const plby2 = lfocy - directrix;
    // parabola in degenerate case where focus is on directrix
    if (!plby2) {
      return lfocx;
    }
    const hl = lfocx - rfocx;
    const aby2 = 1 / pby2 - 1 / plby2;
    const b = hl / plby2;
    if (aby2) {
      return (
        (-b +
          Math.sqrt(
            b * b -
              2 *
                aby2 *
                ((hl * hl) / (-2 * plby2) -
                  lfocy +
                  plby2 / 2 +
                  rfocy -
                  pby2 / 2)
          )) /
          aby2 +
        rfocx
      );
    }
    // both parabolas have same distance to directrix, thus break point is midway
    return (rfocx + lfocx) / 2;
  }

  // calculate the right break point of a particular beach section,
  // given a particular directrix
  static rightBreakPoint(arc, directrix) {
    const rArc = arc.rbNext;
    if (rArc) {
      return Voronoi.leftBreakPoint(rArc, directrix);
    }
    const { site } = arc;
    return site.y === directrix ? site.x : Infinity;
  }

  detachBeachsection(beachsection) {
    this.detachCircleEvent(beachsection); // detach potentially attached circle event
    this.beachline.rbRemoveNode(beachsection); // remove from RB-tree
    this.beachsectionJunkyard.push(beachsection); // mark for reuse
  }

  removeBeachsection(beachsection) {
    const circle = beachsection.circleEvent;
    const { x } = circle;
    const y = circle.ycenter;
    const vertex = this.createVertex(x, y);
    let previous = beachsection.rbPrevious;
    let next = beachsection.rbNext;
    const disappearingTransitions = [beachsection];

    // remove collapsed beachsection from beachline
    this.detachBeachsection(beachsection);

    // there could be more than one empty arc at the deletion point, this
    // happens when more than two edges are linked by the same vertex,
    // so we will collect all those edges by looking up both sides of
    // the deletion point.
    // by the way, there is *always* a predecessor/successor to any collapsed
    // beach section, it's just impossible to have a collapsing first/last
    // beach sections on the beachline, since they obviously are unconstrained
    // on their left/right side.

    // look left
    let lArc = previous;
    while (
      lArc.circleEvent &&
      Math.abs(x - lArc.circleEvent.x) < ε &&
      Math.abs(y - lArc.circleEvent.ycenter) < ε
    ) {
      previous = lArc.rbPrevious;
      disappearingTransitions.unshift(lArc);
      this.detachBeachsection(lArc); // mark for reuse
      lArc = previous;
    }
    // even though it is not disappearing, I will also add the beach section
    // immediately to the left of the left-most collapsed beach section, for
    // convenience, since we need to refer to it later as this beach section
    // is the 'left' site of an edge for which a start point is set.
    disappearingTransitions.unshift(lArc);
    this.detachCircleEvent(lArc);

    // look right
    let rArc = next;
    while (
      rArc.circleEvent &&
      Math.abs(x - rArc.circleEvent.x) < ε &&
      Math.abs(y - rArc.circleEvent.ycenter) < ε
    ) {
      next = rArc.rbNext;
      disappearingTransitions.push(rArc);
      this.detachBeachsection(rArc); // mark for reuse
      rArc = next;
    }
    // we also have to add the beach section immediately to the right of the
    // right-most collapsed beach section, since there is also a disappearing
    // transition representing an edge's start point on its left.
    disappearingTransitions.push(rArc);
    this.detachCircleEvent(rArc);

    // walk through all the disappearing transitions between beach sections and
    // set the start point of their (implied) edge.
    const nArcs = disappearingTransitions.length;
    for (let iArc = 1; iArc < nArcs; iArc++) {
      rArc = disappearingTransitions[iArc];
      lArc = disappearingTransitions[iArc - 1];
      rArc.edge.setEdgeStartpoint(lArc.site, rArc.site, vertex);
    }

    // create a new edge as we have now a new transition between
    // two beach sections which were previously not adjacent.
    // since this edge appears as a new vertex is defined, the vertex
    // actually define an end point of the edge (relative to the site
    // on the left)
    // eslint-disable-next-line prefer-destructuring
    lArc = disappearingTransitions[0];
    rArc = disappearingTransitions[nArcs - 1];
    rArc.edge = this.createEdge(lArc.site, rArc.site, undefined, vertex);

    // create circle events if any for beach sections left in the beachline
    // adjacent to collapsed sections
    this.attachCircleEvent(lArc);
    this.attachCircleEvent(rArc);
  }

  addBeachsection(site) {
    const { x } = site;
    const directrix = site.y;

    // find the left and right beach sections which will surround the newly
    // created beach section.
    // rhill 2011-06-01: This loop is one of the most often executed,
    // hence we expand in-place the comparison-against-epsilon calls.
    let lArc;
    let rArc;
    let dxl;
    let dxr;
    let node = this.beachline.root;

    while (node) {
      dxl = Voronoi.leftBreakPoint(node, directrix) - x;
      // x lessThanWithEpsilon xl => falls somewhere before the left edge of the beachsection
      if (dxl > ε) {
        // this case should never happen
        // if (!node.rbLeft) {
        //    rArc = node.rbLeft;
        //    break;
        //    }
        node = node.rbLeft;
      } else {
        dxr = x - Voronoi.rightBreakPoint(node, directrix);
        // x greaterThanWithEpsilon xr => falls somewhere after the right edge of the beachsection
        if (dxr > ε) {
          if (!node.rbRight) {
            lArc = node;
            break;
          }
          node = node.rbRight;
        } else {
          // x equalWithEpsilon xl => falls exactly on the left edge of the beachsection
          if (dxl > -ε) {
            lArc = node.rbPrevious;
            rArc = node;
          }
          // x equalWithEpsilon xr => falls exactly on the right edge of the beachsection
          else if (dxr > -ε) {
            lArc = node;
            rArc = node.rbNext;
          }
          // falls exactly somewhere in the middle of the beachsection
          else {
            lArc = node;
            rArc = node;
          }
          break;
        }
      }
    }
    // at this point, keep in mind that lArc and/or rArc could be
    // undefined or null.

    // create a new beach section object for the site and add it to RB-tree
    const newArc = this.createBeachsection(site);
    this.beachline.rbInsertSuccessor(lArc, newArc);

    // cases:
    //

    // [null,null]
    // least likely case: new beach section is the first beach section on the
    // beachline.
    // This case means:
    //   no new transition appears
    //   no collapsing beach section
    //   new beachsection become root of the RB-tree
    if (!lArc && !rArc) {
      return;
    }

    // [lArc,rArc] where lArc == rArc
    // most likely case: new beach section split an existing beach
    // section.
    // This case means:
    //   one new transition appears
    //   the left and right beach section might be collapsing as a result
    //   two new nodes added to the RB-tree
    if (lArc === rArc) {
      // invalidate circle event of split beach section
      this.detachCircleEvent(lArc);

      // split the beach section into two separate beach sections
      rArc = this.createBeachsection(lArc.site);
      this.beachline.rbInsertSuccessor(newArc, rArc);

      // since we have a new transition between two beach sections,
      // a new edge is born
      newArc.edge = this.createEdge(lArc.site, newArc.site);
      rArc.edge = newArc.edge;

      // check whether the left and right beach sections are collapsing
      // and if so create circle events, to be notified when the point of
      // collapse is reached.
      this.attachCircleEvent(lArc);
      this.attachCircleEvent(rArc);
      return;
    }

    // [lArc,null]
    // even less likely case: new beach section is the *last* beach section
    // on the beachline -- this can happen *only* if *all* the previous beach
    // sections currently on the beachline share the same y value as
    // the new beach section.
    // This case means:
    //   one new transition appears
    //   no collapsing beach section as a result
    //   new beach section become right-most node of the RB-tree
    if (lArc && !rArc) {
      newArc.edge = this.createEdge(lArc.site, newArc.site);
      return;
    }

    // [null,rArc]
    // impossible case: because sites are strictly processed from top to bottom,
    // and left to right, which guarantees that there will always be a beach section
    // on the left -- except of course when there are no beach section at all on
    // the beach line, which case was handled above.
    // rhill 2011-06-02: No point testing in non-debug version
    // if (!lArc && rArc) {
    //    throw "Voronoi.addBeachsection(): What is this I don't even";
    //    }

    // [lArc,rArc] where lArc != rArc
    // somewhat less likely case: new beach section falls *exactly* in between two
    // existing beach sections
    // This case means:
    //   one transition disappears
    //   two new transitions appear
    //   the left and right beach section might be collapsing as a result
    //   only one new node added to the RB-tree
    if (lArc !== rArc) {
      // invalidate circle events of left and right sites
      this.detachCircleEvent(lArc);
      this.detachCircleEvent(rArc);

      // an existing transition disappears, meaning a vertex is defined at
      // the disappearance point.
      // since the disappearance is caused by the new beachsection, the
      // vertex is at the center of the circumscribed circle of the left,
      // new and right beachsections.
      // http://mathforum.org/library/drmath/view/55002.html
      // Except that I bring the origin at A to simplify
      // calculation
      const lSite = lArc.site;
      const ax = lSite.x;
      const ay = lSite.y;
      const bx = site.x - ax;
      const by = site.y - ay;
      const rSite = rArc.site;
      const cx = rSite.x - ax;
      const cy = rSite.y - ay;
      const d = 2 * (bx * cy - by * cx);
      const hb = bx * bx + by * by;
      const hc = cx * cx + cy * cy;
      const vertex = this.createVertex(
        (cy * hb - by * hc) / d + ax,
        (bx * hc - cx * hb) / d + ay
      );

      // one transition disappear
      rArc.edge.setEdgeStartpoint(lSite, rSite, vertex);

      // two new transitions appear at the new vertex location
      newArc.edge = this.createEdge(lSite, site, undefined, vertex);
      rArc.edge = this.createEdge(site, rSite, undefined, vertex);

      // check whether the left and right beach sections are collapsing
      // and if so create circle events, to handle the point of collapse.
      this.attachCircleEvent(lArc);
      this.attachCircleEvent(rArc);
    }
  }

  attachCircleEvent(arc) {
    const lArc = arc.rbPrevious;
    const rArc = arc.rbNext;
    if (!lArc || !rArc) {
      return;
    } // does that ever happen?
    const lSite = lArc.site;
    const cSite = arc.site;
    const rSite = rArc.site;

    // If site of left beachsection is same as site of
    // right beachsection, there can't be convergence
    if (lSite === rSite) {
      return;
    }

    // Find the circumscribed circle for the three sites associated
    // with the beachsection triplet.
    // rhill 2011-05-26: It is more efficient to calculate in-place
    // rather than getting the resulting circumscribed circle from an
    // object returned by calling Voronoi.circumcircle()
    // http://mathforum.org/library/drmath/view/55002.html
    // Except that I bring the origin at cSite to simplify calculations.
    // The bottom-most part of the circumcircle is our Fortune 'circle
    // event', and its center is a vertex potentially part of the final
    // Voronoi diagram.
    const bx = cSite.x;
    const by = cSite.y;
    const ax = lSite.x - bx;
    const ay = lSite.y - by;
    const cx = rSite.x - bx;
    const cy = rSite.y - by;

    // If points l->c->r are clockwise, then center beach section does not
    // collapse, hence it can't end up as a vertex (we reuse 'd' here, which
    // sign is reverse of the orientation, hence we reverse the test.
    // http://en.wikipedia.org/wiki/Curve_orientation#Orientation_of_a_simple_polygon
    // rhill 2011-05-21: Nasty finite precision error which caused circumcircle() to
    // return infinites: 1e-12 seems to fix the problem.
    const d = 2 * (ax * cy - ay * cx);
    if (d >= -2e-12) {
      return;
    }

    const ha = ax * ax + ay * ay;
    const hc = cx * cx + cy * cy;
    const x = (cy * ha - ay * hc) / d;
    const y = (ax * hc - cx * ha) / d;
    const ycenter = y + by;

    // Important: ybottom should always be under or at sweep, so no need
    // to waste CPU cycles by checking

    // recycle circle event object if possible
    let circleEvent = this.circleEventJunkyard.pop();
    if (!circleEvent) {
      circleEvent = new CircleEvent();
    }
    circleEvent.arc = arc;
    circleEvent.site = cSite;
    circleEvent.x = x + bx;
    circleEvent.y = ycenter + Math.sqrt(x * x + y * y); // y bottom
    circleEvent.ycenter = ycenter;
    arc.circleEvent = circleEvent;

    // find insertion point in RB-tree: circle events are ordered from
    // smallest to largest
    let predecessor = null;
    let node = this.circleEvents.root;
    while (node) {
      if (
        circleEvent.y < node.y ||
        (circleEvent.y === node.y && circleEvent.x <= node.x)
      ) {
        if (node.rbLeft) {
          node = node.rbLeft;
        } else {
          predecessor = node.rbPrevious;
          break;
        }
      } else if (node.rbRight) {
        node = node.rbRight;
      } else {
        predecessor = node;
        break;
      }
    }
    this.circleEvents.rbInsertSuccessor(predecessor, circleEvent);
    if (!predecessor) {
      this.firstCircleEvent = circleEvent;
    }
  }

  detachCircleEvent(arc) {
    const { circleEvent } = arc;
    if (circleEvent) {
      if (!circleEvent.rbPrevious) {
        this.firstCircleEvent = circleEvent.rbNext;
      }
      this.circleEvents.rbRemoveNode(circleEvent); // remove from RB-tree
      this.circleEventJunkyard.push(circleEvent);
      arc.circleEvent = null;
    }
  }

  // ---------------------------------------------------------------------------
  // Diagram completion methods

  // connect dangling edges (not if a cursory test tells us
  // it is not going to be visible.
  // return value:
  //   false: the dangling endpoint couldn't be connected
  //   true: the dangling endpoint could be connected
  connectEdge(edge, bbox) {
    // skip if end point already connected
    let { vb } = edge;
    if (vb) {
      return true;
    }

    // make local copy for performance purpose
    let { va } = edge;
    const { xl } = bbox;
    const { xr } = bbox;
    const { yt } = bbox;
    const { yb } = bbox;
    const { lSite } = edge;
    const { rSite } = edge;
    const lx = lSite.x;
    const ly = lSite.y;
    const rx = rSite.x;
    const ry = rSite.y;
    const fx = (lx + rx) / 2;
    const fy = (ly + ry) / 2;
    let fm;
    let fb;

    // if we reach here, this means cells which use this edge will need
    // to be closed, whether because the edge was removed, or because it
    // was connected to the bounding box.
    this.cells[lSite.voronoiId].closeMe = true;
    this.cells[rSite.voronoiId].closeMe = true;

    // get the line equation of the bisector if line is not vertical
    if (ry !== ly) {
      fm = (lx - rx) / (ry - ly);
      fb = fy - fm * fx;
    }

    // remember, direction of line (relative to left site):
    // upward: left.x < right.x
    // downward: left.x > right.x
    // horizontal: left.x == right.x
    // upward: left.x < right.x
    // rightward: left.y < right.y
    // leftward: left.y > right.y
    // vertical: left.y == right.y

    // depending on the direction, find the best side of the
    // bounding box to use to determine a reasonable start point

    // rhill 2013-12-02:
    // While at it, since we have the values which define the line,
    // clip the end of va if it is outside the bbox.
    // https://github.com/gorhill/Javascript-Voronoi/issues/15
    // TODO: Do all the clipping here rather than rely on Liang-Barsky
    // which does not do well sometimes due to loss of arithmetic
    // precision. The code here doesn't degrade if one of the vertex is
    // at a huge distance.

    // special case: vertical line
    if (fm === undefined) {
      // doesn't intersect with viewport
      if (fx < xl || fx >= xr) {
        return false;
      }
      // downward
      if (lx > rx) {
        if (!va || va.y < yt) {
          va = this.createVertex(fx, yt);
        } else if (va.y >= yb) {
          return false;
        }
        vb = this.createVertex(fx, yb);
      }
      // upward
      else {
        if (!va || va.y > yb) {
          va = this.createVertex(fx, yb);
        } else if (va.y < yt) {
          return false;
        }
        vb = this.createVertex(fx, yt);
      }
    }
    // closer to vertical than horizontal, connect start point to the
    // top or bottom side of the bounding box
    else if (fm < -1 || fm > 1) {
      // downward
      if (lx > rx) {
        if (!va || va.y < yt) {
          va = this.createVertex((yt - fb) / fm, yt);
        } else if (va.y >= yb) {
          return false;
        }
        vb = this.createVertex((yb - fb) / fm, yb);
      }
      // upward
      else {
        if (!va || va.y > yb) {
          va = this.createVertex((yb - fb) / fm, yb);
        } else if (va.y < yt) {
          return false;
        }
        vb = this.createVertex((yt - fb) / fm, yt);
      }
    }
    // closer to horizontal than vertical, connect start point to the
    // left or right side of the bounding box
    else {
      // rightward
      // eslint-disable-next-line no-lonely-if
      if (ly < ry) {
        if (!va || va.x < xl) {
          va = this.createVertex(xl, fm * xl + fb);
        } else if (va.x >= xr) {
          return false;
        }
        vb = this.createVertex(xr, fm * xr + fb);
      }
      // leftward
      else {
        if (!va || va.x > xr) {
          va = this.createVertex(xr, fm * xr + fb);
        } else if (va.x < xl) {
          return false;
        }
        vb = this.createVertex(xl, fm * xl + fb);
      }
    }
    edge.va = va;
    edge.vb = vb;

    return true;
  }

  // line-clipping code taken from:
  //   Liang-Barsky function by Daniel White
  //   http://www.skytopia.com/project/articles/compsci/clipping.html
  // Thanks!
  // A bit modified to minimize code paths
  clipEdge(edge, bbox) {
    const ax = edge.va.x;
    const ay = edge.va.y;
    const bx = edge.vb.x;
    const by = edge.vb.y;
    let t0 = 0;
    let t1 = 1;
    const dx = bx - ax;
    const dy = by - ay;
    // left
    let q = ax - bbox.xl;
    if (dx === 0 && q < 0) {
      return false;
    }
    let r = -q / dx;
    if (dx < 0) {
      if (r < t0) {
        return false;
      }
      if (r < t1) {
        t1 = r;
      }
    } else if (dx > 0) {
      if (r > t1) {
        return false;
      }
      if (r > t0) {
        t0 = r;
      }
    }
    // right
    q = bbox.xr - ax;
    if (dx === 0 && q < 0) {
      return false;
    }
    r = q / dx;
    if (dx < 0) {
      if (r > t1) {
        return false;
      }
      if (r > t0) {
        t0 = r;
      }
    } else if (dx > 0) {
      if (r < t0) {
        return false;
      }
      if (r < t1) {
        t1 = r;
      }
    }
    // top
    q = ay - bbox.yt;
    if (dy === 0 && q < 0) {
      return false;
    }
    r = -q / dy;
    if (dy < 0) {
      if (r < t0) {
        return false;
      }
      if (r < t1) {
        t1 = r;
      }
    } else if (dy > 0) {
      if (r > t1) {
        return false;
      }
      if (r > t0) {
        t0 = r;
      }
    }
    // bottom
    q = bbox.yb - ay;
    if (dy === 0 && q < 0) {
      return false;
    }
    r = q / dy;
    if (dy < 0) {
      if (r > t1) {
        return false;
      }
      if (r > t0) {
        t0 = r;
      }
    } else if (dy > 0) {
      if (r < t0) {
        return false;
      }
      if (r < t1) {
        t1 = r;
      }
    }

    // if we reach this point, Voronoi edge is within bbox

    // if t0 > 0, va needs to change
    // rhill 2011-06-03: we need to create a new vertex rather
    // than modifying the existing one, since the existing
    // one is likely shared with at least another edge
    if (t0 > 0) {
      edge.va = this.createVertex(ax + t0 * dx, ay + t0 * dy);
    }

    // if t1 < 1, vb needs to change
    // rhill 2011-06-03: we need to create a new vertex rather
    // than modifying the existing one, since the existing
    // one is likely shared with at least another edge
    if (t1 < 1) {
      edge.vb = this.createVertex(ax + t1 * dx, ay + t1 * dy);
    }

    // va and/or vb were clipped, thus we will need to close
    // cells which use this edge.
    if (t0 > 0 || t1 < 1) {
      this.cells[edge.lSite.voronoiId].closeMe = true;
      this.cells[edge.rSite.voronoiId].closeMe = true;
    }

    return true;
  }

  // Connect/cut edges at bounding box
  clipEdges(bbox) {
    // connect all dangling edges to bounding box
    // or get rid of them if it can't be done
    const { edges } = this;
    let iEdge = edges.length;
    let edge;

    // iterate backward so we can splice safely
    while (iEdge--) {
      edge = edges[iEdge];
      // edge is removed if:
      //   it is wholly outside the bounding box
      //   it is looking more like a point than a line
      if (
        !this.connectEdge(edge, bbox) ||
        !this.clipEdge(edge, bbox) ||
        (Math.abs(edge.va.x - edge.vb.x) < ε &&
          Math.abs(edge.va.y - edge.vb.y) < ε)
      ) {
        edge.va = null;
        edge.vb = null;
        edges.splice(iEdge, 1);
      }
    }
  }

  // Close the cells.
  // The cells are bound by the supplied bounding box.
  // Each cell refers to its associated site, and a list
  // of halfedges ordered counterclockwise.
  closeCells(bbox) {
    const { xl } = bbox;
    const { xr } = bbox;
    const { yt } = bbox;
    const { yb } = bbox;
    const { cells } = this;
    let iCell = cells.length;
    let cell;
    let iLeft;
    let halfedges;
    let nHalfedges;
    let edge;
    let va;
    let vb;
    let vz;
    let lastBorderSegment;

    while (iCell--) {
      cell = cells[iCell];
      // prune, order halfedges counterclockwise, then add missing ones
      // required to close cells
      if (!cell.prepareHalfedges()) {
        continue;
      }
      if (!cell.closeMe) {
        continue;
      }
      // find first 'unclosed' point.
      // an 'unclosed' point will be the end point of a halfedge which
      // does not match the start point of the following halfedge
      halfedges = cell.halfedges;
      nHalfedges = halfedges.length;
      // special case: only one site, in which case, the viewport is the cell
      // ...

      // all other cases
      iLeft = 0;
      while (iLeft < nHalfedges) {
        va = halfedges[iLeft].getEndpoint();
        vz = halfedges[(iLeft + 1) % nHalfedges].getStartpoint();
        // if end point is not equal to start point, we need to add the missing
        // halfedge(s) up to vz
        if (Math.abs(va.x - vz.x) >= ε || Math.abs(va.y - vz.y) >= ε) {
          // rhill 2013-12-02:
          // "Holes" in the halfedges are not necessarily always adjacent.
          // https://github.com/gorhill/Javascript-Voronoi/issues/16

          // find entry point:
          switch (true) {
            case equalWithEpsilon(va.x, xl) && lessThanWithEpsilon(va.y, yb):
              // walk downward along left side
              lastBorderSegment = equalWithEpsilon(vz.x, xl);
              vb = this.createVertex(xl, lastBorderSegment ? vz.y : yb);
              edge = this.createBorderEdge(cell.site, va, vb);
              iLeft++;
              halfedges.splice(iLeft, 0, new Halfedge(edge, cell.site, null));
              nHalfedges++;
              if (lastBorderSegment) {
                break;
              }
              va = vb;
            // fall through

            case equalWithEpsilon(va.y, yb) && lessThanWithEpsilon(va.x, xr):
              // walk rightward along bottom side
              lastBorderSegment = equalWithEpsilon(vz.y, yb);
              vb = this.createVertex(lastBorderSegment ? vz.x : xr, yb);
              edge = this.createBorderEdge(cell.site, va, vb);
              iLeft++;
              halfedges.splice(iLeft, 0, new Halfedge(edge, cell.site, null));
              nHalfedges++;
              if (lastBorderSegment) {
                break;
              }
              va = vb;
            // fall through

            case equalWithEpsilon(va.x, xr) && greaterThanWithEpsilon(va.y, yt):
              // walk upward along right side
              lastBorderSegment = equalWithEpsilon(vz.x, xr);
              vb = this.createVertex(xr, lastBorderSegment ? vz.y : yt);
              edge = this.createBorderEdge(cell.site, va, vb);
              iLeft++;
              halfedges.splice(iLeft, 0, new Halfedge(edge, cell.site, null));
              nHalfedges++;
              if (lastBorderSegment) {
                break;
              }
              va = vb;
            // fall through

            case equalWithEpsilon(va.y, yt) && greaterThanWithEpsilon(va.x, xl):
              // walk leftward along top side
              lastBorderSegment = equalWithEpsilon(vz.y, yt);
              vb = this.createVertex(lastBorderSegment ? vz.x : xl, yt);
              edge = this.createBorderEdge(cell.site, va, vb);
              iLeft++;
              halfedges.splice(iLeft, 0, new Halfedge(edge, cell.site, null));
              nHalfedges++;
              if (lastBorderSegment) {
                break;
              }
              va = vb;
              // fall through

              // walk downward along left side
              lastBorderSegment = equalWithEpsilon(vz.x, xl);
              vb = this.createVertex(xl, lastBorderSegment ? vz.y : yb);
              edge = this.createBorderEdge(cell.site, va, vb);
              iLeft++;
              halfedges.splice(iLeft, 0, new Halfedge(edge, cell.site, null));
              nHalfedges++;
              if (lastBorderSegment) {
                break;
              }
              va = vb;
              // fall through

              // walk rightward along bottom side
              lastBorderSegment = equalWithEpsilon(vz.y, yb);
              vb = this.createVertex(lastBorderSegment ? vz.x : xr, yb);
              edge = this.createBorderEdge(cell.site, va, vb);
              iLeft++;
              halfedges.splice(iLeft, 0, new Halfedge(edge, cell.site, null));
              nHalfedges++;
              if (lastBorderSegment) {
                break;
              }
              va = vb;
              // fall through

              // walk upward along right side
              lastBorderSegment = equalWithEpsilon(vz.x, xr);
              vb = this.createVertex(xr, lastBorderSegment ? vz.y : yt);
              edge = this.createBorderEdge(cell.site, va, vb);
              iLeft++;
              halfedges.splice(iLeft, 0, new Halfedge(edge, cell.site, null));
              nHalfedges++;
              if (lastBorderSegment) {
                break;
              }
            // fall through

            default:
              throw new Error("Voronoi.closeCells() > this makes no sense!");
          }
        }
        iLeft++;
      }
      cell.closeMe = false;
    }
  }

  /**
   * Debugging helper
   */
  dumpBeachline(y) {
    const { log } = console;
    log("Voronoi.dumpBeachline(%f) > Beachsections, from left to right:", y);
    if (!this.beachline) {
      log("  None");
    } else {
      let bs = this.beachline.getFirst();
      while (bs) {
        log(
          "  site %d: xl: %f, xr: %f",
          bs.site.voronoiId,
          Voronoi.leftBreakPoint(bs, y),
          Voronoi.rightBreakPoint(bs, y)
        );
        bs = bs.rbNext;
      }
    }
  }

  /**
   * Helper: Quantize sites
   *
   * rhill 2013-10-12:
   * This is to solve https://github.com/gorhill/Javascript-Voronoi/issues/15
   * Since not all users will end up using the kind of coordinate values which
   * would cause the issue to arise, I chose to let the user decide whether or
   * not they should sanitize his coord values through this helper. This way,
   * for users who uses coord values which are known to be fine, no overhead is
   * added.
   */
  static quantizeSites(sites) {
    sites.forEach(site => {
      site.x = Math.floor(site.x / ε) * ε;
      site.y = Math.floor(site.y / ε) * ε;
    });
  }

  // ---------------------------------------------------------------------------
  // Helper: Recycle diagram: all vertex, edge and cell objects are
  // "surrendered" to the Voronoi object for reuse.
  // TODO: rhill-voronoi-core v2: more performance to be gained
  // when I change the semantic of what is returned.

  recycle(diagram) {
    if (diagram) {
      if (diagram instanceof Diagram) {
        this.toRecycle = diagram;
      } else {
        throw new Error("Voronoi.recycleDiagram() > Need a Diagram object.");
      }
    }
  }

  // ---------------------------------------------------------------------------
  // Top-level Fortune loop

  // rhill 2011-05-19:
  //   Voronoi sites are kept client-side now, to allow
  //   user to freely modify content. At compute time,
  //   *references* to sites are copied locally.

  compute(sites, bbox) {
    // to measure execution time
    const startTime = new Date();

    // init internal state
    this.reset();

    // any diagram data available for recycling?
    // I do that here so that this is included in execution time
    if (this.toRecycle) {
      this.vertexJunkyard = this.vertexJunkyard.concat(this.toRecycle.vertices);
      this.edgeJunkyard = this.edgeJunkyard.concat(this.toRecycle.edges);
      this.toRecycle = null;
    }

    // Initialize site event queue
    const siteEvents = sites.slice(0);
    siteEvents.sort((a, b) => {
      const r = b.y - a.y;
      if (r) {
        return r;
      }
      return b.x - a.x;
    });

    // process queue
    let site = siteEvents.pop();
    let siteid = 0;
    let xsitex; // to avoid duplicate sites
    let xsitey;
    const { cells } = this;
    let circle;

    // main loop
    for (;;) {
      // we need to figure whether we handle a site or circle event
      // for this we find out if there is a site event and it is
      // 'earlier' than the circle event
      circle = this.firstCircleEvent;

      // add beach section
      if (
        site &&
        (!circle ||
          site.y < circle.y ||
          (site.y === circle.y && site.x < circle.x))
      ) {
        // only if site is not a duplicate
        if (site.x !== xsitex || site.y !== xsitey) {
          // first create cell for new site
          cells[siteid] = new Cell(site);
          site.voronoiId = siteid++;
          // then create a beachsection for that site
          this.addBeachsection(site);
          // remember last site coords to detect duplicate
          xsitey = site.y;
          xsitex = site.x;
        }
        site = siteEvents.pop();
      }

      // remove beach section
      else if (circle) {
        this.removeBeachsection(circle.arc);
      }

      // all done, quit
      else {
        break;
      }
    }

    // wrapping-up:
    //   connect dangling edges to bounding box
    //   cut edges as per bounding box
    //   discard edges completely outside bounding box
    //   discard edges which are point-like
    this.clipEdges(bbox);

    //   add missing edges in order to close opened cells
    this.closeCells(bbox);

    // to measure execution time
    const stopTime = new Date();

    // prepare return values
    const diagram = new Diagram();
    diagram.cells = this.cells;
    diagram.edges = this.edges;
    diagram.vertices = this.vertices;
    diagram.execTime = stopTime.getTime() - startTime.getTime();

    // clean up
    this.reset();

    return diagram;
  }
}

export default Voronoi;
