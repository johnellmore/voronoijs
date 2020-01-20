const Voronoi = require("./Voronoi");

test("Can make a diagram from an empty set of sites", () => {
  const v = new Voronoi();
  const boundingBox = { xl: 0, xr: 100, yt: 0, yb: 100 };
  const diagram = v.compute([], boundingBox);
  expect(diagram.prototype).toBe(Voronoi.Diagram);
});

test("Can make a diagram from two sites", () => {
  const v = new Voronoi();
  const boundingBox = { xl: 0, xr: 100, yt: 0, yb: 100 };
  const sites = [
    { x: 25, y: 25 },
    { x: 75, y: 75 }
  ];
  const diagram = v.compute(sites, boundingBox);
  expect(diagram.prototype).toBe(Voronoi.Diagram);
  expect(diagram.cells).toHaveLength(2);
  const dividingEdges = diagram.edges.filter(edge => edge.lSite && edge.rSite);
  expect(dividingEdges).toHaveLength(1);
  const dividingEdge = dividingEdges[0];
  expect(dividingEdge.va).toEqual({ x: 0, y: 100 });
  expect(dividingEdge.vb).toEqual({ x: 100, y: 0 });
});
