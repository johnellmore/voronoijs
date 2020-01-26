import seedrandom from "seedrandom";
import Voronoi from "./Voronoi";

test("Can make a diagram from an empty set of sites", () => {
  const v = new Voronoi();
  const boundingBox = { xl: 0, xr: 100, yt: 0, yb: 100 };
  const diagram = v.compute([], boundingBox);
  expect(diagram.prototype).toBe(Voronoi.Diagram);
});

test("Makes a correct diagram from two sites", () => {
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

test("Can make a diagram from fifty sites", () => {
  // generate random points
  const width = 600;
  const height = 400;
  const margin = 20;
  const boundingBox = { xl: 0, xr: width, yt: 0, yb: height };
  const effectiveWidth = width - margin * 2;
  const effectiveHeight = height - margin * 2;

  const rand = seedrandom("fity-sites-test");
  const sites = [];
  for (let i = 0; i < 50; i += 1) {
    sites.push({
      x: Math.round((margin + rand() * effectiveWidth) * 10) / 10,
      y: Math.round((margin + rand() * effectiveHeight) * 10) / 10
    });
  }

  const v = new Voronoi();
  const diagram = v.compute(sites, boundingBox);

  expect(diagram.prototype).toBe(Voronoi.Diagram);
  expect(diagram.cells).toHaveLength(50);
});
