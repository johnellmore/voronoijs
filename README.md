# voronoijs

A Javascript implementation of Steven J. Fortune's algorithm to efficiently compute Voronoi diagrams. The Voronoi object's purpose is to solely compute a Voronoi diagram. There are no production dependencies--only dev dependencies for development usage. It also contains no rendering code: that is left to the user of the library.

This implementation is forked from [Raymond Hill's excellent JS library](https://github.com/gorhill/Javascript-Voronoi). The goal with this fork is to support the modified version of the Fortune algorithm which can construct an additively weighted Voronoi diagram. Also, this refactors much of the implementation to use modern JavaScript syntax--ES modules, classes, etc.

## License

This software is released under the MIT License (like the original implementation).

## Further reading

- [Voronoi Diagram](https://en.wikipedia.org/wiki/Voronoi_diagram)
- [Fortune's Algorithm](https://en.wikipedia.org/wiki/Fortune's_algorithm)
- [Weighted Voronoi Algorithm](https://en.wikipedia.org/wiki/Weighted_Voronoi_diagram)
