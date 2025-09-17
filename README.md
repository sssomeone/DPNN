# DPNN: **Efficient Nearest Neighbor Search Using Dynamic Programming**

This repository contains the source code for our research paper **[Efficient Nearest Neighbor Search Using Dynamic Programming](https://pengfei.me/project_page/nns/)**

Abstract: Given a collection of points in R3, KD-Tree and R-Tree are well-known nearest neighbor search (NNS) algorithms that rely on spatial partitioning and indexing techniques. However, when the query point is far from the data points or the data points inherently represent a 2-manifold surface, their query performance may degrade. To address this, we propose a novel dynamic programming technique that precomputes a Directed Acyclic Graph (DAG) to encode the proximity structure between data points. More specifically, the DAG captures how the proximity structure evolves during the incremental construction of the Voronoi diagram of the data points. Experimental results demonstrate that our method achieves a speed increase of 1-10x. Furthermore, our algorithm demonstrates significant practical value in diverse applications. We validated its effectiveness through extensive testing in four key applications: Point-to-Mesh Distance Queries, Iterative Closest Point (ICP) Registration, Density Peak Clustering, and Point-to-Segments Distance Queries. A particularly notable feature of our approach is its unique ability to efficiently identify the nearest neighbor among the first k points in the point cloud, a capability that enables substantial acceleration in low-dimensional applications like Density Peak Clustering. As a natural extension of our incremental construction process, our method can also be readily adapted for farthest-point sampling tasks. These experimental results across multiple domains underscore the broad applicability and practical importance of our approach.

Project Page: https://pengfei.me/project_page/nns/
Doi: https://ieeexplore.ieee.org/document/11165079
## Dependencies

This project requires the following libraries:

- **Boost**
- **Eigen3**
- **CGAL**

## CGAL Kernel Configuration

The code supports two CGAL kernel types:

```cpp
typedef CGAL::Simple_cartesian<double> K;
// typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
```

By default, we use `Simple_cartesian<double>`. For certain datasets that may encounter numerical precision issues, you may need to:
- Switch to `Exact_predicates_inexact_constructions_kernel`, or
- Apply very slight random perturbations to the input data

## Usage

See `src/main.cpp` for usage examples and demonstrations of the algorithm.

## Citation

If you use this code in your research, please cite our paper:

```
@ARTICLE{11165079,
author={Wang, Pengfei and Song, Jiantao and Xin, Shiqing and Chen, Shuangmin and Tu, Changhe and Wang, Wenping and Wang, Jiaye},
journal={IEEE Transactions on Pattern Analysis and Machine Intelligence},
title={Efficient Nearest Neighbor Search Using Dynamic Programming},
year={2025},
volume={},
number={},
pages={1-16},
keywords={Delaunay triangulation;density peak clustering;farthest point sampling;nearest neighbor search;voronoi diagram},
doi={10.1109/TPAMI.2025.3610211}}
```

## License

MIT License