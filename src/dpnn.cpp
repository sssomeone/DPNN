//
// Created by pengfei on 25-7-16.
//
#include <vector>
#include <algorithm>
#include <random>

#include <queue>
#include <map>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Spatial_sort_traits_adapter_3.h>
#include <CGAL/draw_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_with_circumcenter_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/point_generators_3.h>


typedef CGAL::Simple_cartesian<double> K;
// typedef CGAL::Exact_predicates_inexact_constructions_kernel K;


typedef CGAL::Delaunay_triangulation_3<K>                   DT3;
typedef CGAL::Creator_uniform_3<double, K::Point_3>          Creator;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K> Points;


typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K>      Vb;
typedef CGAL::Delaunay_triangulation_cell_base_with_circumcenter_3<K> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb>                            Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds>                    Delaunay;
typedef K::Point_3 Point_3;
typedef typename Delaunay::Concurrency_tag                     Concurrency_tag;

using namespace std;
#ifndef SEGMENTTREE
#define SEGMENTTREE

class segmentTree {
public:
	struct node {
		double v;
		int maxPos;
	};
	vector<node>tree;
	vector<int> PointIndex2NodeIndex;
	vector<vector<int>> CellsPoints;
	vector<Eigen::Vector3d>& pointList;
public:

	segmentTree(vector<Eigen::Vector3d>& _pointlist, bool usage=true) :pointList(_pointlist) {
		if(usage) {
			tree.resize(4 * pointList.size());
			PointIndex2NodeIndex.resize(pointList.size());
			CellsPoints.resize(pointList.size());
			build(1, 0, pointList.size() - 1, pointList[0], pointList, tree);
			CellsPoints[0].resize(pointList.size() - 1);
			for (int i = 1; i < pointList.size(); ++i) {
				CellsPoints[0][i - 1] = i;
			}
		}
	}
	void build(int rot, int l, int r, const Eigen::Vector3d& q, const vector<Eigen::Vector3d>& ps, vector<node>& nodes) {
		if (l == r) {
			nodes[rot].v = (q - ps[l]).squaredNorm();
			nodes[rot].maxPos = l;
			PointIndex2NodeIndex[l] = rot;
			return;
		}
		int mid = (l + r) >> 1;
		build(rot * 2, l, mid, q, ps, nodes);
		build(rot * 2 + 1, mid + 1, r, q, ps, nodes);
		if (nodes[2 * rot].v > nodes[2 * rot + 1].v) {
			nodes[rot].v = nodes[2 * rot].v;
			nodes[rot].maxPos = nodes[2 * rot].maxPos;
		}
		else {
			nodes[rot].v = nodes[2 * rot + 1].v;
			nodes[rot].maxPos = nodes[2 * rot + 1].maxPos;
		}
	}
	int Insert(const int newIndex, vector<int>& nearPoints, int nearPointSize) {
		movePoints(nearPoints, newIndex, nearPointSize);
		return tree[1].maxPos;
	}

	void movePoints(const vector<int>& nearPoints, int newIndex, int nearPointSize) {
		vector<int> startIndex(nearPointSize);
		int size = 0;
		SetNodes(PointIndex2NodeIndex[newIndex], 0, tree);
		for (int i = 0; i < nearPointSize; ++i) {
			int l = 0, r = 0;
			while (r < CellsPoints[nearPoints[i]].size()) {
				double preDis = (pointList[CellsPoints[nearPoints[i]][r]] - pointList[nearPoints[i]]).squaredNorm();
				double newDis = (pointList[CellsPoints[nearPoints[i]][r]] - pointList[newIndex]).squaredNorm();
				if (preDis < newDis) {
					swap(CellsPoints[nearPoints[i]][l], CellsPoints[nearPoints[i]][r]);
					++l;
				}
				++r;
			}
			startIndex[i] = l;
			size += r - l;
		}
		CellsPoints[newIndex].resize(size);
		size = 0;
		for (int i = 0; i < nearPointSize; ++i) {
			for (int j = startIndex[i]; j < CellsPoints[nearPoints[i]].size(); ++j) {
				CellsPoints[newIndex][size++] = CellsPoints[nearPoints[i]][j];
				SetNodes(PointIndex2NodeIndex[CellsPoints[nearPoints[i]][j]], (pointList[CellsPoints[nearPoints[i]][j]] - pointList[newIndex]).squaredNorm(), tree);
			}
			CellsPoints[nearPoints[i]].resize(startIndex[i]);
		}
	}


	void SetNodes(int rot, double v, vector<node>& nodes) {
		nodes[rot].v = min(v, nodes[rot].v);
		while (rot != 1) {
			rot >>= 1;
			if (nodes[2 * rot].v > nodes[2 * rot + 1].v) {
				if (nodes[rot].maxPos == nodes[2 * rot].maxPos && abs(nodes[rot].v - nodes[2 * rot].v) < 1e-9) {
					break;
				}
				nodes[rot].v = nodes[2 * rot].v;
				nodes[rot].maxPos = nodes[2 * rot].maxPos;
			}
			else {
				if (nodes[rot].maxPos == nodes[2 * rot + 1].maxPos && abs(nodes[rot].v - nodes[2 * rot + 1].v) < 1e-9) {
					break;
				}
				nodes[rot].v = nodes[2 * rot + 1].v;
				nodes[rot].maxPos = nodes[2 * rot + 1].maxPos;
			}
		}
	}
};

#endif // !segmentTree

typedef Eigen::Vector3d Point;

using Point_with_info = std::pair<Point_3, unsigned>;
using Accessor      = CGAL::First_of_pair_property_map<Point_with_info>;
using Sort_traits_3 = CGAL::Spatial_sort_traits_adapter_3<K, Accessor>;


class NNTResultSet
{
public:
	using DistanceType = double;
	using IndexType = uint32_t;
	using CountType = uint32_t;

private:
	IndexType* indices;
	DistanceType* dists;
	CountType     capacity;
	CountType     count;

public:
	explicit NNTResultSet(CountType capacity_)
		: indices(nullptr), dists(nullptr), capacity(capacity_), count(0)
	{
	}

	inline void init(IndexType* indices_, DistanceType* dists_)
	{
		indices = indices_;
		dists = dists_;
		count = 0;
		if (capacity)
			dists[capacity - 1] = (std::numeric_limits<DistanceType>::max)();
	}

	CountType size() const { return count; }
	bool      empty() const { return count == 0; }
	bool      full() const { return count == capacity; }

	/**
	 * Called during search to add an element matching the criteria.
	 * @return true if the search should be continued, false if the results are
	 * sufficient
	 */
	int addPoint(DistanceType dist, IndexType index)
	{
		if (dist > dists[capacity - 1])
			return -1;
		int pre_id = indices[capacity - 1];
		CountType i;
		for (i = count; i > 0; --i)
		{
			/** If defined and two points have the same distance, the one with
			 *  the lowest-index will be returned first. */
			if (dists[i - 1] > dist)
			{
				if (i < capacity)
				{
					dists[i] = dists[i - 1];
					indices[i] = indices[i - 1];
				}
			}
			else
				break;
		}
		if (i < capacity)
		{
			dists[i] = dist;
			indices[i] = index;
		}
		else {
			return -1;
		}
		if (count < capacity) {
			count++;
			return -1;
		}
		return pre_id;
	}

	DistanceType worstDist() const { return dists[capacity - 1]; }

	void sort()
	{
		// already sorted
	}
};
class dpnn {
private:
	std::vector<Point> PointList;
	std::vector<std::vector<int>>              SonList;
	std::vector<int>						   ids;
	Delaunay dt3;
	segmentTree tree;
public:
	dpnn(std::vector<Point>& InputPointList, bool using_farthest_sorting = false) :PointList(InputPointList),tree(InputPointList,using_farthest_sorting) {
		ids.resize(InputPointList.size());

		if (using_farthest_sorting) {
			for (int i = 0; i < ids.size(); ++i)
				ids[i] = i;
			GetDelaunaySonListFarthest(PointList);
		}
		else {
			std::vector<Point_with_info> pts;
			pts.reserve(InputPointList.size());
			for (unsigned i = 0; i < InputPointList.size(); ++i){
				const auto& v = InputPointList[i];
				pts.emplace_back(Point_3(v.x(), v.y(), v.z()), i);
			}
			CGAL::spatial_sort< Delaunay::Concurrency_tag >(
					pts.begin(), pts.end(),
					Sort_traits_3()               // ← 告诉算法如何拿 Point_3
			);
			PointList.resize(pts.size());
			for(int i=0;i<InputPointList.size();++i) {
				PointList[i] = Point(pts[i].first.x(), pts[i].first.y(), pts[i].first.z());
			}
			for (int i = 0; i < ids.size(); ++i)
				ids[i] = pts[i].second;
			GetDelaunaySonList(PointList);
		}
		construct_soa();
	}
	struct alignas(32) Node {
		double x, y, z;
		uint32_t firstChild;      // 在 children[] 中的起始偏移
		uint32_t endChild;      // <= 65535 够用了
	};
public:
	uint32_t getNearestNeighbor(const Point& p) {
		auto cellID = getCellId(p);
		return ids[cellID];
	}

private:
	std::vector<Node> SOASon;          // size == PointList.size()
	std::vector<uint32_t> children;   // 扁平孩子索引池
	void construct_soa()
	{
		const std::size_t N = PointList.size();

		SOASon.resize(N);                  // 一次性分配所有节点
		children.clear();

		for (std::size_t i = 0; i < N; ++i)
		{
			Node& n = SOASon[i];
			n.x = PointList[i].x();
			n.y = PointList[i].y();
			n.z = PointList[i].z();
			n.firstChild = static_cast<uint32_t>(children.size());
			n.endChild = n.firstChild + static_cast<uint32_t>(SonList[i].size());
			children.insert(children.end(),
				SonList[i].begin(),
				SonList[i].end());
			SonList[i].clear();
		}
	}
	static inline double sqDist(double ax, double ay, double az,
		double bx, double by, double bz) noexcept {
		return (ax - bx) * (ax - bx) + (ay - by) * (ay - by) + (az - bz) * (az - bz);
	}

	inline uint32_t getCellId(const Point& p)const {
		uint32_t curCellId = 0;
		const double px = p.x(), py = p.y(), pz = p.z();
		double minDis = sqDist(SOASon[curCellId].x, SOASon[curCellId].y, SOASon[curCellId].z, px, py, pz);
		bool goNextLevel = true;
		int end = -1;
		while (goNextLevel) {
			goNextLevel = false;
			end = SOASon[curCellId].endChild;
			for (int i = SOASon[curCellId].firstChild;i < end;++i) {
				const double dist = sqDist(SOASon[children[i]].x, SOASon[children[i]].y, SOASon[children[i]].z, px, py, pz);
				if (dist < minDis) {
					goNextLevel = true;
					minDis = dist;
					curCellId = children[i];
					break;
				}
			}
		}
		return curCellId;
	}
	void GetDelaunaySonListFarthest(const vector<Point>& points) {
		std::vector<Eigen::Vector3d> EigenPoints = points;
		SonList.resize(EigenPoints.size());

		vector<Delaunay::Vertex_handle> iterators;
		iterators.resize(EigenPoints.size());
		auto v = dt3.insert(Delaunay::Point(EigenPoints[0](0), EigenPoints[0](1), EigenPoints[0](2)));
		v->info() = 0;
		iterators[0] = v;
		ids[0] = 0;
		PointList[0] = EigenPoints[0];

		int nxt = tree.tree[1].maxPos;
		for (int i = 1; i < 4; ++i) {
			ids[i] = nxt;
			auto v = dt3.insert(Delaunay::Point(EigenPoints[nxt](0), EigenPoints[nxt](1), EigenPoints[nxt](2)));
			PointList[i] = Point(EigenPoints[nxt](0), EigenPoints[nxt](1), EigenPoints[nxt](2));
			v->info() = i;
			iterators[i] = v;
			vector<int> nearPoints;
			for (int j = 0; j < i; ++j)
				nearPoints.push_back(ids[j]);
			nxt = tree.Insert(nxt, nearPoints, nearPoints.size());
		}

		SonList[0].push_back(1);
		SonList[0].push_back(2);
		SonList[0].push_back(3);
		SonList[1].push_back(2);
		SonList[1].push_back(3);
		SonList[2].push_back(3);

		auto GetNns = [&](const Eigen::Vector3d& p) {
			int curCellId = 0;
			double minDis = (PointList[curCellId] - p).squaredNorm();
			bool goNextLevel = true;
			while (goNextLevel) {
				goNextLevel = false;
				for (int son : SonList[curCellId]) {
					if ((PointList[son] - p).squaredNorm() < minDis) {
						goNextLevel = true;
						minDis = (PointList[son] - p).squaredNorm();
						curCellId = son;
						break;
					}
				}
			}
			return curCellId;
			};

		vector<Delaunay::Vertex_iterator> vss;
		vss.reserve(300);
		vector<int> nearPoints(points.size());
		int nearPointSize = 0;
		for (int i = 4; i < EigenPoints.size(); ++i) {
			PointList[i] = EigenPoints[nxt];
			auto index = GetNns(EigenPoints[nxt]);
			auto v = dt3.insert(Delaunay::Point(EigenPoints[nxt](0), EigenPoints[nxt](1), EigenPoints[nxt](2)), iterators[index]);
			ids[i] = nxt;
			v->info() = i;
			iterators[i] = v;
			vss.clear();
			dt3.finite_adjacent_vertices(v, std::back_inserter(vss));
			nearPointSize = 0;
			for (auto& vv : vss) {
				if (vv->info() != -1 && (SonList[vv->info()].empty() || SonList[vv->info()].back() != i))
				{
					SonList[vv->info()].push_back(i);
					nearPoints[nearPointSize++] = ids[vv->info()];
				}
			}
			if (i % 100000 == 0) {
				cout << "have inserted " << i << " points" << endl;
			}
			int xxx = nxt;
			nxt = tree.Insert(nxt, nearPoints, nearPointSize);//?
		}

	}

	void GetDelaunaySonList(const vector<Point>& points) {

		SonList.resize(points.size());
		for (int i = 0; i < 4; ++i) {
			auto v = dt3.insert(Delaunay::Point(points[i](0), points[i](1), points[i](2)));
			v->info() = i;
		}
		SonList[0].push_back(1);
		SonList[0].push_back(2);
		SonList[0].push_back(3);
		SonList[1].push_back(2);
		SonList[1].push_back(3);
		SonList[2].push_back(3);

		vector<Delaunay::Vertex_iterator> vss;
		vss.reserve(300);

		Delaunay::Vertex_handle pre;
		using Cell_circulator = Delaunay::Cell_circulator;
		using Cell_handle = Delaunay::Cell_handle;
		using Vertex_handle = Delaunay::Vertex_handle;
		std::vector<Cell_handle> inc_cells;
		inc_cells.reserve(16);
		int infos[4];

		for (int i = 4; i < points.size(); ++i) {
			auto v = dt3.insert(Delaunay::Point(points[i](0), points[i](1), points[i](2)), pre);
			v->info() = i;
			pre = v;
			inc_cells.clear();
			vss.clear();
			dt3.finite_adjacent_vertices(v, std::back_inserter(vss));
			for (auto& vv : vss) {
				if (vv->info() != -1 && (SonList[vv->info()].empty() || SonList[vv->info()].back() != i))
				{
					SonList[vv->info()].emplace_back(i);
				}
			}
		}
		cout << "end extract delaunay result" << endl;
	}
};
