#include <iostream>
#include <vector>
using namespace std;
#include <Eigen/Dense>
typedef Eigen::Vector3d Point;
#include <random>
#include <tuple>
#include <fstream>
#include "dpnn.cpp"
#include <chrono>

tuple<pair<double, double>, pair<double, double>, pair<double, double>> get_bounding_box(const vector<Point>& pts) {
    pair<double, double> bounds_x(1e9, -1e9);
    pair<double, double> bounds_y(1e9, -1e9);
    pair<double, double> bounds_z(1e9, -1e9);
    for (int i = 0;i < pts.size();++i) {
        bounds_x.first = min(bounds_x.first, pts[i].x());
        bounds_x.second = max(bounds_x.second, pts[i].x());
        bounds_y.first = min(bounds_y.first, pts[i].y());
        bounds_y.second = max(bounds_y.second, pts[i].y());
        bounds_z.first = min(bounds_z.first, pts[i].z());
        bounds_z.second = max(bounds_z.second, pts[i].z());
    }
    return make_tuple(bounds_x, bounds_y, bounds_z);
}

void perturbPoints(std::vector<Eigen::Vector3d>& points) {
    if (points.empty()) return;

    // 1. 计算包围盒和对角线长度
    auto [bounds_x, bounds_y, bounds_z] = get_bounding_box(points);

    // 2. 计算包围盒的对角线长度
    double diagonal_length = std::sqrt(
        std::pow(bounds_x.second - bounds_x.first, 2) +
        std::pow(bounds_y.second - bounds_y.first, 2) +
        std::pow(bounds_z.second - bounds_z.first, 2)
    );

    // 3. 扰动范围设为对角线长度的 0.1%
    double perturbation_magnitude = diagonal_length * 0.000001; // 0.1% = 0.001

    // 4. 设置随机数生成器
    std::random_device rd;
    std::mt19937 gen(rd());

    // 5. 定义扰动分布
    std::uniform_real_distribution<double> dist(-perturbation_magnitude, perturbation_magnitude);

    // 6. 对每个点应用扰动
    for (Eigen::Vector3d& point : points) {
        point.x() += dist(gen);
        point.y() += dist(gen);
        point.z() += dist(gen);
    }
}

vector<Point> generateRandomPoints(int num_query_points, const vector<Point>& ps, double scale=0.5) {
    vector<Point> query_points;
    Eigen::Vector3d minbox;
    Eigen::Vector3d maxbox;
    auto bouninx_box = get_bounding_box(ps);

    double len_x = get<0>(bouninx_box).second - get<0>(bouninx_box).first;
    double len_y = get<1>(bouninx_box).second - get<1>(bouninx_box).first;
    double len_z = get<2>(bouninx_box).second - get<2>(bouninx_box).first;

    random_device rd;
    mt19937 gen(rd());

    minbox(0) = get<0>(bouninx_box).first - scale * len_x;
    minbox(1) = get<1>(bouninx_box).first - scale * len_y;
    minbox(2) = get<2>(bouninx_box).first - scale * len_z;

    maxbox(0) = get<0>(bouninx_box).second + scale * len_x;
    maxbox(1) = get<1>(bouninx_box).second + scale * len_y;
    maxbox(2) = get<2>(bouninx_box).second + scale * len_z;

    uniform_real_distribution<double> dis_x(get<0>(bouninx_box).first - scale * len_x, get<0>(bouninx_box).second + scale * len_x);
    uniform_real_distribution<double> dis_y(get<1>(bouninx_box).first - scale * len_y, get<1>(bouninx_box).second + scale * len_y);
    uniform_real_distribution<double> dis_z(get<2>(bouninx_box).first - scale * len_z, get<2>(bouninx_box).second + scale * len_z);

    query_points.reserve(num_query_points);
    for (int i = 0;i < num_query_points;++i) {
        query_points.emplace_back(dis_x(gen), dis_y(gen), dis_z(gen));
    }

    return query_points;
}

std::vector<Point> ReadInput(const std::string& filename) {
    std::vector<Point> pts;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file -> " << filename << std::endl;
        return pts; // Return an empty vector
    }
    double x, y, z;
    while (file >> x >> y >> z) {
        pts.push_back(Point(x, y, z));
    }
    file.close();
    return pts;
}
vector<uint32_t> line_search(const vector<Eigen::Vector3d>& points, const vector<Eigen::Vector3d>& query_points) {
    vector<uint32_t> result;
    result.reserve(query_points.size());

    // 对每个查询点进行搜索
    for (const auto& query_point : query_points) {
        uint32_t nearest_index = 0;
        double min_distance_squared = std::numeric_limits<double>::max();

        // 遍历所有点找到最近的
        for (uint32_t i = 0; i < points.size(); ++i) {
            // 使用Eigen的squaredNorm()计算距离平方
            double distance_squared = (points[i] - query_point).squaredNorm();

            if (distance_squared < min_distance_squared) {
                min_distance_squared = distance_squared;
                nearest_index = i;
            }
        }

        result.push_back(nearest_index);
    }

    return result;
}
int main() {

    vector<Point> pts=ReadInput("../data/lucy.xyz");
    std::cout << pts.size() << std::endl;
    auto query_points = generateRandomPoints(1000000,pts,0.5);

    auto start = std::chrono::high_resolution_clock::now();
    dpnn nnt(pts,false);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double, std::milli>(end - start);
    std::cout << "construction Time based on cgal sorting: " << duration.count() << " ms" << std::endl;

    start = std::chrono::high_resolution_clock::now();
    dpnn nnt_fast(pts,true);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration<double, std::milli>(end - start);
    std::cout << "construction Time based on farthest sampling sorting: " << duration.count() << " ms" << std::endl;

    vector<uint32_t> kpIdx1, kpIdx2;
    start = std::chrono::high_resolution_clock::now();
    for (const auto& p : query_points) {
        auto id = nnt.getNearestNeighbor(p);
        kpIdx1.push_back(id);
    }
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration<double, std::milli>(end - start);
    std::cout << "query Time based on Cgal sorting: " << duration.count()*1000/query_points.size() << " ms" << std::endl;


    start = std::chrono::high_resolution_clock::now();
    for (const auto& p : query_points) {
        auto id = nnt_fast.getNearestNeighbor(p);
        kpIdx2.push_back(id);
    }

    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration<double, std::milli>(end - start);
    std::cout << "query Time based on Farthest Point sorting: " << duration.count()*1000/query_points.size() << " ms" << std::endl;

    // auto correct_answer=line_search(pts,query_points);
    // int wrong=0,right=0;
    // for (int i=0;i<query_points.size();++i) {
    //     if (kpIdx1[i] != correct_answer[i]||kpIdx2[i] != correct_answer[i]) {
    //         wrong++;
    //     }
    //     else {
    //         right++;
    //     }
    // }
    // cout<<"right radio: "<<right*1.0/query_points.size()*100.0<<" %"<<endl;
    return 0;
}