#pragma once
#ifndef CP_POLYGON_H
#define CP_POLYGON_H

#include <iostream>
using namespace std;
#include <vector>

extern double ACC;

#define DOUBLE_PI           6.28318530717958647692
#define PI                  3.14159265358979323846
#define HALF_PI             1.57079632679489661923
#define ZERO(x)             (x) < ACC && (x) > -ACC
#define INF                 1e30
#define INF_SMALL           1e-30

enum POINT_STATUS{
	INSIDE,
	BOUNDARY,
	OUTSIDE
};

class CP_Point
{
	public:
		double m_x, m_y;
	public:
		CP_Point(void) :m_x(0.0), m_y(0.0) {};
		CP_Point(double x, double y) : m_x(x), m_y(y) {};
		bool operator == (const CP_Point& p) {
			return abs(m_x - p.m_x) < ACC && abs(m_y - p.m_y) < ACC;
		}
}; // 类CP_Point定义结束

class CP_Segment {
public: 
	CP_Point p1, p2;
public:
	CP_Segment(CP_Point _p1, CP_Point _p2) :p1(_p1), p2(_p2) {};
	bool operator == (const CP_Segment& s) {
		return (p1 == s.p1 && p2 == s.p2);
	}
};

typedef vector<CP_Point> VT_PointArray;

class CP_Polygon;

typedef vector<int> VT_IntArray;
typedef vector<VT_IntArray> VT_IntArray2;

class CP_Loop
{
	public:
		VT_IntArray m_pointIDArray;
		int m_loopIDinRegion;
		int m_regionIDinPolygon;
		CP_Polygon* m_polygon;

	public:
		CP_Loop(void) :m_loopIDinRegion(0), m_regionIDinPolygon(0), m_polygon(NULL) { };
		CP_Loop(int loopId, int regionId, CP_Polygon* polygon)
			:m_loopIDinRegion(loopId), m_regionIDinPolygon(regionId), m_polygon(polygon) {};
}; // 类CP_Loop定义结束

typedef vector<CP_Loop> VT_LoopArray;

class CP_Region
{
	public:
		VT_LoopArray m_loopArray;
		int m_regionIDinPolygon;
		CP_Polygon* m_polygon;

	public:
		CP_Region(void) :m_regionIDinPolygon(0), m_polygon(NULL) { }
		CP_Region(CP_Loop loop) :m_regionIDinPolygon(0), m_polygon(NULL) {
			m_loopArray.push_back(loop);
		}
		CP_Region(int regionId, CP_Polygon* polygon) :m_regionIDinPolygon(regionId), m_polygon(polygon) {}
}; // 类CP_Region定义结束
typedef vector<CP_Region> VT_RegionArray;

class CP_Polygon
{
	public:
		VT_PointArray m_pointArray;
		VT_RegionArray m_regionArray;

	private:
		bool checkLoopInLoop(CP_Loop l1, CP_Loop l2);
		bool checkLoopIntersectLoop(CP_Loop loop1, CP_Loop loop2);
		bool checkLoopsDirection();
		bool checkLoopDirection(CP_Loop loop, bool isOuter);
		bool checkEdgeIntersected();
		bool checkInnerLoopInOuterLoop();
		bool checkInnerLoopInInnerLoop();
		bool checkRegion();
		bool checkRegionIntersectRegion(CP_Region region1, CP_Region region2);
	public:
		CP_Polygon(vector<CP_Segment> segments);
		CP_Polygon() {};
		void mb_clear() { m_pointArray.clear(); m_regionArray.clear(); }
		POINT_STATUS include(CP_Point p);
		bool check(string& message);
		bool check(string& message, string name);
		CP_Polygon Union(const CP_Polygon& b);
		CP_Polygon Intersection(const CP_Polygon& b);
		CP_Polygon Subtract(const CP_Polygon& b);
		CP_Polygon operator = (const CP_Polygon polygon);
}; // 类CP_Polygon定义结束

// 点到多边形所有边的最短距离
extern void     gb_distanceMinPointLoop(double&d, int& idRegion, int& idLoop, CP_Point& pt, CP_Polygon& pn);

// 点到多边形所有点的最短距离
extern void     gb_distanceMinPointPolygon(double&d, int& id, CP_Point& pt, CP_Polygon& pn);

// 点到点的距离
extern double   gb_distancePointPoint(CP_Point& p1, CP_Point& p2);

// 点到线段的最短距离
extern double   gb_distancePointSegment(CP_Point& pt, CP_Point& p1, CP_Point& p2);


// 如果pi到点距离小于eT，vi[i]=i
extern void     gb_getIntArrayPointInPolygon(VT_IntArray& vi, CP_Polygon& pn, CP_Point& p, double eT);

// 找到点出现在多边形的位置
extern bool     gb_findPointInLoop(CP_Polygon& pn, int& idRegion, int& idLoop, int& idPointInLoop, int pointInPolygon);

// 将点插入到多边形指定位置
extern void     gb_insertPointInPolygon(CP_Polygon& pn, int& idRegion, int& idLoop, int& idPointInLoop, CP_Point& newPoint);

// 将vi的所有元素初始化为data
extern void     gb_intArrayInit(VT_IntArray& vi, int data);

extern void     gb_intArrayInitLoop(VT_IntArray& vi, CP_Polygon& pn, int idRgion, int idLoop, double eT);
extern void     gb_intArrayInitPoint(VT_IntArray& vi, CP_Polygon& pn, int v, double eT);
extern void     gb_intArrayInitPointSame(VT_IntArray& vi, CP_Polygon& pn, double eT);
extern void     gb_intArrayInitPolygon(VT_IntArray& vi, CP_Polygon& pn);
extern void     gb_intArrayInitPolygonSamePoint(VT_IntArray& vr, CP_Polygon& pr, VT_IntArray& vs, CP_Polygon& ps, double eT);
extern void     gb_intArrayInitRegion(VT_IntArray& vi, CP_Polygon& pn, int idRegion, double eT);

extern void     gb_moveLoop(CP_Polygon& pn, int idRegion, int idLoop, double vx, double vy);
extern void     gb_movePoint(CP_Polygon& pn, int id, double vx, double vy);
extern void     gb_movePointIntArray(CP_Polygon& pn, VT_IntArray& vi, double vx, double vy);
extern void     gb_movePolygon(CP_Polygon& pn, double vx, double vy);
extern void     gb_moveRegion(CP_Polygon& pn, int idRegion, double vx, double vy);

extern void     gb_pointConvertFromGlobalToScreen(CP_Point& result, CP_Point pointGlobal, double scale, CP_Point translation, int screenX, int screenY);
extern void     gb_pointConvertFromScreenToGlobal(CP_Point& result, CP_Point pointScreen, double scale, CP_Point translation, int screenX, int screenY);
extern bool     gb_polygonNewInLoopRegular(CP_Polygon& p, int idRegion, int n, double r, double cx, double cy);
extern void     gb_polygonNewOutLoopRegular(CP_Polygon& p, int n, double r, double cx, double cy);
extern bool     gb_removeLoop(CP_Polygon& pn, int idRegion, int idLoop);
extern bool     gb_removePoint(CP_Polygon& pn, int id);
extern bool     gb_removeRegion(CP_Polygon& pn, int idRegion);
extern void     gb_subtractOneAboveID(CP_Polygon& pn, int id);

extern double xmult(CP_Point a, CP_Point b, CP_Point c);
extern double distance(CP_Point a, CP_Point b);
extern bool parallel(CP_Point a1, CP_Point a2, CP_Point b1, CP_Point b2);
extern bool parallel(CP_Segment s1, CP_Segment s2);
extern CP_Point middlePoint(CP_Point p1, CP_Point p2);
extern CP_Point middlePoint(CP_Segment s);
extern bool pointInSegment(CP_Point a, CP_Point l1, CP_Point l2, bool include_vertex = true);
extern bool pointInSegment(CP_Point a, CP_Segment s, bool include_vertex = true);
extern POINT_STATUS segmentInPolygon(CP_Segment segment, CP_Polygon polygon);
extern bool inSameSideOfSegment(CP_Point a, CP_Point b, CP_Segment s);
extern bool segmentIntersected(CP_Segment s1, CP_Segment s2);
extern bool segmentIntersected(CP_Point p1, CP_Point p2, CP_Point p3, CP_Point p4);
extern bool inSameLine(CP_Point a, CP_Point b, CP_Point c);
extern CP_Point getIntersection(CP_Point a1, CP_Point a2, CP_Point b1, CP_Point b2);
extern void addIntersectedPoint(CP_Polygon a, CP_Polygon b, CP_Polygon& c, CP_Polygon& d,
	vector<CP_Segment>& sameDir, vector<CP_Segment>& reveDir);
#endif

