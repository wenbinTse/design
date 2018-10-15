#include "stdafx.h"
#include "CP_Polygon.h"
#include<algorithm>

bool order(const CP_Segment& a, const CP_Segment& b) {
	return a.p1.m_x < b.p1.m_x;
}

bool CP_Polygon::check(string& message) {
	if (!checkLoopsDirection()) {
		message = "环的方向不对";
		return false;
	}
	if (!checkEdgeIntersected()) {
		message = "两边要么不相交，要么交点是顶点";
		return false;
	}
	if (!checkInnerLoopInOuterLoop()) {
		message = "内环应该包含在外环中";
		return false;
	}
	if (!checkInnerLoopInInnerLoop()) {
		message = "内环存在交集";
		return  false;
	}
	if (!checkRegion()) {
		message = "区域之间存在交集";
		return false;
	}
	return true;
}

bool CP_Polygon::check(string& message, string name) {
	bool ans = check(message);
	message = "多边形" + name + "不合法：\n" + message;
	return ans;
}

bool CP_Polygon::checkLoopDirection(CP_Loop loop, bool isOuter) {
	int np = loop.m_pointIDArray.size();
	if (np < 3) return false;
	int i, rightest = 0, xMax = INF_SMALL;
	for (i = 1; i < np; i++) {
		CP_Point point = m_pointArray[loop.m_pointIDArray[i]];
		if (point.m_x > xMax) {
			xMax = point.m_x;
			rightest = i;
		}
	}
	CP_Point p0 = m_pointArray[loop.m_pointIDArray[(rightest - 1 + np) % np]];
	CP_Point p1 = m_pointArray[loop.m_pointIDArray[rightest]];
	CP_Point p2 = m_pointArray[loop.m_pointIDArray[(rightest + 1) % np]];
	return isOuter ^ xmult(p0, p1, p2) < 0;
}

// 检查环的方向
bool CP_Polygon::checkLoopsDirection() {
	int nr = m_regionArray.size(), nl, i, j;
	for (i = 0; i < nr; i++) {
		CP_Region region = m_regionArray[i];
		// 检查外环
		CP_Loop loop = region.m_loopArray[0];
		if (!checkLoopDirection(loop, true)) return false;
		nl = region.m_loopArray.size();
		for (j = 1; j < nl; j++)
			if (!checkLoopDirection(region.m_loopArray[j], false)) return false;
	}
	return true;
}

// 不相交返回True
bool CP_Polygon::checkEdgeIntersected() {
	vector<CP_Segment> edges;
	int nr = m_regionArray.size(), nl, nv, i, j, z;
	for (i = 0; i < nr; i++) {
		nl = m_regionArray[i].m_loopArray.size();
		for (j = 0; j < nl; j++) {
			CP_Loop loop = m_regionArray[i].m_loopArray[j];
			nv = loop.m_pointIDArray.size();
			for (z = 0; z < nv; z++) {
				CP_Segment edge = CP_Segment(m_pointArray[loop.m_pointIDArray[z]], m_pointArray[loop.m_pointIDArray[(z + 1) % nv]]);
				edges.push_back(edge);
			}
		}
	}
	int n = edges.size();
	for (i = 0; i < n - 1; i++)
		for (j = i + 1; j < n; j++) {
			if (edges[i] == edges[j]) return false;
			else if (parallel(edges[i], edges[j]) &&
				(pointInSegment(edges[i].p1, edges[j], false) || pointInSegment(edges[i].p2, edges[j], false))
			)
				return false;
			else if (!segmentIntersected(edges[i], edges[j])) {
				// 不相交
				continue;
			}
			else if (edges[i].p1 == edges[j].p1 || edges[i].p1 == edges[j].p2 || edges[i].p2 == edges[j].p1 || edges[i].p2 == edges[j].p2) {
				// 顶点相交
				continue;
			}
			else
				return false;
		}
	return true;
}

// 判断点是否在多边形里
POINT_STATUS CP_Polygon::include(CP_Point p) {
	// x方向做拓展
	CP_Point INF_POINT = CP_Point(INF, p.m_y);
	int nr = m_regionArray.size(), nl, nv, i, j, z;
	int  count = 0;
	for (i = 0; i < nr; i++) {
		nl = m_regionArray[i].m_loopArray.size();
		for (j = 0; j < nl; j++) {
			CP_Loop loop = m_regionArray[i].m_loopArray[j];
			nv = loop.m_pointIDArray.size();
			for (z = 0; z < nv; z++) {
				CP_Point p1 = m_pointArray[loop.m_pointIDArray[z]];
				p1 = CP_Point(p1.m_x, p1.m_y); // 深复制
				CP_Point p2 = m_pointArray[loop.m_pointIDArray[(z + 1) % nv]];
				p2 = CP_Point(p2.m_x, p2.m_y);
				if (p1.m_y > p2.m_y) swap(p1, p2);
				if (pointInSegment(p, p1, p2))
					return BOUNDARY;
				if (p1.m_y == p2.m_y) {
					// 与x轴平行
				}
				else if (pointInSegment(p1, p, INF_POINT)) {
					count++;
				}
				else if (pointInSegment(p2, p, INF_POINT)) {
					// Do nothing
				}
				else if (segmentIntersected(p, INF_POINT, p1, p2)) {
					count++;
				}
			}
		}
	}
	return count % 2 ? INSIDE : OUTSIDE;
}

// 假设已经确保了线段要么不相交，要么交点是端点
// 判断loop2是否被完全包含在loop1, 完全包含返回True
bool CP_Polygon::checkLoopInLoop(CP_Loop loop1, CP_Loop loop2) {
	int i, nv = loop2.m_pointIDArray.size();
	CP_Polygon polygon;
	polygon.m_pointArray = m_pointArray;
	polygon.m_regionArray.push_back(CP_Region(loop1));
	for (i = 0; i < nv; i++) {
		CP_Point a = m_pointArray[loop2.m_pointIDArray[i]];
		CP_Point b = m_pointArray[loop2.m_pointIDArray[(i + 1) % nv]];
		CP_Point m = middlePoint(a, b);
		if (polygon.include(m) == OUTSIDE) return false;
	}
	return true;
}

// 假设已经确保了线段要么不相交，要么交点是端点
// 判断loop1和loop2是否有相交（不含点相交情况）， 有相交返回True
bool CP_Polygon::checkLoopIntersectLoop(CP_Loop loop1, CP_Loop loop2) {
	int i, nv1 = loop1.m_pointIDArray.size(), nv2 = loop2.m_pointIDArray.size();
	CP_Polygon polygon;
	polygon.m_pointArray = m_pointArray;
	polygon.m_regionArray.push_back(CP_Region(loop1)); //  将loop1赋予多边形
	for (i = 0; i < nv2; i++) {
		CP_Point a = m_pointArray[loop2.m_pointIDArray[i]];
		CP_Point b = m_pointArray[loop2.m_pointIDArray[(i + 1) % nv2]];
		CP_Point m = middlePoint(a, b);
		if (polygon.include(m) != OUTSIDE) return false;
	}
	polygon.m_regionArray.clear();
	polygon.m_regionArray.push_back(CP_Region(loop2)); //  将loop2赋予多边形
	for (i = 0; i < nv1; i++) {
		CP_Point a = m_pointArray[loop1.m_pointIDArray[i]];
		CP_Point b = m_pointArray[loop1.m_pointIDArray[(i + 1) % nv1]];
		CP_Point m = middlePoint(a, b);
		if (polygon.include(m) != OUTSIDE) return false;
	}
	return true;
}

// 假设已经确保了线段要么不相交，要么交点是端点
// 判断内环是否被外环包含
bool CP_Polygon::checkInnerLoopInOuterLoop() {
	int nr = m_regionArray.size(), nl, nv, i, j, z;
	for (i = 0; i < nr; i++) {
		nl = m_regionArray[i].m_loopArray.size();
		CP_Loop outerLoop = m_regionArray[i].m_loopArray[0];
		for (j = 1; j < nl; j++) { // 内环遍历
			CP_Loop innerLoop = m_regionArray[i].m_loopArray[j];
			if (!checkLoopInLoop(outerLoop, innerLoop)) return false;
		}
	}
	return true;
}

// 假设已经确保了线段要么不相交，要么交点是端点
// 判断内环是否有交集
bool CP_Polygon::checkInnerLoopInInnerLoop() {
	int nr = m_regionArray.size(), nl, nv, i, j, z;
	for (i = 0; i < nr; i++) {
		CP_Region region = m_regionArray[i];
		int loopsSize = region.m_loopArray.size();
		for (j = 1; j < loopsSize - 1; j++)
			for (z = j + 1; z < loopsSize; z++) {
				CP_Loop loop1 = region.m_loopArray[j];
				CP_Loop loop2 = region.m_loopArray[z];
				if (checkLoopIntersectLoop(loop1, loop2))
					return false;
			}
	}
	return true;
}

// 假设两个Region都是合法的
// 判断这两个Region是否有非法交集，没有则返回True
bool CP_Polygon::checkRegionIntersectRegion(CP_Region region1, CP_Region region2) {
	CP_Polygon polygon;
	polygon.m_pointArray = m_pointArray;
	int nv1 = region1.m_loopArray[0].m_pointIDArray.size(),
		nv2 = region2.m_loopArray[0].m_pointIDArray.size(), i;
	polygon.m_regionArray.push_back(region1); // 将region1赋予多边形
	CP_Loop outerLoop = region2.m_loopArray[0];
	for (i = 0; i < nv2; i++) {
		CP_Point p1 = m_pointArray[outerLoop.m_pointIDArray[i]];
		CP_Point p2 = m_pointArray[outerLoop.m_pointIDArray[(i + 1) % nv2]];
		CP_Point m = middlePoint(p1, p2);
		if (polygon.include(m) != OUTSIDE) return false;
	}

	polygon.m_regionArray.clear();
	polygon.m_regionArray.push_back(region2); // 将region2赋予多边形
	outerLoop = region1.m_loopArray[0];
	for (i = 0; i < nv1; i++) {
		CP_Point p1 = m_pointArray[outerLoop.m_pointIDArray[i]];
		CP_Point p2 = m_pointArray[outerLoop.m_pointIDArray[(i + 1) % nv1]];
		CP_Point m = middlePoint(p1, p2);
		if (polygon.include(m) != OUTSIDE) return false;
	}
	return true;
}

// 假设两个Region都是合法的
// 判断Region中是否有非法交集，有则返回True
bool CP_Polygon::checkRegion() {
	int i, j, nr = m_regionArray.size();
	for (i = 0; i < nr - 1; i++)
		for (j = i + 1; j < nr; j++) {
			if (!checkRegionIntersectRegion(m_regionArray[i], m_regionArray[j]))
				return false;
		}
	return true;
}


CP_Polygon CP_Polygon::Union(CP_Polygon b) {
	CP_Polygon aNew, bNew, ans;
	vector<CP_Segment> segments;
	addIntersectedPoint(*this, b, aNew, bNew);
	int nr = aNew.m_regionArray.size(), nl, nv, i, j, z;
	for (i = 0; i < nr; i++) {
		nl = aNew.m_regionArray[i].m_loopArray.size();
		for (j = 0; j < nl; j++) {
			CP_Loop loop = aNew.m_regionArray[i].m_loopArray[j];
			nv = loop.m_pointIDArray.size();
			for (z = 0; z < nv; z++) {
				CP_Point p1 = aNew.m_pointArray[loop.m_pointIDArray[z]];
				CP_Point p2 = aNew.m_pointArray[loop.m_pointIDArray[(z + 1) % nv]];
				CP_Segment segment = CP_Segment(p1, p2);
				// 属于A不在B
				if (!segmentInPolygon(segment, b)) { // TODO 重合
					segments.push_back(segment);
				}
			}
		}
	}
	nr = bNew.m_regionArray.size();
	for (i = 0; i < nr; i++) {
		nl = bNew.m_regionArray[i].m_loopArray.size();
		for (j = 0; j < nl; j++) {
			CP_Loop loop = bNew.m_regionArray[i].m_loopArray[j];
			nv = loop.m_pointIDArray.size();
			for (z = 0; z < nv; z++) {
				CP_Point p1 = bNew.m_pointArray[loop.m_pointIDArray[z]];
				CP_Point p2 = bNew.m_pointArray[loop.m_pointIDArray[(z + 1) % nv]];
				CP_Segment segment = CP_Segment(p1, p2);
				if (!segmentInPolygon(segment, *this)) { // TODO 重合
					segments.push_back(segment);
				}
			}
		}
	}
	sort(segments.begin(), segments.end(), order);//先对vector中元素按照x-y坐标升序排序
	return CP_Polygon(segments);
}

CP_Polygon::CP_Polygon(vector<CP_Segment> segments) {
	CP_Point leftest = CP_Point(INF_SMALL, 0), now;
	double xMax = INF_SMALL, yMax = INF_SMALL, yMin = INF;
	int rId = -1, lId = -1;
	while (segments.size()) {
		if (rId == -1 || this->include(segments[0].p1) != INSIDE) {
			rId++;
			lId = 0;
			this->m_regionArray.push_back(CP_Region(rId, this));
			this->m_regionArray[rId].m_regionIDinPolygon = rId;
			this->m_regionArray[rId].m_polygon = this;

			this->m_regionArray[rId].m_loopArray.push_back(CP_Loop(lId, rId, this));
			CP_Loop* tmpLoop = & this->m_regionArray[rId].m_loopArray[lId];
			tmpLoop->m_loopIDinRegion = lId;
			tmpLoop->m_polygon = this;
			tmpLoop->m_regionIDinPolygon = rId;

			this->m_pointArray.push_back(segments[0].p1);
			this->m_pointArray.push_back(segments[0].p2);
			int arrSize = this->m_pointArray.size();
			tmpLoop->m_pointIDArray.push_back(arrSize - 2);
			tmpLoop->m_pointIDArray.push_back(arrSize - 1);

			leftest.m_x = segments[0].p1.m_x;
			leftest.m_y = segments[0].p1.m_y;

			xMax = segments[0].p2.m_x;
			yMax = max(segments[0].p1.m_y, segments[0].p2.m_y);
			yMin = min(segments[0].p1.m_y, segments[0].p2.m_y);
			now = segments[0].p2;

			segments.erase(segments.begin());
		}
		else {
			lId++;
			this->m_regionArray[rId].m_loopArray.push_back(CP_Loop());
			CP_Loop* tmpLoop = & this->m_regionArray[rId].m_loopArray[lId];
			tmpLoop->m_loopIDinRegion = lId;
			tmpLoop->m_polygon = this;
			tmpLoop->m_regionIDinPolygon = rId;

			this->m_pointArray.push_back(segments[0].p1);
			this->m_pointArray.push_back(segments[0].p2);
			int arrSize = this->m_pointArray.size();
			tmpLoop->m_pointIDArray.push_back(arrSize - 2);
			tmpLoop->m_pointIDArray.push_back(arrSize - 1);
			
			leftest.m_x = segments[0].p1.m_x;
			leftest.m_y = segments[0].p1.m_y;
			now = segments[0].p2;
			
			segments.erase(segments.begin());
		}

		while (true) {
			vector<CP_Segment>::iterator it;
			for (it = segments.begin(); it != segments.end(); it++) {
				if (it->p1 == now) break;
			}
			if (it->p2 == leftest) {
				segments.erase(it);
				break;
			}
			this->m_pointArray.push_back(it->p2);
			int arrSize = this->m_pointArray.size();
			this->m_regionArray[rId].m_loopArray[lId].m_pointIDArray.push_back(arrSize - 1);
			now = it->p2;
			xMax = max(it->p2.m_x, xMax);
			yMax = max(it->p2.m_y, yMax);
			yMin = min(it->p2.m_y, yMin);
			segments.erase(it);
		}
	}
}

void gb_distanceMinPointLoop(double&d, int& idRegion, int& idLoop,
	CP_Point& pt, CP_Polygon& pn)
{
	d = 0.0;
	idRegion = -1;
	idLoop = -1;
	int nr = pn.m_regionArray.size();
	int i, j, k, nl, nv, v1, v2;
	double dt;
	for (i = 0; i<nr; i++)
	{
		nl = pn.m_regionArray[i].m_loopArray.size();
		for (j = 0; j<nl; j++)
		{
			nv = pn.m_regionArray[i].m_loopArray[j].m_pointIDArray.size();
			for (k = 0; k<nv; k++)
			{
				v1 = pn.m_regionArray[i].m_loopArray[j].m_pointIDArray[k];
				if (k == nv - 1)
					v2 = pn.m_regionArray[i].m_loopArray[j].m_pointIDArray[0];
				else v2 = pn.m_regionArray[i].m_loopArray[j].m_pointIDArray[k + 1];
				dt = gb_distancePointSegment(pt, pn.m_pointArray[v1], pn.m_pointArray[v2]);
				if ((idLoop == -1) || (d>dt))
				{
					d = dt;
					idRegion = i;
					idLoop = j;
				} // if结束
			} // for(k)结束
		} // for(j)结束
	} // for(i)结束
} // 函数gb_distanceMinPointPolygon结束

void gb_distanceMinPointPolygon(double&d, int& id, CP_Point& pt, CP_Polygon& pn)
{
	d = 0.0;
	id = -1;
	int n = pn.m_pointArray.size();
	if (n <= 0)
		return;
	d = gb_distancePointPoint(pt, pn.m_pointArray[0]);
	id = 0;
	int i;
	double dt;
	for (i = 1; i<n; i++)
	{
		dt = gb_distancePointPoint(pt, pn.m_pointArray[i]);
		if (dt<d)
		{
			d = dt;
			id = i;
		} // if结束
	} // for结束
} // 函数gb_distanceMinPointPolygon结束

double gb_distancePointPoint(CP_Point& p1, CP_Point& p2)
{
	double dx = p1.m_x - p2.m_x;
	double dy = p1.m_y - p2.m_y;
	double d2 = dx*dx + dy*dy;
	double d = sqrt(d2);
	return d;
} // 函数gb_distancePointPoint结束

double gb_distancePointSegment(CP_Point& pt, CP_Point& p1, CP_Point& p2)
{
	double dx0 = p2.m_x - p1.m_x;
	double dy0 = p2.m_y - p1.m_y;
	double dx1 = pt.m_x - p1.m_x;
	double dy1 = pt.m_y - p1.m_y;
	double dx2 = pt.m_x - p2.m_x;
	double dy2 = pt.m_y - p2.m_y;
	double d1 = dx1*dx1 + dy1*dy1;
	double d2 = dx2*dx2 + dy2*dy2;
	double d01 = dx1*dx0 + dy1*dy0;
	double d02 = -dx2*dx0 - dy2*dy0;
	double d, d0;
	if ((d01>0) && (d02>0))
	{
		d0 = dx0*dx0 + dy0*dy0;
		d = d01*d01 / d0; // 如果计算溢出，如何处理?
		d = d1 - d;
		d = sqrt(d);
		return d;
	} // if结束
	if (d1>d2)
		d = d2;
	else d = d1;
	d = sqrt(d);
	return d;
} // 函数gb_distancePointPoint结束

void gb_getIntArrayPointInPolygon(VT_IntArray& vi, CP_Polygon& pn, CP_Point& p, double eT)
{
	int i, n;
	double d;
	n = pn.m_pointArray.size();
	for (i = 0; i<n; i++)
	{
		d = gb_distancePointPoint(p, pn.m_pointArray[i]);
		if (d <= eT)
		{
			vi[i] = i;
		} // if结束
	} // for(i)结束
} // 函数gb_getIntArrayPointInPolygon结束

bool gb_findPointInLoop(CP_Polygon& pn, int& idRegion, int& idLoop, int& idPointInLoop, int pointInPolygon)
{
	idRegion = 0;
	idLoop = 0;
	idPointInLoop = 0;
	int nr, nL, nv;
	int i, j, k;
	nr = pn.m_regionArray.size();
	for (i = 0; i<nr; i++)
	{
		nL = pn.m_regionArray[i].m_loopArray.size();
		for (j = 0; j<nL; j++)
		{
			nv = pn.m_regionArray[i].m_loopArray[j].m_pointIDArray.size();
			for (k = 0; k<nv; k++)
			{
				if (pn.m_regionArray[i].m_loopArray[j].m_pointIDArray[k] == pointInPolygon)
				{
					idRegion = i;
					idLoop = j;
					idPointInLoop = k;
					return true;
				} // if结束
			} // for(nv)结束
		} // for(nL)结束
	} // for(nr)结束
	return false;
} // 函数gb_findPointInLoop结束

// 这里假设所有的条件都成立，即函数内部不判断输入的合法性
void gb_insertPointInPolygon(CP_Polygon& pn, int& idRegion, int& idLoop, int& idPointInLoop, CP_Point& newPoint)
{
	int nv = pn.m_pointArray.size();
	pn.m_pointArray.push_back(newPoint);
	pn.m_regionArray[idRegion].m_loopArray[idLoop].m_pointIDArray.insert(
		pn.m_regionArray[idRegion].m_loopArray[idLoop].m_pointIDArray.begin() + idPointInLoop + 1,
		nv);
} // 函数gb_findPointInLoop结束

void gb_intArrayInit(VT_IntArray& vi, int data)
{
	int n = vi.size();
	int i;
	for (i = 0; i<n; i++)
		vi[i] = data;
} // 函数gb_intArrayInit结束

void gb_intArrayInitLoop(VT_IntArray& vi, CP_Polygon& pn, int idRegion, int idLoop, double eT)
{
	int i, v;
	int n = pn.m_pointArray.size();
	vi.resize(n);
	n = pn.m_regionArray[idRegion].m_loopArray[idLoop].m_pointIDArray.size();
	for (i = 0; i<n; i++)
	{
		v = pn.m_regionArray[idRegion].m_loopArray[idLoop].m_pointIDArray[i];
		vi[v] = v;
	} // for结束
	gb_intArrayInitPointSame(vi, pn, eT);
} // 函数gb_intArrayInitLoop结束

void gb_intArrayInitPoint(VT_IntArray& vi, CP_Polygon& pn, int v, double eT)
{
	int n = pn.m_pointArray.size();
	if (n <= 0)
	{
		vi.clear();
		return;
	} // if结束
	vi.resize(n);
	int i;
	double d;
	for (i = 0; i<n; i++)
	{
		if (i == v)
			vi[i] = i;
		else
		{
			d = gb_distancePointPoint(pn.m_pointArray[i], pn.m_pointArray[v]);
			if (d <= eT)
				vi[i] = i;
			else vi[i] = -1;
		} // if/else结束
	} // for结束
} // 函数gb_intArrayInitPoint结束

void gb_intArrayInitPointSame(VT_IntArray& vi, CP_Polygon& pn, double eT)
{
	int i, j, n;
	double d;
	n = vi.size();
	if (n <= 0)
		return;
	for (i = 0; i<n; i++)
	{
		if (vi[i] >= 0)
		{
			for (j = 0; j<n; j++)
			{
				if (vi[j]<0)
				{
					d = gb_distancePointPoint(pn.m_pointArray[i], pn.m_pointArray[j]);
					if (d <= eT)
						vi[j] = j;
				} // if结束
			} // for(j)结束
		} // if结束
	} // for(i)结束
} // 函数gb_intArrayInitPointSame结束

void gb_intArrayInitPolygon(VT_IntArray& vi, CP_Polygon& pn)
{
	int i;
	int n = pn.m_pointArray.size();
	vi.resize(n);
	for (i = 0; i<n; i++)
		vi[i] = i;
} // 函数gb_intArrayInitPolygon结束

void gb_intArrayInitPolygonSamePoint(VT_IntArray& vr, CP_Polygon& pr, VT_IntArray& vs, CP_Polygon& ps, double eT)
{
	int i, j;
	int n0, n1;
	double da;
	n1 = pr.m_pointArray.size();
	if (n1 <= 0)
	{
		vr.clear();
		return;
	} // if结束
	vr.resize(n1);
	gb_intArrayInit(vr, -1);
	n0 = ps.m_pointArray.size();
	for (i = 0; i<n0; i++)
	{
		if (vs[i]<0)
			continue;
		for (j = 0; j<n1; j++)
		{
			if (vr[j]<0)
			{
				da = gb_distancePointPoint(ps.m_pointArray[i], pr.m_pointArray[j]);
				if (da <= eT)
					vr[j] = j;
			} // if结束
		} // for(j)结束
	} // for(i)结束
} // 函数gb_intArrayInitPolygonSamePoint结束

void gb_intArrayInitRegion(VT_IntArray& vi, CP_Polygon& pn, int idRegion, double eT)
{
	int i, j, nr, v;
	int n = pn.m_pointArray.size();
	vi.resize(n);
	nr = pn.m_regionArray[idRegion].m_loopArray.size();
	for (i = 0; i<nr; i++)
	{
		n = pn.m_regionArray[idRegion].m_loopArray[i].m_pointIDArray.size();
		for (j = 0; j<n; j++)
		{
			v = pn.m_regionArray[idRegion].m_loopArray[i].m_pointIDArray[j];
			vi[v] = v;
		} // for(j)结束
	} // for(i)结束
	gb_intArrayInitPointSame(vi, pn, eT);
} // 函数gb_intArrayInitRegion结束

void gb_moveLoop(CP_Polygon& pn, int idRegion, int idLoop, double vx, double vy)
{
	int nr, nL, nv;
	int i, id;
	nr = pn.m_regionArray.size();
	if ((idRegion<0) || (idRegion >= nr))
		return;
	nL = pn.m_regionArray[idRegion].m_loopArray.size();
	if ((idLoop<0) || (idLoop >= nL))
		return;
	nv = pn.m_regionArray[idRegion].m_loopArray[idLoop].m_pointIDArray.size();
	for (i = 0; i<nv; i++)
	{
		id = pn.m_regionArray[idRegion].m_loopArray[idLoop].m_pointIDArray[i];
		pn.m_pointArray[id].m_x += vx;
		pn.m_pointArray[id].m_y += vy;
	} // for结束
} // 函数gb_moveLoop结束

void gb_movePoint(CP_Polygon& pn, int id, double vx, double vy)
{
	int n = pn.m_pointArray.size();
	if ((id<0) || (id >= n))
		return;
	pn.m_pointArray[id].m_x += vx;
	pn.m_pointArray[id].m_y += vy;
} // 函数gb_movePoint结束

void gb_movePointIntArray(CP_Polygon& pn, VT_IntArray& vi, double vx, double vy)
{
	int n = vi.size();
	int i;
	for (i = 0; i<n; i++)
		gb_movePoint(pn, vi[i], vx, vy);
} // 函数gb_movePoint结束

void gb_movePolygon(CP_Polygon& pn, double vx, double vy)
{
	int n = pn.m_pointArray.size();
	int i;
	for (i = 0; i<n; i++)
	{
		pn.m_pointArray[i].m_x += vx;
		pn.m_pointArray[i].m_y += vy;
	} // for结束
} // 函数gb_movePolygon结束

void gb_moveRegion(CP_Polygon& pn, int idRegion, double vx, double vy)
{
	int nr, nL, nv;
	int i, j, k, id;
	nr = pn.m_regionArray.size();
	if ((idRegion<0) || (idRegion >= nr))
		return;
	i = idRegion;
	nL = pn.m_regionArray[i].m_loopArray.size();
	for (j = 0; j<nL; j++)
	{
		nv = pn.m_regionArray[i].m_loopArray[j].m_pointIDArray.size();
		for (k = 0; k<nv; k++)
		{
			id = pn.m_regionArray[i].m_loopArray[j].m_pointIDArray[k];
			pn.m_pointArray[id].m_x += vx;
			pn.m_pointArray[id].m_y += vy;
		} // for结束
	} // for结束
} // 函数gb_moveRegion结束

  // 将在全局坐标系下的点转换成为在屏幕坐标下的点
  // result:      输出的在屏幕坐标下的点;
  // pointGlobal: 输入的在全局坐标系下的点;
  // scale:       输入的比例因子;
  // translation: 输入的平移坐标值。
void gb_pointConvertFromGlobalToScreen(CP_Point& result, CP_Point pointGlobal, double scale, CP_Point translation, int screenX, int screenY)
{
	result.m_x = (pointGlobal.m_x - translation.m_x)*scale;
	result.m_y = (pointGlobal.m_y - translation.m_y)*scale;
	result.m_x += (screenX / 2);
	result.m_y = screenY / 2 - result.m_y;
} // 函数PointConvertFromGlobalToScreen结束

  // 将在屏幕坐标下的点转换成为在全局坐标系下的点
  // result:      输出的在全局坐标系下的点;
  // pointScreen: 输入的在屏幕坐标系下的点;
  // scale:       输入的比例因子;
  // translation: 输入的平移坐标值。
void gb_pointConvertFromScreenToGlobal(CP_Point& result, CP_Point pointScreen, double scale, CP_Point translation, int screenX, int screenY)
{
	result.m_x = pointScreen.m_x - screenX / 2;
	result.m_y = screenY / 2 - pointScreen.m_y;
	result.m_x = result.m_x / scale + translation.m_x;
	result.m_y = result.m_y / scale + translation.m_y;
} // 函数gb_PointConvertFromScreenToGlobal结束

  // 给多边形p增加新的内环，该内环是外接圆半径为r的正n边形。
bool gb_polygonNewInLoopRegular(CP_Polygon& p, int idRegion, int n, double r, double cx, double cy)
{
	if (n<3)
		return false;
	int nr = p.m_regionArray.size();
	if ((idRegion<0) || (idRegion >= nr))
		return false;
	int nL = p.m_regionArray[idRegion].m_loopArray.size();
	if (nL <= 0)
		return false;
	p.m_regionArray[idRegion].m_loopArray.resize(nL + 1);
	int s = p.m_pointArray.size();
	int t = s + n;
	int i, k;
	p.m_pointArray.resize(t);
	double da = DOUBLE_PI / n;
	double d = 0.0;
	for (i = s; i<t; i++, d += da)
	{
		p.m_pointArray[i].m_x = cx + r*cos(d);
		p.m_pointArray[i].m_y = cy + r*sin(d);
	} // for结束
	p.m_regionArray[idRegion].m_loopArray[nL].m_polygon = &p;
	p.m_regionArray[idRegion].m_loopArray[nL].m_regionIDinPolygon = idRegion;
	p.m_regionArray[idRegion].m_loopArray[nL].m_loopIDinRegion = nL;
	p.m_regionArray[idRegion].m_loopArray[nL].m_pointIDArray.resize(n);
	for (i = 0, k = t - 1; i<n; i++, k--)
	{
		p.m_regionArray[idRegion].m_loopArray[nL].m_pointIDArray[i] = k;
	} // for结束
	return true;
} // 函数gb_polygonNewInLoopRegular结束

  // 给多边形p增加新的外环，该外环是外接圆半径为r的正n边形。
void gb_polygonNewOutLoopRegular(CP_Polygon& p, int n, double r, double cx, double cy)
{
	if (n<3)
		return;
	int s = p.m_pointArray.size();
	int t = s + n;
	int i, k;
	p.m_pointArray.resize(t);
	double da = DOUBLE_PI / n;
	double d = 0.0;
	for (i = s; i<t; i++, d += da)
	{
		p.m_pointArray[i].m_x = cx + r*cos(d);
		p.m_pointArray[i].m_y = cy + r*sin(d);
	} // for结束
	int rs = p.m_regionArray.size();
	p.m_regionArray.resize(rs + 1);
	p.m_regionArray[rs].m_polygon = &p;
	p.m_regionArray[rs].m_regionIDinPolygon = rs;
	p.m_regionArray[rs].m_loopArray.resize(1);
	p.m_regionArray[rs].m_loopArray[0].m_polygon = &p;
	p.m_regionArray[rs].m_loopArray[0].m_regionIDinPolygon = rs;
	p.m_regionArray[rs].m_loopArray[0].m_loopIDinRegion = 0;
	p.m_regionArray[rs].m_loopArray[0].m_pointIDArray.resize(n);
	for (i = 0, k = s; i<n; i++, k++)
	{
		p.m_regionArray[rs].m_loopArray[0].m_pointIDArray[i] = k;
	} // for结束
} // 函数gb_polygonNewOutLoopRegular结束

bool gb_removeLoop(CP_Polygon& pn, int idRegion, int idLoop)
{
	int nL, nLv, iLv, v;
	nL = pn.m_regionArray[idRegion].m_loopArray.size();
	if ((idLoop == 0) || (nL<2))
		return(gb_removeRegion(pn, idRegion));
	nLv = pn.m_regionArray[idRegion].m_loopArray[idLoop].m_pointIDArray.size();
	for (iLv = 0; iLv<nLv; iLv++)
	{
		v = pn.m_regionArray[idRegion].m_loopArray[idLoop].m_pointIDArray[iLv];
		pn.m_pointArray.erase(pn.m_pointArray.begin() + v);
		gb_subtractOneAboveID(pn, v);
	} // for(iLv)结束
	pn.m_regionArray[idRegion].m_loopArray.erase(
		pn.m_regionArray[idRegion].m_loopArray.begin() + idLoop);
	return true;
} // 函数gb_removeLoop结束

bool gb_removePoint(CP_Polygon& pn, int id)
{
	int ir, iL, iLv, nLv;
	bool rf = gb_findPointInLoop(pn, ir, iL, iLv, id);
	if (!rf)
		return false;
	nLv = pn.m_regionArray[ir].m_loopArray[iL].m_pointIDArray.size();
	if (nLv<4) // 删除整个环
		return (gb_removeLoop(pn, ir, iL));
	pn.m_regionArray[ir].m_loopArray[iL].m_pointIDArray.erase(
		pn.m_regionArray[ir].m_loopArray[iL].m_pointIDArray.begin() + iLv);
	pn.m_pointArray.erase(pn.m_pointArray.begin() + id);
	gb_subtractOneAboveID(pn, id);
	return true;
} // 函数gb_removePoint结束

bool gb_removeRegion(CP_Polygon& pn, int idRegion)
{
	int nr, nL, nLv, iL, iLv, v;
	nr = pn.m_regionArray.size();
	if (nr<2)
	{
		pn.mb_clear();
		return true;
	} // if结束
	nL = pn.m_regionArray[idRegion].m_loopArray.size();
	for (iL = 0; iL<nL; iL++)
	{
		nLv = pn.m_regionArray[idRegion].m_loopArray[iL].m_pointIDArray.size();
		for (iLv = 0; iLv<nLv; iLv++)
		{
			v = pn.m_regionArray[idRegion].m_loopArray[iL].m_pointIDArray[iLv];
			pn.m_pointArray.erase(pn.m_pointArray.begin() + v);
			gb_subtractOneAboveID(pn, v);
		} // for(iLv)结束
	} // for(iL)结束
	pn.m_regionArray.erase(pn.m_regionArray.begin() + idRegion);
	return true;
} // 函数gb_removeRegion结束

// 所有ID大于等于id的点的ID都减一
void gb_subtractOneAboveID(CP_Polygon& pn, int id)
{
	int nr = pn.m_regionArray.size();
	int nL, nLv, ir, iL, iLv;
	for (ir = 0; ir<nr; ir++)
	{
		nL = pn.m_regionArray[ir].m_loopArray.size();
		for (iL = 0; iL<nL; iL++)
		{
			nLv = pn.m_regionArray[ir].m_loopArray[iL].m_pointIDArray.size();
			for (iLv = 0; iLv<nLv; iLv++)
			{
				if (pn.m_regionArray[ir].m_loopArray[iL].m_pointIDArray[iLv] >= id)
					pn.m_regionArray[ir].m_loopArray[iL].m_pointIDArray[iLv]--;
			} // for(iLv)结束
		} // for(iL)结束
	} // for(ir)结束
} // 函数gb_subtractOneAboveID结束

double xmult(CP_Point a, CP_Point b, CP_Point c) {
	return (a.m_x - c.m_x) * (b.m_y - c.m_y) - (a.m_y - c.m_y) * (b.m_x - c.m_x);
}

bool parallel(CP_Point a1, CP_Point a2, CP_Point b1, CP_Point b2) {
	double help = (a1.m_x - a2.m_x) * (b1.m_y - b2.m_y) - (b1.m_x - b2.m_x) * (a1.m_y - a2.m_y);
	return ZERO(help);
}

bool parallel(CP_Segment s1, CP_Segment s2) {
	return parallel(s1.p1, s1.p2, s2.p1, s2.p2);
}

bool inSameLine(CP_Point a, CP_Point b, CP_Point c) {
	return ZERO(xmult(a, b, c)) || ZERO(xmult(a, c, b));
}

CP_Point middlePoint(CP_Point p1, CP_Point p2) {
	return CP_Point((p1.m_x + p2.m_x) / 2, (p1.m_y + p2.m_y) / 2);
}

CP_Point middlePoint(CP_Segment segment) {
	return middlePoint(segment.p1, segment.p2);
}

bool pointInSegment(CP_Point a, CP_Point l1, CP_Point l2, bool include_vertex) {
	double tmp = xmult(a, l1, l2);
	if (include_vertex) {
		if (a.m_x == l1.m_x)
			return ZERO(tmp) && (
				(a.m_y <= l1.m_y && a.m_y >= l2.m_y) || (a.m_y <= l2.m_y && a.m_y >= l1.m_y)
			);
		else
			return ZERO(tmp) && (
				(a.m_x <= l1.m_x && a.m_x >= l2.m_x) || (a.m_x <= l2.m_x && a.m_x >= l1.m_x)
			);
	}

	else {
		if (a.m_x == l1.m_x)
			return ZERO(tmp) && (
				(a.m_y < l1.m_y && a.m_y > l2.m_y) || (a.m_y < l2.m_y && a.m_y > l1.m_y)
			);
		else
			return ZERO(tmp) && (
				(a.m_x < l1.m_x && a.m_x > l2.m_x) || (a.m_x < l2.m_x && a.m_x > l1.m_x)
			);
	}
}

bool pointInSegment(CP_Point a, CP_Segment s, bool include_vertex) {
	return pointInSegment(a, s.p1, s.p2, include_vertex);
}

// 假设已经确保两线段要么不相交，要么相交于顶点
bool segmentInPolygon(CP_Segment segment, CP_Polygon polygon) {
	return polygon.include(middlePoint(segment)) != OUTSIDE;
}

bool inSameSideOfSegment(CP_Point a, CP_Point b, CP_Segment s) {
	return xmult(a, s.p1, s.p2) * xmult(b, s.p1, s.p2) > 0;
}

bool segmentIntersected(CP_Point p1, CP_Point p2, CP_Point p3, CP_Point p4) {
	return !inSameSideOfSegment(p1, p2, CP_Segment(p3, p4)) &&
		!inSameSideOfSegment(p3, p4, CP_Segment(p1, p2));
}

bool segmentIntersected(CP_Segment s1, CP_Segment s2) {
	return segmentIntersected(s1.p1, s1.p2, s2.p1, s2.p2);
}

bool between(double a, double b, double c) {
	return (a - b) * (a - c) < 0;
}

double distance(CP_Point a, CP_Point b) {
	return gb_distancePointPoint(a, b);
}

CP_Point getIntersection(CP_Point a1, CP_Point a2, CP_Point b1, CP_Point b2) {
	CP_Point ans = a1;
	double t = ((a1.m_x - b1.m_x)*(b1.m_y - b2.m_y) - (a1.m_y - b1.m_y)*(b1.m_x - b2.m_x))
		/ ((a1.m_x - a2.m_x)*(b1.m_y - b2.m_y) - (a1.m_y - a2.m_y)*(b1.m_x - b2.m_x));
	ans.m_x += (a2.m_x - a1.m_x)*t;
	ans.m_y += (a2.m_y - a1.m_y)*t;
	return ans;
}

void insertPoint(CP_Polygon& polygon, CP_Point point, CP_Loop* & loop, int position) {
	if (position > 0 &&
		polygon.m_pointArray[loop->m_pointIDArray[position - 1]] == point) {
		return;
	}
	polygon.m_pointArray.push_back(point);
	loop->m_pointIDArray.insert(loop->m_pointIDArray.begin() + position, polygon.m_pointArray.size() - 1);
}

// 正确插入交点
void addIntersectedPoint(CP_Polygon a_old, CP_Polygon b_old, CP_Polygon& a, CP_Polygon& b) {
	a = a_old;
	b = b_old;
	int anr, anl, anv, bnr, bnl, bnv, ia, ja, za, ib, jb, zb;
	anr = a.m_regionArray.size(), bnr = b.m_regionArray.size();
	vector<CP_Segment> segmentsSame, segmentsReve;
	for (ia = 0; ia < anr; ia++) {
		anl = a.m_regionArray[ia].m_loopArray.size();
		for (ja = 0; ja < anl; ja++) {
			CP_Loop* aLoop = &a.m_regionArray[ia].m_loopArray[ja];
			anv = aLoop->m_pointIDArray.size();
			for (za = 0; za < anv; za++) {
				CP_Point pa1 = a.m_pointArray[aLoop->m_pointIDArray[za]];
				CP_Point pa2 = a.m_pointArray[aLoop->m_pointIDArray[(za + 1) % anv]];

				for (ib = 0; ib < bnr; ib++) {
					bnl = b.m_regionArray[ib].m_loopArray.size();
					for (jb = 0; jb < bnl; jb++) {
						CP_Loop* bLoop = &b.m_regionArray[ib].m_loopArray[jb];
						bnv = bLoop->m_pointIDArray.size();
						for (zb = 0; zb < bnv; zb++) {
							CP_Point pb1 = b.m_pointArray[bLoop->m_pointIDArray[zb]];
							CP_Point pb2 = b.m_pointArray[bLoop->m_pointIDArray[(zb + 1) % bnv]];
						
							if (pa1 == pb1 && pa2 == pb2) {
								segmentsSame.push_back(CP_Segment(pa1, pa2));
							} else if (pa1 == pb2 && pa2 == pb1) {
								segmentsReve.push_back(CP_Segment(pa1, pa2));
							}
							else if (parallel(pa1, pa2, pb1, pb2) && !inSameLine(pa1, pa2, pb1) || !segmentIntersected(pa1, pa2, pb1, pb2)) {
								continue;
							}
							else if (parallel(pa1, pa2, pb1, pb2) && inSameLine(pa1, pa2, pb1)) {
								// 完全重合(a in b)
								if (pa1.m_x != pa2.m_x && between(pa1.m_x, pb1.m_x, pb2.m_x) && between(pa2.m_x, pb1.m_x, pb2.m_x) ||
									(pa1.m_x == pa2.m_x && between(pa1.m_y, pb1.m_y, pb2.m_y) && between(pa2.m_y, pb1.m_y, pb2.m_y))
									) {
									// 同向
									if (distance(pa1, pb1) < distance(pa2, pb1)) { // TODO(考虑pa1 == pb1)
										insertPoint(b, pa1, bLoop, zb + 1);
										insertPoint(b, pa2, bLoop, zb + 2);
										segmentsSame.push_back(CP_Segment(pa1, pa2));
									}
									// 反向
									else {
										insertPoint(b, pa2, bLoop, zb + 1);
										insertPoint(b, pa1, bLoop, zb + 2);
										segmentsReve.push_back(CP_Segment(pa1, pa2));
									}
									zb += 2;
								}
								// 完全重合(b in a)
								else if (pa1.m_x != pa2.m_x && between(pb1.m_x, pa1.m_x, pa2.m_x) && between(pb2.m_x, pa1.m_x, pa2.m_x) ||
									(pa1.m_x == pa2.m_x && between(pb1.m_y, pa1.m_y, pa2.m_y) && between(pb2.m_y, pa1.m_y, pa2.m_y))
									) {
									// 同向
									if (distance(pb1, pa1) < distance(pb2, pa1)) {
										insertPoint(a, pb1, aLoop, za + 1);
										insertPoint(a, pb2, aLoop, za + 2);
										segmentsSame.push_back(CP_Segment(pb1, pb2));
									}
									else {
										insertPoint(a, pb2, aLoop, za + 1);
										insertPoint(a, pb1, aLoop, za + 2);
										segmentsReve.push_back(CP_Segment(pb2, pb1));
									}
									pa2 = a.m_pointArray[aLoop->m_pointIDArray[(za + 1) % anv]];
								}
								// 部分重合 // pa1-pb1-pa2-pb2
								else if (pointInSegment(pb1, pa1, pa2) && pointInSegment(pa2, pb1, pb2)) {
									insertPoint(a, pb1, aLoop, za + 1);
									insertPoint(b, pa2, bLoop, zb + 1);
									segmentsSame.push_back(CP_Segment(pb1, pa2));
								}
								// pb1 pa1 pb2 pa2
								else if (pointInSegment(pa1, pb1, pb2) && pointInSegment(pb2, pa1, pa2)) {
									insertPoint(a, pb2, aLoop, za + 1);
									insertPoint(b, pa1, bLoop, zb + 1);
									segmentsSame.push_back(CP_Segment(pb2, pa1));
								}
								// pa1 pb2 pa2 pb1
								else if (pointInSegment(pb2, pa1, pb2) && pointInSegment(pa2, pb1, pb2)) {
									insertPoint(a, pb2, aLoop, za + 1);
									insertPoint(b, pa2, bLoop, zb + 1);
									segmentsReve.push_back(CP_Segment(pb2, pa2));
								}
								// pb1 pa1 pb1 pa2
								else if (pointInSegment(pb1, pa1, pb2) && pointInSegment(pa1, pb1, pb2)) {
									insertPoint(a, pb1, aLoop, za + 1);
									insertPoint(b, pa1, bLoop, zb + 1);
									segmentsReve.push_back(CP_Segment(pb1, pa1));
								}
							}
							// 普通相交
							else {
								CP_Point intersection = getIntersection(pa1, pa2, pb1, pb2);
								// 端点相交
								if ((intersection == pa1 && (pa1 == pb1 || pa1 == pb2)) ||
									(intersection == pa2 && (pa2 == pb1 || pa2 == pb2)) ||
									(intersection == pb1 && (pb1 == pa1 || pb1 == pa2)) ||
									(intersection == pb2 && (pb2 == pa1 || pb2 == pa2)) ) {
									continue;
								}
								insertPoint(a, intersection, aLoop, za + 1);
								insertPoint(b, intersection, bLoop, zb + 1);
								zb += 1;
							}
							anv = aLoop->m_pointIDArray.size();
							bnv = bLoop->m_pointIDArray.size();
							pa1 = a.m_pointArray[aLoop->m_pointIDArray[za]];
							pa2 = a.m_pointArray[aLoop->m_pointIDArray[(za + 1) % anv]];
						} // zb
					} // jb
				} // ib
				anv = a.m_regionArray[ia].m_loopArray[ja].m_pointIDArray.size();
			} //za
		} // ja
	} // ia
}
