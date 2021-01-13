#include "ofApp.h"
ofPolyline knob;
bool loadedKnobs = false;
typedef std::tuple<int, int, int> i3tuple;
int maxIter = 500;
vector<i3tuple> circles;
vector<ofRectangle> boundingBoxes;
bool drawBoundingBoxes = false;
int scallopR = 350;
ofImage c;
float pi = 3.14159;
//--------------------------------------------------------------
ofPoint segsegintersect(const ofPoint& a, const ofPoint& b, const ofPoint& c, const ofPoint& d) {
	
	float tnum = (a.x - c.x) * (c.y - d.y) - (a.y - c.y) * (c.x - d.x);
	float tdenom = (a.x - b.x) * (c.y - d.y) - (a.y - b.y) * (c.x - d.x);
	float unum = (a.x - b.x) * (a.y - c.y) - (a.y - b.y) * (a.x - c.x);
	float udenom = (a.x - b.x) * (c.y - d.y) - (a.y - b.y) * (c.x - d.x);
	float t = tnum / tdenom;
	return ofPoint(a.x + t * (b.x - a.x), a.y + t * (b.y - a.y));
	
	}
ofPoint segIntersection(
	const ofPoint& p0, 
	const ofPoint& p1, 
	const ofPoint& p2, 
	const ofPoint& p3
	)
{
	ofPoint d0 = p1 - p0;
	ofPoint d1 = p3 - p2;

	float a = d0.dot(d0);
	float b = -d0.dot(d1);
	float c = d0.dot(d1);
	float d = -d1.dot(d1);
	float e = (p2 - p0).dot(d0);
	float f = (p2 - p0).dot(d1);
	float det = a * d - b * c;
	float s, t;
	if (det != 0.0) {
		det = 1.0 / det;
		s = (e * d - b * f) * det;
		t = (a * f - e * c) * det;
	}
	else {	// d0 and d1 parallel
		float s0 = p0.dot(d0);
		float s1 = p1.dot(d0);
		float t0 = p2.dot(d0);
		float t1 = p3.dot(d0);
		bool flip0 = false;
		bool flip1 = false;

		if (s0 > s1) { float f = s0; s0 = s1; s1 = f; flip0 = true; }
		if (t0 > t1) { float f = t0; t0 = t1; t1 = f; flip1 = true; }

		if (s0 >= t1) {
			s = !flip0 ? 0.0 : 1.0;
			t = !flip1 ? 1.0 : 0.0;
		}
		else if (t0 >= s1) {
			s = !flip0 ? 1.0 : 0.0;
			t = !flip1 ? 0.0 : 1.0;
		}
		else {		// overlap
			float mid = (s0 > t0) ? (s0 + t1) * 0.5 : (t0 + s1) * 0.5;
			s = (s0 == s1) ? 0.5 : (mid - s0) / (s1 - s0);
			t = (t0 == t1) ? 0.5 : (mid - t0) / (t1 - t0);
		}
	}
	if (s < 0.0) {
		s = 0.0;
		t = f / d;
	}
	else if (s > 1.0) {
		s = 1.0;
		t = (f + b) / d;
	}
	if (t < 0.0) {
		t = 0.0;
		if (e < 0) {
			s = 0.0;
		}
		else if (e > a) {
			s = 1.0;
		}
		else {
			s = e / a;
		}
	}
	else if (t > 1.0) {
		t = 1.0;
		if (e - b < 0.0) {
			s = 0.0;
		}
		else if ((e - b) > a) {
			s = 1.0;
		}
		else {
			s = (e - b) / a;
		}
	}

	float b0 = 1.0 - s;
	float b1 = s;
	return p0 + s * d0;
}

// Given three colinear points p, q, r, the function checks if 
// point q lies on line segment 'pr' 
bool onSegment(ofPoint p, ofPoint q, ofPoint r)
{
	if (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) &&
		q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y))
		return true;

	return false;
}

// To find orientation of ordered triplet (p, q, r). 
// The function returns following values 
// 0 --> p, q and r are colinear 
// 1 --> Clockwise 
// 2 --> Counterclockwise 
int orientation(ofPoint p, ofPoint q, ofPoint r)
{
	// See https://www.geeksforgeeks.org/orientation-3-ordered-points/ 
	// for details of below formula. 
	int val = (q.y - p.y) * (r.x - q.x) -
		(q.x - p.x) * (r.y - q.y);

	if (val == 0) return 0;  // colinear 

	return (val > 0) ? 1 : 2; // clock or counterclock wise 
}

// The main function that returns true if line segment 'p1q1' 
// and 'p2q2' intersect. 
bool doIntersect(ofPoint p1, ofPoint q1, ofPoint p2, ofPoint q2)
{
	// Find the four orientations needed for general and 
	// special cases 
	int o1 = orientation(p1, q1, p2);
	int o2 = orientation(p1, q1, q2);
	int o3 = orientation(p2, q2, p1);
	int o4 = orientation(p2, q2, q1);
	
	// General case 
	if (o1 != o2 && o3 != o4)
		return true;

	// Special Cases 
	// p1, q1 and p2 are colinear and p2 lies on segment p1q1 
	if (o1 == 0 && onSegment(p1, p2, q1)) return true;

	// p1, q1 and q2 are colinear and q2 lies on segment p1q1 
	if (o2 == 0 && onSegment(p1, q2, q1)) return true;

	// p2, q2 and p1 are colinear and p1 lies on segment p2q2 
	if (o3 == 0 && onSegment(p2, p1, q2)) return true;

	// p2, q2 and q1 are colinear and q1 lies on segment p2q2 
	if (o4 == 0 && onSegment(p2, q1, q2)) return true;

	return false; // Doesn't fall in any of the above cases 
}
void ofApp::concentric() {
	int bigR = 300;
	int center = 400;
	int r = 30;
	for (int i = 0; i < 8; ++i) {

		auto circle = make_tuple(center, center, r);
		addCircle(circle);
		circles.push_back(circle);
		if(r>30){
			int offset = 0;
		for(int theta = 0;theta<360;theta+=50*360./(2*3.15)/r){

			if (theta == 0) {
				offset = ofRandom(0, 360);
			}
		ofPolyline radialEdge;
		float rad = float(theta + offset) * 2 * 3.14159 / 360.;
		ofPoint p1 = { center + (r+3) * cos(rad),center + (r+3) * sin(rad) };
		ofPoint p2 = { center + (r - 43) * cos(rad),center + (r-43) * sin(rad) };
		radialEdge.addVertex(p1);
		radialEdge.addVertex(p2);
		radialEdge = radialEdge.getResampledByCount(30);
		edges.push_back(radialEdge);
		}
		}
		r += 40;
	}
}
void ofApp::makePieces() {
	//vertices : round to 1 pix , list endpoints w/ directed edges and direction
	map<pair<int,int>, vector<pair<int,bool>>> vertices;
	for (int i = 0; i < edges.size();++i) {
		auto e = edges[i];
		auto p0 = e[0];
		auto r0 = make_pair(int(p0.x), int(p0.y));
		auto p1 = e[e.size()-1];
		auto r1 = make_pair(int(p1.x), int(p1.y));
		auto key0 = vertices.find(r0);
		auto key1 = vertices.find(r1);
		if (key0 != vertices.end()) {
			vertices[r0].push_back(make_pair(i, true));
		}
		else {
			vector<pair<int, bool>> a;
			a.push_back(make_pair(i, true));
			vertices[r0] = a;
		}
		if (key1 != vertices.end()) {
			vertices[r1].push_back(make_pair(i, false));
		}
		else {
			vector<pair<int, bool>> a;
			a.push_back(make_pair(i, false));
			vertices[r1] = a;
		}
	}

}
void ofApp::setup() {
	gui.setup();
	srand(time(NULL));
	//addLotsofCircles();
	loadedKnobs = true;
	loadKnobs();
	knobUsage.resize(knobPaths.size());
	//knob = addKnob(edges[0]);
	//intersections();

	vector<ofPolyline> newEdges;

	//edges = newEdges;
	
}
void ofApp::addKnobs() {
	for (int i = 0; i < edges.size(); ++i) {
		auto e = addKnob(i);
		edges[i] = e;
		//newEdges.push_back(e);
	}
	for (int i = 0; i < knobUsage.size(); ++i) {
		cout << knobNames[i] << ": " << knobUsage[i] << endl;
	}
}
//circle in format x,y,R
bool ofApp::floodFillTest(i3tuple circ, vector<i3tuple>& circles){
	cBuffer.allocate(ofGetWidth(), ofGetHeight());
	cBuffer.begin();
	ofClear(0, 0, 0);
	for (auto circle : circles) {
		int x = get<0>(circle);
		int y = get<1>(circle);
		int radius = get<2>(circle);
		int segments = int(radius);
		ofPolyline c;
		for (int i = 0; i < segments; ++i) {
			c.addVertex(ofPoint(x + radius * cos(2. * 3.14159 * float(i) / float(segments)), y + radius * sin(2. * 3.14159 * float(i) / float(segments))));
		}
		c.close();
		c.draw();
	}
	int x = get<0>(circ);
	int y = get<1>(circ);
	int radius = get<2>(circ);
	int segments = int(radius);
	ofPolyline c;
	for (int i = 0; i < segments; ++i) {
		c.addVertex(ofPoint(x + radius * cos(2. * 3.14159 * float(i) / float(segments)), y + radius * sin(2. * 3.14159 * float(i) / float(segments))));
	}
	c.close();
	c.draw();
	cBuffer.end();
	ofPixels pixels; 
	cBuffer.readToPixels(pixels);
	int x0 = x - radius;
	int x1 = x + radius;
	int y0 = y - radius;
	int y1 = y + radius;
	for (int x0 = x0; x0 < x1; x0++) {
		for (int y0 = y0; y0 < y1; y0++) {
			
		}
	}
	return false;
}
float distance(int x1, int x2, int y1, int y2) {
	return sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2));
}
float distance(ofPoint p1, ofPoint p2) {
	return sqrt(pow((p1.x - p2.x), 2) + pow((p1.y - p2.y), 2));
}
int ofApp::overlapped(i3tuple c1, i3tuple c2) {
	float d = sqrt(pow((get<0>(c2) - get<0>(c1)), 2) + pow((get<1>(c2) - get<1>(c1)), 2));
	float r = (float)get<2>(c1);
	float R = (float)get<2>(c2);
	float r1 = min(r, R);
	float r2 = max(r, R);
	r = r1;
	R = r2;
	/*float A = circleOverlap(c1, c2);
	if (A > 0 && (A < 72 * 72 || A > R*R*3-72*72)) {
		return -1;
	}
	*/
	if (abs(r + R - d) < r / 4.) {
		return -1;
	}
	if (abs(r - R - d) < r / 4.) {
		return -1;
	}
	if (abs(R - r - d) < r / 4.) {
		return -1;
	}
	return 0;
}
int ofApp::overlapped1(i3tuple c1, i3tuple c2) {
	float d = sqrt(pow((get<0>(c2) - get<0>(c1)), 2) + pow((get<1>(c2) - get<1>(c1)), 2));
	float r = (float)get<2>(c1);
	float R = (float)get<2>(c2);
	float r1 = min(r, R);
	float r2 = max(r, R);
	r = r1;
	R = r2;
	/*float A = circleOverlap(c1, c2);
	if (A > 0 && (A < 72 * 72 || A > R*R*3-72*72)) {
		return -1;
	}
	*/
	if (abs(r + R - d) < 6) {
		return -1;
	}
	if (abs(r - R - d) < 6) {
		return -1;
	}
	if (abs(R - r - d) < 6) {
		return -1;
	}
	return 0;
}
bool ofApp::testCircleCheap(i3tuple circle, vector<i3tuple>& circles) {
	int overlap = 0;
	bool ok = true;
	for (auto circ : circles) {
		overlap = overlapped1(circ, circle);
		if (overlap == -1) {
			ok = false;
			break;
		}
		if (overlap == 1) {
			ok = true;
		}

	}
	return ok;
}
bool ofApp::testCircle(i3tuple circle, vector<i3tuple>& circles) {
	int overlap = 0;
	bool ok = true;
	for (auto circ : circles) {
		overlap = overlapped(circ, circle);
		if (overlap == -1) {
			ok = false;
			break;
		}
		if (overlap == 1) {
			ok = true;
		}

		}
	
	if(ok){
	auto oldEdges = edges;
	addCircle(circle);
	intersections();
	for (ofPolyline edge : edges) {
		if (edge.getLengthAtIndex(edge.size() - 1) < 24) {
			edges = oldEdges;
			return false;
			
		}
	}
	edges = oldEdges;
	return true;
	}
	return false;
}

bool ofApp::testCircle(i3tuple circle) {
	int overlap = 0;
	bool ok = true;
	ofPoint p0(get<0>(circle), get<1>(circle));
	int r1 = get<2>(circle);
	for (auto circ : circles) {
		ofPoint p1(get<0>(circ), get<1>(circ));
		int r2 = get<2>(circ);
		float d = p0.distance(p1);
		if (abs(d - r1 - r2) < min(r1, r2) / 2 || abs(d + r1 - r2) < min(r1, r2) / 2 || abs(d - r1 + r2) < min(r1, r2) / 2) {
			return false;
		}

	}
	auto oldEdges = edges;
	addCircle(circle);
	intersections();
	for (ofPolyline edge : edges) {
		if (edge.getLengthAtIndex(edge.size() - 1) < 10) {
			edges = oldEdges;
			return false;

		}
	}
	edges = oldEdges;
	return true;
}
float ofApp::circleOverlapChord(const i3tuple c1, const i3tuple c2) {
	float d = sqrt(pow((get<0>(c2) - get<0>(c1)), 2) + pow((get<1>(c2) - get<1>(c1)), 2));
	float r = (float)get<2>(c1);
	float R = (float)get<2>(c2);
	float a = 1. / d * sqrt((-d + r - R) * (-d - r + R) * (-d + r + R) * (d + r + R));
	return a;
}
float ofApp::circleOverlap(const i3tuple c1, const i3tuple c2) {
	float d = sqrt(pow((get<0>(c2) - get<0>(c1)), 2) + pow((get<1>(c2) - get<1>(c1)), 2));
	float r = (float)get<2>(c1);
	float R = (float)get<2>(c2);
	if (d > r && d > R) {
		return 0;
	}
	float t1 = (d * d + r * r - R * R) / (2 * d * r);
	float t2 = (d * d + R * R - r * r) / (2 * d * R);
	float t3 = 1. / 2. * sqrt((-d + r + R) * (d + r - R) * (d - r + R) * (d + r + R));
	float A = r * r * acos(t1) + R * R * acos(t2) - t3;
	return A;
}
bool pointCircle(ofPoint ray) {
	bool intersected = false;
	for (auto circle : circles) {
		float d = sqrt(pow((get<0>(circle) - ray.x), 2) + pow((get<1>(circle) - ray.y), 2));
		if (d < get<2>(circle)) {

			intersected = true;
			break;
		}
	}
	return intersected;

}

int pointCircleIndex(ofPoint ray) {
	bool intersected = false;
	int index = -1;
	for (int i=0;i<circles.size();++i){
		auto circle = circles[i];
		index = i;
		float d = sqrt(pow((get<0>(circle) - ray.x), 2) + pow((get<1>(circle) - ray.y), 2));
		if (d < get<2>(circle)) {

			intersected = true;
			break;
		}
	}
	return index;

}
void ofApp::scatterCircles() {
	int minR = 40;
	int maxR = 41;
	int minD = 40;
	int size = 300;
	int center = 400;
	vector<ofPoint> scatterPoints;
	int origin = 400;
	int maxY = 400;
	int maxX = 400;
	//spiral spacing
	float theta = 0;
	float sr = 0;
	//scatterPoints.push_back(ofPoint(origin, origin));
	/*
	for (int i = 0; i < 1000; ++i) {
		
		
		float x = cos(theta) *sr + origin;
		float y = sin(theta) * sr + origin;
		if (sr == 0) {
			sr = minR;
			theta = 1.25;
		}
		else{
		theta += minR / sr * 1.25;
		sr = theta / (2 * pi) * minR*.9+ minR;
		}
		if(distance(ofPoint(origin,origin),ofPoint(x,y))<size){
		scatterPoints.push_back(ofPoint(x, y));
		}
	}
	*/
	cout << scatterPoints.size();
	 //random points
	/*
	for (int tries = 0; tries < 100000; ++tries) {
		int x = rand() % (size * 2) - size + center;
		int y = rand() % (size * 2) - size + center;
		ofPoint p(x, y);
		bool ok = true;
		for (auto p1 : scatterPoints) {
			if (distance(p, p1) < minD) {
				ok = false;
				break;
			}
		}
		if (ok) {
			scatterPoints.push_back(p);
		}
	}*/
	
	//regular spacing
	
	
	for (int y = origin; y < maxY + origin; y += minR * ySpacing) {
		for (int x = origin; x < maxX + origin; x += minR * xSpacing) {
			scatterPoints.push_back(ofPoint(x, y));
		}
	}
	
	
	//in order 
	/*
	while (scatterPoints.size() > 0) {
		int i = 0;
		ofPoint p = scatterPoints[i];
		scatterPoints.erase(scatterPoints.begin() + i);
		int r = minR + rand() % (maxR - minR);
		i3tuple circle = make_tuple(p.x, p.y, r);
		auto circleSegs = scatterAddCircle(circle);
		if (circleSegs.size() > 0) {
			circles.push_back(circle);
			for (auto e : circleSegs) {
				edges.push_back(e);
			}
		}
	}
	*/
	
	//random ordering
	
	while (scatterPoints.size() > 0) {
		int i = rand() % scatterPoints.size();
		ofPoint p = scatterPoints[i];
		scatterPoints.erase(scatterPoints.begin()+i);
		int r = minR + rand() % (maxR - minR);
		i3tuple circle = make_tuple(p.x, p.y, r);
		auto circleSegs = scatterAddCircle(circle);
		if(circleSegs.size()>0){
		circles.push_back(circle);
		for (auto e : circleSegs) {
			edges.push_back(e);
		}
		}
	}
	
	//greedy ordering
	/*
	int i = rand() % scatterPoints.size();
	ofPoint p = scatterPoints[i];
	scatterPoints.erase(scatterPoints.begin() + i);
	int r = minR + rand() % (maxR - minR);
	i3tuple circle = make_tuple(p.x, p.y, r);
	auto circleSegs = scatterAddCircle(circle);
	circles.push_back(circle);
	for (auto e : circleSegs) {
		edges.push_back(e);
	}
	
	while (scatterPoints.size() > 0) {
		int candidate = 0;
		float minFound = 100000;
		for (int i = 0; i < scatterPoints.size(); ++i) {
			ofPoint p1 = scatterPoints[i];
			float d = distance(p, p1);
			if (d < minFound) {
				ofPoint p = scatterPoints[i];
				i3tuple circle = make_tuple(p.x, p.y, r);
				auto circleSegs = scatterAddCircle(circle);
				if(circleSegs.size()>0){
				minFound = d;
				candidate = i;
				}
			}
		}

		ofPoint p = scatterPoints[candidate];
		scatterPoints.erase(scatterPoints.begin() + candidate);
		/*
		bool ok = true;
		for (int tries = 0; tries < 1000; ++tries) {
			ok = true;
			int r = minR + rand() % (maxR - minR);
			for (auto c : circles) {
				ofPoint p0(get<0>(c), get<1>(c));
				int r1 = get<2>(c);
				float d = distance(p, p0);
				if (abs(d - r - r1) < 9 || abs(d-r+r1)<18 || abs(d+r-r1)<18) {
					ok = false;
				}
			}
			if (ok) {
				break;
			}
		}
		if (!ok) {
			cout << "lost a point" << endl;
			continue;
		}
		*//*
		i3tuple circle = make_tuple(p.x, p.y, r);
		auto circleSegs = scatterAddCircle(circle);
		if (circleSegs.size() > 0) {
			circles.push_back(circle);
			for (auto e : circleSegs) {
				edges.push_back(e);
			}
		}
	}
	*/
}
vector<ofPolyline> ofApp::scatterAddCircle(i3tuple circle) {
	vector<ofPolyline> circleSegs;
	int x = get<0>(circle);
	int y = get<1>(circle);
	ofPoint center(x, y);
	int r = get<2>(circle);
	vector<float> startStop;
	startStop.clear();
	bool cIntersected = true;
	bool intersected = false;
	for (int i = 0; i < 2 * 3.14 * r; ++i) {
		float checkTheta = i*1.0/r;
		ofPoint pp = center + ofPoint(r * cos(checkTheta), r * sin(checkTheta));
		intersected = pointCircle(pp);
		if (cIntersected != intersected) {
			startStop.push_back(checkTheta);
			cIntersected = !cIntersected;
		}
	}
	if (startStop.size() == 1) {
		cout << "no intersection" << endl;
		ofPolyline circ;
		
		for (int i = 0; i < 2 * 3.14 * r; ++i) {
			float theta = i * 1.0 / r;
			ofPoint pp = center + ofPoint(r * cos(theta), r * sin(theta));
			circ.addVertex(pp);
		}
		circ.close();
		circleSegs.push_back(circ);
		return circleSegs;
	}
	if (startStop.size() == 0) {
		cout << "fully intersected";
		return circleSegs;
	}
	if (startStop[0] == 0) {
		startStop[0] = startStop[startStop.size() - 1];
		startStop.pop_back();
	}
	if (startStop.size() % 2 != 0) {
		cout << "uneven";
	}
	bool added = false;
	cout << startStop.size() << endl;
	for (int seg = 0; seg < startStop.size(); seg += 2) {
		float theta1 = startStop[seg];
		float theta2 = startStop[seg + 1];
		ofPolyline circleSeg;
		while (theta2 < theta1) {
			theta2 += 2 * 3.1416;
		}
		theta2 += 3. / float(r);
		theta1 -= 3. / float(r);
		while (theta1 < theta2) {
			ofPoint p1 = center + ofPoint(r * cos(theta1), r * sin(theta1));
			circleSeg.addVertex(p1);
			theta1 += 1.0/r;
		}
		circleSeg.addVertex(center + ofPoint(r * cos(theta2), r * sin(theta2)));
		if (circleSeg.getLengthAtIndex(circleSeg.size() - 1) < 30) {
			continue;
		}
		added = true;
		circleSegs.push_back(circleSeg);
		//circleSeg = addKnob(edges.size() - 1);
		//edges.pop_back();
		//edges.push_back(circleSeg);
	}
	return circleSegs;
}
void ofApp::infinity() {
	cout << " hexagons" << endl;
	float size = 150;
	int minD = 30;
	int minR = 20;
	int maxR = 40;
	int o = 400;
	ofPoint center(o, o);
	vector < ofPoint> scatterPoints;
	for (int tries = 0; tries < 100000; ++tries) {
		ofPoint p(ofRandom(-1200, 1200), ofRandom(-1200, 1200));
		float q = 2./3.*(p.x) / size;
		float r = (-1. / 3. * p.x + sqrt(3) / 3. * p.y) / size;
		if (abs(q)>=1 || abs(r)>=1) {
			continue;
		}
		p = p + center;
		bool ok = true;
		for (auto p1 : scatterPoints) {
			if (distance(p, p1) < minD) {
				ok = false;
				break;
			}
		}
		if (ok) {
			scatterPoints.push_back(p);
		}

	}
	int i = rand() % scatterPoints.size();
	ofPoint p = scatterPoints[i];
	scatterPoints.erase(scatterPoints.begin() + i);
	int r = minR + rand() % (maxR - minR);
	i3tuple circle = make_tuple(p.x, p.y, r);
	auto circleSegs = scatterAddCircle(circle);
	circles.push_back(circle);
	for (auto e : circleSegs) {
		edges.push_back(e);
	}

	while (scatterPoints.size() > 0) {
		int candidate = 0;
		float minFound = 100000;
		for (int i = 0; i < scatterPoints.size(); ++i) {
			ofPoint p1 = scatterPoints[i];
			float d = distance(p, p1);
			if (d < minFound) {
				ofPoint p = scatterPoints[i];
				i3tuple circle = make_tuple(p.x, p.y, r);
				auto circleSegs = scatterAddCircle(circle);
				if (circleSegs.size() > 0) {
					minFound = d;
					candidate = i;
				}
			}
		}

		ofPoint p = scatterPoints[candidate];
		scatterPoints.erase(scatterPoints.begin() + candidate);
		i3tuple circle = make_tuple(p.x, p.y, r);
		auto circleSegs = scatterAddCircle(circle);
		if (circleSegs.size() > 0) {
			circles.push_back(circle);
			for (auto e : circleSegs) {
				edges.push_back(e);
			}
		}
	}
	

}
void ofApp :: regularScallop() {
	int maxY = 72*6;
	int maxX = 72 * 6;
	int origin = 100;
	int r = scallopRadius;
	bool half = false;
	float theta = 3.14159 / 2;
	auto circle = make_tuple(origin, origin, r);
	circles.push_back(circle);
	addCircle(circle);
	for (int y = origin; y < maxY+origin; y += r *ySpacing) {
	for (int x = origin; x < maxX+origin; x += r *xSpacing) {
		if (x == origin && y == origin) { continue; }
			auto ray = ofPoint(x, y);
			bool intersected = false;
			float theta1 = theta, theta2 = theta;
			ofPoint p1, p2, p1old, p2old;
			//right
			while (!intersected) {
				p2old = p2;
				theta2 -= 3. / (2. * 3.14 * r);
				p2 = ray + ofPoint(r * cos(theta2), r * sin(theta2));
				intersected = pointCircle(p2);
			}
			//theta2 -= 7. / (2 * 3.14 * r);;
			p2 = ray + ofPoint(r * cos(theta2), r * sin(theta2));
			intersected = false;
			while (!intersected) {
				p1old = p1;
				theta1 += 3. / (2. * 3.14 * r);
				p1 = ray + ofPoint(r * cos(theta1), r * sin(theta1));
				intersected = pointCircle(p1);
			}
			//theta1 += 7. / (2. * 3.14 * r);
			p1 = ray + ofPoint(r * cos(theta1), r * sin(theta1));




			ofPolyline circleSeg;
			circleSeg.addVertex(p1);
			while (theta1 > theta2) {
				p1 = ray + ofPoint(r * cos(theta1), r * sin(theta1));
				circleSeg.addVertex(p1);
				theta1 -= 2. * 3.14159 / 360.;
			}

			circleSeg.addVertex(p2);
			edges.push_back(circleSeg);
			i3tuple c = make_tuple(x, y, r);
			circles.push_back(c);
		}
	}
}
float overlapRatio = .5;
float minRadius = 20;
float maxRadius = 40;
void ofApp::genGrid() {
	scallopR = 100;
	vector<ofPolyline> holdEdges;
	overlapRatio = 0;
	for(int i=0;i<10;++i){
		overlapRatio += .1;
		minRadius = 20;
		maxRadius = 40;
	for(int j=0;j<10;++j){
		minRadius += 2;
		maxRadius += 2;
		int cx = i * 350;
		int cy = j * 350;
	makeScallop();
	intersections();
	clean();
	addKnobs();
	draw();
	for (auto line : edges) {
		line.translate(ofPoint(cx, cy));
		holdEdges.push_back(line);
	}
	edges.clear();
	circles.clear();
	boundingBoxes.clear();
	drawBoundingBoxes = false;
	
	c.clear();
	}
	}
	edges = holdEdges;
	exportPDF("grid.pdf");
}
void ofApp::makeScallop() {
	int center = 400;

	auto circle = make_tuple(center, center, 30);
	circles.push_back(circle);
	addCircle(circle);
	float cmax = 50;
	float theta = 0;
	int placed = 0;
	int oldplaced = 0;
	for (int i = 0; i < 20000; ++i) {
		theta = 2 * 3.14159 * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		/*if (placed != oldplaced) {
			theta += 2 * 3.14159 * 60. / 360.;
		}
		else {
			theta += float(rand()) / float(RAND_MAX) * .05;
		}
		placed = oldplaced;
		*/
		auto base = ofPoint(cos(theta) , sin(theta));
		auto ray = base*2.*cmax+ofPoint(center,center);
		bool intersected = false;
		int r = rand() % int(maxRadius-minRadius) + minRadius;
		if(i< 5000){
		while (!intersected) {
			intersected = pointCircle(ray);
			ray -= base*4.;
			
		}
		ray += base * (rand() % int(r*overlapRatio) + 4);
		}
		else {
			ray = ofPoint(ofRandom(100, 700), ofRandom(100, 700));
			intersected = pointCircle(ray);
			if (intersected) {
				continue;
			}

		}
		
		
		auto c = make_tuple(int(ray.x), int(ray.y), r);
		if (!testCircleCheap(c, circles)) {
			continue;
		}
		float cmaxold = cmax;
		cmax = max(cmax, abs(ray.x - center));
		cmax = max(cmax, abs(ray.y - center));
		if (cmax > scallopR) {
			cmax = cmaxold;
			continue;
		}
		/*
		intersected = false;
		float theta1=theta,theta2 = theta;
		ofPoint p1,p2,p1old,p2old;
		//right
		while (!intersected) {
			p2old = p2;
			theta2 -= 3. / (2. * 3.14 * r);
			p2 = ray + ofPoint(r * cos(theta2), r * sin(theta2));
			intersected = pointCircle(p2);
		}
		//theta2 -= 7. / (2 * 3.14 * r);;
		p2 = ray + ofPoint(r * cos(theta2), r * sin(theta2));
		intersected = false;
		while (!intersected) {
			p1old = p1;
			theta1 += 3. / (2. * 3.14 * r);
			p1 = ray + ofPoint(r * cos(theta1), r * sin(theta1));
			intersected = pointCircle(p1);
		}
		//theta1 += 7. / (2. * 3.14 * r);
		p1 = ray + ofPoint(r * cos(theta1), r * sin(theta1));
		*/
		vector<float> startStop;
		startStop.clear();
		bool cIntersected = true;
		for (int tDeg = 0; tDeg < 360; ++tDeg) {
			float checkTheta = 2 * 3.14159 * tDeg / 360.;
			ofPoint pp = ray + ofPoint(r * cos(checkTheta), r * sin(checkTheta));
			intersected = pointCircle(pp);
			if (cIntersected != intersected) {
				startStop.push_back(checkTheta);
				cIntersected = !cIntersected;
			}
		}
		if (startStop.size() == 0) {
			cout << "no intersection" << endl;
			continue;
		}
		cout << "startstop" << endl;
		for (auto vv : startStop) {
			cout << vv << endl;
		}
		cout << "k" << endl;
		
		if (startStop[0] == 0) {
			startStop[0] = startStop[startStop.size() - 1];
			startStop.pop_back();
		}
		if (startStop.size() % 2 != 0) {
			cout << "uneven";
			continue;
		}
		circles.push_back(c);
		bool added = false;
		for (int seg = 0; seg < startStop.size(); seg += 2) {
			float theta1 = startStop[seg];
			float theta2 = startStop[seg + 1];
			ofPolyline circleSeg;
			while (theta2 < theta1) {
				theta2 += 2 * 3.1416;
			}
			theta2 += 2. / float(r);
			theta1 -= 2. / float(r);
			//circleSeg.addVertex(ray + ofPoint(r * cos(theta1), r * sin(theta1)));
			while(theta1<theta2){
				ofPoint p1 = ray + ofPoint(r * cos(theta1), r * sin(theta1));
				circleSeg.addVertex(p1);
				theta1 += 2. * 3.14159 / 360.;
			}
			circleSeg.addVertex(ray + ofPoint(r * cos(theta2), r * sin(theta2)));
			if (circleSeg.getLengthAtIndex(circleSeg.size() - 1) < 30) {
				continue;
			}
			added = true;
			edges.push_back(circleSeg);
			//circleSeg = addKnob(edges.size() - 1);
			//edges.pop_back();
			//edges.push_back(circleSeg);
			}
		if (!added) {
			circles.pop_back();
		}
		else{ cout << "added a circle " << i << endl; }
		
		}
		

		
		

	
	

	
}
bool ofApp::noCircleIntersect(i3tuple c1) {
	bool ok = true;
	for (auto c2 : circles) {
		
		float d = sqrt(pow((get<0>(c2) - get<0>(c1)), 2) + pow((get<1>(c2) - get<1>(c1)), 2));
		if (d < get<2>(c2)||d<get<2>(c1)) {
			ok = false;
			break;
		}
	}
	return ok;
}
ofPoint ofApp::segEdgesIntersect(pair<ofPoint, ofPoint> segment) {
	auto width = abs(segment.second.x - segment.first.x)+1;
	auto height = abs(segment.second.y - segment.first.y) + 1;
	auto x0 = min(segment.second.x, segment.first.x);
	auto y0 = min(segment.second.y, segment.first.y);
	ofRectangle bbox1(x0,y0,width,height);
	auto va = segment.first;
	auto vb = segment.second;

	for (int j = 0; j < edges.size(); ++j) {
		auto e2 = &edges[j];
		auto bbox2 = e2->getBoundingBox();
		if (bbox1.intersects(bbox2)) {
			auto v2 = e2->getVertices();
				for (int l = 0; l < v2.size(); ++l) {
					auto vc = v2[l];
					auto vd = v2[(l + 1) % v2.size()];

					
					bool intersects = doIntersect(va, vb, vc, vd);
					if (intersects) {
						ofPoint intersection = segsegintersect(va, vb, vc, vd);
						return intersection;
					}

				}
			}
		}
	return ofPoint(-10000,-10000);
	}
void ofApp::addLotsofCircles() {
	int bigR = 300;
	int center = 400;
	auto circle = make_tuple(center, center, bigR);
	if(testCircle(circle,circles)){
	addCircle(circle);
	circles.push_back(circle);
	}
	
	int numMed = rand() % 3 + 2;
	int added = 0;
	int i = 0;
	int largePacked = 0;
	int medPacked = 0;
	int smallPacked = 0;
	int tinyPacked = 0;
	int medR = 0;
	while (i < maxIter) {
		i += 1;
		if (largePacked < 2) {
			medR = float(float(400 + rand() % 300) / 1000. * float(bigR));
		}
		else {
			if (medPacked < 10) {
				medR = float(float(100 + rand() % 300) / 1000. * float(bigR));
			}
			else {
				if (smallPacked < 80) {
					medR = float(float(100 + rand() % 100) / 1000. * float(bigR));
				}
				else {
					
						medR = float(float(50 + rand() % 50) / 1000. * float(bigR));
					}
				}
			}
		
		int theta = rand() % 360;
		float r = sqrt(static_cast <float> (rand()) / static_cast <float> (RAND_MAX));
		int rOffset = int(r * (float(bigR) * 1));
		int x = int(float(center) + (cos(float(theta) / (2 * 3.14))) * float(rOffset));
		int y = int(float(center) + (sin(float(theta) / (2 * 3.14))) * float(rOffset));
		circle = make_tuple(x, y, medR);
		bool ok = testCircle(circle);
		if (ok) {
					 addCircle(circle);
					 circles.push_back(circle);
					 if (medR > float(bigR) * .5) {
						 largePacked += 1;
					 }
					 if (medR > float(bigR) * .3) {
						 medPacked += 1;
					 }
					 if (medR > float(bigR) * .1) {
						 smallPacked += 1;
					 }
					 else  {
						 tinyPacked += 1;
					 }
					 cout << "edges : " << edges.size() << endl;
					 i = 0;
				 }
			
		

	}
	
}

void ofApp::packCircles() {
	
	int bigR = 300;
	int center = 400;
	auto circle = make_tuple(center, center, bigR);
	circles.push_back(circle);
	addCircle(circle);
	int numMed = rand() % 3 + 2;
	int added = 0;
	for (int i = 0; i < 10000; ++i) {
		int medR = float(float(400 + rand() % 300) / 1000. * float(bigR));
		int theta = rand() % 360;
		float r = sqrt(static_cast <float> (rand()) / static_cast <float> (RAND_MAX));
		int rOffset = int(r * (float(bigR) * 1));
		int x = int(float(center) + (cos(float(theta) / (2 * 3.14))) * float(rOffset));
		int y = int(float(center) + (sin(float(theta) / (2 * 3.14))) * float(rOffset));
		circle = make_tuple(x, y, medR);
		bool ok = testCircle(circle, circles);
		if (ok) {
			addCircle(circle);
			circles.push_back(circle);
			added += 1;
		}
		if (added == numMed) {
			break;
		}
	}
	int numSmall = rand() % 6 + 3;
	added = 0;
	for (int i = 0; i < 1000; ++i) {
		int medR = float(float(200 + rand() % 100) / 1000. * float(bigR));
		int theta = rand() % 360;
		float r = sqrt(static_cast <float> (rand()) / static_cast <float> (RAND_MAX));
		int rOffset = int(r * (float(bigR) * 1.2));
		int x = int(center + cos(float(theta) / (2 * 3.14)) * float(rOffset));
		int y = int(center + sin(float(theta) / (2 * 3.14)) * float(rOffset));
		circle = make_tuple(x, y, medR);
		bool ok = testCircle(circle, circles);
		if (ok) {
			addCircle(circle);
			circles.push_back(circle);
			added += 1;
		}
		if (added == numSmall) {
			break;
		}
	}
	int numTiny = rand() % 7 + 3;
	added = 0;
	for (int i = 0; i < 1000; ++i) {
		int medR = float(float(100 + rand() % 50) / 1000. * float(bigR));
		int theta = rand() % 360;
		float r = sqrt(static_cast <float> (rand()) / static_cast <float> (RAND_MAX));
		int rOffset = int(r * (float(bigR) * 1.2));
		int x = int(center + cos(float(theta) / (2 * 3.14)) * float(rOffset));
		int y = int(center + sin(float(theta) / (2 * 3.14)) * float(rOffset));
		circle = make_tuple(x, y, medR);
		bool ok = testCircle(circle, circles);
		if (ok) {
			addCircle(circle);
			circles.push_back(circle);
			added += 1;
		}
		if (added == numTiny) {
			break;
		}
	}
}
//----------------------------------
void ofApp::randomColor() {
	ofFbo offscreen;
	float scale = 1;
	int w = 1000 * scale, h = 1000 * scale;
	c.allocate(w, h, OF_IMAGE_COLOR);
	offscreen.allocate(w, h);
	offscreen.begin();
	ofBackground(255);
	ofSetColor(255);
	ofDisableSmoothing();
	ofPushMatrix();
	ofColor eColor = { 0,0,0 };
	ofSetColor(eColor);
	for (auto e : edges) {
		auto e1 = e;
		e1.scale(scale, scale);
		e1.draw();
	}

	offscreen.end();
	offscreen.readToPixels(c.getPixels());
	c.update();
	ofColor B = c.getColor(0, 0);
	int count = 0;
	for (int j = 1; j < h - 1; ++j) {
		for (int i = 1; i < w - 1; ++i) {

			if (c.getColor(i, j) == B) {
				 count += 1;
				 cout << "adding color " << count << endl << "circles: " << circles.size();
				float R = ofRandom(0, 255);
				float G = ofRandom(0, 255);
				ofColor rColor = ofColor::fromHsb(R, 255, 255);
				vector<pair<int, int>> neighbors;
				int cX = i;
				int cY = j;
				do {
					if (c.getColor(cX, cY) == B) {
						c.setColor(cX, cY, rColor);
						if (cX < w - 1) {
							neighbors.push_back(make_pair(cX + 1, cY));
						}
						if (cX > 1) {
							neighbors.push_back(make_pair(cX - 1, cY));
						}
						if (cY < h - 1) {
							neighbors.push_back(make_pair(cX, cY + 1));
						}
						if (cY > 1) {
							neighbors.push_back(make_pair(cX, cY - 1));
						}
					}
					cX = neighbors[neighbors.size() - 1].first;
					cY = neighbors[neighbors.size() - 1].second;
					neighbors.pop_back();

				} while (neighbors.size() > 0);
			}
		}
	}
}
//-----------------------------------
void ofApp::addColor() {
	ofFbo offscreen;
	float scale = 300 / 72.;
	int w = 1000*scale, h = 1000*scale;
	c.allocate(w, h, OF_IMAGE_COLOR_ALPHA);
	offscreen.allocate(w, h);
	offscreen.begin();
	ofBackground(255);
	ofSetColor(255);
	ofDisableSmoothing();
	ofPushMatrix();
	ofColor eColor = { 255,0,255 };
	ofSetColor(eColor);
	for (auto e : edges) {
		auto e1 = e;
		e1.scale(scale, scale);
		e1.draw();
	}
	
	offscreen.end();
	offscreen.readToPixels(c.getPixels());
	c.update();
	ofColor B = c.getColor(0, 0);
	unsigned char* pix = c.getPixels().getData();
	int cColor = -1;
	for (int j = 1; j < h - 1; ++j) {
		for (int i = 1; i < w - 1; ++i) {

			if (c.getColor(i, j) == B) {
				cColor += 1;
				if (cColor > 250) {
					cout << "too many pieces" << endl;
				}
				cColor = cColor % 250;
				float R = ofRandom(0, 255);
				float G = ofRandom(0, 255);
				ofColor rColor = { R,G,(255-(R+G)/2) };
				vector<pair<int, int>> neighbors;
				int cX = i;
				int cY = j;
				do {
					if (c.getColor(cX, cY) == B) {
						c.setColor(cX, cY, { 0,float(cColor),0 });
						//c.setColor(cX, cY, rColor);
						if (cX < w - 1) {
							neighbors.push_back(make_pair(cX + 1, cY));
						}
						if (cX > 1) {
							neighbors.push_back(make_pair(cX - 1, cY));
						}
						if (cY < h - 1) {
							neighbors.push_back(make_pair(cX, cY+1));
						}
						if (cY > 1) {
							neighbors.push_back(make_pair(cX, cY-1));
						}
					}
					cX = neighbors[neighbors.size() - 1].first;
					cY = neighbors[neighbors.size() - 1].second;
					neighbors.pop_back();
					
				} while (neighbors.size() > 0);
			}
		}
	}
	bool removedLine = false;
	vector<tuple<int, int, ofColor>> changeColor;
	do {
		removedLine = false;
		for (int j = 1; j < h - 1; ++j) {
			for (int i = 1; i < w - 1; ++i) {
				vector<pair<int, int>> neighbors;
				if (c.getColor(i, j) == ofColor({ 255, 0, 255 })) {
					int cX = i;
					int cY = j;
					if (cX < w - 1) {
						neighbors.push_back(make_pair(cX + 1, cY));
					}
					if (cX > 1) {
						neighbors.push_back(make_pair(cX - 1, cY));
					}
					if (cY < h - 1) {
						neighbors.push_back(make_pair(cX, cY + 1));
					}
					if (cY > 1) {
						neighbors.push_back(make_pair(cX, cY - 1));
					}
					while (neighbors.size() > 0) {
						cX = neighbors[neighbors.size() - 1].first;
						cY = neighbors[neighbors.size() - 1].second;
						neighbors.pop_back();
						if (c.getColor(cX, cY) != eColor) {
							neighbors.clear();
							changeColor.push_back(make_tuple(i, j, c.getColor(cX, cY)));
							removedLine = true;
						}

					}
				}
			}
		}
		while (changeColor.size() > 0) {
			auto pos = changeColor[changeColor.size() - 1];
			c.setColor(get<0>(pos), get<1>(pos), get<2>(pos));
			changeColor.pop_back();
		}
	} while (removedLine == true);
	vector<vector<int>> adjacency(cColor+1, vector<int>(cColor+1));
	
	for (int j = 1; j < h - 1; ++j) {
		for (int i = 1; i < w - 1; ++i) {
			int g1 = c.getColor(i, j).g;
			vector<pair<int, int>> neighbors;
			int cX = i;
			int cY = j;
			if (cX < w - 1) {
				neighbors.push_back(make_pair(cX + 1, cY));
			}
			if (cX > 1) {
				neighbors.push_back(make_pair(cX - 1, cY));
			}
			if (cY < h - 1) {
				neighbors.push_back(make_pair(cX, cY + 1));
			}
			if (cY > 1) {
				neighbors.push_back(make_pair(cX, cY - 1));
			}
			while (neighbors.size() > 0) {
				cX = neighbors[neighbors.size() - 1].first;
				cY = neighbors[neighbors.size() - 1].second;
				neighbors.pop_back();
				int g2 = c.getColor(cX, cY).g;
				if (g2 != g1) {
					adjacency[g1][g2] = 1;
					adjacency[g2][g1] = 1;
				}

			}
		}
	}
	c.update();
	cout << "done with adjacency" << endl;
	auto fColors = fourColor(adjacency);
	for (auto c : fColors) {
		cout << c << endl;
	}
	int maxColor = 0;
	for (int i = 0; i < fColors.size(); ++i) {
		if (fColors[i] > maxColor) {
			maxColor = fColors[i];
		}
	}
	vector<ofColor> tileColors;
	tileColors.clear();
	for (int i = 0; i < maxColor+10; ++i) {
		tileColors.push_back(ofColor::fromHsb((i*60)%255, ofRandom(200,255), 255));
	}
	auto white = ofColor();
	white.r = 255;
	white.g = 255;
	white.b = 255;
	tileColors[0] = white;

	for (int i = 0; i < w; ++i) {
		for (int j = 0; j < h; ++j) {
			int idx = fColors[int(c.getColor(i, j).g)];
			c.setColor(i, j, tileColors[idx]);
		}
	}
	c.update();
	c.save("test.png");
}
vector<int> ofApp::fourColor(vector<vector<int>> adjacency) {
	cout << "four color" << endl;
	int colors = adjacency[0].size();
	vector<pair<int,vector<int>>> removed;
	vector <int> degree;
	vector<int> vertexColor(colors);
	for (int i = 1; i < colors; ++i) {
		int d = 0;
		for (int j = 1; j < colors; ++j) {
			if (adjacency[i][j] > 0) {
				d += 1;
			}
		}
		degree.push_back(d);
	}
	int r = 0;
	vector<int> removedIndices;
	removedIndices.clear();
	bool foundOne = true;
	while (foundOne ==true ) {
		foundOne = false;
		for (int i = 1; i < colors; ++i) {
			if (degree[i] < 5) {
				bool good = true;
				for (int j = 0; j < removedIndices.size(); ++j) {
					if (removedIndices[j] == i) {
						good = false;
						break;
					}
				}
				if (!good) {
					continue;
				}
				removed.push_back(make_pair(i,adjacency[i]));
				removedIndices.push_back(i);
				foundOne = true;
				r += 1;
				for (int j = 1; j < colors; ++j) {
					int con = adjacency[i][j];
					if (con > 0) {
						degree[j] -= 1;
						adjacency[i][j] = 0;
						adjacency[j][i] = 0;
					}
				}
			}
		}
	}
	cout << "removed:" << r <<" circles: " << circles.size() << endl;
	for (auto node : removed) {
		int ind = node.first;
		vector<int> putback = node.second;
		for (int j = 1; j < colors; ++j) {
			int con = putback[j];
			if (con > 0) {
				adjacency[ind][j] = 1;
				adjacency[j][ind] = 1;
			}
		}
		int minColor = 0;
		vector<int> usedColors;
		usedColors.clear();
		for (int j = 0; j < colors; ++j) {
			if (adjacency[ind][j] > 0) {
					usedColors.push_back(vertexColor[j]);
			}
		}
		for(int ii=1;ii<20;++ii){
		if (std::find(usedColors.begin(), usedColors.end(), ii) == usedColors.end()) {
			minColor = ii;
			break;
		}
		}
		vertexColor[ind] = minColor;

	}
	return vertexColor;
}
void ofApp::exportPDF(string filename)
{
	ofFbo offscreen;
	offscreen.allocate(1000, 1000);
	offscreen.begin();
	ofBeginSaveScreenAsPDF(filename + to_string(edges.size())+".pdf");
	c.update();
	c.save(filename+to_string(edges.size())+".png");
	ofSetColor({ 255,0,255 });
	for (auto e : edges) {
		e.draw();
	}
	ofEndSaveScreenAsPDF();
	offscreen.end();
}
//--------------------------------------------------------------
void ofApp::update(){

}

//--------------------------------------------------------------
static char str0[128] = "1";
bool addedColor = false;
bool showColor = false;

void ofApp::draw(){
	if(!addedColor){
	drawEdges();
	}

	gui.begin();
	{
		if (ImGui::Button("Clear")) {
			edges.clear();
			circles.clear();
			boundingBoxes.clear();
			drawBoundingBoxes = false;
			addedColor = false;
			c.clear();
			showColor = false;

		}
		if (ImGui::Button("scatter circle")) {
			scatterCircles();
		}
		if (ImGui::Button("infinity circle")) {
			infinity();
		}
		if (ImGui::Button("gen grid")) {
			genGrid();
		}
		if (ImGui::Button("Concentric")) {
			concentric();
		}
		ImGui::InputText("filename", str0, IM_ARRAYSIZE(str0));
		if (ImGui::Button("Save")) {
			exportPDF(string(str0));
		}
		ImGui::SliderInt("Max iterations", &maxIter, 50, 5000);
		if (ImGui::Button("Add Circles")) {
			addLotsofCircles();
		}
		ImGui::SliderInt("Scallop Size", &scallopR, 100, 400);
		if (ImGui::Button("Scallop")) {
			makeScallop();
		}
		ImGui::InputFloat("overlap", &overlapRatio);
		ImGui::InputFloat("min radius", &minRadius);
		ImGui::InputFloat("max radius", &maxRadius);
		ImGui::InputInt("Radius", & scallopRadius);
		ImGui::InputFloat("x Spacing", & xSpacing);
		ImGui::InputFloat("y Spacing", & ySpacing);
		if (ImGui::Button("Regular Scallop")) {
			regularScallop();
		}
		
		if (ImGui::Button("Calculate Intersections")) {
			intersections();
		}
		if (ImGui::Button("Clean")) {
			clean();
		}
		if (ImGui::Button("draw knob bounding boxes")) {
			drawBoundingBoxes = !drawBoundingBoxes;
		}
		if(drawBoundingBoxes){
			ofSetColor({ 255,0,0 });
			for (auto knobRec : boundingBoxes) {
				ofPolyline testRect;
				testRect.addVertex(knobRec.getTopLeft());
				testRect.addVertex(knobRec.getTopRight());
				testRect.addVertex(knobRec.getBottomRight());
				testRect.addVertex(knobRec.getBottomLeft());
				testRect.addVertex(knobRec.getTopLeft());
				testRect = testRect.getResampledByCount(40);
				testRect.draw();
			}
		}
		if (ImGui::Button("Add knobs")) {
			addKnobs();
		}
		if (ImGui::Button("Add Color")) {
			if(addedColor == false){
				

				cout << "adding color";
			randomColor();
			}
			addedColor = true;
			showColor = !showColor;
		}
		if (showColor) {
			ofSetColor(0xffffff);
			c.draw(0, 0);
		}
	}
	gui.end();
}
void ofApp::drawKnobs() {
	int x = 0;
	int y = 100;
	for (auto path : knobPaths) {
		path.setStrokeWidth(1);
		path.setStrokeColor({ 255,0,0 });
		path.draw(x, y);
		x += 100;
		if (x > 1000) {
			x = 0;
			y += 100;
		}
	}
}
//--------------------------------------------------------------
void ofApp::keyPressed(int key){

}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y){

}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y){

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){

}
bool sortbyx(const pair<pair<ofPoint, ofPoint>, int>& a, const pair<pair<ofPoint, ofPoint>, int>& b) {
	return (a.first.first.x < b.first.first.x);
}

bool sortCutPoints(const pair<int,ofPoint>& a, const pair<int, ofPoint>& b) {
	return (a.first < b.first);
}

bool ofApp::intersections(ofPolyline & e1, int index) 
	{
		
		auto bbox1 = e1.getBoundingBox();

		for (int j =0; j < edges.size(); ++j) {
			if (j == index) {
				continue;
			}
			auto e2 = &edges[j];
			auto bbox2 = e2->getBoundingBox();
			if (bbox1.intersects(bbox2)) {
				auto v1 = e1.getVertices();
				auto v2 = e2->getVertices();
				if (v2.size() < 4 || v1.size()<4) {
					continue;
				}
				int extra1 = 1;
				int extra2 = 1;
				if (e1.isClosed()) {
					extra1 = 0;
				}
				if (e2->isClosed()) {
					extra2 = 0;
				}
				
				for (int k = 3; k < v1.size() - extra1-3; ++k) {
					auto va = v1[k];
					auto vb = v1[(k + 1) % v1.size()];
					for (int l = 3; l < v2.size() - extra1-3; ++l) {
						auto vc = v2[l];
						auto vd = v2[(l + 1) % v2.size()];
						bool inters = doIntersect(va, vb, vc, vd);
						if (inters) {
							
							return true;
						}

					}
				}
			}



		}
		return false;
	}
	
void ofApp::clean() {
	/*
		if (iter->getLengthAtIndex(iter->size() - 1) < 1) {
			iter = edges.erase(iter);
		}
		else {
			++iter;
		}
	}
	*/
for (auto e = edges.begin(); e != edges.end();) {
	bool erased = false;
		if (e->getLengthAtIndex(e->size()-1) < 10) {
			cout << "found an edge to test" << endl;
			ofPoint p0 = e->getVertices()[0];
			ofPoint p1 = e->getVertices()[e->size() - 1];
			float minD0 = 1000;
			float minD1 = 1000;
			for (auto e1 : edges) {
				ofPoint pp = e1.getVertices()[0];
				ofPoint pp1 = e1.getVertices()[e1.size() - 1];
				if ( pp== p0 &&  pp1== p1) {
					continue;
				}
				ofPoint p0a = e1.getClosestPoint(p0);
				ofPoint p1a = e1.getClosestPoint(p1);
				float cMin0 = p0a.distance(p0);
				float cMin1 = p1a.distance(p1);
				minD0 = min(cMin0, minD0);
				minD1 = min(cMin1, minD1);
				if (cMin0 < 1 && cMin1 < 1) {
					e = edges.erase(e);
					//cout << "deleted an edge for completely overlapping" << endl;
					erased = true;
					break;
				}

			}
			if (minD0 > 1 || minD1 > 1 && erased ==false) {
				e = edges.erase(e);
				erased = true;
				//cout << "deleted an edge for dangling" << endl;
			}



		}
		if (erased == false) { ++e; }
	}
}
void ofApp::intersections() {
	map<int, vector<pair<int, ofPoint>>> cutPointMap;
	for (int i = 0; i < edges.size(); ++i) {
		auto e1 = &edges[i];
		auto bbox1 = e1->getBoundingBox();

		for (int j = i + 1; j < edges.size(); ++j) {
			auto e2 = &edges[j];
			auto bbox2 = e2->getBoundingBox();
			if (bbox1.intersects(bbox2)) {
				
				auto v1 = e1->getVertices();
				auto v2 = e2->getVertices();
				int extra1 = 1;
				int extra2 = 1;
				if (e1->isClosed()) {
					extra1 = 0;
				} 
				if (e2->isClosed()) {
					extra2 = 0;
				}

				for (int k = 0; k < v1.size() - extra1; ++k) {
					if (!(e1->isClosed())) {
						if (k < 2 || k>v1.size()-3) {
							continue;
						}
					}
					auto va = v1[k];
					auto vb = v1[(k + 1) % v1.size()];
					for (int l = 0; l < v2.size() - extra1; ++l) {
						if (!(e2->isClosed())) {
							if (l < 2 || l>v2.size() - 3) {
								continue;
							}
						}
						auto vc = v2[l];
						auto vd = v2[(l + 1) % v2.size()];
						bool intersects = doIntersect(va, vb, vc, vd);
						
						if (intersects) {
							ofPoint intersection = segsegintersect(va, vb, vc, vd);
							//intersection = segIntersection(va, vb, vc, vd);
							if(intersection.x>=min(va.x,vb.x)-5 && intersection.x>=min(vc.x,vd.x) - 5 &&
								intersection.x<= max(va.x, vb.x)+5&& intersection.x<= max(vc.x, vd.x) + 5 &&
								intersection.y >=min(va.y, vb.y) - 5 && intersection.y >= min(vc.y, vd.y) - 5 &&
								intersection.y <= max(va.y, vb.y) + 5 && intersection.y<= max(vc.y, vd.y) + 5
								){ 
							if (cutPointMap.count(i)) {
								cutPointMap[i].push_back(make_pair(k, intersection));
							}
							else {
								vector<pair<int, ofPoint>> p;
								p.push_back(make_pair(k, intersection));
								cutPointMap[i] = p;
							}
							if (cutPointMap.count(j)) {
								cutPointMap[j].push_back(make_pair(l, intersection));
							}
							else {
								vector<pair<int, ofPoint>> p;
								p.push_back(make_pair(l, intersection));
								cutPointMap[j] = p;
							}
						}
						}

					}
				}
			}



		}

	}

	std::map<int, vector<pair<int, ofPoint>>>::iterator it = cutPointMap.begin();
	
	while (it != cutPointMap.end()) {
		int idx = it->first;
		auto points = it->second;
		sort(points.begin(), points.end(), sortCutPoints);
		
		
		it++;
	}
	it = cutPointMap.begin();
	while (it != cutPointMap.end()) {
		int idx = it->first;
		auto points = it->second;
		it++;
		sort(points.begin(), points.end(), sortCutPoints);
		auto e = edges[idx];
		auto vertices = e.getVertices();

		if (e.isClosed()) {
			ofPolyline newEdge;
			newEdge.setClosed(false);
			auto p1 = points.back();
			auto p2 = points[0];
			newEdge.addVertex(p1.second);
			for (int i = p1.first; i < vertices.size(); ++i) {
				newEdge.addVertex(vertices[i]);
			}
			for (int i = 0; i <= p2.first; ++i) {
				newEdge.addVertex(vertices[i]);
			}
			newEdge.addVertex(p2.second);
			edges.push_back(newEdge);
		}
		else {
			ofPolyline newEdge;
			newEdge.setClosed(false);
			auto p1 = points[0];

			for (int i = 0; i < points[0].first; ++i) {
				newEdge.addVertex(vertices[i]);
			}
			newEdge.addVertex(p1.second);
			edges.push_back(newEdge);
		}

		for (int s = 0; s < points.size() - 1; ++s) {
ofPolyline newEdge;
auto p1 = points[s];
auto p2 = points[s + 1];
newEdge.addVertex(p1.second);
for (int x = p1.first; x < p2.first; ++x) {
	newEdge.addVertex(vertices[x]);
}
newEdge.addVertex(p2.second);
edges.push_back(newEdge);
		}
		if (!e.isClosed()) {
			ofPolyline newEdge;
			newEdge.addVertex(points.back().second);
			for (int s = points.back().first; s < vertices.size(); ++s) {
				newEdge.addVertex(vertices[s]);
			}
			edges.push_back(newEdge);
		}

	}
	it = cutPointMap.begin();
	vector <int> indices;
	while (it != cutPointMap.end()) {
		int idx = it->first;
		indices.push_back(idx);
		it++;
	}
	sort(indices.begin(), indices.end());
	reverse(indices.begin(), indices.end());
	for (auto i : indices) {
		edges.erase(edges.begin() + i);
	}
	

}
void ofApp::loadKnobs() {
	string inputPath = "svgs";
	ofDirectory dataDirectory(inputPath);
	auto files = dataDirectory.getFiles();
	for (size_t i = 0; i < files.size(); i++)
	{
		if (files[i].getExtension() == "svg") {
			auto svg = ofxSVG();
			knobNames.push_back(files[i].getFileName());
			knobsvgs.push_back(svg);
			svg.load(files[i].getAbsolutePath());
			auto path = svg.getPathAt(0);
			path.setStrokeWidth(1);
			auto p = path.getOutline();
			auto v = p[0].getVertices();


			knobs.push_back(p[0]);
			knobPaths.push_back(path);
		}
	}
}
void ofApp::generateEdge()
{
	ofPolyline edge;
	edge.addVertex(50, 50);
	edge.addVertex(150, 150);
	edges.push_back(edge);

}
void ofApp::drawEdges() {
	int i = 0;
	for (auto edge : edges) {
		i += 1;
		auto c1 = 1234 * i * i % 255;
		auto c2 = 1234 * i * i * i % 255;
		auto c3 = 255. - (c1 + c2) / 2.;
		ofSetColor({ float(c1) ,float(c2),float(c3) });
		edge.draw();
	}
}
void ofApp::addCircle(i3tuple circle) {
	int x = get<0>(circle);
	int y = get<1>(circle);
	int radius = get<2>(circle);
	int segments = max(int(radius), 180);
	ofPolyline c;
	for (int i = 0; i < segments; ++i) {
		c.addVertex(ofPoint(x + radius * cos(2. * 3.14159 * float(i) / float(segments)), y + radius * sin(2. * 3.14159 * float(i) / float(segments))));
	}
	c.close();
	edges.push_back(c);
}

ofPolyline ofApp::addKnob(int index) {
	auto edge1 = edges[index];
	if (edge1.size() < 10) {
		return edge1;
	}
	
	if (edge1.getLengthAtIndex(edge1.size() - 1) < 10){
	return edge1;
}

	std::vector<int> y(knobUsage.size());
	std::size_t n(0);
	std::generate(std::begin(y), std::end(y), [&] { return n++; });
	std::sort(std::begin(y),
		std::end(y),
		[&](int i1, int i2) { return knobUsage[i1] < knobUsage[i2]; });
	ofLog(OF_LOG_NOTICE, "knobs");
	cout << "edge " << index <<" size: " << edge1.size() << endl;
	bool placed = false;
	for (auto i : y) {
		for (int j = 0; j < 100; ++j) {
			//cout << "start of trial index " << i << "try number " << j << endl;
			ofPolyline edge = edge1.getResampledByCount(edge1.getLengthAtIndex(edge1.size() - 1)*5);
			ofPolyline knob = knobs[i];
			auto v = knob.getVertices();
			vector<ofPoint> vertices;
			bool flip = rand() % 2;
			for (auto a : v) {
				if (flip) {
					a.y = -a.y;
				}
				vertices.push_back(a);

			}
			if (vertices.size() < 10) {
				continue;
			}
			if (rand() % 2) {
				reverse(vertices.begin(), vertices.end());
			}
			size_t n = vertices.size();
			float scaleKnobX =  (static_cast <float> (rand()) / static_cast <float> (RAND_MAX) -.5) * .1+1;
			float scaleKnobY = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX) - .5) * .1 + 1;
			for (auto v : vertices) {
				v.x = v.x * scaleKnobX;
				v.y = v.y * scaleKnobY;
			}
			ofPoint a1 = vertices[0];
			ofPoint a2 = vertices[n - 1];
			
			float fuzz = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX) - .5) * .3;
			
			float midIndex = edge.getIndexAtPercent(.5 + fuzz);
			float midLength = edge.getLengthAtIndexInterpolated(midIndex);
			float d = vertices[0].distance(vertices[n - 1]);
			if (d > edge.getLengthAtIndex(edge.size() - 1) / 1.5) {
				continue;
			}
			
			ofPoint p1 = edge.getPointAtLength(midLength - d / 2);
			ofPoint p2 = edge.getPointAtLength(midLength + d / 2);
			
			float i1 = edge.getIndexAtLength(midLength - d / 2);
			float i2 = edge.getIndexAtLength(midLength + d / 2);
			
			//ofLog(OF_LOG_NOTICE, "POINTS");
			//ofLog(OF_LOG_NOTICE, ofToString(p1[0]) + " " + ofToString(p1[1]));
			//ofLog(OF_LOG_NOTICE, ofToString(p2[0]) + " " + ofToString(p2[1]));
			//scale knob to match
			float actualDistance = p1.distance(p2);
			float startDistance = vertices[0].distance(vertices[n - 1]);
			float scaleFactor = actualDistance / startDistance;
			for (int i = 0; i < n; ++i) {
				vertices[i] *= scaleFactor;
			}
			//ofLog(OF_LOG_NOTICE, ofToString(scaleFactor));
			//rotate knob to match 
			float rotAngle = -atan2((p2.x - p1.x) * (a2.y - a1.y) - (p2.y - p1.y) * (a2.x - a1.x),
				(p2.x - p1.x) * (a2.x - a1.x) + (p2.y - p1.y) * (a2.y - a1.y));

			//ofLog(OF_LOG_NOTICE, "Rotation angle " + ofToString(rotAngle));
			for (int i = 0; i < n; ++i) {
				vertices[i] = vertices[i].rotateRad(rotAngle, ofVec3f(0, 0, 1));
			}
			ofVec3f translate = p1 - vertices[0];
			//translate the first point of the knob to our point
			for (int i = 0; i < n; ++i) {
				vertices[i] += translate;
			}
			vector<ofVec3f> newVerts;
			ofPolyline newEdge;
			ofPolyline testKnob;

			auto edgeVerts = edge.getVertices();
			
			for (int i = 0; i <= min(int(i1),int(edgeVerts.size())); ++i) {
				newEdge.addVertex(edgeVerts[i]);
			}
			
			for (int i = 0; i < n; ++i) {
				newEdge.addVertex(vertices[i]);
				testKnob.addVertex(vertices[i]);
				//ofLog(OF_LOG_NOTICE, ofToString(vertices[i][0]) + " " + ofToString(vertices[i][1]));

			}
			
			for (int i = min(int(ceil(i2)), int(edgeVerts.size()-1)); i < edgeVerts.size(); ++i) {
				newEdge.addVertex(edgeVerts[i]);
			}
		
			ofPolyline testRect;
			auto knobRec = testKnob.getBoundingBox();
			knobRec.standardize();
			knobRec.scaleFromCenter(1.5);
			testRect.addVertex(knobRec.getTopLeft());
			testRect.addVertex(knobRec.getTopRight());
			testRect.addVertex(knobRec.getBottomRight());
			testRect.addVertex(knobRec.getBottomLeft());
			testRect.addVertex(knobRec.getTopLeft());
			testRect = testRect.getResampledByCount(300);
			testKnob.clear();
			
			if (!intersections(testRect, index)) {
				
				testRect.clear();
				testRect.addVertex((knobRec.getTopLeft() + knobRec.getTopRight()) / 2);
				testRect.addVertex((knobRec.getBottomLeft() + knobRec.getBottomRight()) / 2);
				testRect = testRect.getResampledByCount(100);
				if (!intersections(testRect, index)) {

					
					knobUsage[i] += 1;
					boundingBoxes.push_back(knobRec);
					return newEdge;

				}
			}
		}
	}
	cout << "failure" << endl;
	return edge1.getResampledByCount(500);
}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}
