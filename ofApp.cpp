#include "ofApp.h"
ofPolyline knob;
bool loadedKnobs = false;
typedef std::tuple<int, int, int> i3tuple;
int maxIter = 500;
vector<i3tuple> circles;
//--------------------------------------------------------------
void ofApp::setup() {
	gui.setup();
	srand(time(NULL));
	//addLotsofCircles();
	loadedKnobs = true;
	loadKnobs();
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
void ofApp::makeScallop() {
	int center = 400;
	auto circle = make_tuple(center, center, 25);
	circles.push_back(circle);
	addCircle(circle);
	float cmax = 50;
	for (int i = 0; i < 100; ++i) {
		float theta = 2 * 3.14159 * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		bool project = true;
		auto ray = make_pair(cos(theta) * cmax, sin(theta) * cmax);
		auto base = make_pair(cos(theta) * cmax*.9, sin(theta) * cmax*.9);
		while (project) {
			ofPolyline 
		}
	}
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
				if (smallPacked < 10) {
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
		bool ok = testCircle(circle, circles);
		if (ok) {
			auto circle1 = make_tuple(x, y, medR + 5);
			 ok = testCircle(circle1, circles);
			if ( ok){
				circle1 = make_tuple(x, y, medR - 5);
				 ok = testCircle(circle1, circles);
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
//-----------------------------------
void ofApp::addColor() {

}
//--------------------------------------------------------------
void ofApp::update(){

}

//--------------------------------------------------------------
void ofApp::draw(){

	drawEdges();
	gui.begin();
	{
		if (ImGui::Button("Clear")) {
			edges.clear();
			circles.clear();

		}
		ImGui::SliderInt("Max iterations", &maxIter, 50, 5000);
		if (ImGui::Button("Add Circles")) {
			addLotsofCircles();
		}
		if (ImGui::Button("Calculate Intersections")) {
			intersections();
		}

		if (ImGui::Button("Add knobs")) {
			addKnobs();
		}
		if (ImGui::Button("Add Color")) {
			addColor();
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
ofPoint segsegintersect(const ofPoint & a, const ofPoint & b, const ofPoint & c, const ofPoint & d) {
	float tnum = (a.x - c.x) * (c.y - d.y) - (a.y - c.y) * (c.x - d.x);
	float tdenom = (a.x - b.x) * (c.y - d.y) - (a.y - b.y) * (c.x - d.x);
	float unum = (a.x - b.x) * (a.y - c.y) - (a.y - b.y) * (a.x - c.x);
	float udenom = (a.x - b.x) * (c.y - d.y) - (a.y - b.y) * (c.x - d.x);
	if ((tnum > 0 && tdenom > 0) || (tnum < 0 && tdenom < 0)) {
		if (abs(tnum) < abs(tdenom)) {
			if ((unum > 0 && udenom > 0) || (unum < 0 && udenom < 0)) {
				if (abs(unum) < abs(udenom)) {
					float t = tnum / tdenom;
					return ofPoint(a.x + t * (b.x - a.x), a.y + t * (b.y - a.y));



				}
			}
		}
	}
	return ofPoint(-10000,-10000);
}
bool sortCutPoints(const pair<int,ofPoint>& a, const pair<int, ofPoint>& b) {
	return (a.first < b.first);
}
bool ofApp::rayIntersection(ofPoint p2) {
	ofPoint p1 = p2;
	ofPoint intersectionPoint;
	bool intersected = true;
	while(intersected == true){
		intersected = false;
	for (auto circle : circles) {
		float d = sqrt(pow((get<0>(circle) - p1.x), 2) + pow((get<1>(circle) - p1.y), 2));
		if (d < get<2>(circle)) {
			intersectionPoint = p1;
			intersected = true;
			break;
		}
	}
	p1 += p2;
	}
	ofPolyline ray;
	ray.addVertex(p1);
	ray.addVertex(p2);

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
						ofPoint intersection = segsegintersect(va, vb, vc, vd);
						if (intersection.x != -10000 && intersection.y != -10000) {
							
							return true;
						}

					}
				}
			}



		}
		return false;
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
						ofPoint intersection = segsegintersect(va, vb, vc, vd);
						if (intersection.x != -10000 && intersection.y != -10000) {
							
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
	string inputPath = "C:\\Users\\jessi\\Nervous System Dropbox\\Nervous System\\puzzles\\puzzle development\\chris yates\\Knobs\\SVGs\\svg resave";
	//inputPath = "C:\\Users\\jules\\Desktop\\svg resave";
	ofDirectory dataDirectory(inputPath);
	auto files = dataDirectory.getFiles();
	for (size_t i = 0; i < files.size(); i++)
	{
		if (files[i].getExtension() == "svg") {
			auto svg = ofxSVG();
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
	ofSetColor({ 255,0,0 });
	for (auto edge : edges) {
		edge.draw();
	}
}
void ofApp::addCircle(i3tuple circle) {
	int x = get<0>(circle);
	int y = get<1>(circle);
	int radius  = get<2>(circle);
	int segments = int(radius );
	ofPolyline c;
	for (int i = 0; i < segments; ++i) {
		c.addVertex(ofPoint(x + radius * cos(2. * 3.14159 * float(i) / float(segments)), y + radius * sin(2. * 3.14159 * float(i) / float(segments))));
	}
	c.close();
	edges.push_back(c);
}

ofPolyline ofApp::addKnob(int index) {
	auto edge1 = edges[index];
	ofLog(OF_LOG_NOTICE, "knobs");
	bool placed = false;
	for (int k = 0; k < 1000; ++k){
		
	ofPolyline edge = edge1.getResampledByCount(300);
	int i = rand() % knobs.size();
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
	int n = vertices.size();
	ofPoint a1 = vertices[0];
	ofPoint a2 = vertices[n-1];
	
	float fuzz = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)-.5)*.3;
	float midIndex = edge.getIndexAtPercent(.5+fuzz);
	float midLength = edge.getLengthAtIndexInterpolated(midIndex);
	float d = vertices[0].distance(vertices[n-1]);
	ofPoint p1 = edge.getPointAtLength(midLength - d / 2);
	ofPoint p2 = edge.getPointAtLength(midLength + d / 2);
	float i1 = edge.getIndexAtLength(midLength - d / 2);
	float i2 = edge.getIndexAtLength(midLength + d / 2);

	//ofLog(OF_LOG_NOTICE, "POINTS");
	//ofLog(OF_LOG_NOTICE, ofToString(p1[0]) + " " + ofToString(p1[1]));
	//ofLog(OF_LOG_NOTICE, ofToString(p2[0]) + " " + ofToString(p2[1]));
	//scale knob to match
	float actualDistance = p1.distance(p2);
	float startDistance = vertices[0].distance(vertices[n-1]);
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
		vertices[i]=vertices[i].rotateRad(rotAngle, ofVec3f(0, 0, 1));
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
	for (int i = 0; i <= i1; ++i) {
		newEdge.addVertex(edgeVerts[i]);
	}
	for (int i = 0; i < n; ++i) {
		newEdge.addVertex(vertices[i]);
		testKnob.addVertex(vertices[i]);
		//ofLog(OF_LOG_NOTICE, ofToString(vertices[i][0]) + " " + ofToString(vertices[i][1]));
	
	}
	
	for (int i = ceil(i2); i < edgeVerts.size(); ++i) {
		newEdge.addVertex(edgeVerts[i]);
	}
	ofPolyline testRect;
	auto knobRec = testKnob.getBoundingBox();
	knobRec.standardize();
	knobRec.scaleFromCenter(1.25);
	testRect.addVertex(knobRec.getTopLeft());
	testRect.addVertex(knobRec.getTopRight());
	testRect.addVertex(knobRec.getBottomRight());
	testRect.addVertex(knobRec.getBottomLeft());
	testRect.close();
	testRect = testRect.getResampledByCount(100);
	testKnob.clear();
	if(!intersections(testRect,index)){
		testRect.clear();
		testRect.addVertex((knobRec.getTopLeft() + knobRec.getTopRight()) / 2);
		testRect.addVertex((knobRec.getBottomLeft() + knobRec.getBottomRight()) / 2);
		testRect = testRect.getResampledByCount(100);
		if (!intersections(testRect, index)) {
			if (!intersections(newEdge, index)) {
				cout << "success" << endl;
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
