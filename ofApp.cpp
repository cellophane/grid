#include "ofApp.h"
ofPolyline knob;
bool loadedKnobs = false;
//--------------------------------------------------------------
void ofApp::setup(){
	generateEdge();
	ofLog(OF_LOG_NOTICE, "Hello");
	
	loadedKnobs = true;
	loadKnobs();
	//knob = addKnob(edges[0]);
	addCircle(400, 425, 100);
	addCircle(550, 350, 95);
	addCircle(200, 350, 80);
	addCircle(450, 450, 300);
	addCircle(600, 250, 75);
	intersections();
	vector<ofPolyline> newEdges;
	for (auto edge : edges) {
		auto e = addKnob(edge);
		newEdges.push_back(e);
	}
	edges = newEdges;
}

//--------------------------------------------------------------
void ofApp::update(){

}

//--------------------------------------------------------------
void ofApp::draw(){

	drawEdges();
	ofSetColor({ 255,255,0 });
	knob.draw();

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
void ofApp::intersections() {
	map<int, vector<pair<int,ofPoint>>> cutPointMap;
	cout << "edges : " << edges.size() << endl;
	for (int i = 0; i < edges.size(); ++i) {
		cout << "current first edge: " << i << endl;
		auto e1 = &edges[i];
		auto bbox1 = e1->getBoundingBox();

		for (int j = i+1; j < edges.size(); ++j) {
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
					auto va = v1[k];
					auto vb = v1[(k + 1) % v1.size()];
					for (int l = 0; l < v2.size() - extra1; ++l) {
						auto vc = v2[l];
						auto vd = v2[(l + 1) % v2.size()];
						ofPoint intersection = segsegintersect(va, vb, vc, vd);
						if (intersection.x != -10000 && intersection.y != -10000) {
							cout << "found intersection! " << intersection.x << ", " << intersection.y << endl;
							if (cutPointMap.count(i)) {
								cutPointMap[i].push_back(make_pair(k,intersection));
							}
							else {
								vector<pair<int,ofPoint>> p;
								p.push_back(make_pair(k,intersection));
								cutPointMap[i] = p;
							}
							if (cutPointMap.count(j)) {
								cutPointMap[j].push_back(make_pair(l,intersection));
							}
							else {
								vector<pair<int,ofPoint>> p;
								p.push_back(make_pair(l,intersection));
								cutPointMap[j] = p;
							}
						}

					}
				}
			}



		}

	}
	
	std::map<int, vector<pair<int,ofPoint>>>::iterator it = cutPointMap.begin();
	cout << "cut point map : " << endl;
	while (it != cutPointMap.end()) {
		int idx = it->first;
		auto points = it->second;
		sort(points.begin(), points.end(), sortCutPoints);
		cout << "index: " << idx << endl;
		cout << "total verts in edge: " << edges[idx].getVertices().size() << endl;
		for (auto p : points) {
			cout << "vertex index: " << p.first << endl;
			cout << "vertex coord: " << p.second.x << ", " << p.second.y << endl;
		}
		it++;
	}
	it = cutPointMap.begin();
	while (it != cutPointMap.end()) {
		int idx = it->first;
		auto points = it->second;
		it++;
		sort(points.begin(), points.end(),sortCutPoints);
		auto e = edges[idx];
		auto vertices = e.getVertices();
		
		if (e.isClosed()) {
			ofPolyline newEdge;
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
			auto p1 = points[0];
			
			for (int i = 0; i < points[0].first; ++i) {
				newEdge.addVertex(vertices[i]);
			}
			newEdge.addVertex(p1.second);
			edges.push_back(newEdge);
		}
		
		for (int s = 0;  s < points.size()-1; ++s) {
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
		if(!e.isClosed()){
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
	inputPath = "C:\\Users\\jules\\Desktop\\svg resave";
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
			cout << v.size() << endl;
			for (auto i : v) {
				cout << i.x << "," << i.y << endl;
			}
			
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
void ofApp::addCircle(float x, float y, float radius) {
	int segments = int(2 * 3.14 * radius * 2);
	ofPolyline c;
	for (int i = 0; i < segments; ++i) {
		c.addVertex(ofPoint(x + radius * cos(2. * 3.14159 * float(i) / float(segments)), y + radius * sin(2. * 3.14159 * float(i) / float(segments))));
	}
	c.close();
	edges.push_back(c);

}
ofPolyline ofApp::addKnob(ofPolyline& edge) {
	ofLog(OF_LOG_NOTICE, "knobs");
	edge = edge.getResampledByCount(1000);
	cout << "knobs";
	int i = rand() % knobs.size();
	ofPolyline knob = knobs[i];
	auto v = knob.getVertices();
	vector<ofPoint> vertices;
	for (auto a : v) {
		vertices.push_back(a);
		ofLog(OF_LOG_NOTICE, ofToString(a[0]) + " " + ofToString(a[1]));
	}
	/*
	ofPoint* vertices = new ofPoint[5];
	edge = edge.getResampledByCount(500);
	vertices[0]=ofPoint(0, 0);
	vertices[1]=ofPoint(5, 0);
	vertices[2]=ofPoint(10, 30);
	vertices[3]=ofPoint(15, 0);
	vertices[4]=ofPoint(20, 0);
		for (int i = 0; i < n;++i) {
	knob.addVertex(vertices[i]);
	}
	*/
	int n = vertices.size();
	ofPoint a1 = vertices[0];
	ofPoint a2 = vertices[n-1];
	ofLog(OF_LOG_NOTICE,ofToString(n));	

	float midIndex = edge.getIndexAtPercent(.5);
	float midLength = edge.getLengthAtIndexInterpolated(midIndex);
	float d = vertices[0].distance(vertices[n-1]);
	ofPoint p1 = edge.getPointAtLength(midLength - d / 2);
	ofPoint p2 = edge.getPointAtLength(midLength + d / 2);
	float i1 = edge.getIndexAtLength(midLength - d / 2);
	float i2 = edge.getIndexAtLength(midLength + d / 2);
	cout << i1 << "first seg idx" << endl;
	cout << i2 << "second seg idx" << endl;
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

	ofLog(OF_LOG_NOTICE, "Rotation angle " + ofToString(rotAngle));
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
	auto edgeVerts = edge.getVertices();
	for (int i = 0; i <= i1; ++i) {
		newEdge.addVertex(edgeVerts[i]);
	}
	for (int i = 0; i < n; ++i) {
		newEdge.addVertex(vertices[i]);
		//ofLog(OF_LOG_NOTICE, ofToString(vertices[i][0]) + " " + ofToString(vertices[i][1]));
	
	}
	cout << "i2 ceil " << i2 << "... newverts size" << newVerts.size() << endl;
	for (int i = ceil(i2); i < edgeVerts.size(); ++i) {
		newEdge.addVertex(edgeVerts[i]);
	}
	return newEdge;
	

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}
