#include "ofApp.h"
ofPolyline knob;
bool loadedKnobs = false;
//--------------------------------------------------------------
void ofApp::setup(){
	generateEdge();
	ofLog(OF_LOG_NOTICE, "Hello");
	
	loadedKnobs = true;
	loadKnobs();
	knob = addKnob(edges[0]);
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
void ofApp::loadKnobs() {
	string inputPath = "C:\\Users\\jessi\\Nervous System Dropbox\\Nervous System\\puzzles\\puzzle development\\chris yates\\Knobs\\SVGs\\svg resave";
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

ofPolyline ofApp::addKnob(ofPolyline& edge) {
	ofLog(OF_LOG_NOTICE, "knobs");
	edge = edge.getResampledByCount(1000);
	cout << "knobs";
	ofPolyline knob = knobs[2];
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
