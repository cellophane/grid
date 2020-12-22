#include "ofApp.h"
ofPolyline knob;
//--------------------------------------------------------------
void ofApp::setup(){
	generateEdge();
	ofLog(OF_LOG_NOTICE, "Hello");
	knob = addKnob(edges[0]);
}

//--------------------------------------------------------------
void ofApp::update(){

}

//--------------------------------------------------------------
void ofApp::draw(){
	drawEdges();
	knob.draw();
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

void ofApp::generateEdge()
{
	ofPolyline edge;
	edge.addVertex(50, 50);
	edge.addVertex(100, 50);
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
	ofPoint* vertices = new ofPoint[5];
	edge = edge.getResampledByCount(500);
	vertices[0]=ofPoint(0, 0);
	vertices[1]=ofPoint(5, 0);
	vertices[2]=ofPoint(10, 30);
	vertices[3]=ofPoint(15, 0);
	vertices[4]=ofPoint(20, 0);
	ofPolyline knob;
	int n = 5;
	ofLog(OF_LOG_NOTICE,ofToString(n));	
	for (int i = 0; i < n;++i) {
	knob.addVertex(vertices[i]);
	}
	float midIndex = edge.getIndexAtPercent(.5);
	float midLength = edge.getLengthAtIndexInterpolated(midIndex);
	float d = vertices[0].distance(vertices[n-1]);
	ofPoint p1 = edge.getPointAtLength(midLength - d / 2);
	ofPoint p2 = edge.getPointAtLength(midLength + d / 2);
	ofLog(OF_LOG_NOTICE, "POINTS");
	ofLog(OF_LOG_NOTICE, ofToString(p1[0]) + " " + ofToString(p1[1]));
	ofLog(OF_LOG_NOTICE, ofToString(p2[0]) + " " + ofToString(p2[1]));
	//scale knob to match
	float actualDistance = p1.distance(p2);
	float startDistance = vertices[0].distance(vertices[n-1]);
	float scaleFactor = actualDistance / startDistance;
	for (int i = 0; i < n; ++i) {
		vertices[i] *= scaleFactor;
	}
	ofLog(OF_LOG_NOTICE, ofToString(scaleFactor));
	//rotate knob to match 
	float finalAngle = p1.angleRad(p2);
	float startAngle = vertices[0].angleRad(vertices[n-1]);
	float rotAngle = finalAngle - startAngle;
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
	for (int i = 0; i < n; ++i) {
		newEdge.addVertex(vertices[i]);
		ofLog(OF_LOG_NOTICE, ofToString(vertices[i][0]) + " " + ofToString(vertices[i][1]));
	
	}

	return newEdge;
	

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}
