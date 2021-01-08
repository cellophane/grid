#pragma once

#include "ofMain.h"
#include "ofxImGui.h"
#include "ofxSVG.h"
typedef std::tuple<int, int, int> i3tuple;
class ofApp : public ofBaseApp{

	public:
		ofxImGui::Gui gui;
		ofFbo cBuffer;
		vector<int> knobUsage;
		void setup();
		bool floodFillTest(i3tuple circle, vector<i3tuple>& circles);
		int overlapped(i3tuple c1, i3tuple c2);
		bool testCircle(i3tuple circle, vector<i3tuple>& circles);
		float circleOverlapChord(const i3tuple c1, const i3tuple c2);
		float circleOverlap(i3tuple c1, i3tuple c2);
		void makeScallop();
		ofPoint segEdgesIntersect(pair<ofPoint, ofPoint> segment);
		void addLotsofCircles();
		void packCircles();
		void update();
		void draw();

		void drawKnobs();

		void keyPressed(int key);
		void keyReleased(int key);
		void mouseMoved(int x, int y );
		void mouseDragged(int x, int y, int button);
		void mousePressed(int x, int y, int button);
		void mouseReleased(int x, int y, int button);
		void mouseEntered(int x, int y);
		void mouseExited(int x, int y);
		void windowResized(int w, int h);
		void dragEvent(ofDragInfo dragInfo);
		void gotMessage(ofMessage msg);
		void generateEdge();
		void drawEdges();
		void loadKnobs();
		vector<ofPolyline> knobs;
		ofPolyline addKnob(int i);
		void addKnobs();
		vector<ofPolyline> edges;
		vector<ofPath> knobPaths;
		vector<ofxSVG> knobsvgs;
		void addCircle(i3tuple circle);
		void intersections();
		bool intersections(ofPolyline & testEdge, int index);
		void addColor();
};
