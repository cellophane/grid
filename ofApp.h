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
		int scallopRadius = 40;
		float xSpacing = 1.8;
		float ySpacing = 1.8;
		void makePieces();
		double threePointAngle(ofPoint P1, ofPoint P2, ofPoint P3);
		void setup();
		void makeTextures();
		void drawPieces();
		bool floodFillTest(i3tuple circle, vector<i3tuple>& circles);
		int overlapped(i3tuple c1, i3tuple c2);
		int overlapped1(i3tuple c1, i3tuple c2);
		bool testCircleCheap(i3tuple circle, vector<i3tuple>& circles);
		bool testCircle(i3tuple circle, vector<i3tuple>& circles);
		bool testCircle(i3tuple circle);
		float circleOverlapChord(const i3tuple c1, const i3tuple c2);
		float circleOverlap(i3tuple c1, i3tuple c2);
		void scatterCircles();
		vector<ofPolyline> scatterAddCircle(i3tuple circle);
		void infinity();
		void regularScallop();
		void genGrid();
		void makeScallop();
		bool noCircleIntersect(i3tuple c1);
		ofPoint segEdgesIntersect(pair<ofPoint, ofPoint> segment);
		void addLotsofCircles();
		void packCircles();
		void randomColor();
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
		vector<string> knobNames;
		vector<vector<pair<int, bool>>> pieces;
		vector<array<float, 2> > pts;
		vector<vector<int> > lines;
		bool doDrawPieces = false;
		void addCircle(i3tuple circle);
		void concentric();
		void intersections();
		bool intersections(ofPolyline & testEdge, int index);
		void addColor();
		vector<int> fourColor(vector<vector<int>> adjacency);
		void exportPDF(string filename);
		void clean();
};
