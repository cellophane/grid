#pragma once

#include "ofMain.h"
#include "ofxImGui.h"
#include "ofxSVG.h"
class ofApp : public ofBaseApp{

	public:
		void setup();
		void update();
		void draw();

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
		ofPolyline addKnob(ofPolyline& edge);
		vector<ofPolyline> edges;
		ofxImGui::Gui gui;
		vector<ofPath> knobPaths;
		vector<ofxSVG> knobsvgs;
};
