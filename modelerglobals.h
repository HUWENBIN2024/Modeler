#ifndef _MODELER_GLOBALS_H
#define _MODELER_GLOBALS_H

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502
#endif

// This is a list of the controls for the SampleModel
// We'll use these constants to access the values 
// of the controls from the user interface.
enum SampleModelControls
{ 
	LEVEL_OF_DETAILS,
	XPOS, YPOS, ZPOS, AXLE_DIR,
		WHEEL_DIR, WHEEL_QUALITY,
		CRANK_DIR,
			//PEDAL_ROTATE,
		TUBE_DIR,
			BEND,
				SADDLE_DIR, SADDLE_QUALITY,
					LABEL_DIR,
	HAPPINESS, 
	FRAME_ALL,
	L_SYSTEM,
	META_BALL_CONTROL,
	META_BALL1_RADIUS,
	META_BALL2_RADIUS,
	META_BALL3_RADIUS,
	META_BALL4_RADIUS,
	META_BALL_SCALE,
	IK_VISUALIZE,
	IK_TARGETX, IK_TARGETY, IK_TARGETZ,
	IK_MODE,
	IK_MOVE,
	NUMCONTROLS
};

//Hyperparameters for Designing Unicycle
const float SCALE_FACTOR = 0.02;

const float AXLE_RADIUS = 1;
const float AXLE_LENGTH_HALF = 30;

const float WHEEL_RADIUS = 80;
const float WHEEL_WIDTH = 20;
const float SPOKE_LENGTH = 75;
const float SPOKE_RADIUS = 0.6;
const float SPOKE_DIST = 4;

const float CRANK_DIST = 27; //distance to the center of the wheel
const float CRANK_LENGTH = 70;
const float CRANK_RADIUS_ROOT = 2;
const float CRANK_RADIUS_TOP = 1.5;

const float PEDAL_LENGTH = 30;
const float PEDAL_WIDTH = 15;
const float PEDAL_HEIGHT = 10;
const float PEDAL_GAP = 2;

const float TUBE_DIST = 18; //distance from center of lower tube to center of wheel
const float TUBE_LENGTH = 20;
const float TUBE_WIDTH = 10;
const float TUBE_HEIGHT_LOWER = 95;
const float TUBE_HEIGHT_MIDDLE = 15;
const float TUBE_HEIGHT_UPPER = 60;
const float TUBE_RADIUS = 10;

const float SEATPOST_RADIUS = 7;
const float SEATPOST_HEIGHT = 40;

const float SADDLE_HEIGHT = 15;
const float SADDLE_LENGTH = 130;
const float SADDLE_WIDTH = 60;
const float SADDLE_DIST = 60; // distance from upper center of seatpost to back,bottom,center of saddle
const float SPRING_X = 17;
const float SPRING_Y = 40; // distance to center of saddle
const float SPRING_CYLINDER_RADIUS = 4;
const float SPRING_CYLINDER_HEIGHT = 50;
const float SPRING_RADIUS_BIG[9] = {8,10,11.6,13,14,13,11.6,10,8};
const float SPRING_RADIUS_SMALL = 2.5;
const float SPRING_DIST_APART = 5.1; //distance between the centers of the springs

const float LABEL_STRING_HEIGHT = 40;
const float LABEL_STRING_RADIUS = 0.6;
const float LABEL_LENGTH = 18;
const float LABEL_HEIGHT = 30;
const float LABEL_TEXTURE_EPSILON = 0.01; //Texture should be "slightly" (i.e. by epsilon) in front of the label box

// Colors
#define COLOR_RED		1.0f, 0.0f, 0.0f
#define COLOR_GREEN		0.0f, 1.0f, 0.0f
#define COLOR_BLUE		0.0f, 0.0f, 1.0f

// We'll be getting the instance of the application a lot; 
// might as well have it as a macro.
#define VAL(name) (ModelerApplication::Instance()->GetControlValue(name))
#define SET(name, val) (ModelerApplication::Instance()->SetControlValue(name, val))

#endif