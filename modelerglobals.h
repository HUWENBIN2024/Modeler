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
		WHEEL_DIR,
		CRANK_DIR,
			//PEDAL_ROTATE,
		TUBE_DIR,
			SEATPOST_HEIGHT,
				SADDLE_DIR,
					LABEL_DIR,
	NUMCONTROLS
};

//Hyperparameters for Designing Unicycle
const float SCALE_FACTOR = 0.02;

const float AXLE_RADIUS = 1;
const float AXLE_LENGTH_HALF = 30;

const float WHEEL_RADIUS = 80;
const float WHEEL_WIDTH = 20;

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
const float SEATPOST_HEIGHT_MIN = 10;
const float SEATPOST_HEIGHT_MAX = 70;

const float SADDLE_HEIGHT = 25;
const float SADDLE_LENGTH = 130;
const float SADDLE_WIDTH = 60;
const float SADDLE_DIST = 60; // distance from upper center of seatpost to back,bottom,center of saddle

const float LABEL_STRING_HEIGHT = 0;

// Colors
#define COLOR_RED		1.0f, 0.0f, 0.0f
#define COLOR_GREEN		0.0f, 1.0f, 0.0f
#define COLOR_BLUE		0.0f, 0.0f, 1.0f

// We'll be getting the instance of the application a lot; 
// might as well have it as a macro.
#define VAL(x) (ModelerApplication::Instance()->GetControlValue(x))

#endif