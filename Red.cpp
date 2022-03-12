// The sample model.  You should build a file
// very similar to this for when you make your model.
#include "modelerview.h"
#include "modelerapp.h"
#include "modelerdraw.h"
#include "modelerui.h"
#include <FL/gl.h>

#include "modelerglobals.h"

//For Animation
int degree_inc = 0;

// To make a Unicycle, we inherit off of ModelerView
class Unicycle : public ModelerView
{
public:
	Unicycle(int x, int y, int w, int h, char* label)
		: ModelerView(x, y, w, h, label) { }

	virtual void draw();
};

// We need to make a creator function, mostly because of
// nasty API stuff that we'd rather stay away from.
ModelerView* createSampleModel(int x, int y, int w, int h, char* label)
{
	return new Unicycle(x, y, w, h, label);
}

// We are going to override (is that the right word?) the draw()
// method of ModelerView to draw out Unicycle
void Unicycle::draw()
{
	// This call takes care of a lot of the nasty projection 
	// matrix stuff.  Unless you want to fudge directly with the 
	// projection matrix, don't bother with this ...
	ModelerView::draw();

	if (ModelerUserInterface::m_controlsAnimOnMenu->value() != 0) degree_inc += 6;
	else degree_inc = 0;

	glPushMatrix();
	glScaled(1, -1, 1);
	glRotated(90, 1, 0, 0);

	// draw the floor
	/*setAmbientColor(1.0f, 1.0f, 1.0f);
	setDiffuseColor(0,0,1);
	glPushMatrix();
	drawCylinder(0.01,7,7);
	setDiffuseColor(1,1,0);
	drawCylinder(0.01,5,5);
	glPopMatrix();*/

	// draw the unicycle
	setAmbientColor(1.0f,1.0f,1.0f);
	setDiffuseColor(0.4,0.4,0.4);
	glTranslated(VAL(XPOS), VAL(YPOS), VAL(ZPOS));
	glScaled(SCALE_FACTOR, SCALE_FACTOR, SCALE_FACTOR);

		//draw the axle
		glPushMatrix();
		glRotated(VAL(AXLE_DIR), 0, 0, 1);
		glRotated(90, 0, 1, 0);
		drawCylinder(AXLE_LENGTH_HALF, AXLE_RADIUS, AXLE_RADIUS);
		glScaled(-1, 1, -1);
		drawCylinder(AXLE_LENGTH_HALF, AXLE_RADIUS, AXLE_RADIUS);
		glScaled(-1, 1, -1);
		glRotated(-90, 0, 1, 0);

		if (VAL(LEVEL_OF_DETAILS) >= 2) {

			//draw the wheel
			setDiffuseColor(0.05, 0.05, 0.05);
			glPushMatrix();
			glRotated(-VAL(WHEEL_DIR), 1, 0, 0);
			glTranslated(-WHEEL_WIDTH / 2, 0, 0);
			glRotated(90, 0, 1, 0);
			drawCylinder(WHEEL_WIDTH, WHEEL_RADIUS, WHEEL_RADIUS);
			glPopMatrix();

			//draw the left crank
			setDiffuseColor(0.4, 0.4, 0.4);
			glPushMatrix();
			glTranslated(-CRANK_DIST, 0, 0);
			glRotated(-VAL(CRANK_DIR) - 90 - degree_inc, 1, 0, 0);
			drawCylinder(CRANK_LENGTH, CRANK_RADIUS_ROOT, CRANK_RADIUS_TOP);
			glTranslated(0, 0, CRANK_LENGTH);
			glRotated(VAL(CRANK_DIR) + 90 + degree_inc, 1, 0, 0);

			if (VAL(LEVEL_OF_DETAILS) >= 3) {

				//draw the left pedal
				setDiffuseColor(0.05, 0.05, 0.05);
				glPushMatrix();
				glTranslated(-PEDAL_LENGTH, -PEDAL_WIDTH - PEDAL_GAP / 2, -PEDAL_HEIGHT / 2);
				drawBox(PEDAL_LENGTH, PEDAL_WIDTH, PEDAL_HEIGHT);
				glTranslated(0, PEDAL_WIDTH + PEDAL_GAP, 0);
				drawBox(PEDAL_LENGTH, PEDAL_WIDTH, PEDAL_HEIGHT);
				glPopMatrix();

			}

			glPopMatrix();

			//draw the right crank
			setDiffuseColor(0.4, 0.4, 0.4);
			glPushMatrix();
			glTranslated(CRANK_DIST, 0, 0);
			glRotated(-VAL(CRANK_DIR) + 90 - degree_inc, 1, 0, 0);
			drawCylinder(CRANK_LENGTH, CRANK_RADIUS_ROOT, CRANK_RADIUS_TOP);
			glTranslated(0, 0, CRANK_LENGTH);
			glRotated(VAL(CRANK_DIR) - 90 + degree_inc, 1, 0, 0);

			if (VAL(LEVEL_OF_DETAILS) >= 3) {

				//draw the right pedal
				setDiffuseColor(0.05, 0.05, 0.05);
				glPushMatrix();
				glTranslated(0, -PEDAL_WIDTH - PEDAL_GAP / 2, -PEDAL_HEIGHT / 2);
				drawBox(PEDAL_LENGTH, PEDAL_WIDTH, PEDAL_HEIGHT);
				glTranslated(0, PEDAL_WIDTH + PEDAL_GAP, 0);
				drawBox(PEDAL_LENGTH, PEDAL_WIDTH, PEDAL_HEIGHT);
				glPopMatrix();

			}

			glPopMatrix();

			//draw the tube
			//left lower tube
			setDiffuseColor(1, 0, 0);
			glPushMatrix();
			glRotated(-VAL(TUBE_DIR), 1, 0, 0);
			glTranslated(-TUBE_DIST - TUBE_WIDTH / 2, -TUBE_LENGTH / 2, -AXLE_RADIUS);
			drawBox(TUBE_WIDTH, TUBE_LENGTH, TUBE_HEIGHT_LOWER);
			glPopMatrix();

			//right lower tube
			glPushMatrix();
			glRotated(-VAL(TUBE_DIR), 1, 0, 0);
			glTranslated(TUBE_DIST - TUBE_WIDTH / 2, -TUBE_LENGTH / 2, -AXLE_RADIUS);
			drawBox(TUBE_WIDTH, TUBE_LENGTH, TUBE_HEIGHT_LOWER);
			glPopMatrix();

			//middle tube
			glPushMatrix();
			glRotated(-VAL(TUBE_DIR), 1, 0, 0);
			glTranslated(-TUBE_DIST - TUBE_WIDTH / 2, -TUBE_LENGTH / 2, -AXLE_RADIUS + TUBE_HEIGHT_LOWER);
			drawBox(TUBE_WIDTH + TUBE_DIST * 2, TUBE_LENGTH, TUBE_HEIGHT_MIDDLE);
			glPopMatrix();

			//upper tube
			glPushMatrix();
			glRotated(-VAL(TUBE_DIR), 1, 0, 0);
			glTranslated(0, 0, TUBE_HEIGHT_LOWER + TUBE_HEIGHT_MIDDLE - AXLE_RADIUS);
			drawCylinder(TUBE_HEIGHT_UPPER, TUBE_RADIUS, TUBE_RADIUS);

			if (VAL(LEVEL_OF_DETAILS) >= 3) {

				//draw the seatpost
				setDiffuseColor(0.4, 0.4, 0.4);
				glPushMatrix();
				glTranslated(0, 0, TUBE_HEIGHT_UPPER);
				drawCylinder(VAL(SEATPOST_HEIGHT), SEATPOST_RADIUS, SEATPOST_RADIUS);

				if (VAL(LEVEL_OF_DETAILS) >= 4) {

					//draw the saddle
					setDiffuseColor(0.05, 0.05, 0.05);
					glPushMatrix();
					glTranslated(0, 0, VAL(SEATPOST_HEIGHT));
					glRotated(-VAL(SADDLE_DIR), 0, 0, 1);
					glTranslated(-SADDLE_WIDTH / 2, -SADDLE_DIST, 0);
					drawTriangle(0,0,0,
								 SADDLE_WIDTH,0,0,
								 SADDLE_WIDTH/2,SADDLE_LENGTH,0);
					drawTriangle(0, 0, SADDLE_HEIGHT,
								 SADDLE_WIDTH / 2, SADDLE_LENGTH, SADDLE_HEIGHT,
								 SADDLE_WIDTH, 0, SADDLE_HEIGHT);
					drawTriangle(0, 0, 0,
								 0, 0, SADDLE_HEIGHT,
								 SADDLE_WIDTH, 0, 0);
					drawTriangle(0, 0, SADDLE_HEIGHT,
								 SADDLE_WIDTH, 0, SADDLE_HEIGHT,
								 SADDLE_WIDTH, 0, 0);
					drawTriangle(0, 0, 0,
								 SADDLE_WIDTH / 2, SADDLE_LENGTH, 0,
								 SADDLE_WIDTH / 2, SADDLE_LENGTH, SADDLE_HEIGHT);
					drawTriangle(0, 0, 0,
								 SADDLE_WIDTH / 2, SADDLE_LENGTH, SADDLE_HEIGHT,
								 0, 0, SADDLE_HEIGHT);
					drawTriangle(SADDLE_WIDTH, 0, 0,
								 SADDLE_WIDTH / 2, SADDLE_LENGTH, SADDLE_HEIGHT,
								 SADDLE_WIDTH / 2, SADDLE_LENGTH, 0);
					drawTriangle(SADDLE_WIDTH, 0, 0,
								 SADDLE_WIDTH, 0, SADDLE_HEIGHT,
								 SADDLE_WIDTH / 2, SADDLE_LENGTH, SADDLE_HEIGHT);
					glPopMatrix();

				}

				glPopMatrix();

			}

			glPopMatrix();

		}

		glPopMatrix();

	glPopMatrix();
}

int main()
{
	// Initialize the controls
	// Constructor is ModelerControl(name, minimumvalue, maximumvalue, 
	// stepsize, defaultvalue)
	ModelerControl controls[NUMCONTROLS];
	controls[LEVEL_OF_DETAILS] = ModelerControl("Level of Details", 1, 4, 1, 4);
	controls[XPOS] = ModelerControl("X Position", -7, 7, 0.1f, 0);
	controls[YPOS] = ModelerControl("Y Position", -7, 7, 0.1f, 0);
	controls[ZPOS] = ModelerControl("Z Position", 0, 5, 0.1f, 0);
	controls[AXLE_DIR] = ModelerControl("Direction of Axle", -180, 180, 1, 0);
	controls[WHEEL_DIR] = ModelerControl("Direction of Wheel", 0, 360, 1, 0);
	controls[CRANK_DIR] = ModelerControl("Direction of Crank", 0, 360, 1, 0);
	//controls[PEDAL_DIR] = ModelerControl("Direction of Paddle", -180, 180, 1, 0);
	controls[TUBE_DIR] = ModelerControl("Direction of Tube", -60, 60, 1, 0);
	controls[SEATPOST_HEIGHT] = ModelerControl("Seatpost Height", SEATPOST_HEIGHT_MIN, SEATPOST_HEIGHT_MAX, 1, (SEATPOST_HEIGHT_MIN+SEATPOST_HEIGHT_MAX)/2);
	controls[SADDLE_DIR] = ModelerControl("Direction of Saddle", -45, 45, 1, 0);
	controls[LABEL_DIR] = ModelerControl("Direction of Label", -180, 180, 1, 0);

	ModelerApplication::Instance()->Init(&createSampleModel, controls, NUMCONTROLS);
	return ModelerApplication::Instance()->Run();
}