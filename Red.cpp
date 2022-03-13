// The sample model.  You should build a file
// very similar to this for when you make your model.
#include "modelerview.h"
#include "modelerapp.h"
#include "modelerdraw.h"
#include "modelerui.h"
#include <FL/gl.h>

#include "modelerglobals.h"
#define PI 3.1415926

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

GLfloat point_cloud[130][60] = {};
GLfloat control_point[5] = {0 , 50, -20, 25, 0};

GLfloat bezier(GLfloat l)
{
	GLfloat t = (l / SADDLE_LENGTH);
	GLfloat b = 0;
	int coef[5] = { 1, 4, 6, 4, 1 };
	for (int i = 0; i < 5; i++)
	{
		b += control_point[i] * coef[i] * pow(t, i) * pow(1-t, 4-i);
	}
	return b;
}

void cal_points()
{

	for (int l = 0; l < SADDLE_LENGTH; l++)
	{
		GLfloat b = bezier(l);
		GLfloat w = (SADDLE_LENGTH - l) / SADDLE_LENGTH * SADDLE_WIDTH;
		GLfloat r = w / 2;
		GLfloat h = (sqrt(2) - 1) * r;
		for (int w_ = -w / 2; w_ < w / 2; w_++)
		{
			GLfloat h_ = sqrt(2 * r * r - w_ * w_) - r;
			GLfloat b_ = (h_ / h) * b;
			point_cloud[l][int(SADDLE_WIDTH) / 2 + w_] = b_ + SADDLE_HEIGHT;
		}
	}
}

void draw_pad()
{
	cal_points();
	for (int l = 0; l < SADDLE_LENGTH; l++)
	{
		GLfloat w = (SADDLE_LENGTH - l) / SADDLE_LENGTH * SADDLE_WIDTH;
		for (int w_ = -w / 2; w_ < w / 2; w_++)
		{
			double x_1 = SADDLE_WIDTH / 2 + w_; double y_1 = l; double z_1 = point_cloud[int(y_1)][int(x_1)];
			double x_2 = SADDLE_WIDTH / 2 + w_ + 1; double y_2 = l; double z_2 = point_cloud[int(y_2)][int(x_2)];
			double x_3 = SADDLE_WIDTH / 2 + w_; double y_3 = l + 1; double z_3 = point_cloud[int(y_3)][int(x_3)];
			double x_4 = SADDLE_WIDTH / 2 + w_ + 1; double y_4 = l + 1; double z_4 = point_cloud[int(y_4)][int(x_4)];

			drawTriangle(x_1, y_1, z_1, 
						x_2, y_2, z_2, 
						x_4, y_4, z_4 ); //draw points 1,2,4 
			drawTriangle(x_1, y_1, z_1,
						x_3, y_3, z_3,
						x_4, y_4, z_4); //draw points 1,3,4 

		}
	};
}

void draw_spring()
{
	float r = 2;
	float R = 8;
	glPushMatrix();
	for (float theta = 0; theta < 360 * 5; theta++)
	{
		double rad = theta * PI / 180;
		glPushMatrix();
		glTranslated(R * cos(rad) + 20, R * sin(rad) + 40, -(1 + 2.5 * theta / 360) * r);
		drawSphere(r);
		glPopMatrix();

		glPushMatrix();		
		glTranslated(R * cos(rad) + 40, R * sin(rad) + 40, -(1 + 2.5 * theta / 360) * r);
		drawSphere(r);
		glPopMatrix();
		// drawSphere(15 * cos(rad)+15, 15*sin(rad)+15, -(1+2*theta/360)*r);
	}
	glPopMatrix();
}

// We are going to override (is that the right word?) the draw()
// method of ModelerView to draw out Unicycle
void Unicycle::draw()
{
	// This call takes care of a lot of the nasty projection 
	// matrix stuff.  Unless you want to fudge directly with the 
	// projection matrix, don't bother with this ...
	ModelerView::draw();

	// draw_surface();

	glEnable(GL_LIGHT2);
	static GLfloat lightDiffuse[] = { 5, 5, 5, 1 };			// strengthen the intencity of the light0 and light1    BW1
	static GLfloat lightAmbient[] = { 0.1, 0.1, 0.1, 1 };   // create an ambient light with low intencity
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiffuse);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, lightDiffuse);
	glLightfv(GL_LIGHT2, GL_AMBIENT, lightAmbient);

	

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

					draw_pad(); // BW5

					draw_spring();

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