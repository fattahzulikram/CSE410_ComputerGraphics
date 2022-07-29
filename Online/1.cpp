#ifdef _WIN32
#include <windows.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <GL/glut.h>

#define MAX_HAND 40.0
#define SPHERE_STACK 20
#define SPHERE_SLICES 24
#define CYLINDER_STACK 50
#define CYLINDER_SLICES 50

// Point Structure
struct point
{
	double x, y, z;

	point operator+(const point &AnotherPoint)
	{
		point RetPoint;
		RetPoint.x = x + AnotherPoint.x;
		RetPoint.y = y + AnotherPoint.y;
		RetPoint.z = z + AnotherPoint.z;
		return RetPoint;
	}

	point operator-(const point &AnotherPoint)
	{
		point RetPoint;
		RetPoint.x = x - AnotherPoint.x;
		RetPoint.y = y - AnotherPoint.y;
		RetPoint.z = z - AnotherPoint.z;
		return RetPoint;
	}

	point &operator=(const point &AnotherPoint)
	{
		x = AnotherPoint.x;
		y = AnotherPoint.y;
		z = AnotherPoint.z;
		return *this;
	}

	point operator*(double amount)
	{
		point RetPoint;
		RetPoint.x = x * amount;
		RetPoint.y = y * amount;
		RetPoint.z = z * amount;
		return RetPoint;
	}
};


// Global Variables
struct point CameraPosition; // Current Position of the Camera
struct point Up;             // Up direction vector, Z Axis here
struct point Right;          // Right direction vector, X Axis here
struct point Look;           // Look Direction Vector, Y Axis here

double MovementAmount = 2.0; // Camera movement multiplier
double RotationAmount = 2.0; // Rotation Angle
double ShapeChange = 2.0;    // Shape Change Amount

double CubeHand = MAX_HAND - 6.0; // Side of the Cube
double SphereRadius = 6.0;        // Radius of the sphere

bool bDrawAxis = true;

// Vector Operations
struct point CrossProduct(const struct point &PointA, const struct point &PointB)
{
	struct point RetPoint;
	// A = a1 * i + a2 * j + a3 * k and B = b1 * i + b2 * j + b3 * k
	// Result = (a2 * b3 – a3 * b2) * i + (a3 * b1 – a1 * b3) * j + (a1 * b2 – a2 * b1) * k
	RetPoint.x = PointA.y * PointB.z - PointA.z * PointB.y;
	RetPoint.y = PointA.z * PointB.x - PointA.x * PointB.z;
	RetPoint.z = PointA.x * PointB.y - PointA.y * PointB.x;
	return RetPoint;
}

// MISC
double Radian(double Angle)
{
	return Angle * M_PI / 180.0; 
}

/*
    DRAWING FUNCTIONS
*/

// Draw the axes
void drawAxes()
{
    if (bDrawAxis)
    {
        glColor3f(1.0, 1.0, 1.0);
        glBegin(GL_LINES);
        {
            glVertex3f(150, 0, 0);
            glVertex3f(-150, 0, 0);

            glVertex3f(0, -150, 0);
            glVertex3f(0, 150, 0);

            glVertex3f(0, 0, 150);
            glVertex3f(0, 0, -150);
        }
        glEnd();
    }
}

// Function for drawing a square
void drawSquare(double a)
{
    glBegin(GL_QUADS);
    {
        glVertex3f(a, a, 0);
        glVertex3f(a, -a, 0);
        glVertex3f(-a, -a, 0);
        glVertex3f(-a, a, 0);
    }
    glEnd();
}

// Function for drawing 1/8th of a sphere
void drawSphere(double radius, int slices, int stacks)
{
    struct point points[stacks + 1][slices + 1];
    int i, j;
    double h, r;
    // generate points
    for (i = 0; i <= stacks; i++)
    {
        h = radius * sin(((double)i / (double)stacks) * (M_PI / 2));
        r = radius * cos(((double)i / (double)stacks) * (M_PI / 2));
        for (j = 0; j <= slices; j++)
        {
            // Pi/2 instead of 2 * Pi in both x and y axes to discard 3 of the 4 parts
            points[i][j].x = r * cos(((double)j / (double)slices) * (M_PI / 2));
            points[i][j].y = r * sin(((double)j / (double)slices) * (M_PI / 2));
            points[i][j].z = h;
        }
    }
    // draw quads using generated points
    for (i = 0; i < stacks; i++)
    {
        for (j = 0; j < slices; j++)
        {
            glBegin(GL_QUADS);
            {
                // Only upper hemisphere to discard lower 4 parts, i.e., discard total 7 parts
                glVertex3f(points[i][j].x, points[i][j].y, points[i][j].z);
                glVertex3f(points[i][j + 1].x, points[i][j + 1].y, points[i][j + 1].z);
                glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, points[i + 1][j + 1].z);
                glVertex3f(points[i + 1][j].x, points[i + 1][j].y, points[i + 1][j].z);
            }
            glEnd();
        }
    }
}

// Function for drawing 1/8th of a cylinder, give 2*height as parameter for 1/4th
void drawCylinder(double height, double radius, int slices, int stacks)
{
    struct point points[stacks + 1][slices + 1];
    int i, j;
    double h, r;
    // generate points
    for (i = 0; i <= stacks; i++)
    {
        h = height * sin(((double)i / (double)stacks) * (M_PI / 2));
        r = radius; // Same radius for all the circles on the stack
        for (j = 0; j <= slices; j++)
        {
            // Pi/2 instead of 2 * Pi in both x and y axes to discard 3 of the 4 parts
            points[i][j].x = r * cos(((double)j / (double)slices) * (M_PI / 2));
            points[i][j].y = r * sin(((double)j / (double)slices) * (M_PI / 2));
            points[i][j].z = h;
        }
    }
    // draw quads using generated points
    for (i = 0; i < stacks; i++)
    {
        for (j = 0; j < slices; j++)
        {
            glBegin(GL_QUADS);
            {
                // Upper Hemisphere
                glVertex3f(points[i][j].x, points[i][j].y, points[i][j].z);
                glVertex3f(points[i][j + 1].x, points[i][j + 1].y, points[i][j + 1].z);
                glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, points[i + 1][j + 1].z);
                glVertex3f(points[i + 1][j].x, points[i + 1][j].y, points[i + 1][j].z);
                // lower hemisphere
                glVertex3f(points[i][j].x, points[i][j].y, -points[i][j].z);
                glVertex3f(points[i][j + 1].x, points[i][j + 1].y, -points[i][j + 1].z);
                glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, -points[i + 1][j + 1].z);
                glVertex3f(points[i + 1][j].x, points[i + 1][j].y, -points[i + 1][j].z);
            }
            glEnd();
        }
    }
}

// Function for drawing a cube
void drawCube()
{
    // All objects drawn after glRotated is called are rotated. Use glPushMatrix and glPopMatrix to save and restore the unrotated coordinate system
    // Color it white
    glColor3f(1.0, 1.0, 1.0);

    // A cube has six sides

    // Bottom Surface
    glPushMatrix();
    glTranslated(0, 0, -MAX_HAND);
    drawSquare(CubeHand);
    glPopMatrix();

    // Top Surface
    glPushMatrix();
    glTranslated(0, 0, MAX_HAND);
    drawSquare(CubeHand);
    glPopMatrix();

    // Side - 3
    glPushMatrix();
    glRotated(90, 1, 0, 0);
    glTranslated(0, 0, MAX_HAND);
    drawSquare(CubeHand);
    glPopMatrix();

    // Side - 4
    glPushMatrix();
    glRotated(-90, 1, 0, 0);
    glTranslated(0, 0, MAX_HAND);
    drawSquare(CubeHand);
    glPopMatrix();

    // Side - 5
    glPushMatrix();
    glRotated(90, 0, 1, 0);
    glTranslated(0, 0, MAX_HAND);
    drawSquare(CubeHand);
    glPopMatrix();

    // Side - 6
    glPushMatrix();
    glRotated(-90, 0, 1, 0);
    glTranslated(0, 0, MAX_HAND);
    drawSquare(CubeHand);
    glPopMatrix();
}

// Function for drawing eight 1/8th of a sphere on the eight corners of the cube
void drawSpheres()
{
    glColor3f(1.0, 0, 0);
    // Upper side, rotate about z to 4 corners before translation
    glPushMatrix();
    glTranslated(CubeHand, CubeHand, CubeHand);
    drawSphere(SphereRadius, SPHERE_SLICES, SPHERE_STACK);
    glPopMatrix();

    glPushMatrix();
    glRotated(90, 0, 0, 1);
    glTranslated(CubeHand, CubeHand, CubeHand);
    drawSphere(SphereRadius, SPHERE_SLICES, SPHERE_STACK);
    glPopMatrix();

    glPushMatrix();
    glRotated(180, 0, 0, 1);
    glTranslated(CubeHand, CubeHand, CubeHand);
    drawSphere(SphereRadius, SPHERE_SLICES, SPHERE_STACK);
    glPopMatrix();

    glPushMatrix();
    glRotated(270, 0, 0, 1);
    glTranslated(CubeHand, CubeHand, CubeHand);
    drawSphere(SphereRadius, SPHERE_SLICES, SPHERE_STACK);
    glPopMatrix();

    /*
        Rotate about (1, 1, 0) by 180 degree to send the sphere down part
        Use -CubeHand as Z to translate to lower half
        Then by same logic as above, distribute them
    */
    glPushMatrix();
    glTranslated(CubeHand, CubeHand, -CubeHand);
    glRotated(180, 1, 1, 0);
    drawSphere(SphereRadius, SPHERE_SLICES, SPHERE_STACK);
    glPopMatrix();

    glPushMatrix();
    glRotated(90, 0, 0, 1);
    glTranslated(CubeHand, CubeHand, -CubeHand);
    glRotated(180, 1, 1, 0);
    drawSphere(SphereRadius, SPHERE_SLICES, SPHERE_STACK);
    glPopMatrix();

    glPushMatrix();
    glRotated(180, 0, 0, 1);
    glTranslated(CubeHand, CubeHand, -CubeHand);
    glRotated(180, 1, 1, 0);
    drawSphere(SphereRadius, SPHERE_SLICES, SPHERE_STACK);
    glPopMatrix();

    glPushMatrix();
    glRotated(270, 0, 0, 1);
    glTranslated(CubeHand, CubeHand, -CubeHand);
    glRotated(180, 1, 1, 0);
    drawSphere(SphereRadius, SPHERE_SLICES, SPHERE_STACK);
    glPopMatrix();
}

// Function for drawing twelve 1/4th of a cylinder on the edges of the cube
void drawCylinders()
{
    glColor3f(0, 1.0, 0);

    // 4 sides, rotate about z to position them
    glPushMatrix();
    glTranslated(CubeHand, CubeHand, 0);
    drawCylinder(CubeHand, SphereRadius, CYLINDER_SLICES, CYLINDER_STACK);
    glPopMatrix();

    glPushMatrix();
    glRotated(90, 0, 0, 1);
    glTranslated(CubeHand, CubeHand, 0);
    drawCylinder(CubeHand, SphereRadius, CYLINDER_SLICES, CYLINDER_STACK);
    glPopMatrix();

    glPushMatrix();
    glRotated(180, 0, 0, 1);
    glTranslated(CubeHand, CubeHand, 0);
    drawCylinder(CubeHand, SphereRadius, CYLINDER_SLICES, CYLINDER_STACK);
    glPopMatrix();

    glPushMatrix();
    glRotated(270, 0, 0, 1);
    glTranslated(CubeHand, CubeHand, 0);
    drawCylinder(CubeHand, SphereRadius, CYLINDER_SLICES, CYLINDER_STACK);
    glPopMatrix();

    /*
        For rotating from vertical to horizontal, rotate 90 degree about (1, 1, 0)
        Then rotate 45 degree about z axis to match the shape
        Finally, rotate about x axis to put it on the upper or lower side
        Then rotate about z by 0, 90, 180 and 270 for four corners
    */

    // Top side, 45 degrees about x to go upper side
    glPushMatrix();
    glRotated(45, 1, 0, 0);
    glRotated(45, 0, 0, 1);
    glTranslated(CubeHand, CubeHand, 0);
    glRotated(90, 1, 1, 0);
    drawCylinder(CubeHand, SphereRadius, CYLINDER_SLICES, CYLINDER_STACK);
    glPopMatrix();

    glPushMatrix();
    glRotated(90, 0, 0, 1);
    glRotated(45, 1, 0, 0);
    glRotated(45, 0, 0, 1);
    glTranslated(CubeHand, CubeHand, 0);
    glRotated(90, 1, 1, 0);
    drawCylinder(CubeHand, SphereRadius, CYLINDER_SLICES, CYLINDER_STACK);
    glPopMatrix();

    glPushMatrix();
    glRotated(180, 0, 0, 1);
    glRotated(45, 1, 0, 0);
    glRotated(45, 0, 0, 1);
    glTranslated(CubeHand, CubeHand, 0);
    glRotated(90, 1, 1, 0);
    drawCylinder(CubeHand, SphereRadius, CYLINDER_SLICES, CYLINDER_STACK);
    glPopMatrix();

    glPushMatrix();
    glRotated(270, 0, 0, 1);
    glRotated(45, 1, 0, 0);
    glRotated(45, 0, 0, 1);
    glTranslated(CubeHand, CubeHand, 0);
    glRotated(90, 1, 1, 0);
    drawCylinder(CubeHand, SphereRadius, CYLINDER_SLICES, CYLINDER_STACK);
    glPopMatrix();

    // Bottom side, -45 degrees about x axis to go lower side
    glPushMatrix();
    glRotated(-45, 1, 0, 0);
    glRotated(45, 0, 0, 1);
    glTranslated(CubeHand, CubeHand, 0);
    glRotated(90, 1, 1, 0);
    drawCylinder(CubeHand, SphereRadius, CYLINDER_SLICES, CYLINDER_STACK);
    glPopMatrix();

    glPushMatrix();
    glRotated(90, 0, 0, 1);
    glRotated(-45, 1, 0, 0);
    glRotated(45, 0, 0, 1);
    glTranslated(CubeHand, CubeHand, 0);
    glRotated(90, 1, 1, 0);
    drawCylinder(CubeHand, SphereRadius, CYLINDER_SLICES, CYLINDER_STACK);
    glPopMatrix();

    glPushMatrix();
    glRotated(180, 0, 0, 1);
    glRotated(-45, 1, 0, 0);
    glRotated(45, 0, 0, 1);
    glTranslated(CubeHand, CubeHand, 0);
    glRotated(90, 1, 1, 0);
    drawCylinder(CubeHand, SphereRadius, CYLINDER_SLICES, CYLINDER_STACK);
    glPopMatrix();

    glPushMatrix();
    glRotated(270, 0, 0, 1);
    glRotated(-45, 1, 0, 0);
    glRotated(45, 0, 0, 1);
    glTranslated(CubeHand, CubeHand, 0);
    glRotated(90, 1, 1, 0);
    drawCylinder(CubeHand, SphereRadius, CYLINDER_SLICES, CYLINDER_STACK);
    glPopMatrix();
}

// Idle animation function
void animate()
{
    // codes for any changes in Models, Camera
    glutPostRedisplay();
}

// Draw on screen
void display()
{

    // clear the display
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0, 0, 0, 0); // color black
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    /********************
    / set-up camera here
    ********************/
    // load the correct matrix -- MODEL-VIEW matrix
    glMatrixMode(GL_MODELVIEW);

    // initialize the matrix
    glLoadIdentity();

    // now give three info
    // 1. where is the camera (viewer)?
    // 2. where is the camera looking?
    // 3. Which direction is the camera's UP direction?
    gluLookAt(CameraPosition.x,
              CameraPosition.y,
              CameraPosition.z,
              CameraPosition.x + Look.x,
              CameraPosition.y + Look.y,
              CameraPosition.z + Look.z,
              Up.x,
              Up.y,
              Up.z);

    // again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);

    /****************************
    / Add your objects from here
    ****************************/
    // add objects
    drawAxes();
    drawCube();
    drawSpheres();
    drawCylinders();

    glutSwapBuffers();
}

/*
    EVENT LISTENER FUNCTIONS
*/

// Handles regular key operations
void keyboardListener(unsigned char key, int x, int y)
{
    switch (key)
    {
    /*
        Rotation about Z Axis
        X = X cosΘ - Y sinΘ
        Y = X sinΘ + Y cosΘ
    */
    // 1 => rotate/look left
    case '1':
    {
        // LOOK LEFT -> Right and Look vectors will change
        Right = Right * cos(Radian(RotationAmount)) + Look * sin(Radian(RotationAmount));
        Look = Look * cos(Radian(RotationAmount)) - Right * sin(Radian(RotationAmount));
        break;
    }
    // 2 => rotate/look right
    case '2':
    {
        // LOOK RIGHT -> Right and Look vectors will change
        Right = Right * cos(Radian(RotationAmount)) - Look * sin(Radian(RotationAmount));
        Look = Right * sin(Radian(RotationAmount)) + Look * cos(Radian(RotationAmount));
        break;
    }
    /*
        Rotation about X axis=>
        Y = Y cosΘ - Z sinΘ
        Z = Y sinΘ + Z cosΘ
    */
    // 3 => rotate/look up
    case '3':
    {
        // LOOK UP -> Up and Look vectors will change
        Look = Look * cos(Radian(RotationAmount)) + Up * sin(Radian(RotationAmount));
        Up = Up * cos(Radian(RotationAmount)) - Look * sin(Radian(RotationAmount));
        break;
    }
    // 4 => rotate/look down
    case '4':
    {
        // LOOK DOWN -> Up and Look vectors will change
        Look = Look * cos(Radian(RotationAmount)) - Up * sin(Radian(RotationAmount));
        Up = Look * sin(Radian(RotationAmount)) + Up * cos(Radian(RotationAmount));
        break;
    }
    /*
        Rotation about Y Axis=>
        X = X cosΘ + Z sinΘ
        Z = Z cosΘ - X sinΘ
    */
    // 5 => tilt clockwise
    case '5':
    {
        // TILT CLOCKWISE -> Right and Up vectors will change
        Right = Right * cos(Radian(RotationAmount)) - Up * sin(Radian(RotationAmount));
        Up = Up * cos(Radian(RotationAmount)) + Right * sin(Radian(RotationAmount));
        break;
    }
    // 6 => tilt counterclockwise
    case '6':
    {
        // TILT COUNTERCLOCKWISE -> Right and Up vectors will change
        Right = Right * cos(Radian(RotationAmount)) + Up * sin(Radian(RotationAmount));
        Up = Up * cos(Radian(RotationAmount)) - Right * sin(Radian(RotationAmount));
        break;
    }
    // 7 => Look up left
    case '7':
    {
        Look = Look * cos(Radian(RotationAmount)) + Up * sin(Radian(RotationAmount));
        Up = Up * cos(Radian(RotationAmount)) - Look * sin(Radian(RotationAmount));
        Right = Right * cos(Radian(RotationAmount)) + Look * sin(Radian(RotationAmount));
        Look = Look * cos(Radian(RotationAmount)) - Right * sin(Radian(RotationAmount));
        break;
    }
    // 8 => Look up right
    case '8':
    {
        Look = Look * cos(Radian(RotationAmount)) + Up * sin(Radian(RotationAmount));
        Up = Up * cos(Radian(RotationAmount)) - Look * sin(Radian(RotationAmount));
        Right = Right * cos(Radian(RotationAmount)) - Look * sin(Radian(RotationAmount));
        Look = Right * sin(Radian(RotationAmount)) + Look * cos(Radian(RotationAmount));
        break;
    }
    // 9 => Look down left
    case '9':
    {
        Look = Look * cos(Radian(RotationAmount)) - Up * sin(Radian(RotationAmount));
        Up = Look * sin(Radian(RotationAmount)) + Up * cos(Radian(RotationAmount));
        Right = Right * cos(Radian(RotationAmount)) + Look * sin(Radian(RotationAmount));
        Look = Look * cos(Radian(RotationAmount)) - Right * sin(Radian(RotationAmount));
        break;
    }
    // 0 => Look down right
    case '0':
    {
        Look = Look * cos(Radian(RotationAmount)) - Up * sin(Radian(RotationAmount));
        Up = Look * sin(Radian(RotationAmount)) + Up * cos(Radian(RotationAmount));
        Right = Right * cos(Radian(RotationAmount)) - Look * sin(Radian(RotationAmount));
        Look = Right * sin(Radian(RotationAmount)) + Look * cos(Radian(RotationAmount));
        break;
    }
    default:
        break;
    }
}

// Handles special key operations
void specialKeyListener(int key, int x, int y)
{
    switch (key)
    {
    // Down Arrow Key => Move Camera Backwards
    case GLUT_KEY_DOWN:
    {
        // Go to the opposite direction of Look direction
        CameraPosition = CameraPosition - Look * MovementAmount;
        break;
    }
    // Down Arrow Key => Move Camera Forward
    case GLUT_KEY_UP:
    {
        // Go forward in the look direction
        CameraPosition = CameraPosition + Look * MovementAmount;
        break;
    }
    // Right Arrow Key => Move Camera Right
    case GLUT_KEY_RIGHT:
    {
        // Go in the right direction
        CameraPosition = CameraPosition + Right * MovementAmount;
        break;
    }
    // Left Arrow Key => Move Camera Left
    case GLUT_KEY_LEFT:
    {
        // Go in the opposite direction of right
        CameraPosition = CameraPosition - Right * MovementAmount;
        break;
    }
    // Page Up Key => Move Camera Up
    case GLUT_KEY_PAGE_UP:
    {
        // Go to up direction
        CameraPosition = CameraPosition + Up * MovementAmount;
        // printf("%f, %f, %f\n", CameraPosition.x, CameraPosition.y, CameraPosition.z);
        break;
    }
    // Page Down Key => Move Camera Down
    case GLUT_KEY_PAGE_DOWN:
    {
        // Do to the opposite of up direction
        CameraPosition = CameraPosition - Up * MovementAmount;
        // printf("%f, %f, %f\n", CameraPosition.x, CameraPosition.y, CameraPosition.z);
        break;
    }
    // Home Key => Cube to Sphere
    case GLUT_KEY_HOME:
    {
        /*
            Reduce cube side and increment Sphere radius
            Lowest cube side is 0
            Highest Sphere radius is MAX_HAND
        */
        if (CubeHand > 0)
        {
            CubeHand -= ShapeChange;
        }
        if (SphereRadius < MAX_HAND)
        {
            SphereRadius += ShapeChange;
        }
        break;
    }
    // End Key => Sphere to Cube
    case GLUT_KEY_END:
    {
        /*
            Reduce sphere radius and increment cube side
            Lowest sphere radius is 0
            Highest cube side is MAX_HAND
        */
        if (SphereRadius > 0)
        {
            SphereRadius -= ShapeChange;
        }
        if (CubeHand < MAX_HAND)
        {
            CubeHand += ShapeChange;
        }
        break;
    }
    default:
        break;
    }
}

// Mouse Click Listener
void mouseListener(int button, int state, int x, int y)
{ // x, y is the x-y of the screen (2D)
    switch (button)
    {
    case GLUT_LEFT_BUTTON:
    {
        if (state == GLUT_DOWN)
        {
            bDrawAxis = !bDrawAxis;
        }
        break;
    }
    default:
        break;
    }
}

/*
    INITIALIZATION FUNCTIONS
*/

// Function that initializes Camera and other directions
void InitializeVectors()
{
    // u = (0, 0, 1)
    Up.x = 0;
    Up.y = 0;
    Up.z = 1.0;

    // r = (-1/√2, 1/√2, 0)
    Right.x = -1 / sqrt(2.0);
    Right.y = 1 / sqrt(2.0);
    Right.z = 0;

    // l = (-1/√2, -1/√2, 0)
    Look.x = -1 / sqrt(2.0);
    Look.y = -1 / sqrt(2.0);
    Look.z = 0;

    // pos = (100, 100, 0)
    CameraPosition.x = 120.0;
    CameraPosition.y = 120.0;
    CameraPosition.z = 0.0;
}

// Initialization Function
void init()
{
    // Initialize Global Variables
    InitializeVectors();

    // clear the screen
    glClearColor(0, 0, 0, 0);

    // load the PROJECTION matrix
    glMatrixMode(GL_PROJECTION);

    // initialize the matrix
    glLoadIdentity();

    // give PERSPECTIVE parameters
    gluPerspective(80, 1, 1, 1000.0);
}

/*
    MAIN FUNCTION
*/

int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB); // Depth, Double buffer, RGB color

    glutCreateWindow("1705058_1_2");

    init();

    glEnable(GL_DEPTH_TEST); // enable Depth Testing

    glutDisplayFunc(display); // display callback function
    glutIdleFunc(animate);    // what you want to do in the idle time (when no drawing is occuring)

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutMouseFunc(mouseListener);

    glutMainLoop(); // The main loop of OpenGL

    return 0;
}
