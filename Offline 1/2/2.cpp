#ifdef _WIN32
#include <windows.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <GL/glut.h>

#define WHEEL_RADIUS 30
#define X_WIDTH_HALF 6.0
#define WHEEL_PROGRESS_ANGLE 3.0
#define WHEEL_PROGRESS_LENGTH (WHEEL_RADIUS * WHEEL_PROGRESS_ANGLE * M_PI) / 180.0; // Arc length = Theta * r
#define WHEEL_ANGLE_CHANGE 3.0

#define WHEEL_SLICES 50

bool bDrawAxis = false;
bool bDrawGrid = true;

double cameraHeight = 100.0;
double cameraAngle = 180;


double Radian(double Angle)
{
    return Angle * M_PI / 180.0;
}

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
} WheelCenter;

// The wheel structure
struct WheelObject
{
private:
    double WheelLookAngle = 0.0;
    double XAngle = 0.0;

    // Function for drawing the wheel body. A modified version of the cylinder code
    void drawWheelBody()
    {
        struct point points[WHEEL_SLICES + 1];

        for (int i = 0; i <= WHEEL_SLICES; i++)
        {
            points[i].x = WHEEL_RADIUS * cos(((double)i / (double)WHEEL_SLICES) * 2 * M_PI);
            points[i].z = WHEEL_RADIUS * sin(((double)i / (double)WHEEL_SLICES) * 2 * M_PI);
        }
        double Color = 0.3;
        for (int i = 0; i < WHEEL_SLICES; i++)
        {
            glBegin(GL_QUADS);
            if(i < WHEEL_SLICES / 2){
                Color += (0.6 / (double) (WHEEL_SLICES / 2.0));
            }
            else{
                Color -= (0.6 / (double) (WHEEL_SLICES / 2.0));
            }
            glColor3f(Color, Color, Color);
            glVertex3f(points[i].x, X_WIDTH_HALF + 0.5, points[i].z + WHEEL_RADIUS);
            glVertex3f(points[i].x, -X_WIDTH_HALF - 0.5, points[i].z + WHEEL_RADIUS);
            glVertex3f(points[i + 1].x, -X_WIDTH_HALF - 0.5, points[i + 1].z + WHEEL_RADIUS);
            glVertex3f(points[i + 1].x, X_WIDTH_HALF + 0.5, points[i + 1].z + WHEEL_RADIUS);
            glEnd();
        }
    }

    // Function for drawing the X shape on the wheel
    void drawX(double height, double width)
    {
        glColor3f(0.7, 0.7, 0.7);

        // A rectangle parallel to Z axis
        glBegin(GL_QUADS);
        glVertex3f(0, width, 0);
        glVertex3f(0, -width, 0);
        glVertex3f(0, -width, height);
        glVertex3f(0, width, height);
        glEnd();

        // A rectangle parallel to X axis
        glBegin(GL_QUADS);
        glVertex3f(-height / 2.0, -width, height / 2.0);
        glVertex3f(-height / 2.0, width, height / 2.0);
        glVertex3f(height / 2.0, width, height / 2.0);
        glVertex3f(height / 2.0, -width, height / 2.0);
        glEnd();
    }

public:
    // Function for drawing Wheel+X
    void drawWheelWithX()
    {
        glPushMatrix();

        // Do transformations
        // Translate to next position
        glTranslated(WheelCenter.x, WheelCenter.y, WheelCenter.z);
        // Rotate Wheel left and right
        glRotated(WheelLookAngle, 0, 0, 1);
        // Rotate the Rim
        // https://stackoverflow.com/a/23500304
        glTranslated(0, 0, WHEEL_RADIUS); // Object Translaation
        glRotated(XAngle, 0, 1, 0); // Rotate Around Point
        glTranslated(0, 0, -WHEEL_RADIUS); // Point Translation
        // Draw the wheel and X
        drawX(WHEEL_RADIUS * 2, X_WIDTH_HALF);
        drawWheelBody();

        glPopMatrix();
    }

    // Function for moving wheel forward
    void moveForward()
    {
        double MoveLen = WHEEL_PROGRESS_LENGTH;

        WheelCenter.x += MoveLen * cos(Radian(WheelLookAngle));
        WheelCenter.y += MoveLen * sin(Radian(WheelLookAngle));
        XAngle += WHEEL_PROGRESS_ANGLE;
    }

    // Function for moving wheel backward
    void moveBackward()
    {
        double MoveLen = WHEEL_PROGRESS_LENGTH;

        WheelCenter.x -= MoveLen * cos(Radian(WheelLookAngle));
        WheelCenter.y -= MoveLen * sin(Radian(WheelLookAngle));
        XAngle -= WHEEL_PROGRESS_ANGLE;
    }

    // Function for Rotating wheel left
    void rotateLeft()
    {
        WheelLookAngle -= WHEEL_ANGLE_CHANGE;
    }

    // Function for rotating wheel right
    void rotateRight()
    {
        WheelLookAngle += WHEEL_ANGLE_CHANGE;
    }
}Wheel;

void drawAxes()
{
    if (bDrawAxis)
    {
        glColor3f(1.0, 1.0, 1.0);
        glBegin(GL_LINES);
        {
            glVertex3f(100, 0, 0);
            glVertex3f(-100, 0, 0);

            glVertex3f(0, -100, 0);
            glVertex3f(0, 100, 0);

            glVertex3f(0, 0, 100);
            glVertex3f(0, 0, -100);
        }
        glEnd();
    }
}

// Function for drawing grid
void drawGrid()
{
    int i;
    if (bDrawGrid)
    {
        glColor3f(0.6, 0.6, 0.6); // grey
        glBegin(GL_LINES);
        {
            for (i = -12; i <= 12; i++)
            {
                // lines parallel to Y-axis
                glVertex3f(i * 10, -130, 0);
                glVertex3f(i * 10, 130, 0);

                // lines parallel to X-axis
                glVertex3f(-130, i * 10, 0);
                glVertex3f(130, i * 10, 0);
            }
        }
        glEnd();
    }
}

void keyboardListener(unsigned char key, int x, int y)
{
    switch (key)
    {

    case 'w':
    {
        Wheel.moveForward();
        break;
    }
    case 's':
    {
        Wheel.moveBackward();
        break;
    }
    case 'a':
    {
        Wheel.rotateLeft();
        break;
    }
    case 'd':
    {
        Wheel.rotateRight();
        break;
    }
    default:
        break;
    }
}

void specialKeyListener(int key, int x, int y)
{
    switch (key)
    {
    case GLUT_KEY_DOWN:
    {
        cameraHeight -= 3.0;
        break;
    }
    case GLUT_KEY_UP:
    {
        cameraHeight += 3.0;
        break;
    }
    case GLUT_KEY_RIGHT:
    {
        cameraAngle += 0.03;
        break;
    }
    case GLUT_KEY_LEFT:
    {
        cameraAngle -= 0.03;
        break;
    }
    default:
        break;
    }
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
    gluLookAt(200 * cos(cameraAngle), 200 * sin(cameraAngle), cameraHeight, 0, 0, 0, 0, 0, 1);

    // again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);

    /****************************
    / Add your objects from here
    ****************************/
    // add objects
    drawAxes();
    drawGrid();
    Wheel.drawWheelWithX();

    // ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
    glutSwapBuffers();
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
    case GLUT_RIGHT_BUTTON:
    {
        if (state == GLUT_DOWN)
        {
            bDrawGrid = !bDrawGrid;
        }
        break;
    }
    default:
        break;
    }
}

// Idle animation function
void animate()
{
    // codes for any changes in Models, Camera
    glutPostRedisplay();
}

void init()
{
    /*
        YOUR CODE HERE
    */

    // clear the screen
    glClearColor(0, 0, 0, 0);

    // load the PROJECTION matrix
    glMatrixMode(GL_PROJECTION);

    // initialize the matrix
    glLoadIdentity();

    // give PERSPECTIVE parameters
    gluPerspective(80, 1, 1, 1000.0);
}

int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB); // Depth, Double buffer, RGB color

    glutCreateWindow("1705058_3");

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
