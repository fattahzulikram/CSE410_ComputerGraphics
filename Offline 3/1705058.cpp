#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "1705058.hpp"
#include "bitmap_image.hpp"

/**********************Global Variables************************/
// Input file name
std::string inputFileName = "scene.txt";

// Window size
int iWindowWidth = 600;
int iWindowHeight = 600;

// Axis and Grid
bool bDrawAxis = true;
bool bDrawGrid = false;
double dAxisLength = 400.0;

// Camera Position
extern Vector vCamera;
Vector vLookAt;
Vector vUp;
Vector vRight;

// Movement
double dMovementAmount = 2.0;
double dRotateAmount = 2.0;

// Image Variables
int iPixelDimension;
int iCurrentImageCount = 1;

// Ray Tracing Variables
extern int iRecursionDepth;

// Object and Light Lists
extern std::vector<Object *> vObjects;
extern std::vector<PointLight *> vPointLights;
extern std::vector<SpotLight *> vSpotLights;

// Misc. Variables
double dFovY = 80.0;

/**********************Axes and Grid Functions************************/
void drawAxes()
{
    if (bDrawAxis)
    {
        glColor3f(1.0, 1.0, 1.0);
        glBegin(GL_LINES);
        {
            glVertex3f(dAxisLength, 0, 0);
            glVertex3f(-dAxisLength, 0, 0);

            glVertex3f(0, -dAxisLength, 0);
            glVertex3f(0, dAxisLength, 0);

            glVertex3f(0, 0, dAxisLength);
            glVertex3f(0, 0, -dAxisLength);
        }
        glEnd();
    }
}

void drawGrid()
{
    int i;
    if (bDrawGrid)
    {
        glColor3f(0.6, 0.6, 0.6); // grey
        glBegin(GL_LINES);
        {
            for (i = -8; i <= 8; i++)
            {

                if (i == 0)
                    continue; // SKIP the MAIN axes

                // lines parallel to Y-axis
                glVertex3f(i * 10, -90, 0);
                glVertex3f(i * 10, 90, 0);

                // lines parallel to X-axis
                glVertex3f(-90, i * 10, 0);
                glVertex3f(90, i * 10, 0);
            }
        }
        glEnd();
    }
}

/**********************Draw Functions************************/
// Draws the objects in the scene
void drawObjects()
{
    for (int i = 0; i < vObjects.size(); i++)
    {
        vObjects[i]->draw();
    }
}

// Draws the point lights in the scene
void drawPointLights()
{
    for (int i = 0; i < vPointLights.size(); i++)
    {
        glBegin(GL_POINTS);
        {
            glColor3f(vPointLights[i]->getLightColor().dRed, vPointLights[i]->getLightColor().dGreen, vPointLights[i]->getLightColor().dBlue);
            glVertex3f(vPointLights[i]->getLightPosition().dX, vPointLights[i]->getLightPosition().dY, vPointLights[i]->getLightPosition().dZ);
        }
        glEnd();
    }
}

// Draws the spot lights in the scene
void drawSpotLights()
{
    for (int i = 0; i < vSpotLights.size(); i++)
    {
        glBegin(GL_POINTS);
        {
            glColor3f(vSpotLights[i]->getLightColor().dRed, vSpotLights[i]->getLightColor().dGreen, vSpotLights[i]->getLightColor().dBlue);
            glVertex3f(vSpotLights[i]->getLightPosition().dX, vSpotLights[i]->getLightPosition().dY, vSpotLights[i]->getLightPosition().dZ);
        }
        glEnd();
    }
}

/**********************Capture Functions************************/
bitmap_image InitializeImage()
{
    bitmap_image image(iPixelDimension, iPixelDimension);

    for (int i = 0; i < iPixelDimension; i++)
    {
        for (int j = 0; j < iPixelDimension; j++)
        {
            image.set_pixel(i, j, 0, 0, 0);
        }
    }

    return image;
}

void Capture()
{
    printf("Capturing Image %d\n", iCurrentImageCount);
    bitmap_image image = InitializeImage();

    double dPlaneDistance = (iWindowHeight / 2.0) / tan(RADIAN(dFovY / 2.0));
    Vector vTopLeft = vCamera + (vLookAt * dPlaneDistance) - (vRight * (iWindowWidth / 2.0)) + (vUp * (iWindowHeight / 2.0));
    double dU = iWindowWidth / (double)iPixelDimension;
    double dV = iWindowHeight / (double)iPixelDimension;

    int iNearest;
    double dT, dTMin;

    // Choose middle of the grid cell
    vTopLeft = vTopLeft + (vRight * dU / 2.0) - (vUp * dV / 2.0);

    for (int iImageWidth = 0; iImageWidth < iPixelDimension; iImageWidth++)
    {
        for (int iImageHeight = 0; iImageHeight < iPixelDimension; iImageHeight++)
        {
            // Get Current Pixel
            Vector vCurrentPixel = vTopLeft + (vRight * (iImageWidth * dU)) - (vUp * (iImageHeight * dV));

            // Cast ray from eye to (curPixel-eye) direction
            Ray ray = Ray(vCamera, vCurrentPixel - vCamera);

            iNearest = -1;
            dTMin = INFINITE;

            for (int Counter = 0; Counter < vObjects.size(); Counter++)
            {
                Color color;
                dT = vObjects[Counter]->Intersect(&ray, color, 0);
                if(dT > 0 && dT < dTMin)
                {
                    dTMin = dT;
                    iNearest = Counter;
                }
            }

            if(iNearest != -1)
            {
                Color cColor;
                dTMin = vObjects[iNearest]->Intersect(&ray, cColor, 1);
                image.set_pixel(iImageWidth, iImageHeight, cColor.dRed * 255, cColor.dGreen * 255, cColor.dBlue * 255);
            }
        }
    }

    image.save_image("Output_1" + std::to_string(iCurrentImageCount) + ".bmp");
    printf("Image %d captured, Check Output_1%d.bmp\n", iCurrentImageCount, iCurrentImageCount);
}

/**********************Keyboard Listeners************************/
void keyboardListener(unsigned char key, int x, int y)
{
    switch (key)
    {
    // 1 => rotate/look left
    case '1':
    {
        vRight = vRight * cos(RADIAN(dMovementAmount)) + vLookAt * sin(RADIAN(dMovementAmount));
        vLookAt = vLookAt * cos(RADIAN(dMovementAmount)) - vRight * sin(RADIAN(dMovementAmount));
        break;
    }
    // 2 => rotate/look right
    case '2':
    {
        vRight = vRight * cos(RADIAN(dMovementAmount)) - vLookAt * sin(RADIAN(dMovementAmount));
        vLookAt = vLookAt * cos(RADIAN(-dMovementAmount)) + vRight * sin(RADIAN(dMovementAmount));
        break;
    }
    // 3 => rotate/look up
    case '3':
    {
        vLookAt = vLookAt * cos(RADIAN(dMovementAmount)) + vUp * sin(RADIAN(dMovementAmount));
        vUp = vUp * cos(RADIAN(dMovementAmount)) - vLookAt * sin(RADIAN(dMovementAmount));
        break;
    }
    // 4 => rotate/look down
    case '4':
    {
        vLookAt = vLookAt * cos(RADIAN(dMovementAmount)) - vUp * sin(RADIAN(dMovementAmount));
        vUp = vUp * cos(RADIAN(dMovementAmount)) + vLookAt * sin(RADIAN(dMovementAmount));
        break;
    }
    // 5 => tilt clockwise
    case '5':
    {
        vRight = vRight * cos(RADIAN(dMovementAmount)) - vUp * sin(RADIAN(dMovementAmount));
        vUp = vUp * cos(RADIAN(dMovementAmount)) + vRight * sin(RADIAN(dMovementAmount));
        break;
    }
    // 6 => tilt counter-clockwise
    case '6':
    {
        vRight = vRight * cos(RADIAN(dMovementAmount)) + vUp * sin(RADIAN(dMovementAmount));
        vUp = vUp * cos(RADIAN(dMovementAmount)) - vRight * sin(RADIAN(dMovementAmount));
        break;
    }
    // 0 => Capture image
    case '0':
    {
        Capture();
        iCurrentImageCount++;
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
    // Down Arrow Key => Move Camera Backwards
    case GLUT_KEY_DOWN:
    {
        vCamera = vCamera - vLookAt * dMovementAmount;
        break;
    }
    // Down Arrow Key => Move Camera Forward
    case GLUT_KEY_UP:
    {
        vCamera = vCamera + vLookAt * dMovementAmount;
        break;
    }
    // Right Arrow Key => Move Camera Right
    case GLUT_KEY_RIGHT:
    {
        vCamera = vCamera + vRight * dMovementAmount;
        break;
    }
    // Left Arrow Key => Move Camera Left
    case GLUT_KEY_LEFT:
    {
        vCamera = vCamera - vRight * dMovementAmount;
        break;
    }
    // Page Up Key => Move Camera Up
    case GLUT_KEY_PAGE_UP:
    {
        vCamera = vCamera + vUp * dMovementAmount;
        break;
    }
    // Page Down Key => Move Camera Down
    case GLUT_KEY_PAGE_DOWN:
    {
        vCamera = vCamera - vUp * dMovementAmount;
        break;
    }
    default:
        break;
    }
}

/**********************Mouse Listener************************/
void mouseListener(int button, int state, int x, int y)
{
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

//**********************Draw On Screen************************/
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
    gluLookAt(vCamera.dX, vCamera.dY, vCamera.dZ, (vCamera.dX + vLookAt.dX), (vCamera.dY + vLookAt.dY), (vCamera.dZ + vLookAt.dZ), vUp.dX, vUp.dY, vUp.dZ);

    // again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);

    /****************************
    / Add your objects from here
    ****************************/
    // add axes
    drawAxes();

    // add objects
    drawObjects();

    // add lights
    if(vPointLights.size() > 0)
    {
        drawPointLights();
    }
    if(vSpotLights.size() > 0)
    {
        drawSpotLights();
    }

    // ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
    glutSwapBuffers();
}

/**********************Idle Animate************************/
void animate()
{
    // codes for any changes in Models, Camera
    glutPostRedisplay();
}

/**********************Misc Functions************************/
// Initialize Camera Vector
void initializeVectors()
{
    vCamera = Vector(100, 100, 0);
    vUp = Vector(0, 0, 1);
    vLookAt = Vector(-1 / sqrt(2.0), -1 / sqrt(2.0), 0);
    vRight = Vector(-1 / sqrt(2.0), 1 / sqrt(2.0), 0);
}

// Read Environment Data
void loadData()
{
    std::ifstream inFile;
    inFile.open(inputFileName);
    if (!inFile)
    {
        printf("Unable to open file scene.txt\n");
        exit(-1); // terminate with error
    }

    int iObjectCount, iPointLightCount, iSpotLightCount;

    // Recursion Level and Pixel Dimension
    inFile >> iRecursionDepth >> iPixelDimension;

    // Object Count
    inFile >> iObjectCount;

    // Object Data
    for (int Counter = 0; Counter < iObjectCount; Counter++)
    {
        std::string sObjectType;
        inFile >> sObjectType;
        Object *pObject = NULL;

        if (sObjectType == "sphere")
        {
            pObject = CreateSphere(inFile);
        }
        else if (sObjectType == "triangle")
        {
            pObject = CreateTriangle(inFile);
        }
        else if (sObjectType == "general")
        {
            pObject = CreateQuadric(inFile);
        }
        else
        {
            printf("Unknown Object Type %s\n", sObjectType.c_str());
            exit(-1); // terminate with error
        }
        vObjects.push_back(pObject);
        pObject = NULL;
    }

    // Point Light Count
    inFile >> iPointLightCount;

    // Point Light Data
    for (int Counter = 0; Counter < iPointLightCount; Counter++)
    {
        PointLight *pPointLight = CreatePointLight(inFile);
        vPointLights.push_back(pPointLight);
        pPointLight = NULL;
    }
    
    // Spot Light Count
    inFile >> iSpotLightCount;

    // Spot Light Data
    for (int Counter = 0; Counter < iSpotLightCount; Counter++)
    {
        SpotLight *pSpotLight = CreateSpotLight(inFile);
        vSpotLights.push_back(pSpotLight);
        pSpotLight = NULL;
    }

    // Add Floor
    vObjects.push_back(new Floor(1000, 20, Color(0, 0, 0), Coefficients(0.5, 0.5, 0.5, 0.5), 15));

    // Close File
    inFile.close();
}

/**********************Initialize OpenGL************************/
void init()
{
    // Initialize Camera
    initializeVectors();

    // clear the screen
    glClearColor(0, 0, 0, 0);

    // load the PROJECTION matrix
    glMatrixMode(GL_PROJECTION);

    // initialize the matrix
    glLoadIdentity();

    // give PERSPECTIVE parameters
    gluPerspective(dFovY, 1, 1, 1000.0);
}

/**********************Main Function************************/
int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitWindowSize(iWindowWidth, iWindowHeight);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB); // Depth, Double buffer, RGB color

    glutCreateWindow("Ray Tracing");

    init();

    glEnable(GL_DEPTH_TEST); // enable Depth Testing

    glutDisplayFunc(display); // display callback function
    glutIdleFunc(animate);    // what you want to do in the idle time (when no drawing is occuring)

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutMouseFunc(mouseListener);

    // Read Environment Data
    loadData();

    glutMainLoop(); // The main loop of OpenGL

    return 0;
}
