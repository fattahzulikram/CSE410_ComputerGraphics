#include <iostream>
#include <fstream>
#include <cmath>
#include <stack>
#include <vector>
#include <iomanip>

#include "bitmap_image.hpp"

using namespace std;
#define PI (2 * acos(0.0))
#define RADIAN(t) (t * PI / 180.0)
#define PRECISION 7
#define INFINITE 9999999

/********************Point Class****************************/
class Point
{
public:
    double x, y, z;
    Point(double x = 0.0, double y = 0.0, double z = 0.0)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }
    // IO stream overloads

    friend ostream &operator<<(ostream &os, const Point &p)
    {
        os << setprecision(PRECISION) << fixed << p.x << " " << p.y << " " << p.z << endl;
        return os;
    }

    friend istream &operator>>(istream &is, Point &p)
    {
        is >> p.x >> p.y >> p.z;
        return is;
    }

    // Overloading Operators

    Point operator-() const
    {
        return Point(-x, -y, -z);
    }

    Point operator+(const Point &p) const
    {
        return Point(x + p.x, y + p.y, z + p.z);
    }

    Point operator-(const Point &p) const
    {
        return Point(x - p.x, y - p.y, z - p.z);
    }

    Point operator*(const Point &p) const
    {
        return Point(x * p.x, y * p.y, z * p.z);
    }

    Point operator*(double d) const
    {
        return Point(x * d, y * d, z * d);
    }

    Point operator/(double d) const
    {
        return Point(x / d, y / d, z / d);
    }

    Point operator=(const Point &p)
    {
        x = p.x;
        y = p.y;
        z = p.z;
        return *this;
    }

    // Normalize the point

    double length() const
    {
        return sqrt(x * x + y * y + z * z);
    }

    void normalize()
    {
        double len = length();
        if (len != 0)
        {
            x /= len;
            y /= len;
            z /= len;
        }
    }
};

/********************Vector operations****************************/
// Dot product
double dot(const Point &p1, const Point &p2)
{
    return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
}

// Cross product
Point cross(const Point &p1, const Point &p2)
{
    return Point(p1.y * p2.z - p1.z * p2.y,
                 p1.z * p2.x - p1.x * p2.z,
                 p1.x * p2.y - p1.y * p2.x);
}

// Rodrigues' Formula
Point RodriguesFormula(Point x, Point a, double angle)
{
    return x * cos(RADIAN(angle)) + cross(a, x) * sin(RADIAN(angle)) + a * dot(a, x) * (1 - cos(RADIAN(angle)));
}

/********************Matrix Struct****************************/
struct Matrix
{
    int row, col;
    vector<vector<double>> matrix;

    Matrix(int row = 4, int col = 4)
    {
        this->row = row;
        this->col = col;
        matrix.resize(row);
        for (int i = 0; i < row; i++)
        {
            matrix[i].resize(col);
        }
    }

    // Functions for setting row and column

    void setRow(int row, const Point &p)
    {
        matrix[row][0] = p.x;
        matrix[row][1] = p.y;
        matrix[row][2] = p.z;
        matrix[row][3] = 0;
    }

    void setColumn(int col, const Point &p)
    {
        matrix[0][col] = p.x;
        matrix[1][col] = p.y;
        matrix[2][col] = p.z;
        matrix[3][col] = 0;
    }

    // IO stream overloads
    friend ostream &operator<<(ostream &os, const Matrix &m)
    {
        // Print Coordinate Matrix
        if (m.col == 1)
        {
            // Print only first 3 rows of matrix
            for (int i = 0; i < m.row - 1; i++)
            {
                os << setprecision(PRECISION) << fixed << m.matrix[i][0] << " ";
            }
            os << endl;
        }
        else
        {
            // Print only 3*3 matrix
            for (int i = 0; i < m.row - 1; i++)
            {
                for (int j = 0; j < m.col - 1; j++)
                {
                    os << setprecision(PRECISION) << fixed << m.matrix[i][j] << " ";
                }
                os << endl;
            }
        }
        return os;
    }
};

/********************Matrix Utilities****************************/
// Matrix Multiplication
Matrix matrixMultiplication(const Matrix &m1, const Matrix &m2)
{
    Matrix m(m1.row, m2.col);
    for (int i = 0; i < m1.row; i++)
    {
        for (int j = 0; j < m2.col; j++)
        {
            for (int k = 0; k < m1.col; k++)
            {
                m.matrix[i][j] += m1.matrix[i][k] * m2.matrix[k][j];
            }
        }
    }
    return m;
}

// Homogenous Coordinate
Matrix CoordinateMatrix(const Point &point)
{
    Matrix m(4, 1);
    m.matrix[0][0] = point.x;
    m.matrix[1][0] = point.y;
    m.matrix[2][0] = point.z;
    m.matrix[3][0] = 1;
    return m;
}

// Identity Matrix
Matrix IdentityMatrix()
{
    Matrix m(4, 4);
    for (int i = 0; i < 4; i++)
    {
        m.matrix[i][i] = 1;
    }
    return m;
}

// Translation Matrix
Matrix TranslationMatrix(const Point &point)
{
    Matrix m = IdentityMatrix();
    m.matrix[0][3] = point.x;
    m.matrix[1][3] = point.y;
    m.matrix[2][3] = point.z;
    m.matrix[3][3] = 1;
    return m;
}

// Scaling Matrix
Matrix ScalingMatrix(const Point &point)
{
    Matrix m = IdentityMatrix();
    m.matrix[0][0] = point.x;
    m.matrix[1][1] = point.y;
    m.matrix[2][2] = point.z;
    m.matrix[3][3] = 1;
    return m;
}

// Rotation Matrix
Matrix RotationMatrix(double angle, const Point &axis)
{
    Point i(1, 0, 0);
    Point j(0, 1, 0);
    Point k(0, 0, 1);

    Point c1 = RodriguesFormula(i, axis, angle);
    Point c2 = RodriguesFormula(j, axis, angle);
    Point c3 = RodriguesFormula(k, axis, angle);

    Matrix RotateMatrix(4, 4);
    RotateMatrix.setColumn(0, c1);
    RotateMatrix.setColumn(1, c2);
    RotateMatrix.setColumn(2, c3);
    RotateMatrix.matrix[3][3] = 1;

    return RotateMatrix;
}

Point ScaleW(const Matrix &m)
{
    double w = m.matrix[m.row - 1][m.col - 1];
    Point p(m.matrix[0][0] / w, m.matrix[1][0] / w, m.matrix[2][0] / w);
    return p;
}

/********************Triangle and Color Structure****************************/
struct Color
{
    double r, g, b;

    Color(double r = 0, double g = 0, double b = 0)
    {
        this->r = r;
        this->g = g;
        this->b = b;
    }

    friend ostream &operator<<(ostream &os, const Color &r)
    {
        os << setprecision(PRECISION) << fixed << r.r << " " << r.g << " " << r.b << endl;
        return os;
    }
};

struct Triangle
{
    Point points[3];
    Color color;
};

Color GetRandomColor()
{
    Color c;
    c.r = rand() % 255;
    c.g = rand() % 255;
    c.b = rand() % 255;
    return c;
}

/********************Global Variables****************************/
Point eye, look, up;                 // gluLookAt parameters
double fov, aspect, zNear, zFar;     // gluPerspective parameters
Matrix viewMatrix, projectionMatrix; // View and Projection Matrices

stack<Matrix> matrixStack;                // Matrix Stack
ifstream inFile, configFile;              // Input file
ofstream stage1, stage2, stage3, zBuffer; // Output files

int Screen_Width, Screen_Height; // Screen size
double xLeft, yBottom;           // Screen coordinates
double zFront, zRear;            // Depth of the screen

vector<Triangle> triangles; // Vector of triangles for stage 4
double dx, dy;              // Scan line distance
double Top_Y, Left_X;       // Top and left pixels of the screen
double Bottom_Y, Right_X;   // Bottom and right pixels of the screen
double **zBufferMatrix;     // zBuffer matrix
Color **colorBufferMatrix;  // colorBuffer matrix

bitmap_image image; // Image

/********************Scan Conversion and Clipping****************************/
// Stage 4 procedure for triangles
void ScanConversion(Triangle triangle)
{
    // Get the maximum and minimum y values
    double maxY = max(triangle.points[0].y, max(triangle.points[1].y, triangle.points[2].y));
    double minY = min(triangle.points[0].y, min(triangle.points[1].y, triangle.points[2].y));
    double maxX = max(triangle.points[0].x, max(triangle.points[1].x, triangle.points[2].x));
    double minX = min(triangle.points[0].x, min(triangle.points[1].x, triangle.points[2].x));

    // Triangle is out of the screen
    if (minY > Top_Y || maxY < Bottom_Y || minX > Right_X || maxX < Left_X)
    {
        return;
    }

    // Triangle is within the screen, clip it if necessary and get the top and bottom scanlines
    int topScanline, bottomScanline;

    if (maxY > Top_Y)
    {
        topScanline = 0;
    }
    else
    {
        topScanline = (int)round((Top_Y - maxY) / dy);
    }
    if (minY < Bottom_Y)
    {
        bottomScanline = Screen_Height - 1;
    }
    else
    {
        bottomScanline = Screen_Height - 1 - (int)round((minY - Bottom_Y) / dy);
    }

    // Iterate through the scanlines
    for (int row = topScanline; row <= bottomScanline; row++)
    {
        // Scanline Y coordinate
        double y = Top_Y - row * dy;

        // Determining intersections unless parallel to X
        Point Intersections[3];
        /*
            Convention:
            Intersections[i][0] = Intersection point
            Intersections[i][1] = One point of the line
            Intersections[i][2] = Another point of the line
        */
        Intersections[0] = Point(INFINITE, 0, 1);
        Intersections[1] = Point(INFINITE, 1, 2);
        Intersections[2] = Point(INFINITE, 2, 0);

        for (int j = 0; j < 3; j++)
        {
            // Three intersection cases - AB, BC, CA
            // j = 0: AB, j = 1: BC, j = 2: CA
            Point PointA, PointB;
            PointA = triangle.points[(int)Intersections[j].y];
            PointB = triangle.points[(int)Intersections[j].z];

            // Unless line is parallel to X, calculate the intersection
            if (PointA.y != PointB.y)
            {
                Intersections[j].x = (y - PointA.y) * (PointB.x - PointA.x) / (PointB.y - PointA.y) + PointA.x;
            }

            // Check if valid or not
            if (Intersections[j].x != INFINITE)
            {
                bool CheckX = (Intersections[j].x < min(PointA.x, PointB.x)) || (Intersections[j].x > max(PointA.x, PointB.x));
                bool CheckY = (y < min(PointA.y, PointB.y)) || (y > max(PointA.y, PointB.y));

                if (CheckX || CheckY)
                {
                    Intersections[j].x = INFINITE;
                }
            }
        }

        // Determine minimum and maximum intersection points
        int minimumIndex = -1;
        int maximumIndex = -1;

        for (int j = 0; j < 3; j++)
        {
            if (Intersections[j].x != INFINITE)
            {
                if (minimumIndex == -1 && maximumIndex == -1)
                {
                    minimumIndex = j;
                    maximumIndex = j;
                }
                else
                {
                    if (Intersections[j].x < Intersections[minimumIndex].x)
                    {
                        minimumIndex = j;
                    }
                    if (Intersections[j].x > Intersections[maximumIndex].x)
                    {
                        maximumIndex = j;
                    }
                }
            }
        }

        // Find Left and Right Intersecting columns after necessary clipping
        int leftColumn, rightColumn;
        if (Intersections[minimumIndex].x <= Left_X)
        {
            leftColumn = 0;
        }
        else
        {
            leftColumn = (int)round((Intersections[minimumIndex].x - Left_X) / dx);
        }
        if (Intersections[maximumIndex].x >= Right_X)
        {
            rightColumn = Screen_Width - 1;
        }
        else
        {
            rightColumn = Screen_Width - 1 - (int)round((Right_X - Intersections[maximumIndex].x) / dx);
        }

        // Determine za and zb
        Point PointA = triangle.points[(int)Intersections[minimumIndex].y];
        Point PointB = triangle.points[(int)Intersections[minimumIndex].z];

        double za = (y - PointA.y) * (PointB.z - PointA.z) / (PointB.y - PointA.y) + PointA.z;

        PointA = triangle.points[(int)Intersections[maximumIndex].y];
        PointB = triangle.points[(int)Intersections[maximumIndex].z];

        double zb = (y - PointA.y) * (PointB.z - PointA.z) / (PointB.y - PointA.y) + PointA.z;

        // Calculate 0th instance of zp
        double zp = za + ((Left_X + leftColumn * dx) - Intersections[minimumIndex].x) * (zb - za) / (Intersections[maximumIndex].x - Intersections[minimumIndex].x);
        double increaseConstant = dx * (zb - za) / (Intersections[maximumIndex].x - Intersections[minimumIndex].x);

        // Iterate from left column to right column
        for (int column = leftColumn; column <= rightColumn; column++)
        {
            // Any instance after 0th instace of zp should add a known constant
            if (column != leftColumn)
            {
                zp = zp + increaseConstant;
            }

            // Update zBuffer if necessary
            if (zp > zFront && zp < zBufferMatrix[row][column])
            {
                zBufferMatrix[row][column] = zp;

                // Update frame buffer
                colorBufferMatrix[row][column] = triangle.color;
            }
        }
    }
}

// Read stage 3 file to get triangles
void ReadTriangles()
{
    ifstream stage3File("stage3.txt");
    // Read the stage 3 file and store the triangles in a vector
    while (true)
    {
        Triangle t;
        for (int i = 0; i < 3; i++)
        {
            stage3File >> t.points[i];
        }
        t.color = GetRandomColor();

        // stage3.txt has endl at the end, so, a triangle with all 0,0,0 might be read. Ignore it
        bool point1IsZero = t.points[0].x == 0 && t.points[0].y == 0 && t.points[0].z == 0;
        bool point2IsZero = t.points[1].x == 0 && t.points[1].y == 0 && t.points[1].z == 0;
        bool point3IsZero = t.points[2].x == 0 && t.points[2].y == 0 && t.points[2].z == 0;
        if (point1IsZero && point2IsZero && point3IsZero)
        {
            break;
        }
        triangles.push_back(t);

        if (stage3File.eof())
        {
            break;
        }
    }
    stage3File.close();
}

// Initialize Necessary Variables
bool init()
{
    inFile.open("scene.txt");
    configFile.open("config.txt");

    if (!inFile || !configFile)
    {
        cout << "File could not be opened" << endl;
        return false;
    }

    // Output files
    stage1.open("stage1.txt");
    stage2.open("stage2.txt");
    stage3.open("stage3.txt");
    zBuffer.open("z_buffer.txt");

    // gluLookAt and gluPerspective parameters
    inFile >> eye >> look >> up >> fov >> aspect >> zNear >> zFar;

    // Initialize view matrix
    Point l = look - eye;
    l.normalize();

    Point r = cross(l, up);
    r.normalize();

    Point u = cross(r, l);

    // Move the eye/camera to origin
    Matrix translation = TranslationMatrix(-eye);

    // Apply the following rotation
    Matrix view = Matrix(4, 4);
    view.setRow(0, r);
    view.setRow(1, u);
    view.setRow(2, -l);
    view.matrix[3][3] = 1;

    // Multiply the translation matrix with the view matrix
    viewMatrix = matrixMultiplication(view, translation);

    // Initialize projection matrix
    double fovX = fov * aspect;
    double tanfov = zNear * tan(RADIAN(fov / 2));
    double tanfovX = zNear * tan(RADIAN(fovX / 2));

    // Apply the following projection
    projectionMatrix = Matrix(4, 4);
    projectionMatrix.matrix[0][0] = zNear / tanfovX;
    projectionMatrix.matrix[1][1] = zNear / tanfov;
    projectionMatrix.matrix[2][2] = -(zFar + zNear) / (zFar - zNear);
    projectionMatrix.matrix[2][3] = -(2 * zFar * zNear) / (zFar - zNear);
    projectionMatrix.matrix[3][2] = -1;

    // Read Config File
    configFile >> Screen_Width >> Screen_Height >> xLeft >> yBottom >> zFront >> zRear;

    return true;
}

// Initialize necessary variables for stage 4, returns the output image file
void stage4init()
{
    double xRight = (-1) * xLeft;
    double yTop = (-1) * yBottom;
    dx = (xRight - xLeft) / Screen_Width;
    dy = (yTop - yBottom) / Screen_Height;

    Top_Y = yTop - dy / 2.0;
    Left_X = xLeft + dx / 2.0;
    Bottom_Y = yBottom + dy / 2.0;
    Right_X = -xLeft - dx / 2.0;

    // Initialize zBuffer matrix
    zBufferMatrix = new double *[Screen_Height];
    for (int i = 0; i < Screen_Height; i++)
    {
        zBufferMatrix[i] = new double[Screen_Width];
        for (int j = 0; j < Screen_Width; j++)
        {
            zBufferMatrix[i][j] = zRear;
        }
    }

    // Initialize framebuffer matrix
    colorBufferMatrix = new Color *[Screen_Height];
    for (int i = 0; i < Screen_Height; i++)
    {
        colorBufferMatrix[i] = new Color[Screen_Width];
        for (int j = 0; j < Screen_Width; j++)
        {
            colorBufferMatrix[i][j] = Color(0, 0, 0);
        }
    }
}

int main()
{
    if (!init())
    {
        cout << "Initialization failed" << endl;
        return -1;
    }
    srand(time(0));
    matrixStack.push(IdentityMatrix()); // Initialize the matrix stack
    while (true)
    {
        string command;
        inFile >> command;

        // End of file
        if (command == "end")
        {
            break;
        }
        // Push the current matrix to the stack
        else if (command == "push")
        {
            matrixStack.push(matrixStack.top());
        }
        // Pop the top matrix from the stack
        else if (command == "pop")
        {
            matrixStack.pop();
        }
        // Translate the current matrix
        else if (command == "translate")
        {
            double tx, ty, tz;
            inFile >> tx >> ty >> tz;
            matrixStack.top() = matrixMultiplication(matrixStack.top(), TranslationMatrix(Point(tx, ty, tz)));
        }
        // Scale the current matrix
        else if (command == "scale")
        {
            double sx, sy, sz;
            inFile >> sx >> sy >> sz;
            matrixStack.top() = matrixMultiplication(matrixStack.top(), ScalingMatrix(Point(sx, sy, sz)));
        }
        // Rotate the current matrix
        else if (command == "rotate")
        {
            double angle, rx, ry, rz;
            inFile >> angle >> rx >> ry >> rz;

            Point axis(rx, ry, rz);
            axis.normalize();
            matrixStack.top() = matrixMultiplication(matrixStack.top(), RotationMatrix(angle, axis));
        }
        // Draw triangle and output necessary values to the files
        else if (command == "triangle")
        {
            Point p[3];
            // Read the triangle points
            for (int i = 0; i < 3; i++)
            {
                inFile >> p[i];
            }

            for (int i = 0; i < 3; i++)
            {
                // Stage 1 Points
                Point model = ScaleW(matrixMultiplication(matrixStack.top(), CoordinateMatrix(p[i])));
                stage1 << model;

                // Stage 2 Points, Uses stage 1 points and view transformation matrix derived during init
                model = ScaleW(matrixMultiplication(viewMatrix, CoordinateMatrix(model)));
                stage2 << model;

                // Stage 3 Points Uses stage 2 points and projection transformation matrix derived during init
                model = ScaleW(matrixMultiplication(projectionMatrix, CoordinateMatrix(model)));
                stage3 << model;
            }
            stage1 << endl;
            stage2 << endl;
            stage3 << endl;
        }
    }

    // Close all files
    inFile.close();
    configFile.close();
    stage1.close();
    stage2.close();
    stage3.close();

    // Begin stage 4
    stage4init();
    ReadTriangles();

    for (Triangle t : triangles)
    {
        ScanConversion(t);
    }

    // Print zBuffer to file
    for (int i = 0; i < Screen_Height; i++)
    {
        for (int j = 0; j < Screen_Width; j++)
        {
            if (zBufferMatrix[i][j] < zRear)
            {
                zBuffer << zBufferMatrix[i][j] << "\t";
            }
        }
        zBuffer << endl;
    }

    // Generate BMP Image
    bitmap_image outputImage(Screen_Width, Screen_Height);
    for (int i = 0; i < Screen_Height; i++)
    {
        for (int j = 0; j < Screen_Width; j++)
        {
            outputImage.set_pixel(j, i, colorBufferMatrix[i][j].r, colorBufferMatrix[i][j].g, colorBufferMatrix[i][j].b);
        }
    }
    outputImage.save_image("out.bmp");

    // Delete Memories
    for (int i = 0; i < Screen_Height; i++)
    {
        delete[] zBufferMatrix[i];
        delete[] colorBufferMatrix[i];
    }
    delete[] zBufferMatrix;
    delete[] colorBufferMatrix;

    zBuffer.close();

    return 0;
}
