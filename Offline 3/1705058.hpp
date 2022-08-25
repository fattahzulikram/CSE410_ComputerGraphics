#ifndef _1705058_HPP
#define _1705058_HPP

#ifdef _WIN32
#include <windows.h>
#endif

#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>

#include <GL/glut.h>

#define PI 2 * acos(0.0)
#define INFINITE 999999999
// MARGIN is used to make sure the check to see whether an object is obscured or not is accurate. If not used, the check will ignore each object
#define MARGIN 0.000001

// Sphere drawing
#define SPHERE_STACK 50
#define SPHERE_SLICES 24

// Degree to radian conversion
#define RADIAN(t) (t * PI / 180.0)

// Recursive Reflection
int iRecursionDepth;

/**************************************Vector Structure*********************************************/
struct Vector
{
    double dX, dY, dZ;

    // Constructors
    Vector()
    {
        dX = dY = dZ = 0;
    }

    Vector(double x, double y, double z)
    {
        this->dX = x;
        this->dY = y;
        this->dZ = z;
    }

    // Operators
    Vector operator+(const Vector &v) const
    {
        return Vector(dX + v.dX, dY + v.dY, dZ + v.dZ);
    }

    Vector operator-(const Vector &v) const
    {
        return Vector(dX - v.dX, dY - v.dY, dZ - v.dZ);
    }

    Vector operator*(const double &d) const
    {
        return Vector(dX * d, dY * d, dZ * d);
    }

    Vector operator/(const double &d) const
    {
        return Vector(dX / d, dY / d, dZ / d);
    }

    Vector operator-() const
    {
        return Vector(-dX, -dY, -dZ);
    }

    Vector &operator=(const Vector &v)
    {
        dX = v.dX;
        dY = v.dY;
        dZ = v.dZ;
        return *this;
    }

    // IO Operators
    friend std::ostream &operator<<(std::ostream &os, const Vector &v)
    {
        os << "(" << v.dX << ", " << v.dY << ", " << v.dZ << ")";
        return os;
    }

    friend std::istream &operator>>(std::istream &is, Vector &v)
    {
        is >> v.dX >> v.dY >> v.dZ;
        return is;
    }

    // Cross product
    Vector operator^(const Vector &v) const
    {
        return Vector(dY * v.dZ - dZ * v.dY, dZ * v.dX - dX * v.dZ, dX * v.dY - dY * v.dX);
    }

    // Dot product
    double operator*(const Vector &v) const
    {
        return dX * v.dX + dY * v.dY + dZ * v.dZ;
    }

    // Normalize
    Vector normalize() const
    {
        return *this / sqrt(dX * dX + dY * dY + dZ * dZ);
    }

    // Magnitude
    double Magnitude() const
    {
        return sqrt(dX * dX + dY * dY + dZ * dZ);
    }
};

double GetDistance(const Vector &v1, const Vector &v2)
{
    return sqrt((v1.dX - v2.dX) * (v1.dX - v2.dX) + (v1.dY - v2.dY) * (v1.dY - v2.dY) + (v1.dZ - v2.dZ) * (v1.dZ - v2.dZ));
}

Vector GetNormalOfSphere(Vector v1, Vector v2)
{
    return (v1 - v2).normalize();
}

Vector GetNormalOfTriangle(Vector v1, Vector v2, Vector v3)
{
    return ((v2 - v1) ^ (v3 - v1)).normalize();
}

Vector GetNormalOfQuadricShape(double dA, double dB, double dC, double dD, double dE, double dF, double dG, double dH, double dI, double dJ, Vector vIntersectingPoint)
{
    // F(x,y,z) = Ax^2+By^2+Cz^2+Dxy+Exz+Fyz+Gx+Hy+Iz+J = 0
    // del x = 2Ax+Dy+Ez+G
    double dDelX = 2 * dA * vIntersectingPoint.dX + dD * vIntersectingPoint.dY + dE * vIntersectingPoint.dZ + dG;

    // del y = 2By+Dx+Fz+H
    double dDelY = 2 * dB * vIntersectingPoint.dY + dD * vIntersectingPoint.dX + dF * vIntersectingPoint.dZ + dH;

    // del z = 2Cz+Ex+Fy+I
    double dDelZ = 2 * dC * vIntersectingPoint.dZ + dE * vIntersectingPoint.dX + dF * vIntersectingPoint.dY + dI;

    return Vector(dDelX, dDelY, dDelZ).normalize();
}

Vector GetReflectedRay(Vector v, Vector n)
{
    return v - n * (v * n) * 2;
}

double GetAngleBetweenVectors(Vector V1, Vector V2)
{
    double dDotProduct = V1 * V2;
    double dMagnitude1 = V1.Magnitude();
    double dMagnitude2 = V2.Magnitude();
    double dCosTheta = dDotProduct / (dMagnitude1 * dMagnitude2);

    return acos(dCosTheta);
}

// Checking purpose
Vector vCamera;

/**************************************Color Structure*********************************************/
struct Color
{
    double dRed, dGreen, dBlue;
    Color()
    {
        dRed = dGreen = dBlue = 0;
    }
    Color(double r, double g, double b)
    {
        this->dRed = r;
        this->dGreen = g;
        this->dBlue = b;
    }

    // Operator Overloading
    Color operator*(const double &d) const
    {
        return Color(dRed * d, dGreen * d, dBlue * d);
    }

    Color operator+(const Color &c) const
    {
        return Color(dRed + c.dRed, dGreen + c.dGreen, dBlue + c.dBlue);
    }

    Color operator*(const Color &c) const
    {
        return Color(dRed * c.dRed, dGreen * c.dGreen, dBlue * c.dBlue);
    }

    Color &operator=(const Color &c)
    {
        dRed = c.dRed;
        dGreen = c.dGreen;
        dBlue = c.dBlue;
        return *this;
    }

    // IO Operator Overloading
    friend std::istream &operator>>(std::istream &is, Color &c)
    {
        is >> c.dRed >> c.dGreen >> c.dBlue;
        return is;
    }

    // Color CLipping
    void clip()
    {
        if (dRed > 1)
            dRed = 1;
        if (dRed < 0)
            dRed = 0;
        if (dGreen > 1)
            dGreen = 1;
        if (dGreen < 0)
            dGreen = 0;
        if (dBlue > 1)
            dBlue = 1;
        if (dBlue < 0)
            dBlue = 0;
    }
};

/**************************************Coefficient Structure*********************************************/
struct Coefficients
{
    double dAmbientCoefficient, dDiffuseCoefficient, dSpecularCoefficient, dReflectionCoefficient;
    Coefficients()
    {
        dAmbientCoefficient = dDiffuseCoefficient = dSpecularCoefficient = dReflectionCoefficient = 0;
    }

    Coefficients(double ambientCoefficient, double diffuseCoefficient, double specularCoefficient, double reflectionCoefficient)
    {
        this->dAmbientCoefficient = ambientCoefficient;
        this->dDiffuseCoefficient = diffuseCoefficient;
        this->dSpecularCoefficient = specularCoefficient;
        this->dReflectionCoefficient = reflectionCoefficient;
    }

    // IO Operator Overloading
    friend std::istream &operator>>(std::istream &is, Coefficients &c)
    {
        is >> c.dAmbientCoefficient >> c.dDiffuseCoefficient >> c.dSpecularCoefficient >> c.dReflectionCoefficient;
        return is;
    }
};

/**************************************Ray Class*********************************************/
class Ray
{
    Vector vOrigin, vDirection;

public:
    // Constructors
    Ray()
    {
        vOrigin = Vector(0, 0, 0);
        vDirection = Vector(0, 0, 0);
    }
    Ray(Vector Origin, Vector Direction)
    {
        vOrigin = Origin;
        vDirection = Direction.normalize();
    }

    // Getters
    inline Vector getOrigin() { return vOrigin; }
    inline Vector getDirection() { return vDirection; }

    // Destructor
    ~Ray()
    {
    }
};

/**************************************Point Light Class*********************************************/
class PointLight
{
    Vector vLightPosition;
    Color cLightColor;

public:
    // Constructors
    PointLight()
    {
        vLightPosition = Vector(0, 0, 0);
        cLightColor = Color(0, 0, 0);
    }
    PointLight(Vector LightPosition, Color LightColor)
    {
        vLightPosition = LightPosition;
        cLightColor = LightColor;
    }

    // Getters
    inline Vector getLightPosition() { return vLightPosition; }
    inline Color getLightColor() { return cLightColor; }

    // Destructor
    ~PointLight()
    {
    }
};

// Extern Variable
std::vector<PointLight *> vPointLights;

/**************************************Spot Light Class*********************************************/
class SpotLight: public PointLight
{
    Vector vLightDirection;
    double dCutoffAngle;

public:
    // Constructors
    SpotLight()
    {
        vLightDirection = Vector(0, 0, 0);
        dCutoffAngle = 0;
    }

    SpotLight(Vector LightPosition, Color LightColor, Vector LightDirection, double CutoffAngle)
    : PointLight(LightPosition, LightColor)
    {
        vLightDirection = LightDirection.normalize();
        dCutoffAngle = CutoffAngle;
    }

    // Getters
    inline Vector getLightDirection() { return vLightDirection; }
    inline double getCutoffAngle() { return dCutoffAngle; }

    // Destructor
    ~SpotLight()
    {
    }
};

// Extern Variable
std::vector<SpotLight *> vSpotLights;

/**************************************Base Object Class*********************************************/
class Object
{
    Color cColor;
    Coefficients coCoefficients;
    int Shine; // exponent term of specular component

public:
    // Constructors
    Object()
    {
        cColor = Color(0, 0, 0);
        coCoefficients = Coefficients(0, 0, 0, 0);
        Shine = 0;
    }
    Object(Color color, Coefficients coefficients, int Shine)
    {
        this->cColor = color;
        this->coCoefficients = coefficients;
        this->Shine = Shine;
    }

    // Destructor
    virtual ~Object()
    {
    }

    // Getters
    inline Color getColor() const
    {
        return cColor;
    }

    inline Coefficients getCoefficients() const
    {
        return coCoefficients;
    }

    inline int getShine() const
    {
        return Shine;
    }

    // Setters
    inline void setColor(Color color)
    {
        this->cColor = color;
    }

    inline void setCoefficients(Coefficients coefficients)
    {
        this->coCoefficients = coefficients;
    }

    inline void setShine(int Shine)
    {
        this->Shine = Shine;
    }

    // Draw Function
    virtual void draw() = 0;

    // Intersect Function
    virtual double Intersect(Ray *r, Color &color, int level) = 0;
};

// Extern Variable
std::vector<Object *> vObjects;

/**************************************Sphere Class*********************************************/
class Sphere : public Object
{
    Vector vCenter;
    double dRadius;

public:
    // Constructors
    Sphere()
    {
        vCenter = Vector(0, 0, 0);
        dRadius = 0;
    }
    Sphere(Vector center, double radius, Color color, Coefficients coefficients, int Shine)
        : Object(color, coefficients, Shine)
    {
        this->vCenter = center;
        this->dRadius = radius;
    }

    // Implementation of draw function
    void draw()
    {
        Vector Points[SPHERE_STACK + 1][SPHERE_SLICES + 1];
        // Generate points
        for (int i = 0; i <= SPHERE_STACK; i++)
        {
            double height = dRadius * sin(((double)i / (double)SPHERE_STACK) * (PI / 2));
            double radius2 = dRadius * cos(((double)i / (double)SPHERE_STACK) * (PI / 2));

            for (int j = 0; j <= SPHERE_SLICES; j++)
            {
                Points[i][j].dX = radius2 * cos(((double)j / (double)SPHERE_SLICES) * 2 * PI);
                Points[i][j].dY = radius2 * sin(((double)j / (double)SPHERE_SLICES) * 2 * PI);
                Points[i][j].dZ = height;
            }
        }

        // draw quads using generated points
        for (int i = 0; i < SPHERE_STACK; i++)
        {
            for (int j = 0; j < SPHERE_SLICES; j++)
            {
                glBegin(GL_QUADS);
                {
                    glColor3f(getColor().dRed, getColor().dGreen, getColor().dBlue);
                    // Upper hemisphere
                    glVertex3f(Points[i][j].dX + vCenter.dX, Points[i][j].dY + vCenter.dY, Points[i][j].dZ + vCenter.dZ);
                    glVertex3f(Points[i][j + 1].dX + vCenter.dX, Points[i][j + 1].dY + vCenter.dY, Points[i][j + 1].dZ + vCenter.dZ);
                    glVertex3f(Points[i + 1][j + 1].dX + vCenter.dX, Points[i + 1][j + 1].dY + vCenter.dY, Points[i + 1][j + 1].dZ + vCenter.dZ);
                    glVertex3f(Points[i + 1][j].dX + vCenter.dX, Points[i + 1][j].dY + vCenter.dY, Points[i + 1][j].dZ + vCenter.dZ);

                    // Lower hemisphere
                    glVertex3f(Points[i][j].dX + vCenter.dX, Points[i][j].dY + vCenter.dY, -Points[i][j].dZ + vCenter.dZ);
                    glVertex3f(Points[i][j + 1].dX + vCenter.dX, Points[i][j + 1].dY + vCenter.dY, -Points[i][j + 1].dZ + vCenter.dZ);
                    glVertex3f(Points[i + 1][j + 1].dX + vCenter.dX, Points[i + 1][j + 1].dY + vCenter.dY, -Points[i + 1][j + 1].dZ + vCenter.dZ);
                    glVertex3f(Points[i + 1][j].dX + vCenter.dX, Points[i + 1][j].dY + vCenter.dY, -Points[i + 1][j].dZ + vCenter.dZ);
                }
                glEnd();
            }
        }
    }

    // Implementation of intersect function
    double Intersect(Ray *r, Color &color, int level)
    {
        // Get Intersecting Point
        double dA = 1.0;
        double dB = (r->getDirection() * (r->getOrigin() - vCenter)) * 2.0;
        double dC = (r->getOrigin() - vCenter) * (r->getOrigin() - vCenter) - dRadius * dRadius;
        double dRoot = -1;

        double dDiscriminant = dB * dB - 4.0 * dA * dC;

        if (dDiscriminant < 0)
        {
            // No intersection
            return -1;
        }
        else if (dDiscriminant == 0)
        {
            // One Intersection
            dRoot = -dB / (2.0 * dA);
        }
        else
        {
            // Two Intersections
            double dDiscriminantSqrt = sqrt(dDiscriminant);
            double dPositiveRoot = (-dB + dDiscriminantSqrt) / (2.0 * dA);
            double dNegativeRoot = (-dB - dDiscriminantSqrt) / (2.0 * dA);
            dRoot = (dNegativeRoot > 0.0) ? dNegativeRoot : dPositiveRoot;
        }

        // Level Check
        if (level <= 0)
        {
            return dRoot;
        }

        // Illumination
        // Ambient Color
        color = getColor() * getCoefficients().dAmbientCoefficient;

        // Get Intersecting Point, Normal and Reflected Ray
        Vector IntersectingPoint = r->getOrigin() + r->getDirection() * dRoot;
        Vector Normal = GetNormalOfSphere(IntersectingPoint, vCenter);
        Normal = (GetDistance(vCamera, vCenter) > dRadius) ? Normal : -Normal;
        Vector vReflectedRay = GetReflectedRay(r->getDirection(), Normal).normalize();

        // Diffuse and Specular Color
        for (int iLightCounter = 0; iLightCounter < vPointLights.size(); iLightCounter++)
        {
            // Get Light Ray
            Ray rLightRay(vPointLights[iLightCounter]->getLightPosition(), IntersectingPoint - vPointLights[iLightCounter]->getLightPosition());

            Vector L = rLightRay.getDirection();
            Vector N = Normal;
            Vector R = (GetReflectedRay(L, N) * -1).normalize();
            Vector V = (IntersectingPoint - r->getOrigin()).normalize();

            // Check if object is in shadow
            double dClosest = INFINITE;
            for (int iObjectCounter = 0; iObjectCounter < vObjects.size(); iObjectCounter++)
            {
                Color cDummyColor;
                double dTemp = vObjects[iObjectCounter]->Intersect(&rLightRay, cDummyColor, 0);

                if (dTemp > 0 && dTemp < dClosest)
                {
                    dClosest = dTemp;
                }
            }

            Vector vClosestPoint = rLightRay.getOrigin() + rLightRay.getDirection() * dClosest;
            if (GetDistance(vClosestPoint, rLightRay.getOrigin()) < GetDistance(rLightRay.getOrigin(), IntersectingPoint) - MARGIN)
            {
                // Object is obscured by another object
                continue;
            }

            // Not in shadow
            double dCosTheta = std::max(0.0, (L * (-1.0)) * N);
            double dCosPhi = std::max(0.0, (R * (-1.0)) * V);

            // Diffuse Color
            color = color + vPointLights[iLightCounter]->getLightColor() * getColor() * getCoefficients().dDiffuseCoefficient * dCosTheta;

            // Specular Color
            color = color + vPointLights[iLightCounter]->getLightColor() * getColor() * getCoefficients().dSpecularCoefficient * pow(dCosPhi, getShine());
        }

        // For Spot Lights
        for(int iSpotLightCounter = 0; iSpotLightCounter < vSpotLights.size(); iSpotLightCounter++)
        {
            // Get Light Ray
            Ray rSpotLightRay(vSpotLights[iSpotLightCounter]->getLightPosition(), IntersectingPoint - vSpotLights[iSpotLightCounter]->getLightPosition());

            Vector L = rSpotLightRay.getDirection();
            Vector N = Normal;
            Vector R = (GetReflectedRay(L, N) * -1).normalize();
            Vector V = (IntersectingPoint - r->getOrigin()).normalize();

            // Check if intersection point is in cutoff angle
            // Get angle between light direction and ray direction
            double dAngleFromSource = GetAngleBetweenVectors(vSpotLights[iSpotLightCounter]->getLightDirection(), rSpotLightRay.getDirection());
            // Since Spot Light sheds light in a cone, both positive and negative angles are allowed
            dAngleFromSource = (dAngleFromSource < 0) ? dAngleFromSource * -1.0 : dAngleFromSource;

            // If the angle is above cutoff angle, this spotlight won't affect the object
            if(dAngleFromSource > vSpotLights[iSpotLightCounter]->getCutoffAngle())
            {
                continue;
            }

            // Check if object is in shadow
            double dClosest = INFINITE;
            for (int iObjectCounter = 0; iObjectCounter < vObjects.size(); iObjectCounter++)
            {
                Color cDummyColor;
                double dTemp = vObjects[iObjectCounter]->Intersect(&rSpotLightRay, cDummyColor, 0);

                if (dTemp > 0 && dTemp < dClosest)
                {
                    dClosest = dTemp;
                }
            }

            Vector vClosestPoint = rSpotLightRay.getOrigin() + rSpotLightRay.getDirection() * dClosest;
            if (GetDistance(vClosestPoint, rSpotLightRay.getOrigin()) < GetDistance(rSpotLightRay.getOrigin(), IntersectingPoint) - MARGIN)
            {
                // Object is obscured by another object
                continue;
            }

            // Not in shadow
            double dCosTheta = std::max(0.0, (L * (-1.0)) * N);
            double dCosPhi = std::max(0.0, (R * (-1.0)) * V);

            // Diffuse Color
            color = color + vSpotLights[iSpotLightCounter]->getLightColor() * getColor() * getCoefficients().dDiffuseCoefficient * dCosTheta;

            // Specular Color
            color = color + vSpotLights[iSpotLightCounter]->getLightColor() * getColor() * getCoefficients().dSpecularCoefficient * pow(dCosPhi, getShine());
        }

        // Recursive Reflection
        if (level >= iRecursionDepth)
        {
            return dRoot;
        }

        // Construct reflected ray from intersection point
        // Moving the intersection point slightly away from the surface
        Ray rReflectedRay(IntersectingPoint + vReflectedRay, vReflectedRay);

        // Find tmin from the nearest intersecting object
        int iNearestObject = -1;
        double dNearestObjectDistance = INFINITE;

        for (int iObjectCounter = 0; iObjectCounter < vObjects.size(); iObjectCounter++)
        {
            Color cDummyColor;
            double dTemp = vObjects[iObjectCounter]->Intersect(&rReflectedRay, cDummyColor, 0);

            if (dTemp > 0 && dTemp < dNearestObjectDistance)
            {
                dNearestObjectDistance = dTemp;
                iNearestObject = iObjectCounter;
            }
        }

        Color cReflectedColor;

        // If found, call intersect(rreflected, colorreflected, level+1)
        if (iNearestObject != -1)
        {
            vObjects[iNearestObject]->Intersect(&rReflectedRay, cReflectedColor, level + 1);
        }

        // update color using the impact of reflection
        color = color + cReflectedColor * getCoefficients().dReflectionCoefficient;
        color.clip();

        return dRoot;
    }
};

/**************************************Triangle Class*********************************************/
class Triangle : public Object
{
    Vector vPointA, vPointB, vPointC;

public:
    // Constructors
    Triangle()
    {
        vPointA = Vector(0, 0, 0);
        vPointB = Vector(0, 0, 0);
        vPointC = Vector(0, 0, 0);
    }
    Triangle(Vector PointA, Vector PointB, Vector PointC, Color color, Coefficients coefficients, int Shine)
        : Object(color, coefficients, Shine)
    {
        this->vPointA = PointA;
        this->vPointB = PointB;
        this->vPointC = PointC;
    }

    // Implementation of draw function
    void draw()
    {
        glBegin(GL_TRIANGLES);
        {
            glColor3f(getColor().dRed, getColor().dGreen, getColor().dBlue);
            glVertex3f(vPointA.dX, vPointA.dY, vPointA.dZ);
            glVertex3f(vPointB.dX, vPointB.dY, vPointB.dZ);
            glVertex3f(vPointC.dX, vPointC.dY, vPointC.dZ);
        }
        glEnd();
    }

    // Implementation of intersect function
    double Intersect(Ray *r, Color &color, int level)
    {
        double dRoot = -1.0;

        //  Use barycentric coordinates
        // |A| = (ax - bx) * ( ((ay - cy) * rz) - ((az- cz) * ry) ) + (ax - cx) * ( ((az-bz) * ry) - ((ay-by) * rz) ) + rx * ( (ay-by) * (az-cz) - (ay-cy) * (az-bz) )
        double dDeterminant = (vPointA.dX - vPointB.dX) * (((vPointA.dY - vPointC.dY) * r->getDirection().dZ) - ((vPointA.dZ - vPointC.dZ) * r->getDirection().dY)) +
                              (vPointA.dX - vPointC.dX) * (((vPointA.dZ - vPointB.dZ) * r->getDirection().dY) - ((vPointA.dY - vPointB.dY) * r->getDirection().dZ)) +
                              r->getDirection().dX * (((vPointA.dY - vPointB.dY) * (vPointA.dZ - vPointC.dZ)) - ((vPointA.dY - vPointC.dY) * (vPointA.dZ - vPointB.dZ)));

        // |B| = (ax - rox) * ( ( (ay-cy) * rz ) - ( (az-cz) * ry ) ) + (ax-cx) * ( ( ( az-roz ) * ry ) - ( (ay-roy) * rz ) ) + rx * ( ( (ay-roy) * (az-cz) ) - ( (ay-cy) * (az-roz) ) )
        double dBeta = (vPointA.dX - r->getOrigin().dX) * (((vPointA.dY - vPointC.dY) * r->getDirection().dZ) - ((vPointA.dZ - vPointC.dZ) * r->getDirection().dY)) +
                       (vPointA.dX - vPointC.dX) * (((vPointA.dZ - r->getOrigin().dZ) * r->getDirection().dY) - ((vPointA.dY - r->getOrigin().dY) * r->getDirection().dZ)) +
                       r->getDirection().dX * (((vPointA.dY - r->getOrigin().dY) * (vPointA.dZ - vPointC.dZ)) - ((vPointA.dY - vPointC.dY) * (vPointA.dZ - r->getOrigin().dZ)));

        // |G| = (ax-bx) * ( ( ( ay-roy ) * rz ) - ( ( az-roz ) * ry ) ) + (ax-rox) * ( ( ( az-bz ) * ry ) - ( ( ay-by ) * rz ) ) + rx * ( ( (ay-by) * (az-roz) ) - ( (ay-roy) * (az-bz) ) )
        double dGamma = (vPointA.dX - vPointB.dX) * (((vPointA.dY - r->getOrigin().dY) * r->getDirection().dZ) - ((vPointA.dZ - r->getOrigin().dZ) * r->getDirection().dY)) +
                        (vPointA.dX - r->getOrigin().dX) * (((vPointA.dZ - vPointB.dZ) * r->getDirection().dY) - ((vPointA.dY - vPointB.dY) * r->getDirection().dZ)) +
                        r->getDirection().dX * (((vPointA.dY - vPointB.dY) * (vPointA.dZ - r->getOrigin().dZ)) - ((vPointA.dY - r->getOrigin().dY) * (vPointA.dZ - vPointB.dZ)));

        // |T| = (ax-bx) * ( ( ( ay-cy ) * ( az-roz ) ) - ( ( az-cz ) * ( ay-roy ) ) ) + (ax-cx) * ( ( ( ay-roy ) * ( az-bz ) ) - ( ( ay-by ) * ( az-roz ) ) ) + (ax-rox) * ( ( ( ay-by ) * ( az-cz ) ) - ( ( ay-cy ) * ( az-bz ) ) )
        double dT = (vPointA.dX - vPointB.dX) * (((vPointA.dY - vPointC.dY) * (vPointA.dZ - r->getOrigin().dZ)) - ((vPointA.dZ - vPointC.dZ) * (vPointA.dY - r->getOrigin().dY))) +
                    (vPointA.dX - vPointC.dX) * (((vPointA.dY - r->getOrigin().dY) * (vPointA.dZ - vPointB.dZ)) - ((vPointA.dY - vPointB.dY) * (vPointA.dZ - r->getOrigin().dZ))) +
                    (vPointA.dX - r->getOrigin().dX) * (((vPointA.dY - vPointB.dY) * (vPointA.dZ - vPointC.dZ)) - ((vPointA.dY - vPointC.dY) * (vPointA.dZ - vPointB.dZ)));

        // If determinant is 0, then the ray is parallel to the plane
        if (dDeterminant != 0.0)
        {
            dBeta = dBeta / dDeterminant;
            dGamma = dGamma / dDeterminant;
            dT = dT / dDeterminant;
            // Intersection if beta + gamma < 1 and beta > 0, gamma > 0, t > 0
            if ((dBeta + dGamma) < 1.0 && dBeta > 0.0 && dGamma > 0.0 && dT > 0.0)
            {
                dRoot = dT;
            }
        }

        // Level Check
        if (level <= 0)
        {
            return dRoot;
        }

        // Illumination
        // Ambient Color
        color = getColor() * getCoefficients().dAmbientCoefficient;

        // Get Intersecting Point, Normal and Reflected Ray
        Vector IntersectingPoint = r->getOrigin() + r->getDirection() * dRoot;
        Vector Normal = GetNormalOfTriangle(vPointA, vPointB, vPointC);
        Normal = ((r->getDirection() * (-1.0)) * Normal > 0) ? Normal : Normal * (-1.0);
        Vector vReflectedRay = GetReflectedRay(r->getDirection(), Normal).normalize();

        // Diffuse and Specular Color
        for (int iLightCounter = 0; iLightCounter < vPointLights.size(); iLightCounter++)
        {
            // Get Light Ray
            Ray rLightRay(vPointLights[iLightCounter]->getLightPosition(), IntersectingPoint - vPointLights[iLightCounter]->getLightPosition());

            Vector L = rLightRay.getDirection();
            Vector N = Normal;
            Vector R = (GetReflectedRay(L, N) * -1).normalize();
            Vector V = (IntersectingPoint - r->getOrigin()).normalize();

            // Check if object is in shadow
            double dClosest = INFINITE;
            for (int iObjectCounter = 0; iObjectCounter < vObjects.size(); iObjectCounter++)
            {
                Color cDummyColor;
                double dTemp = vObjects[iObjectCounter]->Intersect(&rLightRay, cDummyColor, 0);

                if (dTemp > 0 && dTemp < dClosest)
                {
                    dClosest = dTemp;
                }
            }

            Vector vClosestPoint = rLightRay.getOrigin() + rLightRay.getDirection() * dClosest;
            if (GetDistance(vClosestPoint, rLightRay.getOrigin()) < GetDistance(rLightRay.getOrigin(), IntersectingPoint) - MARGIN)
            {
                // Object is obscured by another object
                continue;
            }

            // Not in shadow
            double dCosTheta = std::max(0.0, (L * -1.0) * N);
            double dCosPhi = std::max(0.0, (R * -1.0) * V);

            // Diffuse Color
            color = color + vPointLights[iLightCounter]->getLightColor() * getColor() * getCoefficients().dDiffuseCoefficient * dCosTheta;

            // Specular Color
            color = color + vPointLights[iLightCounter]->getLightColor() * getColor() * getCoefficients().dSpecularCoefficient * pow(dCosPhi, getShine());
        }

        // For Spot Lights
        for(int iSpotLightCounter = 0; iSpotLightCounter < vSpotLights.size(); iSpotLightCounter++)
        {
            // Get Light Ray
            Ray rSpotLightRay(vSpotLights[iSpotLightCounter]->getLightPosition(), IntersectingPoint - vSpotLights[iSpotLightCounter]->getLightPosition());

            Vector L = rSpotLightRay.getDirection();
            Vector N = Normal;
            Vector R = (GetReflectedRay(L, N) * -1).normalize();
            Vector V = (IntersectingPoint - r->getOrigin()).normalize();

            // Check if intersection point is in cutoff angle
            // Get angle between light direction and ray direction
            double dAngleFromSource = GetAngleBetweenVectors(vSpotLights[iSpotLightCounter]->getLightDirection(), rSpotLightRay.getDirection());
            // Since Spot Light sheds light in a cone, both positive and negative angles are allowed
            dAngleFromSource = (dAngleFromSource < 0) ? dAngleFromSource * -1.0 : dAngleFromSource;

            // If the angle is above cutoff angle, this spotlight won't affect the object
            if(dAngleFromSource > vSpotLights[iSpotLightCounter]->getCutoffAngle())
            {
                continue;
            }

            // Check if object is in shadow
            double dClosest = INFINITE;
            for (int iObjectCounter = 0; iObjectCounter < vObjects.size(); iObjectCounter++)
            {
                Color cDummyColor;
                double dTemp = vObjects[iObjectCounter]->Intersect(&rSpotLightRay, cDummyColor, 0);

                if (dTemp > 0 && dTemp < dClosest)
                {
                    dClosest = dTemp;
                }
            }

            Vector vClosestPoint = rSpotLightRay.getOrigin() + rSpotLightRay.getDirection() * dClosest;
            if (GetDistance(vClosestPoint, rSpotLightRay.getOrigin()) < GetDistance(rSpotLightRay.getOrigin(), IntersectingPoint) - MARGIN)
            {
                // Object is obscured by another object
                continue;
            }

            // Not in shadow
            double dCosTheta = std::max(0.0, (L * (-1.0)) * N);
            double dCosPhi = std::max(0.0, (R * (-1.0)) * V);

            // Diffuse Color
            color = color + vSpotLights[iSpotLightCounter]->getLightColor() * getColor() * getCoefficients().dDiffuseCoefficient * dCosTheta;

            // Specular Color
            color = color + vSpotLights[iSpotLightCounter]->getLightColor() * getColor() * getCoefficients().dSpecularCoefficient * pow(dCosPhi, getShine());
        }

        // Recursive Reflection
        if (level >= iRecursionDepth)
        {
            return dRoot;
        }

        // Construct reflected ray from intersection point
        // Moving the intersection point slightly away from the surface
        Ray rReflectedRay(IntersectingPoint + vReflectedRay, vReflectedRay);

        // Find tmin from the nearest intersecting object
        int iNearestObject = -1;
        double dNearestObjectDistance = INFINITE;

        for (int iObjectCounter = 0; iObjectCounter < vObjects.size(); iObjectCounter++)
        {
            Color cDummyColor;
            double dTemp = vObjects[iObjectCounter]->Intersect(&rReflectedRay, cDummyColor, 0);

            if (dTemp > 0 && dTemp < dNearestObjectDistance)
            {
                dNearestObjectDistance = dTemp;
                iNearestObject = iObjectCounter;
            }
        }

        Color cReflectedColor;

        // If found, call intersect(rreflected, colorreflected, level+1)
        if (iNearestObject != -1)
        {
            vObjects[iNearestObject]->Intersect(&rReflectedRay, cReflectedColor, level + 1);
        }

        // update color using the impact of reflection
        color = color + cReflectedColor * getCoefficients().dReflectionCoefficient;
        color.clip();

        return dRoot;
    }
};

/**************************************General Quadric Class*********************************************/
class Quadric : public Object
{
    double dA, dB, dC, dD, dE, dF, dG, dH, dI, dJ;
    Vector vReferencePoint;
    double dLength, dWidth, dHeight;

public:
    // Constructors
    Quadric()
    {
        dA = 0;
        dB = 0;
        dC = 0;
        dD = 0;
        dE = 0;
        dF = 0;
        dG = 0;
        dH = 0;
        dI = 0;
        dJ = 0;
        vReferencePoint = Vector(0, 0, 0);
        dLength = 0;
        dWidth = 0;
        dHeight = 0;
    }

    Quadric(double A, double B, double C, double D, double E, double F, double G, double H, double I, double J, Vector ReferencePoint, double Length, double Width, double Height, Color color, Coefficients coefficients, int Shine)
        : Object(color, coefficients, Shine)
    {
        this->dA = A;
        this->dB = B;
        this->dC = C;
        this->dD = D;
        this->dE = E;
        this->dF = F;
        this->dG = G;
        this->dH = H;
        this->dI = I;
        this->dJ = J;
        this->vReferencePoint = ReferencePoint;
        this->dLength = Length;
        this->dWidth = Width;
        this->dHeight = Height;
    }

    //  However, you do not have to draw the general quadric surfaces
    void draw()
    {
        // Do nothing
    }

    // Implementation of intersect function
    double Intersect(Ray *r, Color &color, int level)
    {
        // Two values of t are possible
        double dRoot = -1.0;
        double dRoot1 = -1.0;
        // Get Intersection Point
        // http://web.cse.ohio-state.edu/~shen.94/681/Site/Slides_files/intersection2.pdf
        // An Introduction to Ray Tracing, Chapter 2, Page 68
        // Aq = A * xd^2 + B * yd^2 + C * zd^2 + D* xd * yd + E * xd * zd + F * yd * zd
        double dAQ = dA * r->getDirection().dX * r->getDirection().dX + dB * r->getDirection().dY * r->getDirection().dY + dC * r->getDirection().dZ * r->getDirection().dZ + dD * r->getDirection().dX * r->getDirection().dY + dE * r->getDirection().dX * r->getDirection().dZ + dF * r->getDirection().dY * r->getDirection().dZ;

        // Bq = 2*A* xr * xd + 2* B * yr * yd + 2 * C * zr * zd + D * (xr * yd + yr * xd) + E * (xr * zd + zr * xd) + F * (yr * zd + yd * zr ) + G * xd + H * yd + I * zd
        double dBQ = 2 * dA * r->getOrigin().dX * r->getDirection().dX + 2 * dB * r->getOrigin().dY * r->getDirection().dY + 2 * dC * r->getOrigin().dZ * r->getDirection().dZ + dD * (r->getOrigin().dX * r->getDirection().dY + r->getOrigin().dY * r->getDirection().dX) + dE * (r->getOrigin().dX * r->getDirection().dZ + r->getOrigin().dZ * r->getDirection().dX) + dF * (r->getOrigin().dY * r->getDirection().dZ + r->getOrigin().dZ * r->getDirection().dY) + dG * r->getDirection().dX + dH * r->getDirection().dY + dI * r->getDirection().dZ;

        // Cq = A * xr^2 + B * yr^2 + C * zr^2 + D * xr * yr + E * xr * zr + F * yr * zr + G * xr + H * yr + I * zr + J
        double dCQ = dA * r->getOrigin().dX * r->getOrigin().dX + dB * r->getOrigin().dY * r->getOrigin().dY + dC * r->getOrigin().dZ * r->getOrigin().dZ + dD * r->getOrigin().dX * r->getOrigin().dY + dE * r->getOrigin().dX * r->getOrigin().dZ + dF * r->getOrigin().dY * r->getOrigin().dZ + dG * r->getOrigin().dX + dH * r->getOrigin().dY + dI * r->getOrigin().dZ + dJ;

        double dDiscriminant = dBQ * dBQ - 4.0 * dAQ * dCQ;

        // Check on dAQ
        if (dAQ != 0.0)
        {
            // Check the discriminant
            // If the discriminant is less than 0, there is no intersection, so no change in dRoot
            if (dDiscriminant == 0)
            {
                // If the discriminant is 0, there is one intersection
                dRoot = -dBQ / (2 * dAQ);
            }
            if (dDiscriminant > 0)
            {
                // If the discriminant is greater than 0, there are two intersections
                dRoot1 = (-dBQ + sqrt(dDiscriminant)) / (2 * dAQ);
                dRoot = (-dBQ - sqrt(dDiscriminant)) / (2 * dAQ);
            }
        }
        // dAQ can be 0
        else
        {
            if (dBQ != 0.0)
            {
                dRoot = -dCQ / dBQ;
            }
        }

        // If two values of t are obtained, check which one (or none or both) falls within the reference cube
        // 0 indicates no clipping along a dimension per specification
        if (dRoot != -1)
        {
            // One root is obtained at least
            Vector vIntersectionPoint = r->getOrigin() + r->getDirection() * dRoot;

            // Check if the intersection point is within the reference cube
            if ((dLength != 0.0 && (vIntersectionPoint.dX < vReferencePoint.dX || vIntersectionPoint.dX > vReferencePoint.dX + dLength)) || (dWidth != 0.0 && (vIntersectionPoint.dY < vReferencePoint.dY || vIntersectionPoint.dY > vReferencePoint.dY + dWidth)) || (dHeight != 0.0 && (vIntersectionPoint.dZ < vReferencePoint.dZ || vIntersectionPoint.dZ > vReferencePoint.dZ + dHeight)))
            {
                dRoot = -1;
            }
            
            if (dRoot1 != -1)
            {
                // Two roots are obtained
                Vector vIntersectionPoint1 = r->getOrigin() + r->getDirection() * dRoot1;
                // Check if the intersection point is within the reference cube
                if ((dLength != 0.0 && (vIntersectionPoint1.dX < vReferencePoint.dX || vIntersectionPoint1.dX > vReferencePoint.dX + dLength)) || (dWidth != 0.0 && (vIntersectionPoint1.dY < vReferencePoint.dY || vIntersectionPoint1.dY > vReferencePoint.dY + dWidth)) || (dHeight != 0.0 && (vIntersectionPoint1.dZ < vReferencePoint.dZ || vIntersectionPoint1.dZ > vReferencePoint.dZ + dHeight)))
                {
                    dRoot1 = -1;
                }
            }
            // Choose the one that is closer to the camera
            if(dRoot > 0 && dRoot1 > 0)
            {
                // Two roots are obtained, check which one is closer to the camera
                if(dRoot > dRoot1)
                {
                    dRoot = dRoot1;
                }
            }
            else
            {
                // One root is obtained, obtain the correct one
                if(dRoot < 0 && dRoot1 > 0)
                {
                    dRoot = dRoot1;
                }
            }
        }

        // Level Check
        if (level <= 0)
        {
            return dRoot;
        }

        // Illumination
        // Ambient Color
        color = getColor() * getCoefficients().dAmbientCoefficient;

        // Get Intersecting Point, Normal and Reflected Ray
        Vector IntersectingPoint = r->getOrigin() + r->getDirection() * dRoot;
        Vector Normal = GetNormalOfQuadricShape(dA, dB, dC, dD, dE, dF, dG, dH, dI, dJ, IntersectingPoint);
        Normal = ((r->getDirection() * (-1.0)) * Normal > 0) ? Normal : Normal * (-1.0);
        Vector vReflectedRay = GetReflectedRay(r->getDirection(), Normal).normalize();

        // Diffuse and Specular Color
        for (int iLightCounter = 0; iLightCounter < vPointLights.size(); iLightCounter++)
        {
            // Get Light Ray
            Ray rLightRay(vPointLights[iLightCounter]->getLightPosition(), IntersectingPoint - vPointLights[iLightCounter]->getLightPosition());

            Vector L = rLightRay.getDirection();
            Vector N = Normal;
            Vector R = (GetReflectedRay(L, N) * -1).normalize();
            Vector V = (IntersectingPoint - r->getOrigin()).normalize();

            // Check if object is in shadow
            double dClosest = INFINITE;
            for (int iObjectCounter = 0; iObjectCounter < vObjects.size(); iObjectCounter++)
            {
                Color cDummyColor;
                double dTemp = vObjects[iObjectCounter]->Intersect(&rLightRay, cDummyColor, 0);

                if (dTemp > 0 && dTemp < dClosest)
                {
                    dClosest = dTemp;
                }
            }

            Vector vClosestPoint = rLightRay.getOrigin() + rLightRay.getDirection() * dClosest;
            if (GetDistance(vClosestPoint, rLightRay.getOrigin()) < GetDistance(rLightRay.getOrigin(), IntersectingPoint) - MARGIN)
            {
                // Object is obscured by another object
                continue;
            }

            // Not in shadow
            double dCosTheta = std::max(0.0, (L * -1.0) * N);
            double dCosPhi = std::max(0.0, (R * -1.0) * V);

            // Diffuse Color
            color = color + vPointLights[iLightCounter]->getLightColor() * getColor() * getCoefficients().dDiffuseCoefficient * dCosTheta;

            // Specular Color
            color = color + vPointLights[iLightCounter]->getLightColor() * getColor() * getCoefficients().dSpecularCoefficient * pow(dCosPhi, getShine());
        }

        // For Spot Lights
        for(int iSpotLightCounter = 0; iSpotLightCounter < vSpotLights.size(); iSpotLightCounter++)
        {
            // Get Light Ray
            Ray rSpotLightRay(vSpotLights[iSpotLightCounter]->getLightPosition(), IntersectingPoint - vSpotLights[iSpotLightCounter]->getLightPosition());

            Vector L = rSpotLightRay.getDirection();
            Vector N = Normal;
            Vector R = (GetReflectedRay(L, N) * -1).normalize();
            Vector V = (IntersectingPoint - r->getOrigin()).normalize();

            // Check if intersection point is in cutoff angle
            // Get angle between light direction and ray direction
            double dAngleFromSource = GetAngleBetweenVectors(vSpotLights[iSpotLightCounter]->getLightDirection(), rSpotLightRay.getDirection());
            // Since Spot Light sheds light in a cone, both positive and negative angles are allowed
            dAngleFromSource = (dAngleFromSource < 0) ? dAngleFromSource * -1.0 : dAngleFromSource;

            // If the angle is above cutoff angle, this spotlight won't affect the object
            if(dAngleFromSource > vSpotLights[iSpotLightCounter]->getCutoffAngle())
            {
                continue;
            }

            // Check if object is in shadow
            double dClosest = INFINITE;
            for (int iObjectCounter = 0; iObjectCounter < vObjects.size(); iObjectCounter++)
            {
                Color cDummyColor;
                double dTemp = vObjects[iObjectCounter]->Intersect(&rSpotLightRay, cDummyColor, 0);

                if (dTemp > 0 && dTemp < dClosest)
                {
                    dClosest = dTemp;
                }
            }

            Vector vClosestPoint = rSpotLightRay.getOrigin() + rSpotLightRay.getDirection() * dClosest;
            if (GetDistance(vClosestPoint, rSpotLightRay.getOrigin()) < GetDistance(rSpotLightRay.getOrigin(), IntersectingPoint) - MARGIN)
            {
                // Object is obscured by another object
                continue;
            }

            // Not in shadow
            double dCosTheta = std::max(0.0, (L * (-1.0)) * N);
            double dCosPhi = std::max(0.0, (R * (-1.0)) * V);

            // Diffuse Color
            color = color + vSpotLights[iSpotLightCounter]->getLightColor() * getColor() * getCoefficients().dDiffuseCoefficient * dCosTheta;

            // Specular Color
            color = color + vSpotLights[iSpotLightCounter]->getLightColor() * getColor() * getCoefficients().dSpecularCoefficient * pow(dCosPhi, getShine());
        }

        // Recursive Reflection
        if (level >= iRecursionDepth)
        {
            return dRoot;
        }

        // Construct reflected ray from intersection point
        // Moving the intersection point slightly away from the surface
        Ray rReflectedRay(IntersectingPoint + vReflectedRay, vReflectedRay);

        // Find tmin from the nearest intersecting object
        int iNearestObject = -1;
        double dNearestObjectDistance = INFINITE;

        for (int iObjectCounter = 0; iObjectCounter < vObjects.size(); iObjectCounter++)
        {
            Color cDummyColor;
            double dTemp = vObjects[iObjectCounter]->Intersect(&rReflectedRay, cDummyColor, 0);

            if (dTemp > 0 && dTemp < dNearestObjectDistance)
            {
                dNearestObjectDistance = dTemp;
                iNearestObject = iObjectCounter;
            }
        }

        Color cReflectedColor;

        // If found, call intersect(rreflected, colorreflected, level+1)
        if (iNearestObject != -1)
        {
            vObjects[iNearestObject]->Intersect(&rReflectedRay, cReflectedColor, level + 1);
        }

        // update color using the impact of reflection
        color = color + cReflectedColor * getCoefficients().dReflectionCoefficient;
        color.clip();

        return dRoot;
    }
};

/**************************************Floor Class*********************************************/
class Floor : public Object
{
    double dFloorWidth, dTileWidth;
    int iTileCount;
    Vector vBottomLeft;

public:
    // Constructors
    Floor()
    {
        dFloorWidth = 0;
        dTileWidth = 0;
        iTileCount = 0;
        vBottomLeft = Vector(0, 0, 0);
    }
    Floor(double FloorWidth, double TileWidth, Color color, Coefficients coefficients, int Shine)
        : Object(color, coefficients, Shine)
    {
        this->dFloorWidth = FloorWidth;
        this->dTileWidth = TileWidth;
        this->iTileCount = FloorWidth / TileWidth;
        this->vBottomLeft = Vector(-FloorWidth / 2, -FloorWidth / 2, 0);
    }

    // Implementation of draw function
    void draw()
    {
        // Draw squares on xy plane alternating colors
        glBegin(GL_QUADS);
        {
            for (int i = 0; i < iTileCount; i++)
            {
                for (int j = 0; j < iTileCount; j++)
                {
                    glColor3f((i + j) % 2, (i + j) % 2, (i + j) % 2);
                    glVertex3f(vBottomLeft.dX + i * dTileWidth, vBottomLeft.dY + j * dTileWidth, vBottomLeft.dZ);
                    glVertex3f(vBottomLeft.dX + i * dTileWidth, vBottomLeft.dY + (j + 1) * dTileWidth, vBottomLeft.dZ);
                    glVertex3f(vBottomLeft.dX + (i + 1) * dTileWidth, vBottomLeft.dY + (j + 1) * dTileWidth, vBottomLeft.dZ);
                    glVertex3f(vBottomLeft.dX + (i + 1) * dTileWidth, vBottomLeft.dY + j * dTileWidth, vBottomLeft.dZ);
                }
            }
        }
    }

    // Implementation of intersect function
    double Intersect(Ray *r, Color &color, int level)
    {
        double dRoot = -1.0;
        // Normal = (0,0,1)
        Vector vNormal = Vector(0, 0, 1);
        vNormal = (vCamera * vNormal > 0.0) ? vNormal : -vNormal;

        // Get Intersecting Point:Ray-Plane Intersection
        // Ray will be parallel to the plane if the ray direction is parallel to the plane normal
        // In this case, the ray will not intersect the plane
        double dDotProduct = r->getDirection() * vNormal;

        if (dDotProduct != 0)
        {
            // Get Intersecting Point: t = -(D + n·R o) / n·Rd
            dRoot = -(vNormal * r->getOrigin()) / dDotProduct;
        }

        //  Checking if the intersection point lies within the span of the floor
        if (dRoot > 0)
        {
            // Get Intersecting Point: P = R o + t·Rd
            Vector vIntersectionPoint = r->getOrigin() + r->getDirection() * dRoot;
            // Check if the intersection point lies within the span of the floor
            // The floor had range of [-FloorWidth/2, FloorWidth/2], [-FloorWidth/2, FloorWidth/2] in x and y
            // So, if the intersection point is outside the span, it did not intersect the plane
            if ((vIntersectionPoint.dX < -dFloorWidth / 2.0 || vIntersectionPoint.dX > dFloorWidth / 2 || vIntersectionPoint.dY < -dFloorWidth / 2 || vIntersectionPoint.dY > dFloorWidth / 2))
            {
                dRoot = -1.0;
            }
        }

        // Level Check
        if (level <= 0)
        {
            return dRoot;
        }

        // Floor has two types of colors: white and black
        // Check which type of color the intersection point is on
        Vector vIntersectionPoint = r->getOrigin() + r->getDirection() * dRoot;
        Vector vFromBottomLeft = vIntersectionPoint - vBottomLeft;
        int iX = vFromBottomLeft.dX / dTileWidth;
        int iY = vFromBottomLeft.dY / dTileWidth;
        Color cIntersectionColor = Color((iX + iY) % 2, (iX + iY) % 2, (iX + iY) % 2);

        // Illumination
        // Ambient Color
        color = cIntersectionColor * getCoefficients().dAmbientCoefficient;

        // Get Intersecting Point, Normal and Reflected Ray
        Vector IntersectingPoint = r->getOrigin() + r->getDirection() * dRoot;
        Vector vReflectedRay = GetReflectedRay(r->getDirection(), vNormal).normalize();

        // Diffuse and Specular Color
        for (int iLightCounter = 0; iLightCounter < vPointLights.size(); iLightCounter++)
        {
            // Get Light Ray
            Ray rLightRay(vPointLights[iLightCounter]->getLightPosition(), IntersectingPoint - vPointLights[iLightCounter]->getLightPosition());

            Vector L = rLightRay.getDirection();
            Vector N = vNormal;
            Vector R = (GetReflectedRay(L, N) * -1).normalize();
            Vector V = (IntersectingPoint - r->getOrigin()).normalize();

            // Check if object is in shadow
            double dClosest = INFINITE;
            for (int iObjectCounter = 0; iObjectCounter < vObjects.size(); iObjectCounter++)
            {
                Color cDummyColor;
                double dTemp = vObjects[iObjectCounter]->Intersect(&rLightRay, cDummyColor, 0);

                if (dTemp > 0 && dTemp < dClosest)
                {
                    dClosest = dTemp;
                }
            }

            Vector vClosestPoint = rLightRay.getOrigin() + rLightRay.getDirection() * dClosest;
            if (GetDistance(vClosestPoint, rLightRay.getOrigin()) < GetDistance(rLightRay.getOrigin(), IntersectingPoint) - MARGIN)
            {
                // Object is obscured by another object
                continue;
            }

            // Not in shadow
            double dCosTheta = std::max(0.0, (L * (-1.0)) * N);
            double dCosPhi = std::max(0.0, (R * (-1.0)) * V);

            // Diffuse Color
            color = color + vPointLights[iLightCounter]->getLightColor() * cIntersectionColor * getCoefficients().dDiffuseCoefficient * dCosTheta;

            // Specular Color
            color = color + vPointLights[iLightCounter]->getLightColor() * cIntersectionColor * getCoefficients().dSpecularCoefficient * pow(dCosPhi, getShine());
        }

        // For Spot Lights
        for(int iSpotLightCounter = 0; iSpotLightCounter < vSpotLights.size(); iSpotLightCounter++)
        {
            // Get Light Ray
            Ray rSpotLightRay(vSpotLights[iSpotLightCounter]->getLightPosition(), IntersectingPoint - vSpotLights[iSpotLightCounter]->getLightPosition());

            Vector L = rSpotLightRay.getDirection();
            Vector N = vNormal;
            Vector R = (GetReflectedRay(L, N) * -1).normalize();
            Vector V = (IntersectingPoint - r->getOrigin()).normalize();

            // Check if intersection point is in cutoff angle
            // Get angle between light direction and ray direction
            double dAngleFromSource = GetAngleBetweenVectors(vSpotLights[iSpotLightCounter]->getLightDirection(), rSpotLightRay.getDirection());
            // Since Spot Light sheds light in a cone, both positive and negative angles are allowed
            dAngleFromSource = (dAngleFromSource < 0) ? dAngleFromSource * -1.0 : dAngleFromSource;

            // If the angle is above cutoff angle, this spotlight won't affect the object
            if(dAngleFromSource > vSpotLights[iSpotLightCounter]->getCutoffAngle())
            {
                continue;
            }

            // Check if object is in shadow
            double dClosest = INFINITE;
            for (int iObjectCounter = 0; iObjectCounter < vObjects.size(); iObjectCounter++)
            {
                Color cDummyColor;
                double dTemp = vObjects[iObjectCounter]->Intersect(&rSpotLightRay, cDummyColor, 0);

                if (dTemp > 0 && dTemp < dClosest)
                {
                    dClosest = dTemp;
                }
            }

            Vector vClosestPoint = rSpotLightRay.getOrigin() + rSpotLightRay.getDirection() * dClosest;
            if (GetDistance(vClosestPoint, rSpotLightRay.getOrigin()) < GetDistance(rSpotLightRay.getOrigin(), IntersectingPoint) - MARGIN)
            {
                // Object is obscured by another object
                continue;
            }

            // Not in shadow
            double dCosTheta = std::max(0.0, (L * (-1.0)) * N);
            double dCosPhi = std::max(0.0, (R * (-1.0)) * V);

            // Diffuse Color
            color = color + vSpotLights[iSpotLightCounter]->getLightColor() * cIntersectionColor * getCoefficients().dDiffuseCoefficient * dCosTheta;

            // Specular Color
            color = color + vSpotLights[iSpotLightCounter]->getLightColor() * cIntersectionColor * getCoefficients().dSpecularCoefficient * pow(dCosPhi, getShine());
        }

        // Recursive Reflection
        if (level >= iRecursionDepth)
        {
            return dRoot;
        }

        // Construct reflected ray from intersection point
        // Moving the intersection point slightly away from the surface
        Ray rReflectedRay(IntersectingPoint + vReflectedRay, vReflectedRay);

        // Find tmin from the nearest intersecting object
        int iNearestObject = -1;
        double dNearestObjectDistance = INFINITE;

        for (int iObjectCounter = 0; iObjectCounter < vObjects.size(); iObjectCounter++)
        {
            Color cDummyColor;
            double dTemp = vObjects[iObjectCounter]->Intersect(&rReflectedRay, cDummyColor, 0);

            if (dTemp > 0 && dTemp < dNearestObjectDistance)
            {
                dNearestObjectDistance = dTemp;
                iNearestObject = iObjectCounter;
            }
        }

        Color cReflectedColor;

        // If found, call intersect(rreflected, colorreflected, level+1)
        if (iNearestObject != -1)
        {
            vObjects[iNearestObject]->Intersect(&rReflectedRay, cReflectedColor, level + 1);
        }

        // update color using the impact of reflection
        color = color + cReflectedColor * getCoefficients().dReflectionCoefficient;
        color.clip();

        return dRoot;
    }
};

// Creator Functions
// Objects
Object *CreateSphere(std::istream &is)
{
    Vector vCenter;
    double dRadius;
    Color cColor;
    Coefficients coCoefficients;
    int iShine;
    is >> vCenter >> dRadius >> cColor >> coCoefficients >> iShine;
    return new Sphere(vCenter, dRadius, cColor, coCoefficients, iShine);
}

Object *CreateTriangle(std::istream &is)
{
    Vector vPointA, vPointB, vPointC;
    Color cColor;
    Coefficients coCoefficients;
    int iShine;
    is >> vPointA >> vPointB >> vPointC >> cColor >> coCoefficients >> iShine;
    return new Triangle(vPointA, vPointB, vPointC, cColor, coCoefficients, iShine);
}

Object *CreateQuadric(std::istream &is)
{
    double dA, dB, dC, dD, dE, dF, dG, dH, dI, dJ;
    Vector vReferencePoint;
    double dLength, dWidth, dHeight;
    Color cColor;
    Coefficients coCoefficients;
    int iShine;
    is >> dA >> dB >> dC >> dD >> dE >> dF >> dG >> dH >> dI >> dJ >> vReferencePoint >> dLength >> dWidth >> dHeight >> cColor >> coCoefficients >> iShine;
    return new Quadric(dA, dB, dC, dD, dE, dF, dG, dH, dI, dJ, vReferencePoint, dLength, dWidth, dHeight, cColor, coCoefficients, iShine);
}

// Lights
PointLight *CreatePointLight(std::istream &is)
{
    Vector vLightPosition;
    Color cLightColor;
    is >> vLightPosition >> cLightColor;
    return new PointLight(vLightPosition, cLightColor);
}

SpotLight *CreateSpotLight(std::istream &is)
{
    Vector vLightPosition;
    Color cLightColor;
    Vector vLightDirection;
    double dLightAngle;
    is >> vLightPosition >> cLightColor >> vLightDirection >> dLightAngle; // dLightAngle is in degrees, convert to radians
    dLightAngle = (dLightAngle < 0) ? -dLightAngle : dLightAngle;   // Save only positive value for the sake of calculation
    return new SpotLight(vLightPosition, cLightColor, vLightDirection, RADIAN(dLightAngle));
}

#endif // !OBJECT_HPP_INCLUDED