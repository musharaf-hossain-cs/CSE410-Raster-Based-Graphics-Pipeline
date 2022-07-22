#include<bits/stdc++.h>
#include "bitmap_image.hpp"
using namespace std;

#define dbg(x) cout<<#x<<": "<<x<<endl;

const double PI = 3.14159265359;

struct Point{
    double x, y, z;
    Point(double x=0, double y=0, double z=0){
        this->x = x;
        this->y = y;
        this->z = z;
    }

    void Normalize(){
        double value = sqrt(x*x + y*y + z*z);
        x = x/value;
        y = y/value;
        z = z/value;
    }

    // for debug
    void Print(){
        dbg(x);
        dbg(y);
        dbg(z);
    }
};

struct Color{
    int r,g,b;
    Color(){
        r=g=b=255;
    }
    Color(int R, int G, int B){
        r = R;
        g = G;
        b = B;
    }
};

struct Column{
    double xs;
    double zs;
    Column(){
        xs=0;
        zs=0;
    }
    Column(double x, double z){
        xs = x;
        zs = z;
    }
};

struct Triangle{
    Point points[3];
    Color color;

    Triangle(){

    }
    Triangle(Point A, Point B, Point C){
        points[0] = A;
        points[1] = B;
        points[2] = C;
    }
    void setPoint(Point A, Point B, Point C){
        points[0] = A;
        points[1] = B;
        points[2] = C;
    }
    void setColor(int r, int g, int b){
        color.r = r;
        color.g = g;
        color.b = b;
    }

    double getMinX(){
        return min(min(points[0].x, points[1].x), points[2].x);
    }

    double getMaxX(){
        return max(max(points[0].x, points[1].x), points[2].x);
    }

    double getMinY(){
        return min(min(points[0].y, points[1].y), points[2].y);
    }

    double getMaxY(){
        return max(max(points[0].y, points[1].y), points[2].y);
    }

    void Print(){
        cout<<"Triangle:"<<endl;
        cout<<"\tP1: "<<fixed<<setprecision(7)<<points[0].x<<' '<<points[0].y<<' '<<points[0].z<<endl;
        cout<<"\tP2: "<<fixed<<setprecision(7)<<points[1].x<<' '<<points[1].y<<' '<<points[1].z<<endl;
        cout<<"\tP3: "<<fixed<<setprecision(7)<<points[2].x<<' '<<points[2].y<<' '<<points[2].z<<endl;
        cout<<"\tColor: rgb("<<color.r<<','<<color.g<<','<<color.b<<")"<<endl;
    }
};


struct Matrix{
    double data[4][4];

    Matrix(){
        for(int i=0; i<4; i++){
            for(int j=0;j<4;j++){
                data[i][j] = 0;
            }
        }
    }

    void MakeIdentityMatrix(){
        for(int i=0; i<4; i++){
            data[i][i] = 1;
        }
        return;
    }

    void Print(){
        for(int i=0;i<4;i++){
            for(int j=0;j<4;j++){
                cout << fixed << setprecision(7) << data[i][j] << ' ';
            }
            cout<<endl;
        }
        cout<<endl;
        return;
    }

    void Print(ofstream &file){
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                file << fixed << setprecision(7) << data[j][i]/data[3][i] << ' ';
            }
            file<<endl;
        }
        file<<endl;
        return;
    }

    Triangle getTriangle(){
        Point p[3];
        for(int i=0;i<3;i++){
            p[i].x = data[0][i] / data[3][i];
            p[i].y = data[1][i] / data[3][i];
            p[i].z = data[2][i] / data[3][i];
        }

        Triangle t(p[0], p[1], p[2]);
        return t;
    }

    Matrix Transpose(){
        Matrix temp;
        for(int i=0; i<4;i++){
            for(int j=0;j<4;j++){
                temp.data[i][j] = data[j][i];
            }
        }
        return temp;
    }
};

Matrix Multiply(Matrix a, Matrix b){
    Matrix temp;
    double entry;
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            entry = 0;
            for(int k=0; k<4; k++){
                entry += a.data[i][k] * b.data[k][j];
            }
            temp.data[i][j] = entry;
            if(abs(entry) < 1e-8) temp.data[i][j] = 0;
        }
    }
    return temp;
}

Point VectorCross(Point a, Point b){
    Point temp;
    temp.x = a.y*b.z - a.z*b.y;
    temp.y = a.z*b.x - a.x*b.z;
    temp.z = a.x*b.y - a.y*b.x;
    return temp;
}


/// **************************** Global variables ***************************************
Point eye, look, up, L, U, R;
double fovY, aspectRatio, near, far;
Matrix P, V;
Matrix transformation;
vector<Triangle> triangles;
stack<Matrix> stk;

void PopMatrix(){
    if(stk.empty()){
        cout<<"Error: The stack is empty. Pop is not possible."<<endl;
        return;
    }
    transformation = stk.top();
    stk.pop();
}

void PushMatrix(){
    stk.push(transformation);
}

void SceneInput(ifstream &scene){
    scene >> eye.x >> eye.y >> eye.z;
    scene >> look.x >> look.y >> look.z;
    scene >> up.x >> up.y >> up.z;
    scene >> fovY >> aspectRatio >> near >> far;

//    eye.Print();
//    look.Print();
//    up.Print();
//    dbg(fovY);
//    dbg(aspectRatio);
//    dbg(near);
//    dbg(far);
}

void Translate(double tx, double ty, double tz){
    Matrix trans;
    trans.MakeIdentityMatrix();
    trans.data[0][3] = tx;
    trans.data[1][3] = ty;
    trans.data[2][3] = tz;
//    cout<<"before"<<endl;
//    transformation.Print();
    transformation = Multiply(transformation, trans);

//    cout<<"after"<<endl;
//    transformation.Print();
}

void Scale(double sx, double sy, double sz){
    Matrix scale;

    scale.data[0][0] = sx;
    scale.data[1][1] = sy;
    scale.data[2][2] = sz;
    scale.data[3][3] = 1;

    transformation = Multiply(transformation, scale);
}

Point RodriguesRotation(Point v, Point k, double angle ){
    ///Ref: https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
    /// Vrot = v.cosA + (k x v) sinA + k(k.v)(1-cosA)

    double sinA = sin(angle * PI / 180);
    double cosA = cos(angle * PI / 180);

    //first term
    Point A(v.x*cosA, v.y*cosA, v.z*cosA);

    //second term
    Point B = VectorCross(k, v);
    B.x = B.x*sinA;
    B.y = B.y*sinA;
    B.z = B.z*sinA;

    //third term
    double dot = (1-cosA) * (v.x*k.x + v.y*k.y + v.z*k.z);
    Point C(k.x * dot, k.y * dot, k.z * dot);

    Point Vrot(A.x+B.x+C.x, A.y+B.y+C.y, A.z+B.z+C.z);
    return Vrot;
}

void Rotate(double angle, Point a){
    a.Normalize();
    Point c1 = RodriguesRotation(Point(1,0,0), a, angle);
    Point c2 = RodriguesRotation(Point(0,1,0), a, angle);
    Point c3 = RodriguesRotation(Point(0,0,1), a, angle);

    Matrix rot;
    rot.data[0][0] = c1.x;
    rot.data[1][0] = c1.y;
    rot.data[2][0] = c1.z;

    rot.data[0][1] = c2.x;
    rot.data[1][1] = c2.y;
    rot.data[2][1] = c2.z;

    rot.data[0][2] = c3.x;
    rot.data[1][2] = c3.y;
    rot.data[2][2] = c3.z;

    rot.data[3][3] = 1;

    transformation = Multiply(transformation, rot);

}

/// not needed now
void Rotate1(double angle, Point a){
    /// Ref: https://www.youtube.com/watch?v=OhgiPknf2mc
    // reverse translation is wrong in the video.

    a.Normalize();

    Matrix rot;
    double sinA = sin(angle * PI / 180);
    double cosA = cos(angle * PI / 180);

    //move rotation axis to center
    Matrix translate2Center;
    translate2Center.MakeIdentityMatrix();
    translate2Center.data[0][3] = -a.x;
    translate2Center.data[1][3] = -a.y;
    translate2Center.data[2][3] = -a.z;

    // Rx, Ry to make the axis align with z-axis
    Matrix Rx;
    double d = sqrt(a.y*a.y + a.z*a.z);
    Rx.data[0][0] = Rx.data[3][3] = 1;
    Rx.data[1][1] = a.z / d;
    Rx.data[1][2] = -a.y / d;
    Rx.data[2][1] = a.y / d;
    Rx.data[2][2] = a.z / d;

    Matrix Ry;
    Ry.data[1][1] = Ry.data[3][3] = 1;
    Ry.data[0][0] = d;
    Ry.data[0][2] = -a.x;
    Ry.data[2][0] = a.x;
    Ry.data[2][2] = d;

    // Now axis is along with z axis
    // rotate along z-axis
    Matrix Rz;
    Rz.data[0][0] = cosA;
    Rz.data[0][1] = -sinA;
    Rz.data[1][0] = sinA;
    Rz.data[1][1] = cosA;
    Rz.data[2][2] = Rz.data[3][3] = 1;

    // rotate axis its original direction by Rx-inverse and Ry-inverse
    Matrix AntiRy = Ry.Transpose();
    Matrix AntiRx = Rx.Transpose();

    // Now move axis its original position
    Matrix translate2Original;
    translate2Original.MakeIdentityMatrix();
    translate2Original.data[0][3] = a.x;
    translate2Original.data[1][3] = a.y;
    translate2Original.data[2][3] = a.z;

    // Now multiply all of them to find rotation matrix
    Matrix rotation;
    rotation.MakeIdentityMatrix();

    rotation = Multiply(rotation, translate2Center);
    rotation = Multiply(rotation, Rx);
    rotation = Multiply(rotation, Ry);
    rotation = Multiply(rotation, Rz);
    rotation = Multiply(rotation, AntiRy);
    rotation = Multiply(rotation, AntiRx);
    rotation = Multiply(rotation, translate2Original);

    transformation = Multiply(transformation, rotation);
}



///*********************************Z-buffer********************************
double leftLimit, rightLimit, topLimit, bottomLimit;
double leftX, rightX;
double topY, bottomY;
double dx, dy;
int scrHeight, scrWidth;
double front, rear;

// value of FindLeftRightScanline...
double z1, z2;

void inputConfig(){
    ifstream in;
    in.open("config.txt");
    in >> scrWidth >> scrHeight;
    in >> leftLimit;
    in >> bottomLimit;
    in >> front >> rear;

    rightLimit = -leftLimit;
    topLimit = -bottomLimit;
}

void setColors(){
    int s = triangles.size();
    int r,g,b;
    for(int i=0; i<s; i++){
        r = rand() % 256;
        g = rand() % 256;
        b = rand() % 256;
        triangles[i].setColor(r,g,b);
    }
}

void printTriangles(){
    for(Triangle t: triangles){
        t.Print();
    }
}

vector<Column> GetColumn(Triangle t, double ys){
    // Ref: given z-buffer calculation
    vector<Column> columns;
    Point a, b;
    double x, z;
    int s;
    int flag1=-1, flag2=-1;
    double minX = t.getMinX();
    double maxX = t.getMaxX();

    a = t.points[0];
    b = t.points[1];
    x = a.x + ((ys-a.y)/(b.y-a.y))*(b.x-a.x);
    z = a.z + ((ys-a.y)/(b.y-a.y))*(b.z-a.z);
    s = (int) round(abs(leftX - x)/dx);
    if( x >= minX && x <= maxX){
        flag1 = s;
        columns.push_back(Column(x,z));
    }

    a = t.points[1];
    b = t.points[2];
    x = a.x + ((ys-a.y)/(b.y-a.y))*(b.x-a.x);
    z = a.z + ((ys-a.y)/(b.y-a.y))*(b.z-a.z);
    s = (int) round(abs(leftX - x)/dx);
    if(x >= minX && x <= maxX && s != flag1){
        columns.push_back(Column(x,z));
        flag2 = s;
    }

    a = t.points[2];
    b = t.points[0];
    x = a.x + ((ys-a.y)/(b.y-a.y))*(b.x-a.x);
    z = a.z + ((ys-a.y)/(b.y-a.y))*(b.z-a.z);
    s = (int) round(abs(leftX - x)/dx);
    if(x >= minX && x <= maxX && s != flag1 && s!= flag2){
        columns.push_back(Column(x,z));
    }

    return columns;
}

void initiateZbuffer(){
    dx = (rightLimit - leftLimit) / (double)scrWidth;
    dy = (topLimit - bottomLimit) / (double)scrHeight;
    topY = topLimit - dy/2;
    bottomY = bottomLimit + dy/2;
    leftX = leftLimit + dx/2;
    rightX = rightLimit - dx/2;
}

pair<int,int> FindLeftRightScanline(Triangle t, double ys){
    vector<Column> columns = GetColumn(t, ys);
    double minX, maxX;
    if(columns.empty()) return make_pair(-1,-1);
    else if(columns.size() == 1){
        minX = columns[0].xs;
        maxX = columns[0].xs;
        z1 = columns[0].zs;
        z2 = columns[0].zs;
    }
    else if(columns[0].xs <= columns[1].xs){
        minX = columns[0].xs;
        maxX = columns[1].xs;
        z1 = columns[0].zs;
        z2 = columns[1].zs;
    }
    else{
        minX = columns[1].xs;
        maxX = columns[0].xs;
        z1 = columns[1].zs;
        z2 = columns[0].zs;
    }
    int left, right;

    if(minX < leftX){
        left = 0;
    }
    else {
        left = (int) round(abs(leftX - minX)/dx);
    }

    if(maxX > rightX){
        right = scrWidth - 1;
    }
    else{
        right = (int) round(abs(leftX-maxX)/dx);
    }

    return make_pair(left,right);
}

pair<int,int> FindTopBottomScanline(Triangle t){
    double maxY = t.getMaxY();
    double minY = t.getMinY();
    int topScanline, bottomScanline;

    if(maxY > topY){
        topScanline = 0;
    }
    else{
        topScanline = (int) round(abs(topY - maxY)/dy);
    }

    if(minY < bottomY){
        bottomScanline = scrHeight - 1;
    }
    else{
        bottomScanline = (int) round(abs(topY-minY)/dy);
    }
    return make_pair(topScanline, bottomScanline);
}

void zBufferAlgo(){
    vector<vector<double>> zBuffer;
    bitmap_image image(scrWidth, scrHeight);

    // background color
    for(int i=0; i<scrHeight; i++){
        vector<double> row;
        for(int j=0; j<scrWidth; j++){
            row.push_back(rear);
            image.set_pixel(i,j,0,0,0);
        }
        zBuffer.push_back(row);
    }

    /// Spec Pseudocode

    // for each object : Triangles
    for(int i=0; i<triangles.size(); i++){
        // Find top_scanline and bottom_scanline after necessary clipping
        Triangle t = triangles[i];
        double minX = t.getMinX();
        double maxX = t.getMaxX();
        pair<int,int> topBottom = FindTopBottomScanline(t);
        int topScanline = topBottom.first;
        int bottomScanline = topBottom.second;

        // for row_no from top_scanline to bottom_scanline
        for(int row = topScanline; row <= bottomScanline; row++){
            // Find left_intersecting_column and right_intersecting_column
            // after necessary clipping
            double ys = topY - row*dy;
            pair<int,int> leftRight = FindLeftRightScanline(t,ys);
            if(leftRight.first < 0) continue;
            int leftScanline = leftRight.first;
            int rightScanline = leftRight.second;

            double increment = (z2-z1)/(rightScanline-leftScanline);
            if(increment < 0) increment = 0;

            if(maxX - minX > 0){
                // for col_no from left_intersecting_column to right_intersecting_column
                double z=z1;
                for(int col=leftScanline; col <= rightScanline; col++, z+=increment){
                    // z < zBuffer[row][col] --> this checks for multiple object. nearest one chosen
                    // Compare with Z-buffer and z_front_limit and update if required
                    if(z >= front && z <= rear && z < zBuffer[row][col]){
                        // z values is calculated in GetColumn function
                        zBuffer[row][col] = z;
                        // Update pixel information if required
                        image.set_pixel(col, row, t.color.r, t.color.g, t.color.b);
                    }
                }
            }
        }
    }

    image.save_image("out.bmp");
    ofstream bufferOut;
    bufferOut.open("z-buffer.txt");

    for(int i=0;i<zBuffer.size();i++){
        for(int j=0; j<zBuffer[i].size(); j++){
            if(zBuffer[i][j] < rear && zBuffer[i][j] >= front){
                bufferOut << fixed <<setprecision(6) << zBuffer[i][j] << "\t";
            }
        }
        bufferOut << endl;
    }

    bufferOut.close();

    // free memory
    for(int i=0; i<zBuffer.size(); i++){
        zBuffer[i].clear();
    }
    zBuffer.clear();
    image.clear();

}

void zBufferUtil(){
    inputConfig();
    initiateZbuffer();
    setColors();
    zBufferAlgo();

//    printTriangles();
}



void Initiate(){
    // calculate L
    L.x = look.x - eye.x;
    L.y = look.y - eye.y;
    L.z = look.z - eye.z;
    L.Normalize();

    // Calculate R
    R = VectorCross(L, up);
    R.Normalize();

    // Calculate U
    U = VectorCross(R, L);
    U.Normalize();


    /// Ref: Projection and view transformation.ppt : Slide-49,50
    Matrix Vtrans;
    Vtrans.MakeIdentityMatrix();
    Vtrans.data[0][3] = -eye.x;
    Vtrans.data[1][3] = -eye.y;
    Vtrans.data[2][3] = -eye.z;

    Matrix Vrotate;
    Vrotate.data[0][0] = R.x;
    Vrotate.data[0][1] = R.y;
    Vrotate.data[0][2] = R.z;

    Vrotate.data[1][0] = U.x;
    Vrotate.data[1][1] = U.y;
    Vrotate.data[1][2] = U.z;

    Vrotate.data[2][0] = -L.x;
    Vrotate.data[2][1] = -L.y;
    Vrotate.data[2][2] = -L.z;

    Vrotate.data[3][3] = 1.0;

    V = Multiply(Vrotate, Vtrans);

    ///Ref: https://unspecified.wordpress.com/2012/06/21/calculating-the-gluperspective-matrix-and-other-opengl-matrix-maths/

    P.data[0][0] = 1.0 / tan(fovY*aspectRatio*PI/360.0); //180*2=360
    P.data[1][1] = 1.0 / tan(fovY*PI/360.0);
    P.data[2][2] = (far + near)/(near - far);
    P.data[3][2] = -1.0;
    P.data[2][3] = 2 * far * near / (near - far);

    // make transformation matrix identity matrix
    transformation.MakeIdentityMatrix();

}


int main(){
    ifstream scene;
    ofstream stage1, stage2, stage3;
    srand(time(0));

    scene.open("scene.txt");
    stage1.open("stage1.txt");
    stage2.open("stage2.txt");
    stage3.open("stage3.txt");

    if(!scene || !stage1 || !stage2 || !stage3){
        cout<<"Files are not opened successfully"<<endl;
    }


    SceneInput(scene);
    Initiate();

    string command;

    while(true){
        scene >> command;
        if(command == "triangle"){
            Matrix triangle;
            scene >> triangle.data[0][0] >> triangle.data[1][0] >> triangle.data[2][0];
            scene >> triangle.data[0][1] >> triangle.data[1][1] >> triangle.data[2][1];
            scene >> triangle.data[0][2] >> triangle.data[1][2] >> triangle.data[2][2];

            triangle.data[0][3] = triangle.data[1][3] = triangle.data[2][3] = 1.0;
            triangle.data[3][0] = triangle.data[3][1] = triangle.data[3][2] = 1.0;

            // stage 1
            triangle = Multiply(transformation, triangle);
            triangle.Print(stage1);

            // stage 2
            triangle = Multiply(V, triangle);
            triangle.Print(stage2);

            // stage 2
            triangle = Multiply(P, triangle);
            triangles.push_back(triangle.getTriangle());
            triangle.Print(stage3);
        }
        else if(command == "translate"){
            double tx, ty, tz;
            scene >> tx >> ty >> tz;
            Translate(tx, ty, tz);
        }
        else if(command == "scale"){
            double sx, sy, sz;
            scene >> sx >> sy >> sz;
            Scale(sx, sy, sz);
        }
        else if(command == "rotate"){
            double angle;
            Point a;
            scene >> angle >> a.x >> a.y >> a.z;
            Rotate(angle, a);
        }
        else if(command == "push"){
            PushMatrix();
        }
        else if(command == "pop"){
            PopMatrix();
        }

        else if(command == "end"){
            break;
        }
    }
    stage1.close();
    stage2.close();
    stage3.close();

    zBufferUtil();
}
