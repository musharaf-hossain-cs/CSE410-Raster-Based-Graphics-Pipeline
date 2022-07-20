#include<bits/stdc++.h>
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

    void Print(ofstream &file){
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                file << fixed << setprecision(7) << data[i][j] << ' ';
            }
            file<<endl;
        }
        file<<endl;
        return;
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
            scene >> triangle.data[0][0] >> triangle.data[0][1] >> triangle.data[0][2];
            scene >> triangle.data[1][0] >> triangle.data[1][1] >> triangle.data[1][2];
            scene >> triangle.data[2][0] >> triangle.data[2][1] >> triangle.data[2][2];

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
            triangle.Print(stage3);
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

}
