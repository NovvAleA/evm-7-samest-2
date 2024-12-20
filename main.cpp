#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;
const double eps       = 0.00001;  // epsilon
const double q         = 0.8;      //
const double x0_        = -0.5;
const double y0_        = -0.85;
const double mu         = 4.8;
double MAX_ITERATIONS   = 1000;

struct vector {
    double x = 0;
    double y = 0;
public:
    vector operator - (vector v) const {
        return vector(x - v.x, y - v.y);
    }
    void operator = (vector v) {
        x = v.x;
        y = v.y;
    }
    double abs() const{
        return sqrt(pow(x,2) + pow(y,2));
    }
};

struct matrix {
    double x11 = 0;
    double x12 = 0;
    double x21 = 0;
    double x22 = 0;
public:
    vector operator * (vector v) const {
        return vector(x11*v.x + x12*v.y, x21*v.x + x22*v.y);
    }
};

//вычисление матрицы D-1
matrix D(vector v){
    matrix D;
    double det_d = fabs(-2*v.y - 2*(v.x)*cos(v.y) + 2*(v.x)*(v.y)*sin(v.y));
    D.x11 = (2*v.y)/det_d;
    D.x12 = (-cos(v.y) + (v.y)*sin(v.y));
    D.x21 = (-2*v.x)/det_d;
    D.x22 = (-1)/det_d;
    return D;
}

vector phi(vector v){
    vector phi = {0,0};
    phi.x = (v.y)*cos(v.y);
    phi.y = - sqrt(1 - pow(v.x, 2));
    return phi;
}

vector f(vector v){
    vector f = {0,0};
    f.x = (v.y)*cos(v.y) - v.x;
    f.y = pow(v.x,2) + pow(v.y,2) - 1;
    return f;
}

vector simple_iteration_method(ostream& res_out){
    int N = 0;
    vector x0(x0_,y0_), x1, temp;
    double delta_x = 1 + eps * (1 - q)/q; //delta_x - разность между n и n+1 результатами итераций
    cout << endl << "simple iteration method" << endl;
    while (delta_x >= eps * (1 - q)/q) {
        x1 = phi(x0);
        temp = x1 - x0;
        delta_x = temp.abs();
        x0 = x1;
        N++;
        cout << "x = " << x1.x << " ,y = " << x1.y << endl;
        if(N > MAX_ITERATIONS){
            cout << "ERROR N = " << MAX_ITERATIONS << " ! Max number of iterations!"<< endl;
            break;
        }
    }
    cout << "N = " << N << endl;
    res_out << "|" << "Простых итераций" << "|" << setw(16) << x1.x << "|"<< setw(16) << x1.y << "|" << setw(15) << f(x1).abs()<< "|" ;
    res_out << "  [" << setw(5) << x0_ << ", " << setw(5) << y0_ << "]"<<"  |" << setw(7) << N << "|" << setw(5) << q<< "|" << setw(5) << mu << "|" << endl;
    res_out << "+----------------+----------------+----------------+---------------+------------------+-------+-----+-----+" << endl;
    return x1;
}

vector newton_method(ostream& res_out){
    int N = 0;
    vector x0(x0_,y0_), x1, temp;
    double delta_x = 1 + eps / mu; //delta_x - разность между n и n+1 результатами итераций
    cout << endl << "Newton method" << endl;
    while (delta_x >= eps / mu) {
        x1 = x0 - D(x0)*f(x0);
        temp = x1 - x0;
        delta_x = temp.abs();
        x0 = x1;
        N++;
        cout << "x = " << x1.x << " ,y = " << x1.y << endl;
        if(N > MAX_ITERATIONS){
            cout << "ERROR N = " << MAX_ITERATIONS << " ! Max number of iterations!"<< endl;
            break;
        }
    }
    cout << "N = " << N << endl;
    res_out << "|" << "Метод Нютона    " << "|" << setw(16) << x1.x << "|"<< setw(16) << x1.y << "|" << setw(15) << f(x1).abs()<< "|" ;
    res_out << "  [" << setw(5) << x0_ << ", " << setw(5) << y0_ << "]"<<"  |" << setw(7) << N << "|" << setw(5) << q<< "|" << setw(5) << mu << "|" << endl;
    res_out << "+----------------+----------------+----------------+---------------+------------------+-------+-----+-----+" << endl;
    return x1;

}

vector zeidel_method(ostream& res_out) {
    int N = 0;
    vector x0(x0_,y0_), x1, temp;
    double delta_x = 1 + eps; //delta_x - разность между n и n+1 результатами итераций
    cout << endl << "Zeidel method" << endl;
    while (delta_x >= eps) {
        x1.x = phi(x0).x;
        x1.y = phi(x1).y;
        temp = x1 - x0;
        delta_x = temp.abs();
        x0 = x1;
        N++;
        cout << "x = " << x1.x << " ,y = " << x1.y << endl;
        if(N > MAX_ITERATIONS){
            cout << "ERROR N = " << MAX_ITERATIONS << " ! Max number of iterations!"<< endl;
            break;
        }
    }
    cout << "N = " << N << endl;
    res_out << "|" << "Метод Зейделя   " << "|" << setw(16) << x1.x << "|"<< setw(16) << x1.y << "|" << setw(15) << f(x1).abs()<< "|" ;
    res_out << "  [" << setw(5) << x0_ << ", " << setw(5) << y0_ << "]"<<"  |" << setw(7) << N << "|" << setw(5) << q<< "|" << setw(5) << mu << "|" << endl;
    res_out << "+----------------+----------------+----------------+---------------+------------------+-------+-----+-----+" << endl;
    return x1;
}

int main() {
    ofstream res_out;
    res_out.open("result.txt");
    res_out << "eps   = " << eps << endl;
    res_out << "+================+================+================+===============+==================+=======+=====+=====+" << endl;
    res_out << "|  Метод         | Корень 1       | Корень 2       | Норма невязки | Начальный вектор | N + 1 | q   | mu  |" << endl;
    res_out << "+================+================+================+===============+==================+=======+=====+=====+" << endl;
    simple_iteration_method(res_out);
    newton_method(res_out);
    zeidel_method(res_out);
    res_out.close();
    return 0;
}

