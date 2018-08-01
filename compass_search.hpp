using namespace std;
typedef pair<char,ulong> pCharLong;

double *compass_search (
    double function_handle ( int m, double x[], const vector<pCharLong> &v, const string &seq, const vector<double> &qua, double lower[], double upper[], double sumBound ),
    int m, double x0[],
    const vector<pCharLong> &v,
    const string &seq, const vector<double> &qua,
    double lower[], double upper[], double sumBound,
    double delta_tol, double delta_init, int k_max, double &fx, int &k );


double r8_abs ( double x );
void r8vec_copy ( int n, double a1[], double a2[] );
void r8vec_print ( int n, double a[], string title );
void timestamp ( );
