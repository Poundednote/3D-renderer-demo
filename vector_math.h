#include <stdexcept>>
#include <vector>
using namespace std;

namespace vector_math {

/* takes two vectors both size n and adds them and returns new resutlting vector
    add([x1, x2...xn-1, xn], [Y1, Y2...yn-1, Yn])
    returns [x1+y1, x2+y2...Xn-1+Yn-1, Xn+Yn]
*/
vector<float> add(vector<float> v1, vector<float> v2) {
    if (!v2.size() == v2.size()) {
        throw invalid_argument("vectors have to be the same size");
    } else {
        vector<float> v3;
        for (int i = 0; i < v1.size(); i++) {
            v3.push_back(v1[i] + v2[i]);
        }
        return v3;
    }
}

/* takes vector [X1, X2...Xn-1, Xn] and multiplies by scalar value such that the
   result is [SX1, SX2,SXn-1, SXn]
    */

vector<float> scalar_mult(vector<float> v1, float s) {
    vector<float> v2;
    for (int i = 0; i < v1.size(); i++) {
        v2.push_back(v1[i] * s);
    }
    return v2;
}

}  // namespace vector_math