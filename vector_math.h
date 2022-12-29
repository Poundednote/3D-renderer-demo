#include <stdexcept>>
#include <vector>
using namespace std;

vector<float> add(vector<float> a, vector<float> b) {
    if (!a.size() == b.size()) {
        throw invalid_argument("vectors have to be the same size");
    } else {
        vector<float> result;
        for (int i = 0; i < a.size(); i++) {
            result.push_back(a[i] + b[i]);
        }
        return result;
    }
}