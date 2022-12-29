#include <vector_math.h>

#include <vector>

using namespace std;

class Particle {
   public:
    float m;          // mass
    vector<float> x;  // position vector
    vector<float> v;  // velocity vector
    vector<float> f;  // force accumulator
};

class ParticleSystem {
   public:
    vector<Particle> p;  // list of particles
    float clock;         // simulation clock

    // returns 6N vector of particle phase space
    vector<float> get_state() {
        vector<float> phase_space;
        for (int i = 0; i < p.size(); i++) {
            phase_space.push_back(p[i].x[0]);
            phase_space.push_back(p[i].x[1]);
            phase_space.push_back(p[i].x[2]);
            phase_space.push_back(p[i].v[0]);
            phase_space.push_back(p[i].v[1]);
            phase_space.push_back(p[i].v[2]);
        }
        return phase_space;
    }

    void set_state(vector<vector<float>> src) {
        for (int i = 0; i < p.size(); i++) {
            p[i].v[0] = src[i][0];
            p[i].v[0] = src[i][1];
            p[i].v[0] = src[i][2];
            p[i].v[0] = src[i][3];
            p[i].v[0] = src[i][4];
            p[i].v[0] = src[i][5];
        }
    }
    vector<float> derivative() {
        vector<float> phase_space_derivative;
        for (int i = 0; i < p.size(); i++) {
            phase_space_derivative.push_back(p[i].v[0]);
            phase_space_derivative.push_back(p[i].v[1]);
            phase_space_derivative.push_back(p[i].v[2]);
            phase_space_derivative.push_back(p[i].f[0] / p[i].m);
            phase_space_derivative.push_back(p[i].f[1] / p[i].m);
            phase_space_derivative.push_back(p[i].f[2] / p[i].m);
        }
        return phase_space_derivative;
    }
};