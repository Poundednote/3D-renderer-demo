#include <vector_math.h>

#include <vector>

using namespace std;

struct Particle {
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

    void set_state(vector<float> v1) {
        int i = 0;
        for (int j = 0; j < p.size(); j++) {
            p[j].x[0] = v1[i++];
            p[j].x[1] = v1[i++];
            p[j].x[2] = v1[i++];
            p[j].v[0] = v1[i++];
            p[j].v[1] = v1[i++];
            p[j].v[2] = v1[i++];
        }
    }

    vector<float> derivative() {
        vector<float> phase_space_derivative;
        // needs to zero the force accummulator first
        // then should apply all relevant forces again
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

void eulerstep(ParticleSystem sys, float deltaT) {
    vector<float> state = sys.get_state();
    vector<float> new_state = vector_math::add(state, vector_math::scalar_mult(sys.derivative(), deltaT));
    sys.set_state(new_state);
    sys.clock += deltaT;
}