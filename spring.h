#include <particlesystem.h>
#include <vector_math.h>

#include <vector>

using namespace std;

class Spring {
   private:
    ParticleSystem sys;
    float kd;
    float ks;
    float r;
    Particle p1;
    Particle p2;

   public:
    Spring(ParticleSystem sys, Particle p1, Particle p2, float ks, float kd) {
    }

    void apply_force() {
        vector<float> deltax = vector_math::subtract(p1.x, p2.x);
        vector<float> deltav = vector_math::subtract(p1.v, p2.v);
        vector<float> unitx = vector_math::scalar_mult(deltax, 1 / vector_math::length(deltax));
        vector<float> f1 = vector_math::scalar_mult(unitx, -(ks * (vector_math::length(deltax) - r) + kd * (vector_math::dot(deltav, deltax), 1 / vector_math::length(deltax))));
        p1.f = vector_math::add(p1.f, f1);
        p2.f = vector_math::subtract(p2.f, f1);
    }
};