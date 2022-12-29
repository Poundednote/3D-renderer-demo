typedef struct
{
    float mass; /* mass */
    float *x;   /* position vector */
    float *v;   /* velocity vector */
    float *f;   /* force accumulator */
} *Particle;