#include <stdint.h>
#include <math.h>
#include <emmintrin.h>
#include <xmmintrin.h>

#include "particle_fun.h"

#define PARTICLE_MASS 1
#define COEFFICIENT_OF_DRAG 0.01f

static V3 gravity; 

static uint32_t parkmiller_rand(uint32_t *state) {
	const uint32_t A = 48271;

	uint32_t low  = (*state & 0x7fff) * A;			// max: 32,767 * 48,271 = 1,581,695,857 = 0x5e46c371
	uint32_t high = (*state >> 15)    * A;			// max: 65,535 * 48,271 = 3,163,439,985 = 0xbc8e4371
	uint32_t x = low + ((high & 0xffff) << 15) + (high >> 16);	// max: 0x5e46c371 + 0x7fff8000 + 0xbc8e = 0xde46ffff

	x = (x & 0x7fffffff) + (x >> 31);
	*state = x;
    return *state;
}

static void create_random_particles(GameState *state, 
                                    int number_of_particles, 
                                    uint32_t seed) {

    for (int i = 0; i < number_of_particles; i++) { 
        //assert(state->particle_count < arraysize(state->particles));
        state->particles.mass[i] = (float)(parkmiller_rand(&seed) % 30) + 10;
        state->particles.radius[i] = state->particles.mass[i];
        state->particles.pos[i].x = (float)(parkmiller_rand(&seed)%WORLD_WIDTH)-WORLD_RIGHT;
        state->particles.pos[i].y = (float)(parkmiller_rand(&seed)%WORLD_HEIGHT)-WORLD_TOP;
        state->particles.pos[i].z = (float)(parkmiller_rand(&seed)%WORLD_DEPTH)-WORLD_FORWARD;
#if 1
        state->particles.vel[i].x = (float)(parkmiller_rand(&seed)%20)-10;
        state->particles.vel[i].y = (float)(parkmiller_rand(&seed)%20)-10;
        state->particles.vel[i].z = 0;
#endif
        state->particles.f_accumulator[i] = {};

        state->particle_count++;
    }
}

inline static void spring_apply_force(Spring *spring) {

    V3 l = spring->p1->pos - spring->p2->pos;
    V3 dl = spring->p1->vel - spring->p2->vel;

    V3 force = -((spring->spring_const*(v3_mag(l) - spring->rest_length)) + 
            spring->damping_const*((v3_dot(l, dl))/v3_mag(l))) * (l/v3_mag(l));

    spring->p1->f_accumulator += force;

    spring->p2->f_accumulator += -force;
}

void game_update_and_render(GameMemory *memory, 
                            OffscreenBuffer *buffer, 
                            GameInput *input) {


    GameState *state = (GameState *)memory->permanent_storage;
    if (!memory->is_initialised) {
        //set GRAVITY
        gravity.y = -9.81f; 

        //init particle system
        int n_particles = 10;
        float pad = WORLD_WIDTH/(float)n_particles;
        V3 def_pos = {};
        def_pos.x = WORLD_LEFT+1; 
        def_pos.y = 0;

#if DEBUG_MODE
        create_random_particles(state, 7000, 101231);
#endif
        
#if 0
        Spring *spring = &state->springs[state->spring_count++];
        spring->p1 = &state->particles[state->particle_count++];
        spring->p2 = &state->particles[state->particle_count++];
        spring->spring_const = 3.0f;
            spring->damping_const = 0.2f;
        spring->rest_length = 3;

        spring->p1->mass = 1;
        spring->p1->radius = 1;
        spring->p1->pos.x = 0;
        spring->p1->pos.y = 0; 
        spring->p1->has_drag = true;
        spring->p1->has_grav = true;
        spring->p1->is_anchor = true;

        spring->p2->mass = 1;
        spring->p2->radius = 1;
        spring->p2->pos.x = 0 + 5;
        spring->p2->pos.y = 0;
        spring->p2->has_drag = true;
        spring->p2->has_grav = true;
        spring->p2->is_anchor = false;
#endif

        state->camera.width = 300;
        state->camera.height = 300;
        state->camera.pos.x = 0;
        state->camera.pos.y = 0;
        state->camera.pos.z = 0;
        state->camera.zfar = 100000; 
        state->camera.znear = 1; 
        state->camera.fov = PI/3.0f;

        memory->is_initialised = true;
    }

#if DEBUG_MODE
    if (state->frame_counter >= GAME_UPDATE_HZ) {
        state->frame_counter = 0;
    }
#endif


    draw_background(buffer, &state->camera);

    state->reference_cubes_v[0].vertices[0].x = _mm_set1_ps(10);
    state->reference_cubes_v[0].vertices[1].x = _mm_set1_ps(10+1);
    state->reference_cubes_v[0].vertices[0].y = _mm_set1_ps(0);
    state->reference_cubes_v[0].vertices[1].y = _mm_set1_ps(0+1);
    state->reference_cubes_v[0].vertices[0].z = _mm_set1_ps(0);
    state->reference_cubes_v[0].vertices[1].z = _mm_set1_ps(0+1);

    state->reference_cubes_v[1].vertices[0].x = _mm_set1_ps(0);
    state->reference_cubes_v[1].vertices[1].x = _mm_set1_ps(0+1);
    state->reference_cubes_v[1].vertices[0].y = _mm_set1_ps(10);
    state->reference_cubes_v[1].vertices[1].y = _mm_set1_ps(10+1);
    state->reference_cubes_v[1].vertices[0].z = _mm_set1_ps(0);
    state->reference_cubes_v[1].vertices[1].z = _mm_set1_ps(0+1);

    state->reference_cubes_v[2].vertices[0].x = _mm_set1_ps(0);
    state->reference_cubes_v[2].vertices[1].x = _mm_set1_ps(0+1);
    state->reference_cubes_v[2].vertices[0].y = _mm_set1_ps(0);
    state->reference_cubes_v[2].vertices[1].y = _mm_set1_ps(0+1);
    state->reference_cubes_v[2].vertices[0].z = _mm_set1_ps(10);
    state->reference_cubes_v[2].vertices[1].z = _mm_set1_ps(10+1);


    Vertex4_to_V2Screen((Vertex4 *)state->reference_cubes_v, 
                        &state->camera, buffer->width, buffer->height, 
                        arraysize(state->reference_cubes_v)*2, 
                        state->reference_cubes_s);

    for (int i = 0; i < arraysize(state->reference_cubes_s); i+=2) {
        for (int j = 0; j < 4; ++j) {
            V2Screen start_pos; 
            start_pos.x = (int)((float *)&state->reference_cubes_s[i].x)[j];
            start_pos.y = (int)((float *)&state->reference_cubes_s[i].y)[j];

            V2Screen end_pos;
            end_pos.x = ((int *)&state->reference_cubes_s[i+1].x)[j];
            end_pos.y = ((int *)&state->reference_cubes_s[i+1].y)[j];
            draw_line(buffer, start_pos, end_pos, 0xFFFFFFFF);
        }
    }

    float timestep = TIME_FOR_FRAME;
    {
        V3 move = {};

        if (input->camleft) {
            state->camera.theta_y -= PI/100.0f;
        }

        if (input->camright) {
            state->camera.theta_y += PI/100.0f;
        }

        if (input->camdown) {
            state->camera.theta_x += PI/100.0f;

        }

        if (input->camup) {
            state->camera.theta_x -= PI/100.0f;
        }

        if (input->left) {
            move -= v3_rotate_on_axis(v3(0,1,0), state->camera.theta_y, v3(1,0,0));
        }

        if (input->right) {
            move += 
                v3_rotate_on_axis(v3(0,1,0), state->camera.theta_y, v3(1,0,0));
        }

        if (input->down) {
            move -= 
                v3_rotate_q4(v3(0,0,1), 
                        (rotation_q4(state->camera.theta_y, v3(0,1,0)) * 
                         rotation_q4(state->camera.theta_x, v3(1,0,0))));

        }

        if (input->up) {
            move += 
                v3_rotate_q4(v3(0,0,1), 
                        (rotation_q4(state->camera.theta_y, v3(0,1,0)) * 
                         rotation_q4(state->camera.theta_x, v3(1,0,0))));

        }

        state->camera.pos += 1*v3_norm(move);



        if (fabs(state->camera.theta_x) > PI/2) {
            state->camera.theta_x = (state->camera.theta_x > 0) ? PI/2 : -PI/2;
        }

        if (fabs(state->camera.theta_y) >= 2*PI) {
            state->camera.theta_y = 0;
        }



        if (state->camera.width > WORLD_WIDTH) {
            state->camera.width = WORLD_WIDTH;
        }

        if (state->camera.height > WORLD_HEIGHT) {
            state->camera.height = WORLD_HEIGHT;
        }

        if ((state->camera.pos.x + state->camera.width/2.0f) > WORLD_RIGHT) {
            state->camera.pos.x = (WORLD_RIGHT - state->camera.width/2.0f);
        }

        if ((state->camera.pos.x - state->camera.width/2.0f) < WORLD_LEFT) {
            state->camera.pos.x = (WORLD_LEFT + state->camera.width/2.0f);
        }

        if ((state->camera.pos.y + state->camera.height/2.0f) > WORLD_TOP) {
            state->camera.pos.y = (WORLD_TOP - state->camera.height/2.0f);
        }

        if ((state->camera.pos.y - state->camera.height/2.0f) < WORLD_BOTTOM) {
            state->camera.pos.y = (WORLD_BOTTOM + state->camera.height/2.0f);
        }
    }

    // transform particles
    int cube_count = 0;
    int screen_count = 0;
    for (int i = 0;i < state->particle_count; i++) {
        uint32_t color = 0xFFFF0000;

        if (i % 2 == 0) {
            color = 0xFFFF00FF;
        }
            
        Vertex4Cube *cube = &state->particle_vert[cube_count++];
        V2Screen4 *screen1 = &state->particle_screen[screen_count++];
        V2Screen4 *screen2 = &state->particle_screen[screen_count++];

        cube->vertices[0].x = _mm_set1_ps(state->particles.pos[i].x);
        cube->vertices[1].x = _mm_set1_ps(state->particles.pos[i].x+1);

        ((float *)&cube->vertices[0].y)[0] = state->particles.pos[i].y;
        ((float *)&cube->vertices[0].y)[1] = state->particles.pos[i].y;
        ((float *)&cube->vertices[0].y)[2] = state->particles.pos[i].y+1;
        ((float *)&cube->vertices[0].y)[3] = state->particles.pos[i].y+1;
        ((float *)&cube->vertices[1].y)[0] = state->particles.pos[i].y;
        ((float *)&cube->vertices[1].y)[1] = state->particles.pos[i].y;
        ((float *)&cube->vertices[1].y)[2] = state->particles.pos[i].y+1;
        ((float *)&cube->vertices[1].y)[3] = state->particles.pos[i].y+1;

        ((float *)&cube->vertices[0].z)[0] = state->particles.pos[i].z;
        ((float *)&cube->vertices[0].z)[1] = state->particles.pos[i].z+1;
        ((float *)&cube->vertices[0].z)[2] = state->particles.pos[i].z;
        ((float *)&cube->vertices[0].z)[3] = state->particles.pos[i].z+1;
        ((float *)&cube->vertices[1].z)[0] = state->particles.pos[i].z;
        ((float *)&cube->vertices[1].z)[1] = state->particles.pos[i].z+1;
        ((float *)&cube->vertices[1].z)[2] = state->particles.pos[i].z;
        ((float *)&cube->vertices[1].z)[3] = state->particles.pos[i].z+1;

        Vertex4_to_V2Screen(cube->vertices, 
                            &state->camera, 
                            buffer->width, 
                            buffer->height, 
                            2, 
                            screen1);

        for (int j = 0; j < 4; ++j) {
            V2Screen start_pos;
            V2Screen end_pos;

            start_pos.x = ((int *)&screen1->x)[j];
            start_pos.y = ((int *)&screen1->y)[j];
            end_pos.x = ((int *)&screen2->x)[j];
            end_pos.y = ((int *)&screen2->y)[j];
            draw_line(buffer, start_pos, end_pos, 0xFFFFFFFF);

            start_pos.x = ((int *)&screen1->x)[j];
            start_pos.y = ((int *)&screen2->y)[j];
            end_pos.x = ((int *)&screen1->x)[j];
            end_pos.y = ((int *)&screen2->y)[j];
            draw_line(buffer, start_pos, end_pos, 0xFFFFFFFF);

            start_pos.x = ((int *)&screen2->x)[j];
            start_pos.y = ((int *)&screen1->y)[j];
            end_pos.x = ((int *)&screen2->x)[j];
            end_pos.y = ((int *)&screen1->y)[j];
            draw_line(buffer, start_pos, end_pos, 0xFFFFFFFF);

        }
    }

    //draw every spring
    for (int i = 0;i < state->spring_count; ++i) {
        Spring *spring = &state->springs[i];
        //render_line(buffer, start_pos, end_pos);
    }

    // TODO: render mouse_spring
            
    // update the simulation time until it syncs with the time after 1 video frame
    float t = state->time + TIME_FOR_FRAME;
    while (state->time < t) {


        // check for collisions and apply global forces

        for (int i = 0; i < state->particle_count; ++i) {

            float particle_left = state->particles.pos[i].x - (state->particles.radius[i]/2.0f);
            float particle_right = state->particles.pos[i].x + (state->particles.radius[i]/2.0f);
            float particle_top = state->particles.pos[i].y + (state->particles.radius[i]/2.0f);
            float particle_bottom = state->particles.pos[i].y - (state->particles.radius[i]/2.0f);

            float particle_x = state->particles.pos[i].x;
            float particle_y = state->particles.pos[i].y;

            float particle_dx = state->particles.pos[i].x;
            float particle_dy = state->particles.vel[i].y;

            if (particle_right > WORLD_RIGHT) {
                particle_dx = -particle_dx;
            }

            if (particle_left < WORLD_LEFT) {
                particle_dx = -particle_dx;
            }


            if (particle_top > WORLD_TOP) {
                particle_dy = -particle_dy;
            }

            if (particle_bottom < WORLD_BOTTOM) {
                particle_dy = -particle_dy;
            }


            //update spatial mask
            //
#if 0
            int left_grid_pos = (int)floorf((particles.pos.x - 
                        particles.radius) / 
                    (WORLD_WIDTH / 16.0f) + 8);

            int right_grid_pos = (int)floorf((particle->pos.x + 
                        particle->radius) / 
                    (WORLD_WIDTH / 16.0f) + 8);

            particle->spatial_mask = (1 << (left_grid_pos-1));
            particle->spatial_mask = particle->spatial_mask | (1 << (right_grid_pos-1));

            int top_grid_pos = (int)floorf(((particle->pos.y - particle->radius) / (WORLD_HEIGHT / 16.0f))+8);
            int bottom_grid_pos = (int)floorf(((particle->pos.y + particle->radius) / (WORLD_HEIGHT / 16.0f))+8);
            particle->spatial_mask = particle->spatial_mask | (1 << (16+top_grid_pos-1));
            particle->spatial_mask = particle->spatial_mask | (1 << (16+bottom_grid_pos));
#endif

#if 0 
            for (int j = 0; j< state->particle_count; ++j) {
                Particle *potential_collider = &state->particles[j];
                if (!potential_collider->active) {
                    continue;
                }
                
                // particles cant collide with themselves
                if (potential_collider == particle) {
                    continue;
                }
                
                // if not in the same grid then move onto the next
                int particle_x_mask = (particle->spatial_mask & 0x0000FFFF); // low 16
                int particle_y_mask = (particle->spatial_mask & 0xFFFF0000); // high 16
                                                                            //
                int collider_x_mask = (potential_collider->spatial_mask & 0x0000FFFF); // low 16
                int collider_y_mask = (potential_collider->spatial_mask & 0xFFFF0000);// high 16
                // if not colliding

                if (!((particle_x_mask & collider_x_mask) && 
                        (particle_y_mask & collider_y_mask))) {

                    continue;
                }

                
                if ((v3_mag(potential_collider->pos - particle->pos)) >= 
                    (particle->radius + potential_collider->radius)) {
                    continue;
                }

                if(particle == player) {
                    if (player->mass > potential_collider->mass) {
                        player->mass += 0.5f;    
                        potential_collider->active = 0;
                        continue;
                    }
                    
                    else {
                        player->mass *= 0.5f;
                        if (player->mass < 0.1f) {
                            player->mass = 0.1f;
                        }
                    }

                }

                V3 unit_normal = (potential_collider->pos - particle->pos) / 
                                 (v3_mag(potential_collider->pos - particle->pos));
                
                V3 unit_tangent = {};
                unit_tangent.x = -unit_normal.y;
                unit_tangent.y = unit_normal.x;

                float p1_dot = v3_dot((particle->vel-potential_collider->vel),
                                           (particle->pos-potential_collider->pos));

                float p2_dot = v3_dot((potential_collider->vel-particle->vel),
                                           (potential_collider->pos-particle->pos));


                
                V3 p1_after_vel =  ((2.0f*potential_collider->mass)/
                                  (particle->mass+potential_collider->mass)) *
                                 (p1_dot/((v3_mag(particle->pos-potential_collider->pos)) *
                                          (v3_mag(particle->pos-potential_collider->pos)))) *
                                 (particle->pos - potential_collider->pos);

                V3 p2_after_vel =  ((2.0f*particle->mass)/
                                  (particle->mass+potential_collider->mass)) *
                                 (p2_dot/((v3_mag(potential_collider->pos-particle->pos)) * 
                                          (v3_mag(potential_collider->pos-particle->pos)))) *
                                 (potential_collider->pos - particle->pos);

                //rollback for colliding particles
                while ((v3_mag(potential_collider->pos - particle->pos)) <
                        (particle->radius + potential_collider->radius)) {

                    particle->pos += -(timestep) * particle->vel;
                    potential_collider->pos += -(timestep) * potential_collider->vel;
                }

                particle->vel = particle->vel - p1_after_vel;
                potential_collider->vel = potential_collider->vel - p2_after_vel;
            }
#endif

            //zero forces
            state->particles.f_accumulator[i] = {}; 
            state->particles.f_accumulator[i] += -(COEFFICIENT_OF_DRAG*state->particles.vel[i]); 
        }

        // apply springs forces
        for (int i = 0; i < state->spring_count; ++i) {
            Spring *spring = &state->springs[i];
            spring_apply_force(spring);

        }

        if (state->mouse_spring.p1) {
            Spring *spring = &state->mouse_spring; 
            spring_apply_force(spring);
        }

        //update particles
        
        int n_particles = state->particle_count;
#if 1
        //update particles;
        for (int i = 0; i < arraysize(state->particles.pos); ++i) {
            state->particles.vel[i] += timestep * (state->particles.f_accumulator[i]/state->particles.mass[i]);
            state->particles.pos[i] += timestep * state->particles.vel[i];
        }

#endif

        
        state->time += timestep;
    }
}
