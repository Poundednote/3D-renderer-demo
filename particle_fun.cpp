#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <emmintrin.h>

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

inline static void string_inc_til_char(const char *string, char character, uint32_t *i) {
    for (;string[*i] != character; ++(*i)) {
        continue;
    }
}

static void load_mesh_from_file_right(char *filename, Mesh *out) {
    ReadFileResult mesh_obj = PlatformReadFile(filename);
    char *string = (char *)mesh_obj.file;
    for (uint32_t i = 0; i < mesh_obj.size; ++i) {
        if (string[i] == '#') {
            string_inc_til_char(string, '\n', &i);
        }
        else if (string[i] == 'v') {
            ++i;
            if (string[i] == 'n') {
                float x, y, z;
                sscanf_s(string+(++i), "%f %f %f\n", &x, &y, &z); 
                out->vertexn[out->vertexn_count++] = v3_norm((v3(x, y, -z)));
            }
            else if (string[i] == 't') {
                continue;
            }

            else {
                float x, y, z;
                sscanf_s(string+(i), "%f %f %f\n", &x, &y, &z); 
                out->vertices[out->vert_count++] = v3(x, y, -z);
            }
        }
        else if (string[i] == 'f') {
            Triangle polygon;
            sscanf_s(string+(++i), "%d/", &polygon.v3);
            string_inc_til_char(string, '/', &i);
            ++i;
            string_inc_til_char(string, '/', &i);
            sscanf_s(string+(++i), "%d", &polygon.vn3);

            string_inc_til_char(string, ' ', &i);
            sscanf_s(string+(++i), "%d/", &polygon.v2);
            string_inc_til_char(string, '/', &i);
            ++i;
            string_inc_til_char(string, '/', &i);
            sscanf_s(string+(++i), "%d", &polygon.vn2);

            string_inc_til_char(string, ' ', &i);
            sscanf_s(string+(++i), "%d/", &polygon.v1);
            string_inc_til_char(string, '/', &i);
            ++i;
            string_inc_til_char(string, '/', &i);
            sscanf_s(string+(++i), "%d", &polygon.vn1);

            polygon.v1 -= 1;
            polygon.v2 -= 1;
            polygon.v3 -= 1;
            polygon.vn1 -= 1;
            polygon.vn2 -= 1;
            polygon.vn3 -= 1;

            out->polygons[out->poly_count++] = polygon;
        }

        else {
            string_inc_til_char(string, '\n', &i);
        }
    }
    PlatformFreeFile(mesh_obj.file);
}
static V3 randomly_distribute_around_object(ParticleSystem *particles, int object_id, float close_scale, float far_scale, uint32_t *seed) {
    float x_offset = ((float)(parkmiller_rand(seed)%200)-100)/100.0f;
    float y_offset = ((float)(parkmiller_rand(seed)%200)-100)/100.0f;
    float z_offset = ((float)(parkmiller_rand(seed)%200)-100)/100.0f;
    V3 result = {};
    result.x = (particles->pos[object_id].x+((particles->mass[object_id]/2)*x_offset*close_scale) + 
                (x_offset)*(far_scale));

    result.y = (particles->pos[object_id].y+((particles->mass[object_id]/2)*y_offset*close_scale) + 
                (y_offset)*(far_scale));

    result.z = (particles->pos[object_id].z+((particles->mass[object_id]/2)*z_offset*close_scale) + 
                (z_offset)*(far_scale));

    return result;
}



#if DEBUG_MODE

static void generate_world(ParticleSystem *particles, 
                           RendererState *render_state,
                           int light_source_id,
                           int number_of_particles, 
                           Mesh *mesh,
                           uint32_t seed) {

    for (int i = 0; i < number_of_particles; i++) { 
        //assert(state->particle_count < arraysize(state->particles));
        particles->id[particles->particle_count] = particles->particle_count;
        particles->mass[particles->particle_count] = (float)(parkmiller_rand(&seed) % 100);
        if (particles->mass[particles->particle_count] < 90) {
            particles->mass[particles->particle_count] = (float)(parkmiller_rand(&seed) % 10) + 1;
        }

        else {
            particles->mass[particles->particle_count] = (float)(parkmiller_rand(&seed) % 5)+10;
        }

        particles->radius[particles->particle_count] = particles->mass[particles->particle_count];
        particles->pos[particles->particle_count] = randomly_distribute_around_object(particles, 
                                                                                      light_source_id, 
                                                                                      100, 100, &seed);
#if 1
        particles->vel[particles->particle_count].x = (float)(parkmiller_rand(&seed)%20)-10;
        particles->vel[particles->particle_count].y = (float)(parkmiller_rand(&seed)%20)-10;
        particles->vel[particles->particle_count].z = 0;
#endif
        particles->f_accumulator[particles->particle_count] = {};
        float r = (float)(parkmiller_rand(&seed) % 100) / 100;
        float g = (float)(parkmiller_rand(&seed) % 100) / 100;
        float b = (float)(parkmiller_rand(&seed) % 100) / 100;

        particles->render_obj[particles->particle_count] = create_render_obj(render_state, mesh, v3(r,g,b));
        float mass = particles->mass[particles->particle_count];
        int parent_id = particles->id[particles->particle_count];
        ++particles->particle_count;
        if (mass > 15) {
            for (uint32_t moon = 0; moon < parkmiller_rand(&seed) % ((int)mass-20); ++moon) {
                particles->id[particles->particle_count] = particles->particle_count;
                particles->mass[particles->particle_count] = (float)(parkmiller_rand(&seed) % ((int)particles->mass[parent_id])-20)*0.5f;

                particles->radius[particles->particle_count] = particles->mass[particles->particle_count];
                particles->pos[particles->particle_count] = randomly_distribute_around_object(particles, 
                                                                                              parent_id, 
                                                                                              10, 10, &seed);
                particles->vel[particles->particle_count] = particles->vel[parent_id];
                particles->f_accumulator[particles->particle_count] = {};
                particles->render_obj[particles->particle_count] = create_render_obj(render_state, mesh, v3(r,g,b));
                ++particles->particle_count;

            }
        }
    }
}

static void create_side_by_side_particles(ParticleSystem *particles, 
                                          RendererState *render_state,
                                          int number_of_particles, 
                                          V3 pad,
                                          V3 start,
                                          Mesh *mesh) {

    for (int i = 0; i < number_of_particles; ++i) { 
        //assert(state->particle_count < arraysize(state->particles));
        particles->id[particles->particle_count] = particles->particle_count;
        particles->mass[particles->particle_count] = PARTICLE_MASS;
        particles->radius[particles->particle_count] = particles->mass[particles->particle_count];
        particles->pos[particles->particle_count] = ((float)i*pad)+start;
        particles->f_accumulator[particles->particle_count] = {};
        particles->render_obj[particles->particle_count] = create_render_obj(render_state, mesh, v3(1.0f,1.0f,1.0f));

        particles->particle_count++;
    }
}
#endif

static int create_particle(ParticleSystem *particles, 
                           RendererState *render_state,
                           float mass,
                           V3 position,
                           V3 velocity,
                           Mesh *mesh) {

    //assert(state->particle_count < arraysize(state->particles));
    int id = particles->particle_count;
    particles->id[particles->particle_count] = id;
    particles->mass[particles->particle_count] = mass;
    particles->radius[particles->particle_count] = particles->mass[particles->particle_count];
    particles->pos[particles->particle_count] = position;
    particles->f_accumulator[particles->particle_count] = {};
    particles->render_obj[particles->particle_count] = create_render_obj(render_state, mesh, v3(1,1,1));
    particles->particle_count++;

    return id;
}

inline static void spring_apply_force(GameState *state, Spring *spring) {

    V3 l = state->particles.pos[spring->p1index] - state->particles.pos[spring->p2index];
    V3 dl = state->particles.vel[spring->p1index] - state->particles.vel[spring->p1index];

    V3 force = -((spring->spring_const*(v3_mag(l) - spring->rest_length)) + 
            spring->damping_const*((v3_dot(l, dl))/v3_mag(l))) * (l/v3_mag(l));

    state->particles.f_accumulator[spring->p1index] += force;

    state->particles.f_accumulator[spring->p2index] += -force;
}

void game_update_and_render(GameMemory *memory, 
                            OffscreenBuffer *buffer, 
                            OffscreenBuffer *zbuffer,
                            GameInput *input) {


    GameState *state = (GameState *)memory->permanent_storage;
    RendererState *render_state = (RendererState *)memory->transient_storage;
    Assets *assets = (Assets *)((char*)memory->transient_storage+sizeof(RendererState));
    Mesh *sphere_mesh = &assets->meshes[0];
    Mesh *spring_mesh = &assets->meshes[1];
    if (!memory->is_initialised) {
        render_state->polygon_count = 0; 
        
       load_mesh_from_file_right("sphere.obj", sphere_mesh);
       load_mesh_from_file_right("spring.obj", spring_mesh);
        //set GRAVITY
        gravity.y = -9.81f; 
#if DEBUG_MODE
        int light_source = create_particle(&state->particles, render_state, 100, v3(0, 10, -30), v3_zero(), sphere_mesh);
        generate_world(&state->particles, render_state, light_source, 50, sphere_mesh, 1);
        create_side_by_side_particles(&state->particles, render_state, 10, v3(0,10,0), v3(0,0,0), sphere_mesh);
        create_side_by_side_particles(&state->particles, render_state, 10, v3(10,0,0), v3(0,0,0), sphere_mesh);

        render_state->light_sources[render_state->light_sources_count].position = state->particles.pos[light_source];
        render_state->light_sources[render_state->light_sources_count].color = v3(1,1,1);
        render_state->light_sources[render_state->light_sources_count].obj = &state->particles.render_obj[light_source];
        render_state->light_sources[render_state->light_sources_count].color = v3(1,1,1);
        ++render_state->light_sources_count;
#endif

        Spring *spring = &state->springs[state->spring_count++];
        spring->p1index = state->particles.particle_count++;
        spring->p2index = state->particles.particle_count++;
        spring->spring_const = 3.0f;
        spring->damping_const = 0.2f;
        spring->rest_length = 3;


        state->particles.mass[spring->p1index] = 1;
        state->particles.radius[spring->p1index] = 1;
        state->particles.pos[spring->p1index] = v3(0,-5,0);
        state->particles.render_obj[spring->p1index] = create_render_obj(render_state, sphere_mesh, v3(1,1,1));

        state->particles.mass[spring->p2index] = 1;
        state->particles.radius[spring->p2index] = 1;
        state->particles.pos[spring->p2index] = v3(0,0,0);
        state->particles.render_obj[spring->p2index] = create_render_obj(render_state, sphere_mesh, v3(1,1,1));
        
        spring->render_obj = create_render_obj(render_state, spring_mesh, v3(1,0,0));

        state->camera.width = 300;
        state->camera.height = 300;
        state->camera.pos.x = 0;
        state->camera.pos.y = 0;
        state->camera.pos.z = 0;
        state->camera.zfar = 10000; 
        state->camera.znear = 0.01f; 
        state->camera.fov = PI/3.0f;
        state->move_speed = 89.4f*TIME_FOR_FRAME;
        memory->is_initialised = true;
    }


#if DEBUG_MODE
    if (state->frame_counter >= GAME_UPDATE_HZ) {
        state->frame_counter = 0;
    }
#endif

    renderer_draw_background(buffer, 0xFF030303); 

    //clear z buffer every frame
    uint32_t zbuffer_size = zbuffer->width * zbuffer->height;
    float *depth_value = ((float *)zbuffer->memory);
    for (uint32_t i = 0; i < zbuffer_size; ++i) {
        *depth_value++ = FLT_MAX;
    }


    float timestep = TIME_FOR_FRAME;


    {
        V3 move = {};
        if (input->action) {
            state->move_speed *= 10.0f;
            if (state->move_speed >= 50) {
                state->move_speed = 0.05f;
            }
        }
            
        float mouse_x = (float)input->mouse_pos.x*((float)buffer->width/(float)buffer->height);
        float mouse_y = (float)input->mouse_pos.y;
        mouse_y = mouse_y; 
        float mouse_mag = (sqrtf(mouse_x*mouse_x + mouse_y*mouse_y));

        if (input->mouse_pos.x) {
            state->camera.theta_y += ((mouse_x/mouse_mag)*(PI/20.0f));
        }

        if (input->mouse_pos.y) {
            state->camera.theta_x += ((mouse_y/mouse_mag)*(PI/20.0f));
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

        state->camera.pos += state->move_speed*v3_norm(move);


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

    int cube_count = 0;
    int screen_count = 0;
    //draw particles
    for (int i = 0;i < state->particles.particle_count; i++) {
        draw_render_obj(render_state, 
                          &state->particles.render_obj[i], 
                          state->particles.mass[i]*v3(1,1,1), q4(1, v3(0,0,0)), 
                          state->particles.pos[i]);
    }
    //draw_springs
    for (int i = 0;i < state->spring_count; ++i) {
        draw_render_obj(render_state, 
                        &state->springs[i].render_obj,
                        0.1f*v3(1,
                                (1/5.904f*(state->particles.pos[state->springs[i].p2index].y-state->particles.pos[state->springs[i].p1index].y)),
                                 1),
                          q4(1, v3(0,0,0)),
                          v3(0,0,0));

    }

    // update the simulation time until it syncs with the time after 1 video frame
    float t = state->time + TIME_FOR_FRAME;
    while (state->time < t) {

        // check for collisions and apply global forces

        for (int i = 0; i < state->particles.particle_count; ++i) {


            float particle_x = state->particles.pos[i].x;
            float particle_y = state->particles.pos[i].y;
            float particle_z = state->particles.pos[i].z;

            float particle_dx = state->particles.vel[i].x;
            float particle_dy = state->particles.vel[i].y;
            float particle_dz = state->particles.vel[i].z;

            if (particle_x > WORLD_RIGHT) {
                particle_dx = -particle_dx;
            }

            if (particle_x < WORLD_LEFT) {
                particle_dx = -particle_dx;
            }


            if (particle_y > WORLD_TOP) {
                particle_dy = -particle_dy;
            }

            if (particle_y < WORLD_BOTTOM) {
                particle_dy = -particle_dy;
            }

            if (particle_z > WORLD_FORWARD) {
                particle_dy = -particle_dy;
            }

            if (particle_z < WORLD_BACK) {
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
            //state->particles.f_accumulator[i] += -(COEFFICIENT_OF_DRAG*state->particles.vel[i]); 
#if GRAVITY
            state->particles.f_accumulator[i] += gravity; 
#endif
        }

        // apply springs forces
        for (int i = 0; i < state->spring_count; ++i) {
            Spring *spring = &state->springs[i];
            spring_apply_force(state, spring);

        }

        //update particles
        //update particles;
        for (int i = 0; i < state->particles.particle_count; ++i) {
            state->particles.vel[i] += timestep * (state->particles.f_accumulator[i]/state->particles.mass[i]);
            state->particles.pos[i] += timestep * state->particles.vel[i];
        }

       state->time += timestep;

    }

    renderer_transform_light_and_cull(render_state,
                                      &state->camera, 
                                      buffer->width, 
                                      buffer->height);


    renderer_draw_triangles_filled(buffer,
                                   zbuffer,
                                   render_state->vertex_list,
                                   render_state->vertex_colors,
                                   render_state->polygons_to_draw, 
                                   render_state->draw_count);

#if DEBUG_MODE
       state->frame_counter++;
#endif

} 
