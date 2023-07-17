#include <cstdint>
#include <stdint.h>
#include <stdio.h>
#include <vcruntime.h>

#include "particle_fun.h"

#define PARTICLE_MASS 1
#define GRAVCONST 0.001f

static V3 gravity;

static uint32_t parkmiller_rand(uint32_t *state) {
const uint32_t A = 48271;

    uint32_t low  = (*state & 0x7fff) * A; // max: 32,767 * 48,271 = 1,581,695,857 = 0x5e46c371
    uint32_t high = (*state >> 15)    * A; // max: 65,535 * 48,271 = 3,163,439,985 = 0xbc8e4371
    uint32_t x = low + ((high & 0xffff) << 15) + (high >> 16); // max: 0x5e46c371 + 0x7fff8000 + 0xbc8e = 0xde46ffff

    x = (x & 0x7fffffff) + (x >> 31);
    *state = x;
    return *state;
}

inline static void string_inc_til_char(const char *string, char character, uint32_t *i) {
    for (;string[*i] != character; ++(*i)) {
        continue;
    }
}

// TODO Cleanup Messy Code
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
                out->vertexn[out->vertexn_count++] = v3(x, y, -z);
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

    while (out->vert_count % 4 != 0) {
        out->vertices[out->vert_count++] = v3(0, 0, 0);
    } 

    while (out->vertexn_count % 4 != 0) {
        out->vertexn[out->vertexn_count++] = v3(0, 0, 0);
    } 
    PlatformFreeFile(mesh_obj.file);
}

static V3 randomly_distribute_around_object(ParticleSystem *particles, int object_id, float close_scale, float far_scale, uint32_t *seed) {
    float x_offset = ((float)(parkmiller_rand(seed)%200)-100)/100.0f;
    float y_offset = ((float)(parkmiller_rand(seed)%200)-100)/100.0f;
    float z_offset = ((float)(parkmiller_rand(seed)%200)-100)/100.0f;
    V3 result = {};
    result.x = (particles->pos[object_id].x+((particles->radius[object_id]/2)*x_offset*close_scale) + 
                (x_offset)*(far_scale));

    result.y = (particles->pos[object_id].y+((particles->radius[object_id]/2)*y_offset*close_scale) + 
                (y_offset)*(far_scale));

    result.z = (particles->pos[object_id].z+((particles->radius[object_id]/2)*z_offset*close_scale) + 
                (z_offset)*(far_scale));

    return result;
}


static inline V3 calculate_chunk_position_offset(WorldChunk current_chunk, WorldChunk chunk) { 
    return v3((float)((chunk.x-current_chunk.x)*WORLD_WIDTH), 0, (float)((chunk.z-current_chunk.z)*WORLD_HEIGHT));
}

static inline V3 calculate_orbital_velocity(float center_body_mass,
                                            V3 center_body_position, 
                                            V3 orbiting_body_position,
                                            uint32_t *random_state) {

    V3 center_to_orbiting_vector = orbiting_body_position-center_body_position;

    if (center_to_orbiting_vector.x == 0 &&
        center_to_orbiting_vector.y == 0 &&
        center_to_orbiting_vector.z == 0) {

        return v3_zero();
    }

    // TODO investigate why the laws of physics in this game are broken
    // for some reason you have to multiply by 1000 so that the orbit is somewhat circular
    float magnitude_velocity = sqrtf(((1000*GRAVCONST*(center_body_mass))));
                                     //(v3_mag(center_to_orbiting_vector)));

    V3 velocity_vector = {};
    if (center_to_orbiting_vector.z != 0) {
        velocity_vector.x = (float)(parkmiller_rand(random_state) % 100);
        velocity_vector.y = (float)(parkmiller_rand(random_state) % 100);
        velocity_vector.z = (-(center_to_orbiting_vector.x*velocity_vector.x) -
                             (center_to_orbiting_vector.y*velocity_vector.y)) /
                            center_to_orbiting_vector.z;

    }

    else {
        velocity_vector.x = (float)(parkmiller_rand(random_state) % 100);
        velocity_vector.z = (float)(parkmiller_rand(random_state) % 100);
        velocity_vector.y = (-(center_to_orbiting_vector.x*velocity_vector.x) -
                             (center_to_orbiting_vector.z*velocity_vector.z)) /
                            center_to_orbiting_vector.y;
    }

    return magnitude_velocity * v3_norm(velocity_vector);
}

static void generate_chunk(ParticleSystem *particles, 
                           RendererState *render_state,
                           int number_of_particles, 
                           Mesh *mesh,
                           WorldChunk current_chunk,
                           WorldChunk chunk,
                           bool lighting) {

    uint32_t seed = (chunk.x << 16) | (chunk.z);
    int light_source_id = particles->particle_count;
    {
        float mass_mul = 100000.0f;
        particles->id[light_source_id] = light_source_id;
        particles->mass[light_source_id] = ((((float)(parkmiller_rand(&seed)%100)/100)*(300-100))+100)*mass_mul;
        particles->radius[light_source_id] = particles->mass[particles->particle_count]/mass_mul;

        particles->pos[light_source_id].x = (float)((parkmiller_rand(&seed)%WORLD_WIDTH)-WORLD_RIGHT);
        particles->pos[light_source_id].y = (float)((parkmiller_rand(&seed)%WORLD_HEIGHT)-WORLD_TOP);
        particles->pos[light_source_id].z = (float)((parkmiller_rand(&seed)%WORLD_DEPTH)-WORLD_FORWARD);

        particles->pos[light_source_id] += calculate_chunk_position_offset(current_chunk, chunk);
        particles->f_accumulator[light_source_id] = {};
        particles->vel[light_source_id] = {};
        particles->orbiting_body_id[light_source_id] = 0xFFFFFFFF;

        float sun_r = (float)(parkmiller_rand(&seed) % 100) / 100;
        float sun_g = (float)(parkmiller_rand(&seed) % 100) / 100;
        float sun_b = (float)(parkmiller_rand(&seed) % 100) / 100;

        particles->render_obj[particles->particle_count] =
            renderer_render_obj_create(render_state, mesh, v3(sun_r,sun_g,sun_b),
                                       particles->radius[particles->particle_count]*v3(1,1,1), 
                                       q4_identity(),
                                       particles->pos[particles->particle_count]);

        if (lighting) {
            renderer_render_obj_make_light_source(render_state,
                                                  &particles->render_obj[light_source_id],
                                                  particles->pos[light_source_id],
                                                  v3(sun_r,sun_g,sun_b),
                                                  (1/particles->radius[light_source_id]*0.01f));
        }
    }


    ++particles->particle_count;
    for (int i = 0; i < number_of_particles; i++) {
        //assert(state->particle_count < arraysize(state->particles));
        uint32_t particle_id = particles->particle_count;

        particles->id[particle_id] = particle_id;
        float mass_mul = 0.0001f;
        particles->mass[particle_id] = (float)(parkmiller_rand(&seed) % 200);
        if (particles->mass[particle_id] < 70) {
            particles->mass[particle_id] = ((float)(parkmiller_rand(&seed) % 10) + 1)*mass_mul;
        }

        else {
            particles->mass[particle_id] = ((float)(parkmiller_rand(&seed) % 5)+10)*mass_mul;
        }

        particles->radius[particle_id] = particles->mass[particle_id]/mass_mul;
        particles->pos[particle_id] = randomly_distribute_around_object(particles, 
                light_source_id, 
                300, 1000, &seed);

#if 1
       particles->vel[particle_id] = calculate_orbital_velocity(particles->mass[light_source_id],
                                                                 particles->pos[light_source_id],
                                                                 particles->pos[particle_id],
                                                                 &seed);

        particles->f_accumulator[particle_id] = {};
        float r = (float)(parkmiller_rand(&seed) % 100) / 100;
        float g = (float)(parkmiller_rand(&seed) % 100) / 100;
        float b = (float)(parkmiller_rand(&seed) % 100) / 100;

        particles->orbiting_body_id[particle_id] = light_source_id;

        particles->render_obj[particle_id] = renderer_render_obj_create(render_state, mesh, v3(1,1,1),
                                                                                      particles->radius[particle_id]*v3(1,1,1), 
                                                                                      q4_identity(),
                                                                                      particles->pos[particle_id]);
        float mass = particles->mass[particle_id];
        uint32_t parent_id = particles->id[particle_id];
        ++particles->particle_count;
#if 1
        if (mass > 10/mass_mul) {
            for (uint32_t moon = 0; moon < parkmiller_rand(&seed) % ((int)mass); ++moon) {
                uint32_t moon_id = particles->particle_count;
                particles->id[moon_id] = particles->particle_count;
                particles->mass[moon_id] = ((float)(parkmiller_rand(&seed) % ((int)particles->mass[parent_id])-3)*0.5f)*mass_mul;

                particles->radius[moon_id] = particles->mass[moon_id]/mass_mul;

                particles->pos[moon_id] = randomly_distribute_around_object(particles, 
                        parent_id, 
                        10, 100, &seed);

#if 1
                particles->vel[moon_id] = calculate_orbital_velocity(particles->mass[parent_id],
                                                                     particles->pos[parent_id],
                                                                     particles->pos[moon_id],
                                                                     &seed);
#endif

                particles->f_accumulator[moon_id] = {};
                particles->orbiting_body_id[moon_id] = parent_id;
                particles->render_obj[moon_id] = renderer_render_obj_create(render_state, 
                                                                            mesh, 
                                                                            v3(1,1,1),
                                                                            particles->radius[moon_id]*v3(1,1,1),
                                                                            q4_identity(),
                                                                            particles->pos[moon_id]);
                ++particles->particle_count;

            }
        }
#endif
    }
#endif
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
        particles->render_obj[particles->particle_count] = 
            renderer_render_obj_create(render_state, mesh, v3(1.0f,1.0f,1.0f),
                                       particles->mass[particles->particle_count]*v3(1,1,1),
                                       q4_identity(),
                                       particles->pos[particles->particle_count]);

        particles->particle_count++;
    }
}

static int create_particle(ParticleSystem *particles, 
                           RendererState *render_state,
                           float mass,
                           V3 position,
                           V3 velocity,
                           Mesh *mesh,
                           WorldChunk current_chunk,
                           WorldChunk chunk) {

    //assert(state->particle_count < arraysize(state->particles));
    int id = particles->particle_count;
    particles->id[particles->particle_count] = id;
    particles->mass[particles->particle_count] = mass;
    particles->radius[particles->particle_count] = particles->mass[particles->particle_count];
    particles->pos[particles->particle_count] = position;
    particles->pos[particles->particle_count] += calculate_chunk_position_offset(current_chunk, chunk);
    particles->f_accumulator[particles->particle_count] = {};
    particles->render_obj[particles->particle_count] = 
        renderer_render_obj_create(render_state, mesh, v3(1.0f,1.0f,1.0f),
                                   particles->mass[particles->particle_count]*v3(1,1,1),
                                   q4_identity(),
                                   particles->pos[particles->particle_count]);
    particles->particle_count++;

    return id;
}

inline static void spring_apply_force(GameState *state, Spring *spring) {

    V3 l = state->particles.pos[spring->p1_id] - state->particles.pos[spring->p2_id];
    V3 dl = state->particles.vel[spring->p1_id] - state->particles.vel[spring->p2_id];

    V3 force = -((spring->spring_const*(v3_mag(l) - spring->rest_length)) + 
            spring->damping_const*((v3_dot(l, dl))/v3_mag(l))) * (l/v3_mag(l));

    state->particles.f_accumulator[spring->p1_id] += force;

    state->particles.f_accumulator[spring->p2_id] += -force;
}

static void zero_buffer(void *buffer, size_t buffer_size) {
    for (int i = 0; i < buffer_size; ++i) {
        uint8_t *byte = (uint8_t *)buffer;
        *byte++ = 0;
    }
}

static void generate_world(GameState *state,
                           RendererState *render_state,
                           Mesh *sphere_mesh) {

    zero_buffer(render_state->vertex_buffer, render_state->vertex_count);
    zero_buffer(render_state->normal_buffer, render_state->normal_count);
    zero_buffer(render_state->polygons, render_state->polygon_count);
    zero_buffer(render_state->light_sources, render_state->light_sources_count);
    //reset ParticleSystem
    {
        zero_buffer(state->particles.pos, state->particles.particle_count);
        zero_buffer(state->particles.mass, state->particles.particle_count);
        zero_buffer(state->particles.id, state->particles.particle_count);
        zero_buffer(state->particles.vel, state->particles.particle_count);
        zero_buffer(state->particles.spatial_mask, state->particles.particle_count);
        zero_buffer(state->particles.render_obj, state->particles.particle_count);
        zero_buffer(state->particles.orbiting_body_id, state->particles.particle_count);
        zero_buffer(state->particles.radius, state->particles.particle_count);

    }


    render_state->vertex_count = 0;
    render_state->normal_count = 0;
    render_state->polygon_count = 0;
    render_state->light_sources_count = 0;
    state->particles.particle_count = 0;
    state->spring_count = 0;

    int offset_z = state->render_distance;
    int offset_x = state->render_distance;

    if (state->render_distance == 0) {
            generate_chunk(&state->particles, 
                           render_state, 10, 
                           sphere_mesh, 
                           state->current_chunk, 
                           state->current_chunk,
                           true);
    }

    for (int z = -offset_z; z <= offset_z; ++z) {
        for (int x = -offset_x; x <= offset_x; ++x) {
            WorldChunk chunk = {};
            chunk.x = state->current_chunk.x + x;
            chunk.z = state->current_chunk.z + z;
            generate_chunk(&state->particles, 
                           render_state, 10, 
                           sphere_mesh, 
                           state->current_chunk, 
                           chunk,
                           true);
        }
    }
} 
                                        
void game_update_and_render(GameMemory *memory, 
                            OffscreenBuffer *buffer, 
                            OffscreenBuffer *zbuffer,
                            OffscreenBuffer *normal_buffer,
                            OffscreenBuffer *postfx_buffer,
                            GameInput *input) {


    GameState *state = (GameState *)memory->permanent_storage;
    RendererState *render_state = (RendererState *)memory->transient_storage;
    int size_render = sizeof(RendererState);
    int size_assets = sizeof(Assets);
    assert((sizeof(RendererState)+sizeof(Assets)) < memory->transient_storage_size);
    Assets *assets = (Assets *)((char*)memory->transient_storage+sizeof(RendererState));
    Mesh *sphere_mesh = &assets->meshes[0];
    Mesh *spring_mesh = &assets->meshes[1];
    float timestep = 0.001f;
    
    if (!memory->is_initialised) {
        render_state->polygon_count = 0; 
        load_mesh_from_file_right("sphere.obj", sphere_mesh);
        load_mesh_from_file_right("spring.obj", spring_mesh);
        //set GRAVITY
        gravity.y = -9.81f; 
        state->current_chunk = {1000, 1000};
        state->chunk_id = 1;
        state->render_distance = 2;
        generate_world(state, render_state, sphere_mesh);
        //create_side_by_side_particles(&state->particles, render_state, 10, v3(0,10,0), v3(0,0,0), sphere_mesh);
        //create_side_by_side_particles(&state->particles, render_state, 10, v3(10,0,0), v3(0,0,0), sphere_mesh);
#if 0
        Spring *spring = &state->springs[state->spring_count++];
        spring->p1_id = state->particles.particle_count++;
        spring->p2_id = state->particles.particle_count++;
        spring->spring_const = 3.0f;
        spring->damping_const = 0.2f;
        spring->rest_length = 3;


        state->particles.mass[spring->p1_id] = 1;
        state->particles.radius[spring->p1_id] = 1;
        state->particles.pos[spring->p1_id] = v3(0,-20,0);
        state->particles.render_obj[spring->p1_id] = 
            renderer_render_obj_create(render_state, sphere_mesh, v3(1,1,1),
                                       state->particles.mass[spring->p1_id]*v3(1,1,1),
                                       q4_identity(),
                                       state->particles.pos[spring->p1_id]);

        state->particles.mass[spring->p2_id] = 1;
        state->particles.radius[spring->p2_id] = 1;
        state->particles.pos[spring->p2_id] = v3(0,0,0);
        state->particles.render_obj[spring->p2_id] = 
            renderer_render_obj_create(render_state, sphere_mesh, v3(1,1,1),
                                       state->particles.mass[spring->p1_id]*v3(1,1,1),
                                       q4_identity(),
                                       state->particles.pos[spring->p1_id]);
        
        spring->render_obj = 
        renderer_render_obj_create(render_state, 
                                   spring_mesh,
                                   v3(1,1,1),
                                   0.1f*v3(1,
                                           (1/5.904f*(state->particles.pos[spring->p2_id].y-state->particles.pos[spring->p2_id].y)),
                                            1),
                                     q4_identity(),
                                     v3(0,0,0));
#endif
        state->camera.pos.x = 0;
        state->camera.pos.y = 0;
        state->camera.pos.z = 0;
        state->camera.zfar = 500000;
        state->camera.znear = 0.01f; 
        state->camera.fov = PI/3.0f;
        state->camera.theta_y = PI;
        state->move_speed = 0.1f;
        memory->is_initialised = true;
    }


#if DEBUG_MODE
    if (state->frame_counter >= GAME_UPDATE_HZ) {
        state->frame_counter = 0;
    }
#endif



    {

        V3 move = {};
        float mouse_x = (float)input->mouse_pos.x*((float)buffer->width/(float)buffer->height);
        float mouse_y = (float)input->mouse_pos.y;
        mouse_y = mouse_y; 
        float mouse_mag = (sqrtf(mouse_x*mouse_x + mouse_y*mouse_y));

        if (input->mouse_pos.x) {
            state->camera.theta_y += ((mouse_x/mouse_mag)*(PI/50.0f));
        }

        if (input->mouse_pos.y) {
            state->camera.theta_x += ((mouse_y/mouse_mag)*(PI/50.0f));
        }

        if (input->camdown) {
            state->camera.theta_x += PI/100.0f;

        }

        if (input->camup) {
            state->camera.theta_x -= PI/100.0f;
        }
        
        if (input->camright) {
            state->camera.theta_y += PI/100.0f;

        }

        if (input->camleft) {
            state->camera.theta_y -= PI/100.0f;
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

#if 0
        if (state->camera.pos.z > WORLD_FORWARD) {
            state->camera.pos.z = WORLD_BACK;
            state->current_chunk = {state->current_chunk.x, state->current_chunk.z+1};
            generate_world(state, render_state, sphere_mesh);
        }

        if (state->camera.pos.z < WORLD_BACK) {
            state->camera.pos.z = WORLD_FORWARD;
            state->current_chunk = {state->current_chunk.x, state->current_chunk.z-1};
            generate_world(state, render_state, sphere_mesh);
        }

        if (state->camera.pos.x > WORLD_RIGHT) {
            state->camera.pos.x = WORLD_LEFT;
            state->current_chunk = {state->current_chunk.x+1, state->current_chunk.z};
            generate_world(state, render_state, sphere_mesh);
        }

        if (state->camera.pos.x < WORLD_LEFT) {
            state->camera.pos.x = WORLD_RIGHT;
            state->current_chunk = {state->current_chunk.x-1, state->current_chunk.z};
            generate_world(state, render_state, sphere_mesh);
        }
#endif



        if (fabs(state->camera.theta_x) > PI/2) {
            state->camera.theta_x = (state->camera.theta_x > 0) ? PI/2 : -PI/2;
        }

        if (fabs(state->camera.theta_y) >= 2*PI) {
            state->camera.theta_y = 0;
        }
    }

    state->camera.pos = state->particles.pos[state->locked_planet];

        if (input->action) {
            state->locked_planet ++;
            if (state->move_speed > 50000) {

                state->move_speed = 0.05f;
            }
        }

    float target_time = state->time+TIME_FOR_FRAME;
    while (state->time < target_time) {

        //draw_springs

        // check for collisions and apply global forces
        for (int i = 0; i < state->particles.particle_count; ++i) {
            //zero forces
            state->particles.f_accumulator[i] = {}; 
            //state->particles.f_accumulator[i] += -(COEFFICIENT_OF_DRAG*state->particles.vel[i]); 
            //
            uint32_t orbiting_body_id = state->particles.orbiting_body_id[i];
            if (!(orbiting_body_id == 0xFFFFFFFF)) {
                V3 planet_to_sun_vector =  state->particles.pos[i] - state->particles.pos[orbiting_body_id];

                float radius = v3_mag(planet_to_sun_vector);
                state->particles.f_accumulator[i] -= 
                    (((GRAVCONST*state->particles.mass[orbiting_body_id]) / (radius)) *
                     v3_norm(planet_to_sun_vector));
            }
#if GRAVITY
            state->particles.f_accumulator[i] += gravity; 
#endif
        }
        // apply springs forces
#if 0
        for (int i = 0; i < state->spring_count; ++i) {
            Spring *spring = &state->springs[i];
            spring_apply_force(state, spring);

        }
#endif

        for (int i = 0; i < state->particles.particle_count; ++i) {
            state->particles.vel[i] += (timestep * (state->particles.f_accumulator[i]/state->particles.mass[i]));
            state->particles.pos[i] += (timestep * state->particles.vel[i]);
        }
        state->time += timestep;
    }

    for (int i = 0;i < state->spring_count; ++i) {
        renderer_render_obj_update(render_state, 
                &state->springs[i].render_obj,
                0.1f*v3(1,
                    (1*(state->particles.pos[state->springs[i].p2_id].y-state->particles.pos[state->springs[i].p1_id].y)),
                    1),
                q4_identity(),
                v3(0,0,0));

    }
        //draw particles
        for (int i = 0;i < state->particles.particle_count; i++) {
#if 1
            renderer_render_obj_update(render_state,
                    &state->particles.render_obj[i], 
                    state->particles.render_obj[i].scale, 
                    state->particles.render_obj[i].rotation, 
                    state->particles.pos[i]);
#endif
        }



    renderer_draw(render_state, &state->camera, buffer, zbuffer, normal_buffer, postfx_buffer);
    printf("vertices: %d\n polygonsin: %d", render_state->vertex_count, render_state->polygon_count);

#if DEBUG_MODE
       state->frame_counter++;
#endif

} 
