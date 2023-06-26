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
    result.x = (particles->pos[object_id].x+((particles->mass[object_id]/2)*x_offset*close_scale) + 
                (x_offset)*(far_scale));

    result.y = (particles->pos[object_id].y+((particles->mass[object_id]/2)*y_offset*close_scale) + 
                (y_offset)*(far_scale));

    result.z = (particles->pos[object_id].z+((particles->mass[object_id]/2)*z_offset*close_scale) + 
                (z_offset)*(far_scale));

    return result;
}


static inline V3 calculate_chunk_position_offset(WorldChunk current_chunk, WorldChunk chunk) { 
    return v3((float)((chunk.x-current_chunk.x)*WORLD_WIDTH), 0, (float)((chunk.z-current_chunk.z)*WORLD_HEIGHT));
}

static void generate_chunk(ParticleSystem *particles, 
                           RendererState *render_state,
                           int number_of_particles, 
                           Mesh *mesh,
                           WorldChunk current_chunk,
                           WorldChunk chunk,
                           bool lighting) {

    uint32_t seed = (chunk.x << 16) | (chunk.z);
    int light_source = particles->particle_count;
    particles->id[particles->particle_count] = light_source;
    particles->mass[particles->particle_count] = (((float)(parkmiller_rand(&seed)%100)/100)*(300-100))+100;
    particles->radius[particles->particle_count] = particles->mass[particles->particle_count];
    particles->pos[particles->particle_count].x = (float)((parkmiller_rand(&seed)%WORLD_WIDTH)-WORLD_RIGHT);
    particles->pos[particles->particle_count].y = (float)((parkmiller_rand(&seed)%WORLD_HEIGHT)-WORLD_TOP);
    particles->pos[particles->particle_count].z = (float)((parkmiller_rand(&seed)%WORLD_DEPTH)-WORLD_FORWARD);
    particles->pos[particles->particle_count] += calculate_chunk_position_offset(current_chunk, chunk);
    particles->f_accumulator[particles->particle_count] = {};
    float sun_r = (float)(parkmiller_rand(&seed) % 100) / 100;
    float sun_g = (float)(parkmiller_rand(&seed) % 100) / 100;
    float sun_b = (float)(parkmiller_rand(&seed) % 100) / 100;
    particles->render_obj[particles->particle_count] = renderer_render_obj_create(render_state, mesh, v3(sun_r,sun_g,sun_b),
                                                                                  particles->mass[particles->particle_count]*v3(1,1,1), 
                                                                                  q4_identity(),
                                                                                  particles->pos[particles->particle_count]);


    if (lighting) {
        renderer_light_source_create(render_state,
                &particles->render_obj[light_source],
                particles->pos[light_source],
                    v3(sun_r,sun_g,sun_b));
    }


    ++particles->particle_count;

    for (int i = 0; i < number_of_particles; i++) { 
        //assert(state->particle_count < arraysize(state->particles));
        particles->id[particles->particle_count] = particles->particle_count;
        particles->mass[particles->particle_count] = (float)(parkmiller_rand(&seed) % 100);
        if (particles->mass[particles->particle_count] < 70) {
            particles->mass[particles->particle_count] = (float)(parkmiller_rand(&seed) % 10) + 1;
        }

        else {
            particles->mass[particles->particle_count] = (float)(parkmiller_rand(&seed) % 5)+10;
        }

        particles->radius[particles->particle_count] = particles->mass[particles->particle_count];
        particles->pos[particles->particle_count] = randomly_distribute_around_object(particles, 
                light_source, 
                300, 1000, &seed);

#if 0
        particles->vel[particles->particle_count].x = (float)(parkmiller_rand(&seed)%20)-10;
        particles->vel[particles->particle_count].y = (float)(parkmiller_rand(&seed)%20)-10;
        particles->vel[particles->particle_count].z = 0;
#endif
        particles->f_accumulator[particles->particle_count] = {};
        float r = (float)(parkmiller_rand(&seed) % 100) / 100;
        float g = (float)(parkmiller_rand(&seed) % 100) / 100;
        float b = (float)(parkmiller_rand(&seed) % 100) / 100;

        particles->render_obj[particles->particle_count] = renderer_render_obj_create(render_state, mesh, v3(1,1,1),
                                                                                      particles->mass[particles->particle_count]*v3(1,1,1), 
                                                                                      q4_identity(),
                                                                                      particles->pos[particles->particle_count]);
        float mass = particles->mass[particles->particle_count];
        int parent_id = particles->id[particles->particle_count];
        ++particles->particle_count;
        if (mass > 10) {
            for (uint32_t moon = 0; moon < parkmiller_rand(&seed) % ((int)mass); ++moon) {
                particles->id[particles->particle_count] = particles->particle_count;
                particles->mass[particles->particle_count] = (float)(parkmiller_rand(&seed) % ((int)particles->mass[parent_id])-3)*0.5f;

                particles->radius[particles->particle_count] = particles->mass[particles->particle_count];

                particles->pos[particles->particle_count] = randomly_distribute_around_object(particles, 
                        parent_id, 
                        10, 100, &seed);

                particles->vel[particles->particle_count] = particles->vel[parent_id];
                particles->f_accumulator[particles->particle_count] = {};
                particles->render_obj[particles->particle_count] = renderer_render_obj_create(render_state, 
                                                                                              mesh, 
                                                                                              v3(1,1,1),
                                                                                              particles->mass[particles->particle_count]*v3(1,1,1),
                                                                                              q4_identity(),
                                                                                              particles->pos[particles->particle_count]);
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

static void generate_world(GameState *state,
                           RendererState *render_state,
                           Mesh *sphere_mesh) {

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
        create_side_by_side_particles(&state->particles, render_state, 10, v3(0,10,0), v3(0,0,0), sphere_mesh);
        create_side_by_side_particles(&state->particles, render_state, 10, v3(10,0,0), v3(0,0,0), sphere_mesh);
        Spring *spring = &state->springs[state->spring_count++];
        spring->p1_id = state->particles.particle_count++;
        spring->p2_id = state->particles.particle_count++;
        spring->spring_const = 3.0f;
        spring->damping_const = 0.2f;
        spring->rest_length = 3;


        state->particles.mass[spring->p1_id] = 1;
        state->particles.radius[spring->p1_id] = 1;
        state->particles.pos[spring->p1_id] = v3(0,-5,0);
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
        state->camera.pos.x = 0;
        state->camera.pos.y = 0;
        state->camera.pos.z = -20;
        state->camera.zfar = 100000; 
        state->camera.znear = 0.01f; 
        state->camera.fov = PI/3.0f;
        state->move_speed = 300.0f;
        memory->is_initialised = true;
    }


#if DEBUG_MODE
    if (state->frame_counter >= GAME_UPDATE_HZ) {
        state->frame_counter = 0;
    }
#endif

    renderer_draw_background(buffer, 0); 

    //clear z buffer every frame
    {
        uint32_t zbuffer_size = zbuffer->width * zbuffer->height;
        float *depth_value = ((float *)zbuffer->memory);
        for (uint32_t i = 0; i < zbuffer_size; ++i) {
            *depth_value++ = FLT_MAX;
        }
    }

    //clear normal buffer every frame
    {
        uint32_t normal_buffer_size = normal_buffer->width * normal_buffer->height;
        float *normal_value = ((float *)normal_buffer->memory);
        for (uint32_t i = 0; i < normal_buffer_size; ++i) {
            *normal_value++ = 0;
        }
    }


    float timestep = TIME_FOR_FRAME;


    {
        V3 move = {};
        if (input->action) {
            state->move_speed *= 100;
            if (state->move_speed > 50000) {

                state->move_speed = 0.05f;
            }
        }
            
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



        if (fabs(state->camera.theta_x) > PI/2) {
            state->camera.theta_x = (state->camera.theta_x > 0) ? PI/2 : -PI/2;
        }

        if (fabs(state->camera.theta_y) >= 2*PI) {
            state->camera.theta_y = 0;
        }
    }

    //draw particles
    for (int i = 0;i < state->particles.particle_count; i++) {
#if 1
        renderer_render_obj_update(render_state,
                                   &state->particles.render_obj[i], 
                                   state->particles.mass[i]*v3(1,1,1), q4_identity(),
                                   state->particles.pos[i]);
#endif
    }

    //draw_springs
    for (int i = 0;i < state->spring_count; ++i) {
        renderer_render_obj_update(render_state, 
                                   &state->springs[i].render_obj,
                                   0.1f*v3(1,
                                           (1/5.904f*(state->particles.pos[state->springs[i].p2_id].y-state->particles.pos[state->springs[i].p2_id].y)),
                                            1),
                                     q4_identity(),
                                     v3(0,0,0));

    }

    // update the simulation time until it syncs with the time after 1 video frame
    float t = state->time + TIME_FOR_FRAME;
    while (state->time < t) {

        // check for collisions and apply global forces

        for (int i = 0; i < state->particles.particle_count; ++i) {
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

        for (int i = 0; i < state->particles.particle_count; ++i) {
            state->particles.vel[i] += timestep * (state->particles.f_accumulator[i]/state->particles.mass[i]);
            state->particles.pos[i] += timestep * state->particles.vel[i];
        }
       state->time += timestep;
    }

    renderer_draw(render_state, &state->camera, buffer, zbuffer, normal_buffer, postfx_buffer);

#if DEBUG_MODE
       state->frame_counter++;
#endif

} 
