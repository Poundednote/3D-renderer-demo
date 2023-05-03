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


#if DEBUG_MODE

static void create_random_particles(GameState *state, 
                                    RendererState *render_state,
                                    int number_of_particles, 
                                    Mesh *mesh,
                                    uint32_t seed) {

    for (int i = 0; i < number_of_particles; i++) { 
        //assert(state->particle_count < arraysize(state->particles));
        state->particles.mass[state->particle_count] = (float)(parkmiller_rand(&seed) % 30) + 10;
        state->particles.radius[state->particle_count] = state->particles.mass[state->particle_count];
        state->particles.pos[state->particle_count].x = (float)(parkmiller_rand(&seed)%WORLD_WIDTH)-WORLD_RIGHT;
        state->particles.pos[state->particle_count].y = (float)(parkmiller_rand(&seed)%WORLD_HEIGHT)-WORLD_TOP;
        state->particles.pos[state->particle_count].z = (float)(parkmiller_rand(&seed)%WORLD_DEPTH)-WORLD_FORWARD;
#if 1
        state->particles.vel[state->particle_count].x = (float)(parkmiller_rand(&seed)%20)-10;
        state->particles.vel[state->particle_count].y = (float)(parkmiller_rand(&seed)%20)-10;
        state->particles.vel[state->particle_count].z = 0;
#endif
        state->particles.f_accumulator[state->particle_count] = {};

        state->particle_count++;
    }

    for (int i = 0; i < number_of_particles; ++i) {
        int vertex_start = render_state->vertex_count;
        for (int vert = 0; vert < mesh->vert_count; ++vert) {
            render_state->vertex_list[render_state->vertex_count++] = 
                (0.5*mesh->vertices[vert])+state->particles.pos[i];
        }

        for (int vert = 0; vert < mesh->vertexn_count; ++vert) {
            render_state->vertexn_list[render_state->vertexn_count++] =
                mesh->vertexn[vert];
        }

        for (int polygon = 0; polygon < mesh->poly_count; ++polygon) {
            render_state->polygons[render_state->polygon_count] = mesh->polygons[polygon];
            render_state->polygons[render_state->polygon_count].v1 += vertex_start;
            render_state->polygons[render_state->polygon_count].v2 += vertex_start;
            render_state->polygons[render_state->polygon_count].v3 += vertex_start;
            render_state->polygons[render_state->polygon_count].vn1 += vertex_start;
            render_state->polygons[render_state->polygon_count].vn2 += vertex_start;
            render_state->polygons[render_state->polygon_count].vn3 += vertex_start;
            render_state->polygon_count++;
        }

    }

}

static void create_side_by_side_particles(GameState *state, 
                                          RendererState *render_state,
                                          int number_of_particles, 
                                          V3 pad,
                                          V3 start,
                                          Mesh *mesh) {

    for (int i = 0; i < number_of_particles; ++i) { 
        //assert(state->particle_count < arraysize(state->particles));
        state->particles.mass[state->particle_count] = PARTICLE_MASS;
        state->particles.radius[state->particle_count] = state->particles.mass[state->particle_count];
        state->particles.pos[state->particle_count] = ((float)i*pad)+start;
        state->particles.f_accumulator[state->particle_count] = {};

        state->particle_count++;
    }

    for (int i = 0; i < number_of_particles; ++i) {
        int vertex_start = render_state->vertex_count;
        for (int vert = 0; vert < mesh->vert_count; ++vert) {
            render_state->vertex_list[render_state->vertex_count++] 
                = (0.1f*mesh->vertices[vert])+state->particles.pos[i];
        }

        int vertexnorm_start = render_state->vertexn_count;
        for (int normal = 0; normal < mesh->vertexn_count; ++normal) {
            render_state->vertexn_list[render_state->vertexn_count++] =
                mesh->vertexn[normal];
        }

        for (int polygon = 0; polygon < mesh->poly_count; ++polygon) {
            render_state->polygons[render_state->polygon_count] = mesh->polygons[polygon];
            render_state->polygons[render_state->polygon_count].v1 += vertex_start;
            render_state->polygons[render_state->polygon_count].v2 += vertex_start;
            render_state->polygons[render_state->polygon_count].v3 += vertex_start;
            render_state->polygons[render_state->polygon_count].vn1 += vertexnorm_start;
            render_state->polygons[render_state->polygon_count].vn2 += vertexnorm_start;
            render_state->polygons[render_state->polygon_count].vn3 += vertexnorm_start;
            render_state->polygon_count++;
        }

    }
}

#endif

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
                            OffscreenBuffer *zbuffer,
                            GameInput *input) {


    GameState *state = (GameState *)memory->permanent_storage;
    RendererState *render_state = (RendererState *)memory->transient_storage;
    Mesh *sphere_mesh = (Mesh *)((char *)memory->transient_storage+sizeof(RendererState));
    if (!memory->is_initialised) {
        ReadFileResult sphere_obj = PlatformReadFile("sphere.obj");
        render_state->polygon_count = 0;
        for (uint32_t i = 0; i < sphere_obj.size; ++i) {
            //find first character begging with v
            char *string = (char *)sphere_obj.file;
            if (string[i] == 'v') {
                if (string[++i] == 'n') {
                    float x, y, z;
                    sscanf_s(string+(++i), "%f %f %f\n", &x, &y, &z); 
                    sphere_mesh->vertexn[sphere_mesh->vertexn_count++] = v3_norm((v3(x, y, -z)));
                }

                else {
                    float x, y, z;
                    sscanf_s(string+(i), "%f %f %f\n", &x, &y, &z); 
                    sphere_mesh->vertices[sphere_mesh->vert_count++] = v3(x, y, -z);
                }
            }
            if (string[i] == 'f') {
                Triangle polygon;
                sscanf_s(string+(++i), "%d/%*d/%d %d/%*d/%d %d/%*d/%d\n", 
                         &polygon.v3, &polygon.vn3, &polygon.v2, &polygon.vn2, &polygon.v1, &polygon.vn1); 
                polygon.v1 -= 1;
                polygon.v2 -= 1;
                polygon.v3 -= 1;
                polygon.vn1 -= 1;
                polygon.vn2 -= 1;
                polygon.vn3 -= 1;
                
                sphere_mesh->polygons[sphere_mesh->poly_count++] = polygon;
            }
        }
        PlatformFreeFile(sphere_obj.file);
        
        //set GRAVITY
        gravity.y = -9.81f; 
#if DEBUG_MODE
        create_side_by_side_particles(state, render_state, 10, v3(10,0,0), v3(0,0,0), sphere_mesh);
        create_side_by_side_particles(state, render_state, 10, v3(0,0,10), v3(0,0,0), sphere_mesh);
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
        state->camera.zfar = 5000; 
        state->camera.znear = 1; 
        state->camera.fov = PI/3.0f;
        state->move_speed = 1.0f;
        memory->is_initialised = true;
    }


#if DEBUG_MODE
    if (state->frame_counter >= GAME_UPDATE_HZ) {
        state->frame_counter = 0;
    }
#endif


    render_state->vertex_count = 0;
    render_state->vertexn_count = 0;
    render_state->screen_vertex_count = 0;
    renderer_draw_background(buffer, 0xFFFF00FF); 

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
    for (int i = 0;i < state->particle_count; i++) {

#if SSE
        V2Screen4 screen[2];
        Vertex4Cube *cube = &state->particle_vert[cube_count++];
#else
        for (int vert = 0; vert < sphere_mesh->vert_count; ++vert) {
            render_state->vertex_list[render_state->vertex_count++] = (0.1f*sphere_mesh->vertices[vert]) + state->particles.pos[i];
        }
        for (int normal = 0; normal < sphere_mesh->vertexn_count; ++normal) {
            render_state->vertexn_list[render_state->vertexn_count++] = sphere_mesh->vertexn[normal]; 
        }
#endif
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

       renderer_transform_light_and_cull(render_state,
                                         &state->camera, 
                                         buffer->width, 
                                         buffer->height);


       renderer_draw_triangles_filled(buffer,
                                      zbuffer,
                                      render_state->screen_vertices,
                                      render_state->polygons_to_draw, 
                                      render_state->draw_count);

       state->time += timestep;
    }

} 
