#include <emmintrin.h>
#include <stdint.h>

#include "renderer.h"

static void renderer_vertex4_to_v2screen(Vertex4 *in, 
                         GameCamera *camera, 
                         int screen_width,
                         int screen_height,
                         int count, 
                         V3Screen4 *out) {

    Quaternion rotation = q4_norm(rotation_q4(-camera->theta_x, v3(1,0,0))*
                           rotation_q4(-camera->theta_y, v3(0,1,0)));

    float aspect_ratio = (float)screen_width / (float)screen_height;
    
    __m128 _qw = _mm_set1_ps(rotation.scalar); 
    __m128 _qx = _mm_set1_ps(rotation.vector.x); 
    __m128 _qy = _mm_set1_ps(rotation.vector.y); 
    __m128 _qz = _mm_set1_ps(rotation.vector.z); 

    __m128 _cam_x = _mm_set1_ps(camera->pos.x);
    __m128 _cam_y = _mm_set1_ps(camera->pos.y);
    __m128 _cam_z = _mm_set1_ps(camera->pos.z);

    for (int i = 0; i < count; ++i) {
        Vertex4 _vertex = in[i];

        // translate position in world relative to camera
        _vertex.x = _mm_sub_ps(_vertex.x, _cam_x);
        _vertex.y = _mm_sub_ps(_vertex.y, _cam_y);
        _vertex.z = _mm_sub_ps(_vertex.z, _cam_z);

        //rotate world relative to comara
        __m128 _cross_x = _mm_sub_ps(_mm_mul_ps(_qy, _vertex.z), _mm_mul_ps(_qz, _vertex.y));
        __m128 _cross_y = _mm_sub_ps(_mm_mul_ps(_qz, _vertex.x), _mm_mul_ps(_qx, _vertex.z));
        __m128 _cross_z = _mm_sub_ps(_mm_mul_ps(_qx, _vertex.y), _mm_mul_ps(_qy, _vertex.x));

        __m128 _scaled_x = _mm_mul_ps(_vertex.x, _qw);
        __m128 _scaled_y = _mm_mul_ps(_vertex.y, _qw);
        __m128 _scaled_z = _mm_mul_ps(_vertex.z, _qw);
        

        __m128 _dot = _mm_add_ps(_mm_add_ps(_mm_mul_ps(_vertex.x, _qx), _mm_mul_ps(_vertex.y, _qy)),
                                 _mm_mul_ps(_vertex.z, _qz));

        __m128 _tmp_w = _mm_mul_ps(_mm_set1_ps(-1.0f), _dot);
        __m128 _tmp_x = _mm_add_ps(_cross_x, _scaled_x);
        __m128 _tmp_y = _mm_add_ps(_cross_y, _scaled_y);
        __m128 _tmp_z = _mm_add_ps(_cross_z, _scaled_z);


        __m128 _conj_w = _qw;
        __m128 _conj_x = _mm_mul_ps(_mm_set1_ps(-1.0f), _qx);
        __m128 _conj_y = _mm_mul_ps(_mm_set1_ps(-1.0f), _qy);
        __m128 _conj_z = _mm_mul_ps(_mm_set1_ps(-1.0f), _qz);

        __m128 _tmp2_x = _mm_add_ps(_mm_mul_ps(_tmp_w, _conj_x), _mm_mul_ps(_conj_w, _tmp_x));
        __m128 _tmp2_y = _mm_add_ps(_mm_mul_ps(_tmp_w, _conj_y), _mm_mul_ps(_conj_w, _tmp_y));
        __m128 _tmp2_z = _mm_add_ps(_mm_mul_ps(_tmp_w, _conj_z), _mm_mul_ps(_conj_w, _tmp_z));

        __m128 _tmp_cross_x = _mm_sub_ps(_mm_mul_ps(_tmp_y, _conj_z), _mm_mul_ps(_tmp_z, _conj_y));
        __m128 _tmp_cross_y = _mm_sub_ps(_mm_mul_ps(_tmp_z, _conj_x), _mm_mul_ps(_tmp_x, _conj_z));
        __m128 _tmp_cross_z = _mm_sub_ps(_mm_mul_ps(_tmp_x, _conj_y), _mm_mul_ps(_tmp_y, _conj_x));

        __m128 _rotated_x = _mm_add_ps(_tmp_cross_x, _tmp2_x);
        __m128 _rotated_y = _mm_add_ps(_tmp_cross_y, _tmp2_y);
        __m128 _rotated_z = _mm_add_ps(_tmp_cross_z, _tmp2_z);


        __m128 _cam_znear = _mm_set1_ps(camera->znear);
        __m128 _cam_zfar = _mm_set1_ps(camera->zfar);

        __m128 _length_x = _mm_set1_ps(tanf(camera->fov/2.0f)*camera->znear*aspect_ratio);
        __m128 _length_y = _mm_set1_ps(tanf(camera->fov/2.0f)*camera->znear);

        __m128 _mask = _mm_cmplt_ps(_rotated_z, _cam_znear);
        _mask = _mm_or_ps(_mask, _mm_cmpgt_ps(_rotated_z, _cam_zfar));

        __m128 _normal_x = _mm_div_ps(_mm_mul_ps(_rotated_x, _mm_div_ps(_cam_znear, _rotated_z)),
                                      _length_x);

        __m128 _normal_y = _mm_div_ps(_mm_mul_ps(_rotated_y, _mm_div_ps(_cam_znear, _rotated_z)),
                                      _length_y);

        __m128 _buffer_x = _mm_set1_ps(screen_width/2.0f);
        __m128 _buffer_y = _mm_set1_ps(screen_height/2.0f);

        __m128 _float_res_x = _mm_add_ps(_mm_mul_ps(_normal_x, _buffer_x), _buffer_x);
        __m128 _float_res_y = _mm_add_ps(_mm_mul_ps(_mm_mul_ps(_normal_y, _mm_set1_ps(-1.0f)), _buffer_y), _buffer_y);

        out[i].x = _mm_or_si128(_mm_cvtps_epi32(_float_res_x), *(__m128i *)&_mask);
        out[i].y = _mm_or_si128(_mm_cvtps_epi32(_float_res_y), *(__m128i *)&_mask);
    }
}

static V3 renderer_world_vertex_to_view(V3 world_pos, GameCamera *camera) {
    V3 result = v3_rotate_q4(world_pos-camera->pos, 
                          (rotation_q4(-camera->theta_x, v3(1,0,0))*
                          rotation_q4(-camera->theta_y, v3(0,1,0))));

    return result;
}

static bool renderer_v3_should_clip(V3 pos, 
                                    GameCamera *camera, 
                                    float aspect_ratio) {

    float max_x = tanf(camera->fov/2.0f)*camera->znear * (pos.z/camera->znear) * 
        aspect_ratio;
    float max_y = tanf(camera->fov/2.0f)*camera->znear * (pos.z/camera->znear);

    if (pos.z <= camera->znear) return true;
    if (pos.z > camera->zfar) return true;
    if (fabs(pos.x) > max_x) return true;
    if (fabs(pos.y) > max_y) return true;

    return false;
}

static V3Screen renderer_world_vertex_to_screen(V3 world_pos, 
                                                GameCamera *camera,
                                                int buffer_width, 
                                                int buffer_height) {

    V3Screen result;

    V3 relative_pos = v3_rotate_q4(world_pos-camera->pos, 
            (rotation_q4(-camera->theta_x, v3(1,0,0))*
             rotation_q4(-camera->theta_y, v3(0,1,0))));

    float aspect_ratio = ((float)buffer_width/(float)buffer_height);

    if (renderer_v3_should_clip(relative_pos, camera, aspect_ratio)) {
        result.x = 0xFFFFFFFF;
        result.y = 0xFFFFFFFF;
        return result;
    }

    float normal_x = (relative_pos.x / relative_pos.z*tanf(camera->fov/2.0f)*aspect_ratio*camera->znear); 

    float normal_y = (relative_pos.y / relative_pos.z*tanf(camera->fov/2.0f)*camera->znear); 

    result.x = (int)(normal_x * (buffer_width/2.0f) + buffer_width/2.0f);
    result.y = (int)(-normal_y * (buffer_height/2.0f) + buffer_height/2.0f);
    return result;
}

static void renderer_transform_light_and_cull(RendererState *render_state,
                                              GameCamera *camera, 
                                              int buffer_width, 
                                              int buffer_height) {

    V3 light_source = v3(0, 10, -30);
    for (int triangle = 0; triangle < render_state->polygon_count; ++triangle) {
        Triangle *current = &render_state->polygons[triangle];
#ifndef DEBUG_SHADING

        float light_inten = v3_dot(v3_norm(light_source-(render_state->vertex_list[current->v1])), 
                render_state->vertexn_list[current->vn1]);

        if (light_inten < 0) {
            current->v1_color = {};
        }

        else {
            current->v1_color = v3(light_inten, light_inten, light_inten);
        }

        light_inten = v3_dot(v3_norm(light_source-(render_state->vertex_list[current->v2])), 
                render_state->vertexn_list[current->vn2]);

        if (light_inten < 0) {
            current->v2_color = {};
        }

        else {
            current->v2_color = v3(light_inten, light_inten, light_inten);
        }

        light_inten = v3_dot(v3_norm(light_source-(render_state->vertex_list[current->v3])), 
                render_state->vertexn_list[current->vn3]);

        if (light_inten < 0) {
            current->v3_color = {};
        }

        else {
            current->v3_color = v3(light_inten, light_inten, light_inten);
        }
#else 
        current->v1_color = v3(1,0,0);
        current->v2_color = v3(0,1,0);
        current->v3_color = v3(0,0,1);
#endif // !DEBUG_SHADING
    }

    for (int vertex = 0; vertex < render_state->vertex_count; ++vertex) {
        render_state->vertex_list[vertex] = 
            renderer_world_vertex_to_view(render_state->vertex_list[vertex], camera);
    }

    render_state->draw_count = 0;
    for (int triangle = 0; triangle < render_state->polygon_count; ++triangle) {
        Triangle current = render_state->polygons[triangle];
        if (v3_dot(render_state->vertex_list[current.v1], 
                v3_cross(render_state->vertex_list[current.v2]-render_state->vertex_list[current.v1], 
                    render_state->vertex_list[current.v3]-render_state->vertex_list[current.v1])) >= 0) {
            continue;
        }

        render_state->polygons_to_draw[render_state->draw_count++] = current;
    }


    float aspect_ratio = ((float)buffer_width/(float)buffer_height);
    for (int vertex = 0; vertex < render_state->vertex_count; ++vertex) {
        V3 relative_pos = render_state->vertex_list[vertex];

        if (renderer_v3_should_clip(relative_pos, camera, aspect_ratio)) {
            render_state->screen_vertices[vertex].x = 0xFFFFFFFF;
            render_state->screen_vertices[vertex].y = 0xFFFFFFFF;
            render_state->screen_vertices[vertex].z = FLT_MAX;
            continue;
        }
        
        float normal_x = (relative_pos.x/relative_pos.z) / 
                         (tanf(camera->fov/2.0f)*aspect_ratio);

        float normal_y = (relative_pos.y/relative_pos.z) / 
                         (tanf(camera->fov/2.0f));

        int x = (int)(normal_x * (buffer_width/2.0f) + buffer_width/2.0f);
        int y = (int)(-normal_y * (buffer_height/2.0f) + buffer_height/2.0f);
        render_state->screen_vertices[vertex].x = x;
        render_state->screen_vertices[vertex].y = y;
        render_state->screen_vertices[vertex].z = relative_pos.z;
    }
}

static inline bool renderer_v2screen_isinvalid(V3Screen screen,
                                          int screen_width,
                                          int screen_height) {

    if (screen.x == 0xFFFFFFFF || screen.y == 0xFFFFFFFF) return true;
    if (screen.x < 0 || screen.x > screen_width) return true;
    if (screen.y < 0 || screen.y > screen_height) return true;

    return false;
}

static void renderer_draw_background(OffscreenBuffer *buffer, uint32_t color) {
    uint8_t *row = (uint8_t *)buffer->memory;
    for (int y = 0;y < buffer->height; ++y) {
        uint32_t *pixel = (uint32_t *)row;

        for (int x =0;x < buffer->width; ++x) {
            *pixel++ = color;
        }

        row += buffer->pitch;
    }
}

static void renderer_draw_line(OffscreenBuffer *buffer,
                               V3Screen start,
                               V3Screen end,
                               uint32_t color) {
            

    if (renderer_v2screen_isinvalid(start, buffer->width, buffer->height)) {
        return;
    }

    if (renderer_v2screen_isinvalid(end, buffer->width, buffer->height)) {
        return;
    }

    if (start.x > end.x) {
        V3Screen temp = end;
        end = start;
        start = temp;
    }

    int dx = end.x - start.x;
    int dy = end.y - start.y;
    float dz = end.z - start.z;
    uint32_t *endofbuffer = (uint32_t *)((uint8_t *)buffer->memory + (buffer->width*buffer->height*buffer->bytes_per_pixel));

    // case where the gradient is below 1;
    int step = 0;
    if (abs(dx) > abs(dy)) {
        step = abs(dx);
    }

    else {
        step = abs(dy);
    }

    float x = (float)start.x;
    float y = (float)start.y;

    float x_inc = (float)dx/(float)step;
    float y_inc = (float)dy/(float)step;

    for (int draw = 0; draw < step; ++draw) {
        int offset = (int)(x)*buffer->bytes_per_pixel + 
                     (int)(y)*buffer->pitch;

            uint32_t *pixel = (uint32_t *)((uint8_t *)buffer->memory + offset);
            if ((pixel < endofbuffer) && (pixel >= buffer->memory)) {
                *pixel = color;
        }

        x += x_inc;
        y += y_inc;
    }
}

static void renderer_draw_line_zcull(OffscreenBuffer *buffer,
                                     OffscreenBuffer *zbuffer,
                                     V3Screen start,
                                     V3Screen end,
                                     V3 start_color,
                                     V3 end_color) {
            

    if (renderer_v2screen_isinvalid(start, buffer->width, buffer->height)) {
        return;
    }

    if (renderer_v2screen_isinvalid(end, buffer->width, buffer->height)) {
        return;
    }

    if (start.x > end.x) {
        V3Screen temp = end;
        end = start;
        start = temp;
    }

    int dx = end.x - start.x;
    int dy = end.y - start.y;
    float dz = end.z - start.z;
    uint32_t *endofbuffer = (uint32_t *)((uint8_t *)buffer->memory + (buffer->width*buffer->height*buffer->bytes_per_pixel));

    // case where the gradient is below 1;
    int step = 0;
    if (abs(dx) > abs(dy)) {
        step = abs(dx);
    }

    else {
        step = abs(dy);
    }

    float x = (float)start.x;
    float y = (float)start.y;
    float z = start.z;
    float sr = start_color.x;
    float sg = start_color.y;
    float sb = start_color.z;

    float er = end_color.x;
    float eg = end_color.y;
    float eb = end_color.z;

    float dr = er - sr;
    float dg = eg - sg;
    float db = eb - sb;

    float r_inc = dr/(float)step;
    float g_inc = dg/(float)step;
    float b_inc = db/(float)step;

    float x_inc = (float)dx/(float)step;
    float y_inc = (float)dy/(float)step;
    float z_inc = (float)dz/(float)step;

    for (int draw = 0; draw < step; ++draw) {
        int offset = (int)x*buffer->bytes_per_pixel + 
                     (int)y*buffer->pitch;

        float *depth_value = (float *)((uint8_t *)zbuffer->memory + offset);
    
            if ((z < *depth_value)) {
                uint32_t color = (((uint8_t)(0xFF*sr)) << 16) | 
                                 (((uint8_t)(0xFF*sg))) << 8 | 
                                 (uint8_t)(0xFF*sb);

                uint32_t *pixel = (uint32_t *)((uint8_t *)buffer->memory + offset);
                *pixel = color;
                *depth_value = z;
        }

        x += x_inc;
        y += y_inc;
        z += z_inc;

        sr = sr + r_inc;
        sg = sg + g_inc;
        sb = sb + b_inc;
    }
}


static void renderer_draw_triangle_wireframe(OffscreenBuffer *buffer,
                           V3Screen v1, 
                           V3Screen v2, 
                           V3Screen v3,
                           uint32_t color) {

    renderer_draw_line(buffer, v1, v2, color);
    renderer_draw_line(buffer, v2, v3, color);
    renderer_draw_line(buffer, v3, v1, color);
}

static void renderer_v2screen4_draw_triangles_wireframe(OffscreenBuffer *buffer,
                                                      V3Screen4 *vertices,
                                                      Triangle *triangles,
                                                      int count) {

    for (int i = 0; i < count; ++i) {
        int index1 = triangles->v1 / 4;
        int index2 = triangles->v3 / 4;
        int index3 = triangles->v2 / 4;

        int offset1 = triangles->v1 % 4;
        int offset2 = triangles->v3 % 4;
        int offset3 = triangles->v2 % 4;
        V3Screen vertex1 = {((int *)&vertices[index1].x)[offset1], ((int* )&vertices[index1].y)[offset1]};
        V3Screen vertex2 = {((int *)&vertices[index2].x)[offset2], ((int *)&vertices[index2].y)[offset2]};
        V3Screen vertex3 = {((int *)&vertices[index3].x)[offset3], ((int *)&vertices[index3].y)[offset3]};
//        renderer_draw_triangle_wireframe(buffer, vertex1, vertex2, vertex3, triangles->v1_color);
        triangles++;
    }
}

static void renderer_draw_flat_top_triangle(OffscreenBuffer *buffer,
                                            OffscreenBuffer *zbuffer,
                                            V3Screen vert1,
                                            V3Screen vert2,
                                            V3Screen vert3,
                                            V3 v1_color,
                                            V3 v2_color,
                                            V3 v3_color) {
    // vertices passed in ascending order
        
    float steps = (float)(vert1.y - vert3.y);
    float x_inc_left = (float)(vert1.x-vert2.x)/steps;
    float x_inc_right = (float)(vert1.x-vert3.x)/steps;

    float x_left = (float)vert2.x;
    float x_right = (float)vert3.x;
    float z_left = vert2.z;
    float z_right = vert3.z;

    float rl = v2_color.x;
    float gl = v2_color.y;
    float bl = v2_color.z;
    float rr = v3_color.x;
    float gr = v3_color.y;
    float br = v3_color.z; 

    float rt = v1_color.x;
    float gt = v1_color.y;
    float bt = v1_color.z; 

    float drl = rt - rl;
    float dgl = gt - gl;
    float dbl = bt - bl;
    float drr = rt - rr;
    float dgr = gt - gr;
    float dbr = bt - br;

    float rl_inc = drl/steps;
    float gl_inc = dgl/steps;
    float bl_inc = dbl/steps;

    float rr_inc = drr/steps;
    float gr_inc = dgr/steps;
    float br_inc = dbr/steps;

    float dzl = vert1.z - vert2.z;
    float dzr = vert1.z - vert3.z;
    float zl_inc = dzl / steps;
    float zr_inc = dzr / steps;

    for (int y = vert2.y; y <= vert1.y; ++y) {
        V3Screen left;
        left.x = (int)x_left;
        left.y = y;
        left.z = z_left;

        V3Screen right;
        right.x = (int)x_right;
        right.y = y;
        right.z = z_right;

        V3 color_left = v3(rl, gl, bl);
        V3 color_right = v3(rr, gr, br);

        renderer_draw_line_zcull(buffer, zbuffer, left, right, color_left, color_right);

        x_left += x_inc_left;
        x_right += x_inc_right;
        z_left += zl_inc;
        z_right += zr_inc;
        rl = rl + rl_inc;
        gl = gl + gl_inc;
        bl = bl + bl_inc;

        rr = rr + rr_inc;
        gr = gr + gr_inc;
        br = br + br_inc;

    }
}

static void renderer_draw_flat_bottom_triangle(OffscreenBuffer *buffer,
                                               OffscreenBuffer *zbuffer,
                                               V3Screen vert1,
                                               V3Screen vert2,
                                               V3Screen vert3,
                                               V3 v1_color,
                                               V3 v2_color,
                                               V3 v3_color) {
    // vertices passed in ascending order

    float steps = (float)(vert3.y - vert2.y);
    float x_inc_left = (float)(vert3.x-vert2.x)/(steps);
    float x_inc_right = (float)(vert3.x-vert1.x)/(steps);

    float x_left = (float)vert2.x;
    float x_right = (float)vert1.x;
    float z_left = vert2.z;
    float z_right = vert1.z;

    float rl = v2_color.x;
    float gl = v2_color.y;
    float bl = v2_color.z;
    float rr = v1_color.x;
    float gr = v1_color.y;
    float br = v1_color.z; 

    float rt = v3_color.x;
    float gt = v3_color.y;
    float bt = v3_color.z; 

    float drl = rt - rl;
    float dgl = gt - gl;
    float dbl = bt - bl;
    float drr = rt - rr;
    float dgr = gt - gr;
    float dbr = bt - br;

    float dzl = vert3.z - vert2.z;
    float dzr = vert3.z - vert1.z;
    float zl_inc = dzl / steps;
    float zr_inc = dzr / steps;

    float rl_inc = drl/steps;
    float gl_inc = dgl/steps;
    float bl_inc = dbl/steps;

    float rr_inc = drr/steps;
    float gr_inc = dgr/steps;
    float br_inc = dbr/steps;
    
    for (int y = vert2.y; y >= vert3.y; --y) {
        V3Screen left;
        left.x = (int)x_left;
        left.y = y;
        left.z = z_left;

        V3Screen right;
        right.x = (int)x_right;
        right.y = y;
        right.z = z_right;

        V3 color_left = v3(rl, gl, bl);
        V3 color_right = v3(rr, gr, br);

        renderer_draw_line_zcull(buffer, zbuffer, left, right, color_left, color_right);

        x_left -= x_inc_left;
        x_right -= x_inc_right;
        z_left -= zl_inc;
        z_right -= zr_inc;

        rl = rl - rl_inc;
        gl = gl - gl_inc;
        bl = bl - bl_inc;

        rr = rr - rr_inc;
        gr = gr - gr_inc;
        br = br - br_inc;
    }
}


static void renderer_draw_triangles_filled(OffscreenBuffer *buffer,
                                           OffscreenBuffer *zbuffer,
                                           V3Screen *vertices,
                                           Triangle *triangles,
                                           int count) {

    for (int i = 0; i < count; ++i) {
        V3Screen vert1 = vertices[triangles[i].v1];
        V3Screen vert2 = vertices[triangles[i].v2];
        V3Screen vert3 = vertices[triangles[i].v3];
        V3 v1_color = triangles[i].v1_color;
        V3 v2_color = triangles[i].v2_color;
        V3 v3_color = triangles[i].v3_color;

        if (renderer_v2screen_isinvalid(vert1, buffer->width, buffer->height)) {
            continue;
        }

        if (renderer_v2screen_isinvalid(vert2, buffer->width, buffer->height)) {
            continue;
        }

        if (renderer_v2screen_isinvalid(vert3, buffer->width, buffer->height)) {
            continue;
        }

        //sort vertexes v3 is the biggest
        if (vert3.y > vert2.y) {
            V3Screen temp = vert2;
            vert2 = vert3;
            vert3 = temp;
            V3 color_temp = v2_color;
            v2_color = v3_color;
            v3_color = color_temp;

        }

        if (vert2.y > vert1.y) {
            V3Screen temp = vert1;
            vert1 = vert2;
            vert2 = temp;
            V3 color_temp = v1_color;
            v1_color = v2_color;
            v2_color = color_temp;
        }

        if (vert3.y > vert2.y) {
            V3Screen temp = vert2;
            vert2 = vert3;
            vert3 = temp;
            V3 color_temp = v2_color;
            v2_color = v3_color;
            v3_color = color_temp;
        }

        if (vert3.y == vert2.y) {
            renderer_draw_flat_top_triangle(buffer, zbuffer, vert1, vert2, vert3, v1_color, v2_color, v3_color);
            continue;
        }

        if (vert1.y == vert2.y) {
            renderer_draw_flat_bottom_triangle(buffer, zbuffer, vert1, vert2, vert3, v1_color, v2_color, v3_color);
            continue;
        }

        V3Screen vert4;
        vert4.x = (int)(vert3.x - (((float)(vert3.y-vert2.y)/(float)(vert3.y-vert1.y))*(float)(vert3.x-vert1.x)));
        vert4.y = vert2.y;
        vert4.z = ((vert1.z - vert2.z)/2) + vert2.z; 

        float v1r = v1_color.x;
        float v1g = v1_color.y;
        float v1b = v1_color.z; 


        float v3r = v3_color.x;
        float v3g = v3_color.y;
        float v3b = v3_color.z; 

        float dr = v3r-v1r;
        float dg = v3g-v1g;
        float db = v3b-v1b;

        float deltax = (float)vert3.y-vert1.y;
        float dx = (float)vert2.y-vert1.y;

        float v4r = v1r + (dr/deltax)*(dx);
        float v4g = v1g + (dg/deltax)*(dx);
        float v4b = v1b + (db/deltax)*(dx);


        V3 v4_color = v3(v4r, v4g, v4b);
        //renderer_draw_flat_top_triangle(buffer, zbuffer, vert1, vert2, vert4, v1_color, v2_color, v4_color);
        //renderer_draw_flat_bottom_triangle(buffer, zbuffer, vert2, vert4, vert3, v2_color, v4_color, v3_color);
    }

}

#if 0
static void renderer_v2screen4_draw_triangles_filled(OffscreenBuffer *buffer,
                                                    V3Screen4 *vertices,
                                                    Triangle *triangles,
                                                    uint32_t *colors,
                                                    int colors_size,
                                                    int tricount) {

    int color_index = 0;
    for (int i = 0; i < tricount; ++i) {
        int index1 = triangles->v1 / 4;
        int index2 = triangles->v3 / 4;
        int index3 = triangles->v2 / 4;

        int offset1 = triangles->v1 % 4;
        int offset2 = triangles->v3 % 4;
        int offset3 = triangles->v2 % 4;
        V3Screen vertex1 = {((int *)&vertices[index1].x)[offset1], ((int *)&vertices[index1].y)[offset1]};
        V3Screen vertex2 = {((int *)&vertices[index2].x)[offset2], ((int *)&vertices[index2].y)[offset2]};
        V3Screen vertex3 = {((int *)&vertices[index3].x)[offset3], ((int *)&vertices[index3].y)[offset3]};
        renderer_draw_triangles_filled(buffer, vertex1, vertex2, vertex3, colors[color_index]);
        if (i % 2) {
            color_index++;
        }
        if (color_index >= colors_size) {
            color_index = 0;
        }
        ++triangles;
    }
}
#endif
