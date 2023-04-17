#include "renderer.h"
#include <emmintrin.h>
#include <stdint.h>

static void renderer_vertex4_to_v2screen(Vertex4 *in, 
                         GameCamera *camera, 
                         int screen_width,
                         int screen_height,
                         int count, 
                         V2Screen4 *out) {

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

static V3 renderer_world_vertex_to_view(V3 world_pos, GameCamera *camera, 
                                  int buffer_width,
                                  int buffer_height) {
    V3 result;
    result = v3_rotate_q4(world_pos-camera->pos, 
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

    if (pos.z < camera->znear) return true;
    if (pos.z > camera->zfar) return true;

    return false;
}

static V2Screen renderer_world_vertex_to_screen(V3 world_pos, 
                                                GameCamera *camera,
                                                int buffer_width, 
                                                int buffer_height) {

    V2Screen result;

    V3 relative_pos = v3_rotate_q4(world_pos-camera->pos, 
            (rotation_q4(-camera->theta_x, v3(1,0,0))*
             rotation_q4(-camera->theta_y, v3(0,1,0))));

    float aspect_ratio = ((float)buffer_width/(float)buffer_height);

    if (renderer_v3_should_clip(relative_pos, camera, aspect_ratio)) {
        result.x = 0xFFFFFFFF;
        result.y = 0xFFFFFFFF;
        return result;
    }

    float normal_x = (relative_pos.x*(camera->znear/relative_pos.z)) / 
        ((tanf(camera->fov/2.0f)*camera->znear)*aspect_ratio);

    float normal_y = (relative_pos.y*(camera->znear/relative_pos.z)) / 
        (tanf(camera->fov/2.0f)*camera->znear);

    result.x = (int)(normal_x * (buffer_width/2.0f) + buffer_width/2.0f);
    result.y = (int)(-normal_y * (buffer_height/2.0f) + buffer_height/2.0f);
    return result;
}

static void renderer_world_vertices_to_screen_and_cull(V3 *in_vertices, 
                                                       int v_count, 
                                                       Triangle *in_triangles,
                                                       int tri_count,
                                                       GameCamera *camera, 
                                                       int buffer_width, 
                                                       int buffer_height,
                                                       Triangle *out_triangles,
                                                       int *tri_out_count,
                                                       V2Screen *out_vertices,
                                                       int *out_vertices_count) {



    for (int vertex = 0; vertex < v_count; ++vertex) {
        in_vertices[vertex] = v3_rotate_q4(in_vertices[vertex]-camera->pos, 
                (rotation_q4(-camera->theta_x, v3(1,0,0))*
                 rotation_q4(-camera->theta_y, v3(0,1,0))));
    }

    *tri_out_count = 0;
    for (int triangle = 0; triangle < tri_count; triangle++) {
        Triangle current = in_triangles[triangle];

        if (v3_dot(in_vertices[current.v1]-camera->pos, 
                    v3_cross(in_vertices[current.v2]-in_vertices[current.v1], 
                        in_vertices[current.v3]-in_vertices[current.v1])) >= 0) {
            continue;
        }
        out_triangles[(*tri_out_count)++] = current;
    }

    *out_vertices_count = 0;
    float aspect_ratio = ((float)buffer_width/(float)buffer_height);
    for (int i = 0; i < v_count; ++i) {
        V3 relative_pos = in_vertices[i];

        if (renderer_v3_should_clip(relative_pos, camera, aspect_ratio)) {
            continue;
        }
        
        float normal_x = (relative_pos.x*(camera->znear/relative_pos.z)) / 
            ((tanf(camera->fov/2.0f)*camera->znear)*aspect_ratio);

        float normal_y = (relative_pos.y*(camera->znear/relative_pos.z)) / 
            (tanf(camera->fov/2.0f)*camera->znear);

        out_vertices[*out_vertices_count].x = (int)(normal_x * (buffer_width/2.0f) + buffer_width/2.0f);
        out_vertices[*out_vertices_count].y = (int)(-normal_y * (buffer_height/2.0f) + buffer_height/2.0f);
        (*out_vertices_count)++;

    }
}

static inline bool renderer_v2screen_isinvalid(V2Screen screen,
                                          int screen_width,
                                          int screen_height) {

    if (screen.x == 0xFFFFFFFF || screen.y == 0xFFFFFFFF) return true;
    if (screen.x < 0 || screen.x > screen_width) return true;
    if (screen.y < 0 || screen.y > screen_height) return true;

    return false;
}

static void renderer_draw_background(OffscreenBuffer *buffer) {
    uint8_t *row = (uint8_t *)buffer->memory;
    for (int y = 0;y < buffer->height; ++y) {
        uint32_t *pixel = (uint32_t *)row;

        for (int x =0;x < buffer->width; ++x) {
            *pixel++ = 0xFF000000;
        }

        row += buffer->pitch;
    }
}

static void renderer_draw_line(OffscreenBuffer *buffer,
                 V2Screen start,
                 V2Screen end,
                 uint32_t color) {
    

    if (renderer_v2screen_isinvalid(start, buffer->width, buffer->height)) {
        return;
    }

    if (renderer_v2screen_isinvalid(end, buffer->width, buffer->height)) {
        return;
    }

    if (start.x > end.x) {
        V2Screen temp = end;
        end = start;
        start = temp;
    }

    int dx = end.x - start.x;
    int dy = end.y - start.y;
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
        uint32_t *pixel = (uint32_t *)((uint8_t *)buffer->memory + 
                (int)roundf(x)*buffer->bytes_per_pixel +
                (int)roundf(y)*buffer->pitch);
        if ((pixel < endofbuffer) && (pixel >= buffer->memory)) {
            *pixel = color;
        }

        x += x_inc;
        y += y_inc;
    }
}

static void renderer_transform_and_draw_line(OffscreenBuffer *buffer,
                 V3 v_start,
                 V3 v_end,
                 GameCamera *camera,
                 uint32_t color) {

    V2Screen start = renderer_world_vertex_to_screen(v_start, 
                                                     camera, 
                                                     buffer->width, 
                                                     buffer->height);

    V2Screen end = renderer_world_vertex_to_screen(v_end, 
                                                   camera, 
                                                   buffer->width, 
                                                   buffer->height);

    if (renderer_v2screen_isinvalid(start, buffer->width, buffer->height)) {
        return;
    }

    if (renderer_v2screen_isinvalid(end, buffer->width, buffer->height)) {
        return;
    }

    int dx = end.x - start.x;
    int dy = end.y - start.y;
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
        uint32_t *pixel = (uint32_t *)((uint8_t *)buffer->memory + 
                (int)roundf(x)*buffer->bytes_per_pixel +
                (int)roundf(y)*buffer->pitch);
        if ((pixel < endofbuffer) && (pixel >= buffer->memory)) {
            *pixel = color;
        }

        x += x_inc;
        y += y_inc;
    }
}

static void renderer_draw_triangle_wireframe(OffscreenBuffer *buffer,
                           V2Screen v1, 
                           V2Screen v2, 
                           V2Screen v3,
                           uint32_t color) {

    renderer_draw_line(buffer, v1, v2, color);
    renderer_draw_line(buffer, v2, v3, color);
    renderer_draw_line(buffer, v3, v1, color);
}

static void renderer_v2screen4_draw_triangles_wireframe(OffscreenBuffer *buffer,
                      V2Screen4 *vertices,
                      Triangle *triangles,
                      uint32_t color,
                      int count) {

    for (int i = 0; i < count; ++i) {
        int index1 = triangles->v1 / 4;
        int index2 = triangles->v3 / 4;
        int index3 = triangles->v2 / 4;

        int offset1 = triangles->v1 % 4;
        int offset2 = triangles->v3 % 4;
        int offset3 = triangles->v2 % 4;
        V2Screen vertex1 = {((int *)&vertices[index1].x)[offset1], ((int* )&vertices[index1].y)[offset1]};
        V2Screen vertex2 = {((int *)&vertices[index2].x)[offset2], ((int *)&vertices[index2].y)[offset2]};
        V2Screen vertex3 = {((int *)&vertices[index3].x)[offset3], ((int *)&vertices[index3].y)[offset3]};
        renderer_draw_triangle_wireframe(buffer, vertex1, vertex2, vertex3, color);
        triangles++;
    }
}

static void renderer_draw_flat_top_triangle(OffscreenBuffer *buffer,
                                      V2Screen v1,
                                      V2Screen v2,
                                      V2Screen v3,
                                      uint32_t color) {
        
    float x_inc_left = (float)(v1.x-v2.x)/(float)(v1.y-v2.y);
    float x_inc_right = (float)(v1.x-v3.x)/(float)(v1.y-v3.y);

    float x_left = (float)v2.x;
    float x_right = (float)v3.x;

    for (int y = v2.y; y <= v1.y; ++y) {
        V2Screen left;
        left.x = (int)x_left;
        left.y = y;

        V2Screen right;
        right.x = (int)x_right;
        right.y = y;
        renderer_draw_line(buffer, left, right, color);

        x_left += x_inc_left;
        x_right += x_inc_right;

    }
}

static void renderer_draw_flat_bottom_triangle(OffscreenBuffer *buffer,
                                   V2Screen v1,
                                   V2Screen v2,
                                   V2Screen v3,
                                   uint32_t color) {
    float x_inc_left = (float)(v3.x-v2.x)/(float)(v3.y-v2.y);
    float x_inc_right = (float)(v3.x-v1.x)/(float)(v3.y-v1.y);

    float x_left = (float)v2.x;
    float x_right = (float)v1.x;
    
    for (int y = v2.y; y >= v3.y; --y) {
        V2Screen left;
        left.x = (int)x_left;
        left.y = y;

        V2Screen right;
        right.x = (int)x_right;
        right.y = y;
        renderer_draw_line(buffer, left, right, color);

        x_left -= x_inc_left;
        x_right -= x_inc_right;
    }
}


static void renderer_draw_triangles_filled(OffscreenBuffer *buffer,
                                           V2Screen *vertices,
                                           Triangle *triangles,
                                           int count) {

    for (int i = 0; i < count; ++i) {
        V2Screen v1 = vertices[triangles[i].v1];
        V2Screen v2 = vertices[triangles[i].v2];
        V2Screen v3 = vertices[triangles[i].v3];
        uint32_t color = triangles[i].color;

        if (renderer_v2screen_isinvalid(v1, buffer->width, buffer->height)) {
            continue;
        }

        if (renderer_v2screen_isinvalid(v2, buffer->width, buffer->height)) {
            continue;
        }

        if (renderer_v2screen_isinvalid(v3, buffer->width, buffer->height)) {
            continue;
        }

        //sort vertexes 
        if (v3.y > v2.y) {
            V2Screen temp = v2;
            v2 = v3;
            v3 = temp;
        }

        if (v2.y > v1.y) {
            V2Screen temp = v1;
            v1 = v2;
            v2 = temp;
        }

        if (v3.y > v2.y) {
            V2Screen temp = v2;
            v2 = v3;
            v3 = temp;
        }

        if (v3.y == v2.y) {
            renderer_draw_flat_top_triangle(buffer, v1, v2, v3, color);
            continue;
        }

        if (v1.y == v2.y) {
            renderer_draw_flat_bottom_triangle(buffer, v1, v2, v3, color);
            continue;
        }

        V2Screen v4;
        v4.x = (int)(v3.x - (((float)(v3.y-v2.y)/(float)(v3.y-v1.y))*(float)(v3.x-v1.x)));
        v4.y = v2.y;
        renderer_draw_flat_top_triangle(buffer, v1, v2, v4, color);
        renderer_draw_flat_bottom_triangle(buffer, v2, v4, v3, color);
    }

}

#if 0
static void renderer_v2screen4_draw_triangles_filled(OffscreenBuffer *buffer,
                                                    V2Screen4 *vertices,
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
        V2Screen vertex1 = {((int *)&vertices[index1].x)[offset1], ((int *)&vertices[index1].y)[offset1]};
        V2Screen vertex2 = {((int *)&vertices[index2].x)[offset2], ((int *)&vertices[index2].y)[offset2]};
        V2Screen vertex3 = {((int *)&vertices[index3].x)[offset3], ((int *)&vertices[index3].y)[offset3]};
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
