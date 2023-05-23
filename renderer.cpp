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

    result.x = (int)roundf((normal_x * (buffer_width/2.0f) + buffer_width/2.0f));
    result.y = (int)roundf((-normal_y * (buffer_height/2.0f) + buffer_height/2.0f));
    return result;
}

static V3 compute_light_intensity(LightSource source, V3 vertex, V3 normal) {
    V3 dl = source.position-vertex; 
    float intensity = v3_dot(v3_norm(dl), normal);
    if (intensity > 0) {
        return intensity*source.color;
    }
    return source.color = {};
}

static bool vertex_is_light_source(uint32_t vertex_index, LightSource *source, int count) {
    for (int i = 0; i < count; ++i) {
        if (vertex_index >= source->obj->vstart && 
                vertex_index <= source->obj->vend) {
            return true;
        }
        ++source;
    }
    return false;
}

static inline float rgb_to_luminance(V3 rgb) {
    return v3_dot(rgb, v3(0.2126f, 0.7152f, 0.0722f));
}

static V3 tone_map_light(V3 light, V3 white_point) {
    float luminance = rgb_to_luminance(light);
    float white_l = rgb_to_luminance(white_point);
    float new_luminance = luminance / (1.0f + luminance);
    V3 result = (new_luminance/luminance)*light;
    return result;
}

static void renderer_transform_light_and_cull(RendererState *render_state,
                                              GameCamera *camera, 
                                              int buffer_width, 
                                              int buffer_height) {

    render_state->draw_count = 0;
    for (uint32_t triangle = 0; triangle < render_state->polygon_count; ++triangle) {
        Triangle *current = &render_state->polygons[triangle];
        Triangle out = *current;
        V3 vertex1 = render_state->vertex_list[current->v1];
        V3 vertex2 = render_state->vertex_list[current->v2];
        V3 vertex3 = render_state->vertex_list[current->v3];
        V3 *v1_color = &render_state->vertex_colors[current->v1];
        V3 *v2_color = &render_state->vertex_colors[current->v2];
        V3 *v3_color = &render_state->vertex_colors[current->v3];
        V3 vn1 = render_state->vertexn_list[current->vn1];
        V3 vn2 = render_state->vertexn_list[current->vn2];
        V3 vn3 = render_state->vertexn_list[current->vn3];

#if SHADING == 1
        V3 white_point = 10*v3(1,1,1);
        // ambient light
        V3 v1_light = {};
        V3 v2_light = {};
        V3 v3_light = {};

        if (!vertex_is_light_source(current->v1, 
                                    render_state->light_sources, 
                                    render_state->light_sources_count)) {
            for (uint32_t light = 0; light < render_state->light_sources_count; ++light) {
                LightSource *source = &render_state->light_sources[light];
                v1_light += compute_light_intensity(*source, vertex1, vn1);
            }
            //v1_light = v1_light / (float)render_state->light_sources_count;
            v1_light = tone_map_light(v1_light, white_point);
            *v1_color = v3_pariwise_mul(*v1_color, v1_light);

        }

        else {
            *v1_color = tone_map_light(*v1_color, white_point);
        }

        if (!vertex_is_light_source(current->v2, 
                                    render_state->light_sources, 
                                    render_state->light_sources_count)) {

            for (uint32_t light = 0; light < render_state->light_sources_count; ++light) {
                LightSource *source = &render_state->light_sources[light];
                v2_light += compute_light_intensity(*source, vertex2, vn2);
            }
            //v2_light = v2_light / (float)render_state->light_sources_count;
            v2_light = tone_map_light(v2_light, white_point);
            *v2_color = v3_pariwise_mul(*v2_color, v2_light);
        }

        else {
            *v2_color = tone_map_light(*v2_color, white_point);
        }


        if (!vertex_is_light_source(current->v3, 
                                    render_state->light_sources, 
                                    render_state->light_sources_count)) {

            for (uint32_t light = 0; light < render_state->light_sources_count; ++light) {
                LightSource *source = &render_state->light_sources[light];
                v3_light += compute_light_intensity(*source, vertex3, vn3);
            }
            //v3_light = v3_light / (float)render_state->light_sources_count;
            v3_light = tone_map_light(v3_light, white_point);
            *v3_color = v3_pariwise_mul(v3_light,*v3_color);
        }

        else {
            *v3_color = tone_map_light(*v3_color, white_point);
        }

#endif // SHADING
    }

    for (uint32_t vertex = 0; vertex < render_state->vertex_count; ++vertex) {
        render_state->vertex_list[vertex] = 
            renderer_world_vertex_to_view(render_state->vertex_list[vertex], camera);
    }

    for (uint32_t triangle = 0; triangle < render_state->polygon_count; ++triangle) {
        Triangle current = render_state->polygons[triangle];

        if (v3_dot(render_state->vertex_list[current.v1], 
                v3_cross(render_state->vertex_list[current.v2]-render_state->vertex_list[current.v1], 
                    render_state->vertex_list[current.v3]-render_state->vertex_list[current.v1])) >= 0) {
            continue;
        }

        render_state->polygons_to_draw[render_state->draw_count++] = current;
    }


    float aspect_ratio = ((float)buffer_width/(float)buffer_height);
    for (uint32_t vertex = 0; vertex < render_state->vertex_count; ++vertex) {
        V3 relative_pos = render_state->vertex_list[vertex];

        if (renderer_v3_should_clip(relative_pos, camera, aspect_ratio)) {
            render_state->vertex_list[vertex].x = FLT_MAX;
            render_state->vertex_list[vertex].y = FLT_MAX;
            render_state->vertex_list[vertex].z = FLT_MAX;
            continue;
        }
        
        float normal_x = (relative_pos.x/relative_pos.z) / 
                         (tanf(camera->fov/2.0f)*aspect_ratio);

        float normal_y = (relative_pos.y/relative_pos.z) / 

                         (tanf(camera->fov/2.0f));

        if (normal_x <= -1) {
            render_state->vertex_list[vertex].x = FLT_MAX;
            render_state->vertex_list[vertex].y = FLT_MAX;
            render_state->vertex_list[vertex].z = FLT_MAX;
            continue;
        }

        if (normal_y <= -1) {
            render_state->vertex_list[vertex].x = FLT_MAX;
            render_state->vertex_list[vertex].y = FLT_MAX;
            render_state->vertex_list[vertex].z = FLT_MAX;
            continue;
        }

        if (normal_x >= 1) {
            render_state->vertex_list[vertex].x = FLT_MAX;
            render_state->vertex_list[vertex].y = FLT_MAX;
            render_state->vertex_list[vertex].z = FLT_MAX;
            continue;
        }

        if (normal_y >= 1) {
            render_state->vertex_list[vertex].x = FLT_MAX;
            render_state->vertex_list[vertex].y = FLT_MAX;
            render_state->vertex_list[vertex].z = FLT_MAX;
            continue;
        }

        float x = floorf(((normal_x * (buffer_width/2.0f) + buffer_width/2.0f)));
        float y = floorf(((-normal_y * (buffer_height/2.0f) + buffer_height/2.0f)));

        render_state->vertex_list[vertex].x = x;
        render_state->vertex_list[vertex].y = y;
        render_state->vertex_list[vertex].z = relative_pos.z;
    }
}

static inline bool renderer_v3_clipped(V3 screen) {
    if (screen.x == FLT_MAX || screen.y == FLT_MAX) return true;
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
                               V3 start,
                               V3 end,
                               uint32_t color) {
            
    if (start.x > end.x) {
        V3 temp = end;
        end = start;
        start = temp;
    }

    float dx = end.x - start.x;
    float dy = end.y - start.y;
    float dz = end.z - start.z;
    uint32_t *endofbuffer = (uint32_t *)((uint8_t *)buffer->memory + (buffer->width*buffer->height*buffer->bytes_per_pixel));

    // case where the gradient is below 1;
    float step = 0;
    if (fabsf(dx) > fabsf(dy)) {
        step = fabsf(dx);
    }

    else {
        step = fabsf(dy);
    }

    float x = (float)start.x;
    float y = (float)start.y;

    float x_inc = (float)dx/(float)step;
    float y_inc = (float)dy/(float)step;

    for (int draw = 0; draw < step; ++draw) {
        int offset = (int)(roundf(x))*buffer->bytes_per_pixel + 
                     ((int)roundf(y))*buffer->pitch;

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
                                     V3 start,
                                     V3 end,
                                     V3 start_color,
                                     V3 end_color) {

    if (start.x < 0) {
        start.x = 0;
    }

    else if (start.x >= 1280) {
        start.x = 1279;
    }

    if (start.y < 0) {
        start.y = 0;
    }

    else if (start.y >= 720) {
        start.y = 719;
    }

    if (end.x < 0) {
        end.x = 0;
    }

    else if (end.x >= 1280) {
        end.x = 1279;
    }

    if (end.y < 0) {
        end.y = 0;
    }

    else if (end.y >= 720) {
        end.y = 719;
    }


    if (start.x > end.x) {
        V3 temp = end;
        end = start;
        start = temp;
    }

    float dx = (float)(end.x - start.x);
    float dy = (float)(end.y - start.y);
    float dz = (float)(end.z - start.z);
    uint32_t *endofbuffer = (uint32_t *)((uint8_t *)buffer->memory + (buffer->width*buffer->height*buffer->bytes_per_pixel));

    // case where the gradient is below 1;
    float step = 0;
    if (fabsf(dx) > fabsf(dy)) {
        step = (float)fabsf(dx);
    }

    else {
        step = (float)fabsf(dy);
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

    float r_inc = dr/step;
    float g_inc = dg/step;
    float b_inc = db/step;

    float x_inc = (float)dx/step;
    float y_inc = (float)dy/step;
    float z_inc = (float)dz/step;

    for (int draw = 0; draw < step; ++draw) {
        int offset = (int)(x)*buffer->bytes_per_pixel + 
                     (int)(y)*buffer->pitch;
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
                                            V3 left,
                                            V3 right,
                                            V3 bottom,
                                            V3 left_color,
                                            V3 right_color,
                                            V3 bottom_color) {
        
    float steps = (float)(bottom.y - right.y);
    float x_inc_left = (float)(bottom.x-left.x)/steps;
    float x_inc_right = (float)(bottom.x-right.x)/steps;

    float x_left = (float)left.x;
    float x_right = (float)right.x;
    float z_left = left.z;
    float z_right = right.z;

    float rl = left_color.x;
    float gl = left_color.y;
    float bl = left_color.z; 

    float rr = right_color.x;
    float gr = right_color.y;
    float br = right_color.z; 

    float rbot = bottom_color.x;
    float gbot = bottom_color.y;
    float bbot = bottom_color.z;

    float drl = rbot - rl;
    float dgl = gbot - gl;
    float dbl = bbot - bl;
    float drr = rbot - rr;
    float dgr = gbot - gr;
    float dbr = bbot - br;

    float rl_inc = drl/steps;
    float gl_inc = dgl/steps;
    float bl_inc = dbl/steps;

    float rr_inc = drr/steps;
    float gr_inc = dgr/steps;
    float br_inc = dbr/steps;

    float dzl = left.z - bottom.z;
    float dzr = left.z - right.z;
    float zl_inc = dzl / steps;
    float zr_inc = dzr / steps;

    for (float y = right.y; y <= bottom.y; ++y) {
        V3 left_v = v3(x_left, y, z_left);
        V3 right_v = v3(x_right, y, z_right);
        V3 color_left = v3(rl, gl, bl);
        V3 color_right = v3(rr, gr, br);

        renderer_draw_line_zcull(buffer, zbuffer, left_v, right_v, color_left, color_right);

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
                                               V3 left,
                                               V3 right,
                                               V3 top,
                                               V3 left_color,
                                               V3 right_color,
                                               V3 top_color) {

    float steps = (float)(top.y - right.y);
    float x_inc_left = (float)(top.x-left.x)/(steps);
    float x_inc_right = (float)(top.x-right.x)/(steps);

    float x_left = (float)left.x;
    float x_right = (float)right.x;
    float z_left = left.z;
    float z_right = right.z;

    float rl = left_color.x;
    float gl = left_color.y;
    float bl = left_color.z;
    float rr = right_color.x;
    float gr = right_color.y;
    float br = right_color.z; 

    float rtop = top_color.x;
    float gtop = top_color.y;
    float btop = top_color.z; 

    float drl = rtop - rl;
    float dgl = gtop - gl;
    float dbl = btop - bl;
    float drr = rtop - rr;
    float dgr = gtop - gr;
    float dbr = btop - br;

    float dzl = top.z - left.z;
    float dzr = top.z - right.z;
    float zl_inc = dzl / steps;
    float zr_inc = dzr / steps;

    float rl_inc = drl/steps;
    float gl_inc = dgl/steps;
    float bl_inc = dbl/steps;

    float rr_inc = drr/steps;
    float gr_inc = dgr/steps;
    float br_inc = dbr/steps;
    
    for (float y = right.y; y >= top.y; --y) {
        V3 left_v;
        left_v.x = x_left;
        left_v.y = y;
        left_v.z = z_left;

        V3 right_v;
        right_v.x = x_right;
        right_v.y = y;
        right_v.z = z_right;

        V3 color_left = v3(rl, gl, bl);
        V3 color_right = v3(rr, gr, br);

        renderer_draw_line_zcull(buffer, zbuffer, left_v, right_v, color_left, color_right);

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
                                           V3 *vertices,
                                           V3 *colors,
                                           Triangle *triangles,
                                           int count) {

    for (int i = 0; i < count; ++i) {
        V3 vert1 = vertices[triangles[i].v1];
        V3 vert2 = vertices[triangles[i].v2];
        V3 vert3 = vertices[triangles[i].v3];

        V3 v1_color = colors[triangles[i].v1];
        V3 v2_color = colors[triangles[i].v2];
        V3 v3_color = colors[triangles[i].v3];

        if (renderer_v3_clipped(vert1)) continue;
        if (renderer_v3_clipped(vert2)) continue;
        if (renderer_v3_clipped(vert3)) continue;

        //sort vertexes v3 is the biggest
        if (vert3.y > vert2.y) {
            V3 temp = vert2;
            vert2 = vert3;
            vert3 = temp;
            V3 color_temp = v2_color;
            v2_color = v3_color;
            v3_color = color_temp;

        }

        if (vert2.y > vert1.y) {
            V3 temp = vert1;
            vert1 = vert2;
            vert2 = temp;
            V3 color_temp = v1_color;
            v1_color = v2_color;
            v2_color = color_temp;
        }

        if (vert3.y > vert2.y) {
            V3 temp = vert2;
            vert2 = vert3;
            vert3 = temp;
            V3 color_temp = v2_color;
            v2_color = v3_color;
            v3_color = color_temp;
        }

        if (vert3.y == vert2.y) {
            if (vert3.x > vert2.x) {
                renderer_draw_flat_top_triangle(buffer, zbuffer, vert2, vert3, vert1, v2_color, v3_color, v1_color);
            }
            else {
                renderer_draw_flat_top_triangle(buffer, zbuffer, vert3, vert2, vert1, v3_color, v2_color, v1_color);
            }

            continue;
        }

        if (vert1.y == vert2.y) {
            if (vert1.x > vert2.x) {
                renderer_draw_flat_bottom_triangle(buffer, zbuffer, vert2, vert1, vert3, v2_color, v1_color, v3_color);
            }
            else {
                renderer_draw_flat_bottom_triangle(buffer, zbuffer, vert1, vert2, vert3, v1_color, v2_color, v3_color);
            }
            continue;
        }

        V3 vert4;
        vert4.x = (vert3.x - (((vert3.y-vert2.y)/(vert3.y-vert1.y))*(vert3.x-vert1.x)));
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

        float deltax = vert3.y-vert1.y;
        float dx = vert2.y-vert1.y;

        float v4r = v1r + (dr/deltax)*(dx);
        float v4g = v1g + (dg/deltax)*(dx);
        float v4b = v1b + (db/deltax)*(dx);


        V3 v4_color = v3(v4r, v4g, v4b);
        if (vert4.x > vert2.x) {
            renderer_draw_flat_top_triangle(buffer, zbuffer, vert2, vert4, vert1, v2_color, v4_color, v1_color);
            renderer_draw_flat_bottom_triangle(buffer, zbuffer, vert2, vert4, vert3, v2_color, v4_color, v3_color);
        }
        else {
            renderer_draw_flat_top_triangle(buffer, zbuffer, vert4, vert2, vert1, v4_color, v2_color, v1_color);
            renderer_draw_flat_bottom_triangle(buffer, zbuffer, vert4, vert2, vert3, v4_color, v2_color, v3_color);
        }
    }
    // TODO: post fx
    //
}
