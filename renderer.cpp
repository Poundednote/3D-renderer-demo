#include <emmintrin.h>
#include <stdint.h>
#include "renderer.h"

#if DEBUG_MODE
#define assert(expression) if(!(expression)) {*(int *)0 = 0;}
#else
#define assert(exppression)
#endif

#define arraysize(array) (sizeof(array) / sizeof((array)[0]))

static V3 renderer_world_vertex_to_view(V3 world_pos, GameCamera *camera) {
    V3 result = v3_rotate_q4(world_pos-camera->pos, 
                          (rotation_q4(-camera->theta_x, v3(1,0,0))*
                          rotation_q4(-camera->theta_y, v3(0,1,0))));

    return result;
}
static void vertex4_buffer_to_view(Vertex4 *buffer, 
                                   GameCamera *camera, 
                                   Quaternion rotation, 
                                   uint32_t count) {

        __m128 _qw = _mm_set1_ps(rotation.scalar); 
        __m128 _qx = _mm_set1_ps(rotation.vector.x); 
        __m128 _qy = _mm_set1_ps(rotation.vector.y); 
        __m128 _qz = _mm_set1_ps(rotation.vector.z); 

        __m128 _cam_x = _mm_set1_ps(camera->pos.x);
        __m128 _cam_y = _mm_set1_ps(camera->pos.y);
        __m128 _cam_z = _mm_set1_ps(camera->pos.z);


    for (uint32_t i = 0; i < count; ++i) {
        Vertex4 _vertex = buffer[i];
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

        buffer[i].x = _rotated_x;
        buffer[i].y = _rotated_y;
        buffer[i].z = _rotated_z;
    }
}

static bool v3_should_clip(V3 pos, GameCamera *camera, float aspect_ratio) {
    float max_x = tanf(camera->fov/2.0f)*camera->znear * (pos.z/camera->znear) * 
        aspect_ratio;
    float max_y = tanf(camera->fov/2.0f)*camera->znear * (pos.z/camera->znear);

    if (pos.z <= camera->znear) return true;
    if (pos.z > camera->zfar) return true;
    if (fabs(pos.x) > max_x) return true;
    if (fabs(pos.y) > max_y) return true;

    return false;
}

static V3 compute_light_intensity(LightSource source, V3 vertex, V3 normal) {
    V3 dl = source.position-vertex; 
    float intensity = v3_dot(v3_norm(dl), normal);
    if (intensity > 0) {
        return intensity*source.color;
    }
    return source.color = {};
}

static inline float rgb_to_luminance(V3 rgb) {
    return v3_dot(rgb, v3(0.2126f, 0.7152f, 0.0722f));
}

static inline V3 triangle_get_attribute_from_buffer(int index, Vertex4 *buffer) {
        V3 result = v3(((float *)&buffer[index/4].x)[index%4],
                      ((float *)&buffer[index/4].y)[index%4],
                      ((float *)&buffer[index/4].z)[index%4]);

        return result;
}

static inline void triangle_set_attribute_to_buffer(V3 value, int index, Vertex4 *buffer) {
    ((float * )&buffer[index/4].x)[index%4] = value.x;
    ((float * )&buffer[index/4].y)[index%4] = value.y;
    ((float * )&buffer[index/4].z)[index%4] = value.z;
}

static void renderer_transform_light_and_cull(RendererState *render_state,
                                              GameCamera *camera, 
                                              int buffer_width, 
                                              int buffer_height) {
    render_state->draw_count = 0;
    // transform vertices to view space
    { 
        Quaternion rotation = rotation_q4(-camera->theta_x, v3(1,0,0))*
            rotation_q4(-camera->theta_y, v3(0,1,0));

        vertex4_buffer_to_view(render_state->vertex_buffer,
                               camera,
                               rotation,
                               render_state->vertex_count);

        vertex4_buffer_to_view(render_state->normal_buffer,
                               camera,
                               rotation,
                               render_state->normal_count);
    }

    // Backface culling
    for (uint32_t triangle = 0; triangle < render_state->polygon_count; ++triangle) {
        Triangle current = render_state->polygons[triangle];
        V3 vert1 = triangle_get_attribute_from_buffer(current.v1, render_state->vertex_buffer);
        V3 vert2 = triangle_get_attribute_from_buffer(current.v2, render_state->vertex_buffer);
        V3 vert3 = triangle_get_attribute_from_buffer(current.v3, render_state->vertex_buffer);


#if 1
        if (v3_dot(vert1, v3_cross(vert2-vert1, vert3-vert1)) >= 0) {
            continue;
        }
#endif

        render_state->polygons_to_draw[render_state->draw_count++] = current;
    }

#if SHADING
    // gourad shading
    for (uint32_t triangle = 0; triangle < render_state->draw_count; ++triangle) {
        Triangle current = render_state->polygons_to_draw[triangle];

        V3 vertex1 = triangle_get_attribute_from_buffer(current.v1, 
                                                        render_state->vertex_buffer);
        V3 vertex2 = triangle_get_attribute_from_buffer(current.v2, 
                                                        render_state->vertex_buffer);
        V3 vertex3 = triangle_get_attribute_from_buffer(current.v3, 
                                                        render_state->vertex_buffer);

        V3 v1_color = triangle_get_attribute_from_buffer(current.v1, 
                                                         render_state->vertex_colors);
        V3 v2_color = triangle_get_attribute_from_buffer(current.v2, 
                                                         render_state->vertex_colors);
        V3 v3_color = triangle_get_attribute_from_buffer(current.v3, 
                                                         render_state->vertex_colors);

        V3 vn1 = triangle_get_attribute_from_buffer(current.vn1, 
                                                    render_state->vertex_buffer);
        V3 vn2 = triangle_get_attribute_from_buffer(current.vn2, 
                                                    render_state->vertex_buffer);
        V3 vn3 = triangle_get_attribute_from_buffer(current.vn3, 
                                                    render_state->vertex_buffer);

        // ambient light
        V3 v1_light = {};
        V3 v2_light = {};
        V3 v3_light = {};

        for (uint32_t light = 0; light < render_state->light_sources_count; ++light) {
            LightSource source = render_state->light_sources[light];
            source.position = renderer_world_vertex_to_view(source.position, 
                                                            camera);
            v1_light += compute_light_intensity(source, vertex1, vn1);
            v2_light += compute_light_intensity(source, vertex2, vn2);
            v3_light += compute_light_intensity(source, vertex3, vn3);
        }

        v1_color = v3_pairwise_mul(v1_color, v1_light);
        v2_color = v3_pairwise_mul(v2_color, v2_light);
        v3_color = v3_pairwise_mul(v3_color, v3_light);

        triangle_set_attribute_to_buffer(v1_color, 
                                         current.v1, 
                                         render_state->vertex_colors);
        triangle_set_attribute_to_buffer(v2_color, 
                                         current.v2, 
                                         render_state->vertex_colors);
        triangle_set_attribute_to_buffer(v3_color, 
                                         current.v3, 
                                         render_state->vertex_colors);
    }
#endif // SHADING

    // change light sources to just be the original color
    for (uint32_t i = 0; i < render_state->light_sources_count; ++i) {
        LightSource light_source = render_state->light_sources[i]; 
        for (uint32_t j = light_source.obj->vstart; j < light_source.obj->vend; j++) {
            Vertex4 _vertex = render_state->vertex_buffer[j];
            Vertex4 _color = render_state->vertex_colors[j];

            _color.x = _mm_add_ps(_mm_set1_ps(light_source.color.x), _mm_set1_ps(100));
            _color.y = _mm_add_ps(_mm_set1_ps(light_source.color.y), _mm_set1_ps(100));
            _color.z = _mm_add_ps(_mm_set1_ps(light_source.color.z), _mm_set1_ps(100));

            render_state->vertex_colors[j] = _color;
        }
    }


    // perspective projection & clipping
    {
        float aspect_ratio = (float)buffer_width / (float)buffer_height;
        for (uint32_t vertex = 0; vertex < render_state->vertex_count; ++vertex) {

            __m128 _vertex_x = render_state->vertex_buffer[vertex].x;
            __m128 _vertex_y = render_state->vertex_buffer[vertex].y;
            __m128 _vertex_z = render_state->vertex_buffer[vertex].z;

            __m128 _cam_znear = _mm_set1_ps(camera->znear);
            __m128 _cam_zfar = _mm_set1_ps(camera->zfar);

            __m128 _length_x = _mm_set1_ps(tanf(camera->fov/2.0f)*camera->znear*aspect_ratio);
            __m128 _length_y = _mm_set1_ps(tanf(camera->fov/2.0f)*camera->znear);

            __m128 _mask = _mm_cmplt_ps(_vertex_z, _cam_znear);
            _mask = _mm_or_ps(_mask, _mm_cmpgt_ps(_vertex_z, _cam_zfar));

            __m128 _normal_x = _mm_div_ps(_mm_mul_ps(_vertex_x, _mm_div_ps(_cam_znear, _vertex_z)),
                    _length_x);

            __m128 _normal_y = _mm_div_ps(_mm_mul_ps(_vertex_y, _mm_div_ps(_cam_znear, _vertex_z)),
                    _length_y);

            _mask = _mm_or_ps(_mask, _mm_cmpgt_ps(_normal_x, _mm_set1_ps(1.0f)));
            _mask = _mm_or_ps(_mask, _mm_cmplt_ps(_normal_x, _mm_set1_ps(-1.0f)));

            _mask = _mm_or_ps(_mask, _mm_cmpgt_ps(_normal_y, _mm_set1_ps(1.0f)));
            _mask = _mm_or_ps(_mask, _mm_cmplt_ps(_normal_y, _mm_set1_ps(-1.0f)));

            __m128 _buffer_x = _mm_set1_ps(buffer_width/2.0f);
            __m128 _buffer_y = _mm_set1_ps(buffer_height/2.0f);

            __m128 _float_res_x = _mm_add_ps(_mm_mul_ps(_normal_x, _buffer_x), _buffer_x);
            __m128 _float_res_y = _mm_add_ps(_mm_mul_ps(_mm_mul_ps(_normal_y, _mm_set1_ps(-1.0f)), _buffer_y), _buffer_y);
            _float_res_x = _mm_cvtepi32_ps(_mm_cvttps_epi32(_float_res_x));
            _float_res_y = _mm_cvtepi32_ps(_mm_cvttps_epi32(_float_res_y));

            render_state->vertex_buffer[vertex].x = _mm_or_ps(_float_res_x, _mask);
            render_state->vertex_buffer[vertex].y = _mm_or_ps(_float_res_y, _mask);
            render_state->vertex_buffer[vertex].z = _mm_or_ps(_vertex_z, _mask);
        }
    }
}


static inline bool renderer_v3_isclipped(V3 screen) {
    if (*((uint32_t *)&screen.x) == 0xFFFFFFFF || *((uint32_t *)&screen.y) == 0xFFFFFFFF) return true;
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

static void clamp_to_screen(V3 *vertex, int buffer_width, int buffer_height) {
    if (vertex->x < 0) {
        vertex->x = 0;
    }

    else if (vertex->x >= buffer_width) {
        vertex->x = (float)(buffer_width - 1);
    }

    if (vertex->y < 0) {
        vertex->y = 0;
    }

    else if (vertex->y >= buffer_height) {
        vertex->y = (float)(buffer_height - 1);
    }
}
static void renderer_draw_line(OffscreenBuffer *buffer,
                               OffscreenBuffer *zbuffer,
                               V3 start,
                               V3 end,
                               V3 start_color,
                               V3 end_color) {


    clamp_to_screen(&start, buffer->width, buffer->height);
    clamp_to_screen(&end, buffer->width, buffer->height);
    if (start.x > end.x) {
        v3_swap(&start, &end);
    }

    float dx = (float)(end.x - start.x);
    float dy = (float)(end.y - start.y);
    //
    // case where the gradient is below 1;

    float step = 0;
    if (fabsf(dx) > fabsf(dy)) {
        step = (float)fabsf(dx);
    }

    else {
        step = (float)fabsf(dy);
    }

    /* I am using the padding byte in the windows bitmap buffer to seperate 
     * the bright objects for post processing all bright pixels go through the draw line call
     * with a color of > 1 this has to be then corrected and then return either 0xFF or 0 for 
     * bright and not bright respectively
     */
    uint8_t bright = 0x0;
    if (end_color.x > 1 && end_color.y > 1 && end_color.z > 1) {
        bright = 0xFF;
        start_color -= v3(1,1,1);
        end_color -= v3(1,1,1);
    }

    V3 color_inc = (end_color - start_color) / step;
    V3 start_inc = (end-start)/step;

    for (int draw = 0; draw <= step; ++draw) {
        int offset = (int)(roundf(start.x-0.5f))*buffer->bytes_per_pixel + 
                     (int)(roundf(start.y-0.5f))*buffer->pitch;
        float *depth_value = (float *)((uint8_t *)zbuffer->memory + offset);
    
           if ((start.z < *depth_value)) {
               uint32_t color = (bright << 24) |
                                (((uint8_t)(0xFF*start_color.x)) << 16) | 
                                (((uint8_t)(0xFF*start_color.y))) << 8 | 
                                (uint8_t)(0xFF*start_color.z);

                uint32_t *pixel = (uint32_t *)((uint8_t *)buffer->memory + offset);
                *pixel = color;
                *depth_value = start.z;
        }

        start += start_inc;
        start_color += color_inc;
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

    if (!steps) {
        return;
    }

    V3 right_inc = v3((bottom.x-right.x)/steps, 1, -(bottom.z-right.z)/steps);
    V3 left_inc = v3((bottom.x-left.x)/steps, 1, -(bottom.z-left.z)/steps);
    V3 left_color_inc = (bottom_color - left_color)/steps;
    V3 right_color_inc = (bottom_color - right_color)/steps;

    for (float y = right.y; y <= bottom.y; ++y) {
        renderer_draw_line(buffer, zbuffer, left, right, left_color, right_color);

        left += left_inc;
        right += right_inc;
        left_color += left_color_inc;
        right_color += right_color_inc;
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

    if (!steps) {
        return;
    }

    V3 left_inc = v3((top.x-left.x)/steps, 1, (top.z-left.z)/steps);
    V3 right_inc = v3((top.x-right.x)/steps, 1, (top.z-right.z)/steps);
    V3 right_color_inc = (top_color - right_color)/steps;
    V3 left_color_inc = (top_color - left_color)/steps;

    for (float y = right.y; y >= top.y; --y) {
        renderer_draw_line(buffer, zbuffer, left, right, left_color, right_color);

        left -= left_inc;
        right -= right_inc;
        left_color -= left_color_inc;
        right_color -= right_color_inc;
    }
}

static void tone_map_color(V3 *color, V3 white_point) {
    float luminance = rgb_to_luminance(*color);
    float white_l = rgb_to_luminance(white_point);
    float numerator = luminance * (1.0f + (luminance/(white_l * white_l)));
    float new_luminance = numerator / (1.0f + luminance);
    *color = (new_luminance/luminance)*(*color);
    (color->x > 1) ? color->x = 1 : color->x;
    (color->y > 1) ? color->y = 1 : color->y;
    (color->z > 1) ? color->z = 1 : color->z;
}

static void color_correct_brightness(V3 *color, V3 white_point) {
        if (color->x > 100 && 
            color->y > 100 && 
            color->z > 100) {

            *color -= v3(100,100,100);
            *color = 100*(*color);
            tone_map_color(color, white_point);
            *color += v3(1,1,1);
        }

        else {
            tone_map_color(color, white_point);
        }
}

static void renderer_draw_triangles_filled(OffscreenBuffer *buffer,
                                           OffscreenBuffer *zbuffer,
                                           OffscreenBuffer *postbuffer,
                                           Vertex4 *vertices,
                                           Vertex4 *colors,
                                           Triangle *triangles,
                                           int count) {

    // Draw Triangles
    {
        V3 white_point = v3(100,100,100);
        for (int i = 0; i < count; ++i) {


            V3 vert1 = triangle_get_attribute_from_buffer(triangles[i].v1, vertices);
            V3 vert2 = triangle_get_attribute_from_buffer(triangles[i].v2, vertices);
            V3 vert3 = triangle_get_attribute_from_buffer(triangles[i].v3, vertices);

            V3 v1_color = triangle_get_attribute_from_buffer(triangles[i].v1, colors);
            V3 v2_color = triangle_get_attribute_from_buffer(triangles[i].v2, colors);
            V3 v3_color = triangle_get_attribute_from_buffer(triangles[i].v3, colors);

            if (renderer_v3_isclipped(vert1)) {continue;}
            if (renderer_v3_isclipped(vert2)) {continue;}
            if (renderer_v3_isclipped(vert3)) {continue;}

            color_correct_brightness(&v1_color, white_point);
            color_correct_brightness(&v2_color, white_point);
            color_correct_brightness(&v3_color, white_point);

            //sort vertexes v3 is the biggest
            if (vert3.y > vert2.y) {
                v3_swap(&vert3, &vert2);
                v3_swap(&v3_color, &v2_color);

            }

            if (vert2.y > vert1.y) {
                v3_swap(&vert2, &vert1);
                v3_swap(&v2_color, &v1_color);
            }

            if (vert3.y > vert2.y) {
                v3_swap(&vert3, &vert2);
                v3_swap(&v3_color, &v2_color);
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

            float lerp_factor = vert2.y-vert1.y;
            float delta_x = vert3.y-vert1.y;

            V3 vert4;
            vert4.x = f_lerp(delta_x, lerp_factor, vert1.x, vert3.x);
            vert4.y = vert2.y;
            vert4.z = f_lerp(delta_x, lerp_factor, vert1.z, vert3.z);

            V3 v4_color = v3_lerp(delta_x, lerp_factor, v1_color, v3_color);

            if (vert4.x > vert2.x) {
                renderer_draw_flat_top_triangle(buffer, zbuffer, vert2, vert4, vert1, v2_color, v4_color, v1_color);
                renderer_draw_flat_bottom_triangle(buffer, zbuffer, vert2, vert4, vert3, v2_color, v4_color, v3_color);
            }
            else {
                renderer_draw_flat_top_triangle(buffer, zbuffer, vert4, vert2, vert1, v4_color, v2_color, v1_color);
                renderer_draw_flat_bottom_triangle(buffer, zbuffer, vert4, vert2, vert3, v4_color, v2_color, v3_color);
            }
        }
    }

#if 1
    // POST PROCESSING BLOOM EFFECT
    
    // zero the buffer
    {
        uint32_t *pixels = (uint32_t *)postbuffer->memory;
        for (int i = 0; i < postbuffer->width*postbuffer->height; ++i) {
            pixels[i] = 0; 
        }
    }

    // box sample buffer into the post processing buffer 
    int ratio = buffer->width / postbuffer->width;
    {
        uint32_t *in_pixels = (uint32_t *)buffer->memory;
        uint32_t *out_pixels = (uint32_t *)postbuffer->memory;
        int radius = 2;
        for (int y = 0; y < postbuffer->height; y++) {
            if (y < radius || y + radius == buffer->height) {
                continue;
            }

            for (int x = 0; x < postbuffer->width; x++) {
                if (x < radius || x + radius == buffer->width) {
                    continue;
                }

                int sumbright = 0;
                for (int i = -radius; i <radius+1; ++i) {
                    for (int j = -radius; j < radius+1; ++j) {
                        sumbright += ((in_pixels[(y*ratio+i)*buffer->width+x*ratio+j] & 0xFF000000) >> 24);
                    }
                }

                if (sumbright < 255) {
                    out_pixels[(y)*postbuffer->width+x] = 0;
                    continue;
                }

                int sumr = 0;
                int sumg = 0;
                int sumb = 0;
                for (int i = -radius; i < radius+1; ++i) {
                    for (int j = -radius; j < radius+1; ++j) {
                        sumr += ((in_pixels[(y*ratio+i)*buffer->width+x*ratio+j] & 0x00FF0000) >> 16);
                        sumg += ((in_pixels[(y*ratio+i)*buffer->width+x*ratio+j] & 0x0000FF00) >> 8);
                        sumb += (in_pixels[(y*ratio+i)*buffer->width+x*ratio+j] & 0x000000FF);
                    }
                }

                uint8_t red = (uint8_t)(sumr/(radius*2+1)*(radius*2+1)); 
                uint8_t green = (uint8_t)(sumg/(radius*2+1)*(radius*2+1)); 
                uint8_t blue = (uint8_t)(sumb/(radius*2+1)*(radius*2+1)); 
                out_pixels[(y)*postbuffer->width+x] = red << 16 | green << 8 | blue;
            }
        }
    }

    // Apply box blur
    {
        uint32_t *pixels = (uint32_t *)postbuffer->memory;
        for (int i = 0; i < ratio; ++i) {
            for (int y = 0; y < postbuffer->height; ++y) {
                if (y < 1 || y + 1 == postbuffer->height) {
                    continue;
                }

                int sumr = (((pixels[(y+1)*postbuffer->width]) & 0x00FF0000) >> 16) +
                    (((pixels[(y)*postbuffer->width]) & 0x00FF0000) >> 16) +
                    (((pixels[(y-1)*postbuffer->width]) & 0x00FF0000) >> 16); 

                int sumg = (((pixels[(y+1)*postbuffer->width]) & 0x0000FF00) >> 8) +
                    (((pixels[(y)*postbuffer->width]) & 0x0000FF00) >> 8) +
                    (((pixels[(y-1)*postbuffer->width]) & 0x0000FF00) >> 8); 

                int sumb = ((pixels[(y+1)*postbuffer->width]) & 0x000000FF) +
                    ((pixels[(y)*postbuffer->width]) & 0x000000FF) +
                    ((pixels[(y-1)*postbuffer->width]) & 0x000000FF); 

                int counter = 3;
                int leftr = 0;
                int leftg = 0;
                int leftb = 0;

                for (int x = 0; x < postbuffer->width; ++x) {
                    if ((x + 1) < postbuffer->width) {

                        sumr += (((pixels[(y-1)*postbuffer->width+x+1] & 0x00FF0000) >> 16) + 
                                ((pixels[(y)*postbuffer->width+x+1] & 0x00FF0000) >> 16) +
                                ((pixels[(y+1)*postbuffer->width+x+1] & 0x00FF0000) >> 16));

                        sumg += (((pixels[(y-1)*postbuffer->width+x+1] & 0x0000FF00) >> 8) + 
                                ((pixels[(y)*postbuffer->width+x+1] & 0x0000FF00) >> 8) +
                                ((pixels[(y+1)*postbuffer->width+x+1] & 0x0000FF00) >> 8));

                        sumb += ((pixels[(y-1)*postbuffer->width+x+1] & 0x000000FF) + 
                                (pixels[(y)*postbuffer->width+x+1] & 0x000000FF) +
                                (pixels[(y+1)*postbuffer->width+x+1] & 0x000000FF));

                        counter += 3;
                    }

                    uint8_t red = (uint8_t)(sumr/counter); 
                    uint8_t green = (uint8_t)(sumg/counter); 
                    uint8_t blue = (uint8_t)(sumb/counter); 

                    sumr -= leftr;
                    sumg -= leftg;
                    sumb -= leftb;

                    leftr = (((pixels[(y-1)*postbuffer->width+x] & 0x00FF0000) >> 16) + 
                            ((pixels[(y)*postbuffer->width+x] & 0x00FF0000) >> 16) +
                            ((pixels[(y+1)*postbuffer->width+x] & 0x00FF0000) >> 16));

                    leftg = (((pixels[(y-1)*postbuffer->width+x] & 0x0000FF00) >> 8) + 
                            ((pixels[(y)*postbuffer->width+x] & 0x0000FF00) >> 8) +
                            ((pixels[(y+1)*postbuffer->width+x] & 0x0000FF00) >> 8));

                    leftb = ((pixels[(y-1)*postbuffer->width+x] & 0x000000FF) + 
                            (pixels[(y)*postbuffer->width+x] & 0x000000FF) +
                            (pixels[(y+1)*postbuffer->width+x] & 0x000000FF));

                    if (x >= 1) {
                        counter -= 3;
                    }

                    pixels[(y)*postbuffer->width+x] = red << 16 | green << 8 | blue;
                }
            }
        }
    }


    // Upscale buffer using bilinear interpolation and blend into the display buffer
    for (int y = 0; y < postbuffer->height; ++y) {
        if (y + 1 == postbuffer->height) {
            continue;
        }

        for (int x = 0; x < postbuffer->width; ++x) {
            if (x + 1 == postbuffer->width) {
                continue;
            }

            uint32_t *color_start = (uint32_t *)postbuffer->memory +
                y*postbuffer->width+x;

            uint32_t *ycolor_end = (uint32_t *)postbuffer->memory +
                (y+1)*postbuffer->width+x;

            uint32_t *xcolor_end = (uint32_t *)postbuffer->memory +
                (y)*postbuffer->width+x+1;

            if (*color_start == 0 && *ycolor_end == 0 && *xcolor_end == 0) {
                continue;
            }

            float yr_inc = ((float)((*ycolor_end & 0x00FF0000) >> 16) -
                    (float)((*color_start & 0x00FF0000) >> 16)) / (float)ratio; 

            float yg_inc = ((float)((*ycolor_end & 0x0000FF00) >> 8) -
                    (float)((*color_start & 0x0000FF00) >> 8)) / (float)ratio; 

            float yb_inc = ((float)(*ycolor_end & 0x000000FF) -
                    (float)(*color_start & 0x000000FF)) / (float)ratio;

            float xr_inc = ((float)((*xcolor_end & 0x00FF0000) >> 16) -
                    (float)((*color_start & 0x00FF0000) >> 16)) / (float)ratio;

            float xg_inc = ((float)((*xcolor_end & 0x0000FF00) >> 8) -
                    (float)((*color_start & 0x0000FF00) >> 8)) / (float)ratio;

            float xb_inc = ((float)(*xcolor_end & 0x000000FF) -
                    (float)(*color_start & 0x000000FF)) / (float)ratio;

            for (int i = 0; i < ratio; ++i) {
                float r = (float)((*color_start & 0x00FF0000) >> 16) + i*yr_inc;
                float g = (float)((*color_start & 0x0000FF00) >> 8) + i*yg_inc;
                float b = (float)(*color_start & 0x000000FF) + i*yb_inc;
                for (int j = 0; j < ratio; ++j) {
                    uint32_t *out_pixel = (uint32_t *)buffer->memory +
                        (y*ratio+i)*buffer->width+(x*ratio)+j;

                    float out_r = (uint8_t)((*out_pixel & 0x00FF0000) >> 16);
                    float out_g = (uint8_t)((*out_pixel & 0x0000FF00) >> 8);
                    float out_b = (uint8_t)(*out_pixel & 0x000000FF);

                    r + out_r > 255 ? out_r = 255 : out_r += r;
                    g + out_g > 255 ? out_g = 255 : out_g += g;
                    b + out_b > 255 ? out_b = 255 : out_b += b;
                    out_r < 0 ? out_r = 0: out_r;
                    out_g < 0 ? out_g = 0: out_g;
                    out_b < 0 ? out_b = 0: out_b;

                    *out_pixel = ((uint8_t)out_r << 16) | ((uint8_t)out_g << 8) | (uint8_t)out_b;

                    r += xr_inc;
                    g += xg_inc;
                    b += xb_inc;
                }
            }
        }
    }
#endif
}


static RenderObj renderer_render_obj_create(RendererState *render_state, Mesh *mesh, V3 color) {
    RenderObj result = {};
    uint32_t vertex_start = render_state->vertex_count;
    assert(mesh->vert_count % 4 == 0);
    for (uint32_t vert = 0; vert < mesh->vert_count; vert+=4) {
    assert(vert+3 < mesh->vert_count);

        V3 vert1 = mesh->vertices[vert];
        V3 vert2 = mesh->vertices[vert+1];
        V3 vert3 = mesh->vertices[vert+2];
        V3 vert4 = mesh->vertices[vert+3];

        ((float * )&render_state->vertex_buffer[render_state->vertex_count].x)[0] = vert1.x; 
        ((float * )&render_state->vertex_buffer[render_state->vertex_count].x)[1] = vert2.x; 
        ((float * )&render_state->vertex_buffer[render_state->vertex_count].x)[2] = vert3.x; 
        ((float * )&render_state->vertex_buffer[render_state->vertex_count].x)[3] = vert4.x; 

        ((float * )&render_state->vertex_buffer[render_state->vertex_count].y)[0] = vert1.y; 
        ((float * )&render_state->vertex_buffer[render_state->vertex_count].y)[1] = vert2.y; 
        ((float * )&render_state->vertex_buffer[render_state->vertex_count].y)[2] = vert3.y; 
        ((float * )&render_state->vertex_buffer[render_state->vertex_count].y)[3] = vert4.y; 

        ((float * )&render_state->vertex_buffer[render_state->vertex_count].z)[0] = vert1.z; 
        ((float * )&render_state->vertex_buffer[render_state->vertex_count].z)[1] = vert2.z; 
        ((float * )&render_state->vertex_buffer[render_state->vertex_count].z)[2] = vert3.z; 
        ((float * )&render_state->vertex_buffer[render_state->vertex_count].z)[3] = vert4.z; 

        __m128 _color_r = _mm_set1_ps(color.x);
        __m128 _color_g = _mm_set1_ps(color.y);
        __m128 _color_b = _mm_set1_ps(color.z);

        render_state->vertex_colors[render_state->vertex_count].x = _color_r;
        render_state->vertex_colors[render_state->vertex_count].y = _color_g;
        render_state->vertex_colors[render_state->vertex_count].z = _color_b;

        ++render_state->vertex_count;
    }


    int vertexnorm_start = render_state->normal_count;
    for (uint32_t normal = 0; normal < mesh->vertexn_count; ++normal) {

        V3 vert1 = mesh->vertexn[normal];
        V3 vert2 = mesh->vertexn[normal+1];
        V3 vert3 = mesh->vertexn[normal+2];
        V3 vert4 = mesh->vertexn[normal+3];
        ((float * )&render_state->normal_buffer[render_state->normal_count].x)[0] = vert1.x; 
        ((float * )&render_state->normal_buffer[render_state->normal_count].x)[1] = vert2.x; 
        ((float * )&render_state->normal_buffer[render_state->normal_count].x)[2] = vert3.x; 
        ((float * )&render_state->normal_buffer[render_state->normal_count].x)[3] = vert4.x; 

        ((float * )&render_state->normal_buffer[render_state->normal_count].y)[0] = vert1.y; 
        ((float * )&render_state->normal_buffer[render_state->normal_count].y)[1] = vert2.y; 
        ((float * )&render_state->normal_buffer[render_state->normal_count].y)[2] = vert3.y; 
        ((float * )&render_state->normal_buffer[render_state->normal_count].y)[3] = vert4.y; 

        ((float * )&render_state->normal_buffer[render_state->normal_count].z)[0] = vert1.z; 
        ((float * )&render_state->normal_buffer[render_state->normal_count].z)[1] = vert2.z; 
        ((float * )&render_state->normal_buffer[render_state->normal_count].z)[2] = vert3.z; 
        ((float * )&render_state->normal_buffer[render_state->normal_count].z)[3] = vert4.z; 

        ++render_state->normal_count;
    }

    int index_start = render_state->polygon_count;
    for (uint32_t polygon = 0; polygon < mesh->poly_count; ++polygon) {
        render_state->polygons[render_state->polygon_count] = mesh->polygons[polygon];
        render_state->polygons[render_state->polygon_count].v1 += vertex_start*4;
        render_state->polygons[render_state->polygon_count].v2 += vertex_start*4;
        render_state->polygons[render_state->polygon_count].v3 += vertex_start*4;
        render_state->polygons[render_state->polygon_count].vn1 += vertexnorm_start*4;
        render_state->polygons[render_state->polygon_count].vn2 += vertexnorm_start*4;
        render_state->polygons[render_state->polygon_count].vn3 += vertexnorm_start*4;
        render_state->polygon_count++;
    }


    result.vstart = vertex_start;
    result.vend = render_state->vertex_count;
    result.index_start = index_start;
    result.index_end = index_start+mesh->poly_count;
    result.mesh = mesh;
    result.color = color;
    assert(render_state->vertex_count < arraysize(render_state->vertex_buffer));

    return result;
}


static void renderer_render_obj_update(RendererState *render_state, 
                                       RenderObj *obj, 
                                       V3 scale, 
                                       Quaternion rotation, 
                                       V3 translation) {

    assert(render_state->vertex_count < arraysize(render_state->vertex_buffer));
    assert(obj->mesh->vert_count % 4 == 0);
    for (uint32_t vertex = 0; vertex < obj->mesh->vert_count; vertex+=4) {
        assert(vertex+3 < obj->mesh->vert_count);
        V3 vert1 = v3_rotate_q4(v3_pairwise_mul(scale,obj->mesh->vertices[vertex]), rotation); 
        V3 vert2 = v3_rotate_q4(v3_pairwise_mul(scale,obj->mesh->vertices[vertex+1]), rotation); 
        V3 vert3 = v3_rotate_q4(v3_pairwise_mul(scale,obj->mesh->vertices[vertex+2]), rotation); 
        V3 vert4 = v3_rotate_q4(v3_pairwise_mul(scale,obj->mesh->vertices[vertex+3]), rotation);

        vert1 += translation;
        vert2 += translation;
        vert3 += translation;
        vert4 += translation;

        ((float * )&render_state->vertex_buffer[obj->vstart+vertex/4].x)[0] = vert1.x; 
        ((float * )&render_state->vertex_buffer[obj->vstart+vertex/4].x)[1] = vert2.x; 
        ((float * )&render_state->vertex_buffer[obj->vstart+vertex/4].x)[2] = vert3.x; 
        ((float * )&render_state->vertex_buffer[obj->vstart+vertex/4].x)[3] = vert4.x; 

        ((float * )&render_state->vertex_buffer[obj->vstart+vertex/4].y)[0] = vert1.y; 
        ((float * )&render_state->vertex_buffer[obj->vstart+vertex/4].y)[1] = vert2.y; 
        ((float * )&render_state->vertex_buffer[obj->vstart+vertex/4].y)[2] = vert3.y; 
        ((float * )&render_state->vertex_buffer[obj->vstart+vertex/4].y)[3] = vert4.y; 

        ((float * )&render_state->vertex_buffer[obj->vstart+vertex/4].z)[0] = vert1.z; 
        ((float * )&render_state->vertex_buffer[obj->vstart+vertex/4].z)[1] = vert2.z; 
        ((float * )&render_state->vertex_buffer[obj->vstart+vertex/4].z)[2] = vert3.z; 
        ((float * )&render_state->vertex_buffer[obj->vstart+vertex/4].z)[3] = vert4.z; 

        __m128 _color_r = _mm_set1_ps(obj->color.x);
        __m128 _color_g = _mm_set1_ps(obj->color.y);
        __m128 _color_b = _mm_set1_ps(obj->color.z);

        render_state->vertex_colors[obj->vstart+vertex/4].x = _color_r;
        render_state->vertex_colors[obj->vstart+vertex/4].y = _color_g;
        render_state->vertex_colors[obj->vstart+vertex/4].z = _color_b;
    }
}

/* TODO: ADD Dynamic Loading and deloading of chunks in 
 * renderer instead of just overwriting the whole buffer evertime */
#if 0
static void destroy_render_obj(RendererState *render_state,
                               RenderObj *obj) {

    for (uint32_t vertex = obj->vstart; vertex < obj->vend; ++vertex) {
        render_state->vertex_buffer[vertex] = v3_zero();
        render_state->vertex_colors[vertex] = v3_zero();
    }

    for (uint32_t index = obj->index_start; index < obj->index_end; ++index) {
        render_state->polygons[index] = {};
    }
}
#endif

static void renderer_draw(RendererState *render_state, 
                          GameCamera *camera, 
                          OffscreenBuffer *buffer, 
                          OffscreenBuffer *zbuffer,
                          OffscreenBuffer *postbuffer) {

    renderer_transform_light_and_cull(render_state, camera, buffer->width, buffer->height);
    renderer_draw_triangles_filled(buffer, zbuffer, postbuffer, 
                                   render_state->vertex_buffer, 
                                   render_state->vertex_colors, 
                                   render_state->polygons_to_draw, 
                                   render_state->draw_count);
}
