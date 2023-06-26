#include <emmintrin.h>
#include <stdint.h>
#include <xmmintrin.h>
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

static inline V3 renderer_view_vertex_to_screen(V3 view_pos, OffscreenBuffer *buffer, GameCamera *camera) {

    float aspect_ratio = (float)buffer->width / (float)buffer->height;
    float length_x = tanf(camera->fov/2.0f)*camera->znear*aspect_ratio;
    float length_y = tanf(camera->fov/2.0f)*camera->znear;
    float normal_x = ((view_pos.x*camera->znear) / view_pos.z*length_x);
    float normal_y = ((view_pos.y*camera->znear) / view_pos.z*length_y);

    float x = floorf((normal_x * (buffer->width/2.0f)) + buffer->width/2.0f);
    float y = floorf((normal_y * (buffer->height/2.0f)) + buffer->height/2.0f);
    return v3(x, y, view_pos.z);
}

static Vertex4 vertex4_rotate(Vertex4 _vertex, Quaternion4 _q) {
        Vertex4 result = {};

        Vertex4 _cross = vertex4_cross(_q.vector, _vertex);
        __m128 _dot = vertex4_dot(_vertex, _q.vector);
        Vertex4 _scaled = vertex4_scale(_vertex, _q.w);

        Quaternion4 _tmp;
        _tmp.w = _mm_mul_ps(_mm_set1_ps(-1.0f), _dot);
        _tmp.vector = vertex4_add(_cross, _scaled);

        Quaternion4 _conj;
        _conj.w = _q.w;
        _conj.vector = vertex4_scale(_q.vector, _mm_set1_ps(-1.0f));

        Vertex4 _tmp2 = vertex4_add(vertex4_scale(_conj.vector, _tmp.w), 
                                    vertex4_scale(_tmp.vector, _conj.w));

        Vertex4 _tmp_cross = vertex4_cross(_tmp.vector, _conj.vector);

        result = vertex4_add(_tmp2, _tmp_cross);
        return result;
}

static void vertex4_view_to_screen(Vertex4 *buffer,
                                   GameCamera *camera,
                                   int buffer_width,
                                   int buffer_height,
                                   uint32_t count) {

    float aspect_ratio = (float)buffer_width / (float)buffer_height;
    for (uint32_t i = 0; i < count; ++i) {

        __m128 _vertex_x = buffer[i].x;
        __m128 _vertex_y = buffer[i].y;
        __m128 _vertex_z = buffer[i].z;

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

        buffer[i].x = _mm_or_ps(_float_res_x, _mask);
        buffer[i].y = _mm_or_ps(_float_res_y, _mask);
        buffer[i].z = _mm_or_ps(_vertex_z, _mask);
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

static inline V3 compute_light_intensity(V3 light_position, V3 light_color, V3 vertex, V3 normal) {
    V3 dl = light_position-vertex; 
    float intensity = (v3_dot(v3_norm(dl), v3_norm(normal))) /
                      (0.0005f*v3_mag(dl));
    if (intensity > 0) {
        return intensity*light_color;
    }
    return light_color = {};
}

static inline float rgb_to_luminance(V3 rgb) {
    return v3_dot(rgb, v3(0.2126f, 0.7152f, 0.0722f));
}

static inline V3 triangle_get_attribute_from_buffer(uint32_t index, Vertex4 *buffer, uint32_t buffer_size) {
        assert(index/4 < buffer_size);
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

    // sample the radius of the light source at different x values to create a sphere to rder
    {

    }
    render_state->draw_count = 0;
    // transform vertices to view space
    { 
        Quaternion rotation = rotation_q4(-camera->theta_x, v3(1,0,0))*
                              rotation_q4(-camera->theta_y, v3(0,1,0));

        Vertex4 _cam = {};
        _cam.x = _mm_set1_ps(camera->pos.x);
        _cam.y = _mm_set1_ps(camera->pos.y);
        _cam.z = _mm_set1_ps(camera->pos.z);

        Quaternion4 _q4 = {};
        _q4.w = _mm_set1_ps(rotation.scalar);
        _q4.vector.x = _mm_set1_ps(rotation.vector.x);
        _q4.vector.y = _mm_set1_ps(rotation.vector.y);
        _q4.vector.z = _mm_set1_ps(rotation.vector.z);

        for (uint32_t vertex = 0; vertex < render_state->vertex_count; ++vertex) {
            Vertex4 _v4 = render_state->vertex_buffer[vertex];
            Vertex4 *_out = &render_state->vertex_out_buffer[vertex];
            // translate position in world relative to camera
            _v4 = vertex4_sub(_v4, _cam);

            //rotate world relative to camera
            *_out = vertex4_rotate(_v4, _q4);
        }
#if 1

        for (uint32_t normal = 0; normal < render_state->normal_count; ++normal) {
            Vertex4 *_n4 = &render_state->normal_buffer[normal];
            Vertex4 *_out = &render_state->normal_out_buffer[normal];
            *_out = vertex4_rotate(*_n4, _q4);
        }
#endif

    //
    }

    // backface culling
    for (uint32_t triangle = 0; triangle < render_state->polygon_count; ++triangle) {
        Triangle current = render_state->polygons[triangle];
        V3 vertex1 = triangle_get_attribute_from_buffer(current.v1, render_state->vertex_out_buffer, 
                                                       render_state->vertex_count);
        V3 vertex2 = triangle_get_attribute_from_buffer(current.v2, render_state->vertex_out_buffer, 
                                                       render_state->vertex_count);
        V3 vertex3 = triangle_get_attribute_from_buffer(current.v3, render_state->vertex_out_buffer, 
                                                       render_state->vertex_count);

#if 1
        if (v3_dot(vertex1, v3_cross(vertex2-vertex1, vertex3-vertex1)) >= 0) {
        }
#endif
            render_state->polygons_to_draw[triangle] = current;
    }

    // change light sources to be 100 times the brightness for post processing
        for (uint32_t i = 0; i < render_state->light_sources_count; ++i) {
            LightSource light_source = render_state->light_sources[i]; 
            render_state->light_colors[i] = light_source.color;
            render_state->light_positions[i] = renderer_world_vertex_to_view(light_source.position, camera);


            for (uint32_t j = light_source.obj->index_start; j < light_source.obj->index_end; ++j) {
                Triangle *polygon = &render_state->polygons_to_draw[j];
                polygon->color = v3(100,100,100)+polygon->color;
            }
        }

    // perspective projection & clipping
    vertex4_view_to_screen(render_state->vertex_out_buffer, 
                           camera, 
                           buffer_width, buffer_height, 
                           render_state->vertex_count);
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
            *color = 100.0f*(*color);
            tone_map_color(color, white_point);
            *color += v3(1,1,1);
        }
}

static inline uint32_t convert_v3_to_RGB(V3 vector) {
    //clamp to 1
    vector.x > 1 ? vector.x = 1: vector.x;
    vector.y > 1 ? vector.y = 1: vector.y;
    vector.z > 1 ? vector.z = 1: vector.z;
    uint32_t result = (((uint8_t)(0xFF*vector.x)) << 16) | 
                      (((uint8_t)(0xFF*vector.y))) << 8 | 
                      (uint8_t)(0xFF*vector.z);
    return result;
}

static void renderer_draw_line(OffscreenBuffer *buffer,
                               OffscreenBuffer *zbuffer,
                               OffscreenBuffer *normal_buffer,
                               V3 start,
                               V3 end,
                               V3 start_normal,
                               V3 end_normal,
                               V3 color) {


    clamp_to_screen(&start, buffer->width, buffer->height);
    clamp_to_screen(&end, buffer->width, buffer->height);
    if (start.x > end.x) {
        v3_swap(&start, &end);
        v3_swap(&start_normal, &end_normal);
    }

    float dx = (floorf(end.x) - floorf(start.x));
    float dy = (floorf(end.y) - floorf(start.y));
    //
    // case where the gradient is below 1;

    float step = 0;
    if (fabsf(dx) > fabsf(dy)) {
        step = (float)fabsf(dx);
    }

    else {
        step = (float)fabsf(dy);
    }

    if (step == 0) {return;};
    uint8_t bright = 0x0;
    if (color.x > 1 && color.y > 1 && color.z > 1) {
        bright = 0x01;
        color -= v3(1,1,1);
    }

    /* I am using the padding byte in the windows bitmap buffer to seperate 
     * the bright objects for post processing all bright pixels go through the draw line call
     * with a color of > 1 this has to be then corrected and then return either 0xFF or 0 for 
     * bright and not bright respectively
     */

    V3 normal_inc = (end_normal - start_normal) / step;
    V3 start_inc = (end-start)/step;

    for (int draw = 0; draw <= step; ++draw) {
        int offset = (int)(floorf(start.x))*buffer->bytes_per_pixel + 
            (int)(floorf(start.y))*buffer->pitch;
        float *depth_value = (float *)((uint8_t *)zbuffer->memory + offset);

        if ((start.z < *depth_value)) {
            uint32_t pixel_color = convert_v3_to_RGB(color);
            pixel_color = ((bright) << 24) | pixel_color;

            // normalise between 0 and 1
            //
            V3 normalised = v3_norm(start_normal);
            normalised.x = roundf((start_normal.x*0.5f+0.5f)*65535.0f);
            normalised.y = roundf((start_normal.y*0.5f+0.5f)*65535.0f);
            uint32_t normal_value = (((uint16_t)(normalised.x)) << 16) | 
                                    (((uint16_t)(normalised.y)));


            uint32_t *pixel = (uint32_t *)((uint8_t *)buffer->memory + offset);
            uint32_t *normal = (uint32_t *)((uint8_t *)normal_buffer->memory + offset);

            *pixel = pixel_color;
            *normal = normal_value;
            *depth_value = start.z;
        }

        start += start_inc;
        start_normal += normal_inc;
    }
}

static void renderer_draw_flat_top_triangle(OffscreenBuffer *buffer,
                                            OffscreenBuffer *zbuffer,
                                            OffscreenBuffer *normal_buffer,
                                            V3 left,
                                            V3 right,
                                            V3 bottom,
                                            V3 left_normal,
                                            V3 right_normal,
                                            V3 bottom_normal,
                                            V3 color) {
        
    float dy = floorf(bottom.y)-floorf(right.y);

    V3 right_inc = v3((bottom.x-right.x)/dy, 1, -(bottom.z-right.z)/dy);
    V3 left_inc = v3((bottom.x-left.x)/dy, 1, -(bottom.z-left.z)/dy);
    V3 left_normal_inc = (bottom_normal - left_normal)/dy;
    V3 right_normal_inc = (bottom_normal - right_normal)/dy;

    for (int y = 0; y <= fabsf(dy); ++y) {
        renderer_draw_line(buffer, zbuffer, normal_buffer, left, right, 
                           left_normal, right_normal, color);

        left += left_inc;
        right += right_inc;
        left_normal += left_normal_inc;
        right_normal += right_normal_inc;
    }
}

static void renderer_draw_flat_bottom_triangle(OffscreenBuffer *buffer,
                                               OffscreenBuffer *zbuffer,
                                               OffscreenBuffer *normal_buffer,
                                               V3 left,
                                               V3 right,
                                               V3 top,
                                               V3 left_normal,
                                               V3 right_normal,
                                               V3 top_normal,
                                               V3 color) {

    float dy = floorf(top.y)-floorf(right.y);

    V3 left_inc = v3((top.x-left.x)/dy, 1, (top.z-left.z)/dy);
    V3 right_inc = v3((top.x-right.x)/dy, 1, (top.z-right.z)/dy);
    V3 right_normal_inc = (top_normal - right_normal)/dy;
    V3 left_normal_inc = (top_normal - left_normal)/dy;

    for (int y = 0; y <= fabsf(dy); ++y) {
        renderer_draw_line(buffer, zbuffer, normal_buffer, left, right, 
                           left_normal, right_normal, color);

        left -= left_inc;
        right -= right_inc;
        left_normal -= left_normal_inc;
        right_normal -= right_normal_inc;
    }
}

static void renderer_draw_triangles(OffscreenBuffer *buffer,
                                    OffscreenBuffer *zbuffer,
                                    OffscreenBuffer *normal_buffer,
                                    OffscreenBuffer *postbuffer,
                                    GameCamera *camera,
                                    RendererState *render_state) {

    V3 white_point = v3(100,100,100);
    // Draw Triangles
    {
        Triangle *triangles = render_state->polygons_to_draw;
        Vertex4 *vertices = render_state->vertex_out_buffer;
        Vertex4 *normals = render_state->normal_out_buffer;

        for (uint32_t i = 0; i < render_state->polygon_count; ++i) {

            V3 vert1 = triangle_get_attribute_from_buffer(triangles[i].v1, vertices, render_state->vertex_count);
            V3 vert2 = triangle_get_attribute_from_buffer(triangles[i].v2, vertices, render_state->vertex_count);
            V3 vert3 = triangle_get_attribute_from_buffer(triangles[i].v3, vertices, render_state->vertex_count);

            V3 normal1 = triangle_get_attribute_from_buffer(triangles[i].vn1, normals, render_state->normal_count);
            V3 normal2 = triangle_get_attribute_from_buffer(triangles[i].vn2, normals, render_state->normal_count);
            V3 normal3 = triangle_get_attribute_from_buffer(triangles[i].vn3, normals, render_state->normal_count);

            V3 color = triangles[i].color;
            color_correct_brightness(&color, white_point);

            if (renderer_v3_isclipped(vert1)) {continue;}
            if (renderer_v3_isclipped(vert2)) {continue;}
            if (renderer_v3_isclipped(vert3)) {continue;}


            //sort vertexes v3 is the biggest
            if (vert3.y > vert2.y) {
                v3_swap(&vert3, &vert2);
                v3_swap(&normal3, &normal2);

            }

            if (vert2.y > vert1.y) {
                v3_swap(&vert2, &vert1);
                v3_swap(&normal2, &normal1);
            }

            if (vert3.y > vert2.y) {
                v3_swap(&vert3, &vert2);
                v3_swap(&normal3, &normal2);
            }

            if (vert3.y == vert2.y) {
                if (vert3.x > vert2.x) {
                    renderer_draw_flat_top_triangle(buffer, zbuffer, normal_buffer, vert2, vert3, vert1, 
                                                    normal2, normal3, normal1, color);
                }
                else {
                    renderer_draw_flat_top_triangle(buffer, zbuffer, normal_buffer, vert3, vert2, vert1, 
                                                    normal3, normal2, normal1, color);
                }

                continue;
            }

            if (vert1.y == vert2.y) {
                if (vert1.x > vert2.x) {
                    renderer_draw_flat_bottom_triangle(buffer, zbuffer, normal_buffer, vert2, vert1, vert3, 
                                                       normal2, normal1, normal3, color);
                }
                else {
                    renderer_draw_flat_bottom_triangle(buffer, zbuffer, normal_buffer, vert1, vert2, vert3, 
                                                       normal1, normal2, normal3, color);
                }
                continue;
            }

            float lerp_factor = vert2.y-vert1.y;
            float delta_x = vert3.y-vert1.y;

            V3 vert4;
            vert4.x = f_lerp(delta_x, lerp_factor, vert1.x, vert3.x);
            vert4.y = vert2.y;
            vert4.z = f_lerp(delta_x, lerp_factor, vert1.z, vert3.z);

            V3 normal4 = v3_lerp(delta_x, lerp_factor, normal1, normal3);

            if (vert4.x > vert2.x) {
                renderer_draw_flat_top_triangle(buffer, zbuffer, normal_buffer, vert2, vert4, vert1, 
                                                normal2, normal4, normal1, color);
                renderer_draw_flat_bottom_triangle(buffer, zbuffer, normal_buffer, vert2, vert4, vert3, 
                                                   normal2, normal4, normal3, color);
            }

            else {
                renderer_draw_flat_top_triangle(buffer, zbuffer, normal_buffer, vert4, vert2, vert1, 
                                                normal4, normal2, normal1, color);

                renderer_draw_flat_bottom_triangle(buffer, zbuffer, normal_buffer, vert4, vert2, vert3, 
                                                   normal4, normal2, normal3, color);
            }
        }
    }

    // Lighting Pass
#if 1
    {
        uint32_t *colors = (uint32_t *)buffer->memory;
        uint32_t *normals = (uint32_t *)normal_buffer->memory;
        float *depth = (float *)zbuffer->memory;
        V3 *light_positions = render_state->light_positions;
        V3 *light_colors = render_state->light_colors;

        for (int pixel = 0; pixel < buffer->width*buffer->height; pixel+=4) {
            if (colors[pixel] & 0x01000000) {
                continue;
            }

            if (depth[pixel] == FLT_MAX) {
                continue;
            }

            Vertex4 _extracted_norm = {};
            for (int i = 0; i < 4; ++i) {
                ((float *)&_extracted_norm.x)[i] = (((uint16_t)((normals[pixel+i] & 0xFFFF0000) >> 16)/65535.0f)*2) - 1;
                ((float *)&_extracted_norm.y)[i] = (((uint16_t)(normals[pixel+i] & 0x0000FFFF)/65535.0f)*2) - 1;

            }
            __m128 _sum_sqrs = _mm_sub_ps(_mm_set1_ps(1.0f), _mm_add_ps(_mm_mul_ps(_extracted_norm.x, _extracted_norm.x),   
                        _mm_mul_ps(_extracted_norm.y, _extracted_norm.y)));

            _extracted_norm.z = _mm_mul_ps(_mm_set1_ps(-1.0f), _mm_mul_ps(_mm_rsqrt_ps(_sum_sqrs), _sum_sqrs));

            Vertex4 _extracted_pos = {};
            for (int i = 0; i < 4; ++i) {
                ((float *)&_extracted_pos.x)[i] = (float)((pixel+i) % 1280); 
                ((float *)&_extracted_pos.y)[i] = (float)((pixel+i) / 1280);
                ((float *)&_extracted_pos.z)[i] = (depth[pixel+i]);
            }

            // transform back to view space using the z component
            {
                __m128 _length_x = _mm_set1_ps(tanf(camera->fov/2.0f)*camera->znear *
                        ((float)buffer->width/(float)buffer->height));

                __m128 _length_y = _mm_set1_ps(tanf(camera->fov/2.0f)*camera->znear);

                __m128 _buffer_width = _mm_set1_ps(buffer->width/2.0f);
                __m128 _buffer_height = _mm_set1_ps(buffer->height/2.0f);

                __m128 _znear = _mm_set1_ps(camera->znear);

                _extracted_pos.x = _mm_div_ps(_mm_sub_ps(_extracted_pos.x, _buffer_width), (_buffer_width));
                _extracted_pos.x = _mm_div_ps(_mm_sub_ps(_extracted_pos.y, _buffer_height), (_buffer_height));
                _extracted_pos.y = _mm_mul_ps(_extracted_pos.y, _mm_set1_ps(-1.0f)); 

                _extracted_pos.x = _mm_mul_ps(_mm_mul_ps(_extracted_pos.x,_length_x), 
                        _mm_div_ps(_extracted_pos.z, _znear));

                _extracted_pos.x = _mm_mul_ps(_mm_mul_ps(_extracted_pos.y,_length_y), 
                        _mm_div_ps(_extracted_pos.z, _znear));
            }


            Vertex4 _extracted_color {};
            for (int i = 0; i < 4; ++i) {
                ((float *)&_extracted_color.x)[i] = ((colors[pixel+i] & 0x00FF0000) >> 16)/255.0f;
                ((float *)&_extracted_color.y)[i] = ((colors[pixel+i] & 0x0000FF00) >> 8)/255.0f;
                ((float *)&_extracted_color.z)[i] = (colors[pixel+i] & 0x000000FF)/255.0f;
            }
            Vertex4 _color = {};

            // compute the light intensity on every pixel for each light and blend additively 
            for (uint32_t light = 0; light < render_state->light_sources_count; ++light) {
                Vertex4 _light_position = {};
                _light_position.x = _mm_set1_ps(light_positions[light].x);
                _light_position.y = _mm_set1_ps(light_positions[light].y);
                _light_position.z = _mm_set1_ps(light_positions[light].z);

                Vertex4 _light_color = {};
                _light_color.x = _mm_set1_ps(light_colors[light].x);
                _light_color.y = _mm_set1_ps(light_colors[light].y);
                _light_color.z = _mm_set1_ps(light_colors[light].z);

                Vertex4 _dl = vertex4_sub(_light_position, _extracted_pos);
                __m128 _intensity = _mm_div_ps(vertex4_dot(vertex4_norm(_dl), _extracted_norm),
                        _mm_mul_ps(_mm_set1_ps(0.0005f), 
                            _mm_sqrt_ps(vertex4_dot(_dl, _dl))));
                __m128 _mask = _mm_cmpgt_ps(_intensity, _mm_setzero_ps());
                _intensity = _mm_and_ps(_intensity, _mask);
                _color = vertex4_add(_color, vertex4_scale(_light_color, _intensity));
            }

            // tone map colors
            {
                _color = vertex4_pairwise_mul(_color, _extracted_color);
                Vertex4 _luminance_const = {};
                _luminance_const.x = _mm_set1_ps(0.2126f);
                _luminance_const.y = _mm_set1_ps(0.7152f);
                _luminance_const.z = _mm_set1_ps(0.0722f);
                Vertex4 _white_point = {};
                _white_point.x = _mm_set1_ps(white_point.x);
                _white_point.y = _mm_set1_ps(white_point.y);
                _white_point.z = _mm_set1_ps(white_point.z);

                __m128 _luminance = vertex4_dot(_color, _luminance_const);
                __m128 _white_l = vertex4_dot(_white_point, _luminance_const);
                __m128 _numerator = _mm_mul_ps(_luminance, 
                        _mm_add_ps(_mm_set1_ps(1.0f), 
                            _mm_div_ps(_luminance, _mm_mul_ps(_white_l, _white_l))));

                __m128 _new_luminance = _mm_div_ps(_numerator, _mm_add_ps(_mm_set1_ps(1.0f), _luminance));
                _color = vertex4_scale(_color ,_mm_div_ps(_new_luminance, _luminance));
            }

            // unpack into the main frambuffer
            for (int i = 0; i < 4; ++i) {
                V3 color = {};
                color.x = ((float *)&_color.x)[i];
                color.y = ((float *)&_color.y)[i];
                color.z = ((float *)&_color.z)[i];
                colors[pixel+i] = convert_v3_to_RGB(color);
            }
        }
    }
#endif

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

                if (sumbright < 0x01) {
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
        for (int i = 0; i < 5; ++i) {
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


static RenderObj renderer_render_obj_create(RendererState *render_state, Mesh *mesh, V3 color, V3 scale, Quaternion rotation, V3 translation) {
    RenderObj result = {};
    uint32_t vertex_start = render_state->vertex_count;
    assert(mesh->vert_count % 4 == 0);

    for (uint32_t vertex = 0; vertex < mesh->vert_count; vertex+=4) {
    assert(vertex+3 < mesh->vert_count);

        V3 vert1 = v3_rotate_q4(v3_pairwise_mul(scale,mesh->vertices[vertex]), rotation); 
        V3 vert2 = v3_rotate_q4(v3_pairwise_mul(scale,mesh->vertices[vertex+1]), rotation); 
        V3 vert3 = v3_rotate_q4(v3_pairwise_mul(scale,mesh->vertices[vertex+2]), rotation); 
        V3 vert4 = v3_rotate_q4(v3_pairwise_mul(scale,mesh->vertices[vertex+3]), rotation);


        vert1 += translation;
        vert2 += translation;
        vert3 += translation;
        vert4 += translation;

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

        ++render_state->vertex_count;
    }


    int normal_start = render_state->normal_count;
    assert(mesh->vertexn_count % 4 == 0);
    for (uint32_t normal = 0; normal < mesh->vertexn_count; normal+=4) {
        assert(normal+3 < mesh->vertexn_count);

        V3 vert1 = v3_rotate_q4(mesh->vertexn[normal], rotation);
        V3 vert2 = v3_rotate_q4(mesh->vertexn[normal+1], rotation);
        V3 vert3 = v3_rotate_q4(mesh->vertexn[normal+2], rotation);
        V3 vert4 = v3_rotate_q4(mesh->vertexn[normal+3], rotation);
 
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
        render_state->polygons[render_state->polygon_count].vn1 += normal_start*4;
        render_state->polygons[render_state->polygon_count].vn2 += normal_start*4;
        render_state->polygons[render_state->polygon_count].vn3 += normal_start*4;
        render_state->polygons[render_state->polygon_count].color = color;
        render_state->polygon_count++;
    }


    result.vstart = vertex_start;
    result.vend = render_state->vertex_count;
    result.nstart = normal_start;
    result.nend = render_state->normal_count;
    result.index_start = index_start;
    result.index_end = index_start+mesh->poly_count;
    result.mesh = mesh;
    result.color = color;
    assert(render_state->vertex_count < arraysize(render_state->vertex_buffer));
    assert(render_state->normal_count < arraysize(render_state->normal_buffer));

    return result;
}


static void renderer_render_obj_update(RendererState *render_state, 
                                       RenderObj *obj, 
                                       V3 scale, 
                                       Quaternion rotation, 
                                       V3 translation) {

    assert(render_state->vertex_count < arraysize(render_state->vertex_buffer));
    assert(obj->mesh->vert_count % 4 == 0);
    assert(obj->mesh->vertexn_count % 4 == 0);

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
        
    }

    for (uint32_t normal = 0; normal < obj->mesh->vertexn_count; normal+=4) {
        assert(normal+3 < obj->mesh->vertexn_count);

        V3 normal1 = v3_rotate_q4(obj->mesh->vertexn[normal], rotation);
        V3 normal2 = v3_rotate_q4(obj->mesh->vertexn[normal+1], rotation);
        V3 normal3 = v3_rotate_q4(obj->mesh->vertexn[normal+2], rotation);
        V3 normal4 = v3_rotate_q4(obj->mesh->vertexn[normal+3], rotation);

        ((float * )&render_state->normal_buffer[obj->nstart+normal/4].x)[0] = normal1.x; 
        ((float * )&render_state->normal_buffer[obj->nstart+normal/4].x)[1] = normal2.x; 
        ((float * )&render_state->normal_buffer[obj->nstart+normal/4].x)[2] = normal3.x; 
        ((float * )&render_state->normal_buffer[obj->nstart+normal/4].x)[3] = normal4.x; 

        ((float * )&render_state->normal_buffer[obj->nstart+normal/4].y)[0] = normal1.y; 
        ((float * )&render_state->normal_buffer[obj->nstart+normal/4].y)[1] = normal2.y; 
        ((float * )&render_state->normal_buffer[obj->nstart+normal/4].y)[2] = normal3.y; 
        ((float * )&render_state->normal_buffer[obj->nstart+normal/4].y)[3] = normal4.y; 

        ((float * )&render_state->normal_buffer[obj->nstart+normal/4].z)[0] = normal1.z; 
        ((float * )&render_state->normal_buffer[obj->nstart+normal/4].z)[1] = normal2.z; 
        ((float * )&render_state->normal_buffer[obj->nstart+normal/4].z)[2] = normal3.z; 
        ((float * )&render_state->normal_buffer[obj->nstart+normal/4].z)[3] = normal4.z; 
    }

    for (uint32_t polygon = 0; polygon < obj->mesh->poly_count; ++polygon) {
        render_state->polygons[obj->index_start+polygon].color = obj->color;
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

static void renderer_light_source_create(RendererState *render_state, RenderObj *obj, V3 position, V3 color) {
        render_state->light_sources[render_state->light_sources_count].position = position;
        render_state->light_sources[render_state->light_sources_count].color = color;
        render_state->light_sources[render_state->light_sources_count].obj = obj;
        ++render_state->light_sources_count;
}


static void renderer_draw(RendererState *render_state, 
                          GameCamera *camera, 
                          OffscreenBuffer *buffer, 
                          OffscreenBuffer *zbuffer,
                          OffscreenBuffer *normal_buffer,
                          OffscreenBuffer *postbuffer) {

    renderer_transform_light_and_cull(render_state, camera, buffer->width, buffer->height);
    renderer_draw_triangles(buffer, zbuffer, normal_buffer, postbuffer, camera, render_state);
}
