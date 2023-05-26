#include <emmintrin.h>
#include <stdint.h>
#include "renderer.h"

static V3 renderer_world_vertex_to_view(V3 world_pos, GameCamera *camera) {
    V3 result = v3_rotate_q4(world_pos-camera->pos, 
                          (rotation_q4(-camera->theta_x, v3(1,0,0))*
                          rotation_q4(-camera->theta_y, v3(0,1,0))));

    return result;
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
        // ambient light
        V3 v1_light = {};
        V3 v2_light = {};
        V3 v3_light = {};

        for (uint32_t light = 0; light < render_state->light_sources_count; ++light) {
            LightSource source = render_state->light_sources[light];
            v1_light += compute_light_intensity(source, vertex1, vn1);
            v2_light += compute_light_intensity(source, vertex2, vn2);
            v3_light += compute_light_intensity(source, vertex3, vn3);
        }

        *v1_color = v3_pairwise_mul(*v1_color, v1_light);
        *v2_color = v3_pairwise_mul(*v2_color, v2_light);
        *v3_color = v3_pairwise_mul(*v3_color, v3_light);
#endif // SHADING
    }

    // change light sources to just be their original color
    for (uint32_t i = 0; i < render_state->light_sources_count; ++i) {
        LightSource light_source = render_state->light_sources[i]; 
        for (uint32_t j = light_source.obj->vstart; j <= light_source.obj->vend; ++j) {
            V3 vertex = render_state->vertex_list[j];
            V3 *color = &render_state->vertex_colors[j];
            *color = v3_pairwise_mul(v3(1,1,1), light_source.color) + 100*v3(1,1,1);
        }
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

        if (v3_should_clip(relative_pos, camera, aspect_ratio)) {
            render_state->vertex_list[vertex] = v3(FLT_MAX,FLT_MAX,FLT_MAX);
            continue;
        }
        
        float normal_x = (relative_pos.x/relative_pos.z) / 
                         (tanf(camera->fov/2.0f)*aspect_ratio);

        float normal_y = (relative_pos.y/relative_pos.z) / 
                         (tanf(camera->fov/2.0f));

        if (normal_x <= -1) {
            render_state->vertex_list[vertex] = v3(FLT_MAX,FLT_MAX,FLT_MAX);
            continue;
        }

        if (normal_y <= -1) {
            render_state->vertex_list[vertex] = v3(FLT_MAX,FLT_MAX,FLT_MAX);
            continue;
        }

        if (normal_x >= 1) {
            render_state->vertex_list[vertex] = v3(FLT_MAX,FLT_MAX,FLT_MAX);
            continue;
        }

        if (normal_y >= 1) {
            render_state->vertex_list[vertex] = v3(FLT_MAX,FLT_MAX,FLT_MAX);
            continue;
        }

        float x = floorf(((normal_x * (buffer_width/2.0f) + buffer_width/2.0f)));
        float y = floorf(((-normal_y * (buffer_height/2.0f) + buffer_height/2.0f)));

        render_state->vertex_list[vertex].x = x;
        render_state->vertex_list[vertex].y = y;
        render_state->vertex_list[vertex].z = relative_pos.z;
    }
}

static inline bool renderer_v3_isclipped(V3 screen) {
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

static void clamp_to_screen(V3 *vertex) {
    if (vertex->x < 0) {
        vertex->x = 0;
    }

    else if (vertex->x >= 1280) {
        vertex->x = 1279;
    }

    if (vertex->y < 0) {
        vertex->y = 0;
    }

    else if (vertex->y >= 720) {
        vertex->y = 719;
    }
}
static void renderer_draw_line(OffscreenBuffer *buffer,
                               OffscreenBuffer *zbuffer,
                               V3 start,
                               V3 end,
                               V3 start_color,
                               V3 end_color) {


    clamp_to_screen(&start);
    clamp_to_screen(&end);
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
                                           V3 *vertices,
                                           V3 *colors,
                                           Triangle *triangles,
                                           int count) {

    V3 white_point = v3(100,100,100);
    for (int i = 0; i < count; ++i) {

        V3 vert1 = vertices[triangles[i].v1];
        V3 vert2 = vertices[triangles[i].v2];
        V3 vert3 = vertices[triangles[i].v3];
        V3 v1_color = colors[triangles[i].v1];
        V3 v2_color = colors[triangles[i].v2];
        V3 v3_color = colors[triangles[i].v3];


        color_correct_brightness(&v1_color, white_point);
        color_correct_brightness(&v2_color, white_point);
        color_correct_brightness(&v3_color, white_point);

        if (renderer_v3_isclipped(vert1)) continue;
        if (renderer_v3_isclipped(vert2)) continue;
        if (renderer_v3_isclipped(vert3)) continue;

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


    // POST PROCESSING BLOOM EFFECT
    for (int y = 0; y < buffer->height; ++y) {
        for(int x = 0; x < buffer->width; ++x) {
            uint32_t *zpixel = (uint32_t *)zbuffer->memory+(y*buffer->width)+x;
            uint32_t *bpixel = (uint32_t *)buffer->memory+y*buffer->width+x;
            uint8_t bright = (*bpixel >> 24);
            if (bright == 0xFF) {
                *zpixel = *bpixel;
            }

            else {
                *(zpixel) = 0;
            }
        }
    }

    for (int i = 0; i < 5; ++i) {
        uint32_t *pixels = (uint32_t *)zbuffer->memory;
        for (int y = 0; y < zbuffer->height; ++y) {
            if (y < 1 || y + 1 == zbuffer->height) {
                continue;
            }

            int sumr = (((pixels[(y+1)*zbuffer->width]) & 0x00FF0000) >> 16) +
                (((pixels[(y)*zbuffer->width]) & 0x00FF0000) >> 16) +
                (((pixels[(y-1)*zbuffer->width]) & 0x00FF0000) >> 16); 

            int sumg = (((pixels[(y+1)*zbuffer->width]) & 0x0000FF00) >> 8) +
                (((pixels[(y)*zbuffer->width]) & 0x0000FF00) >> 8) +
                (((pixels[(y-1)*zbuffer->width]) & 0x0000FF00) >> 8); 

            int sumb = ((pixels[(y+1)*zbuffer->width]) & 0x000000FF) +
                ((pixels[(y)*zbuffer->width]) & 0x000000FF) +
                ((pixels[(y-1)*zbuffer->width]) & 0x000000FF); 

            int counter = 3;
            int leftr = 0;
            int leftg = 0;
            int leftb = 0;

            for (int x = 0; x < zbuffer->width; ++x) {
                if ((x + 1) < zbuffer->width) {

                    sumr += (((pixels[(y-1)*zbuffer->width+x+1] & 0x00FF0000) >> 16) + 
                            ((pixels[(y)*zbuffer->width+x+1] & 0x00FF0000) >> 16) +
                            ((pixels[(y+1)*zbuffer->width+x+1] & 0x00FF0000) >> 16));

                    sumg += (((pixels[(y-1)*zbuffer->width+x+1] & 0x0000FF00) >> 8) + 
                            ((pixels[(y)*zbuffer->width+x+1] & 0x0000FF00) >> 8) +
                            ((pixels[(y+1)*zbuffer->width+x+1] & 0x0000FF00) >> 8));

                    sumb += ((pixels[(y-1)*zbuffer->width+x+1] & 0x000000FF) + 
                            (pixels[(y)*zbuffer->width+x+1] & 0x000000FF) +
                            (pixels[(y+1)*zbuffer->width+x+1] & 0x000000FF));

                    counter += 3;
                }

                uint8_t red = (uint8_t)(sumr/counter); 
                uint8_t green = (uint8_t)(sumg/counter); 
                uint8_t blue = (uint8_t)(sumb/counter); 

                sumr -= leftr;
                sumg -= leftg;
                sumb -= leftb;

                leftr = (((pixels[(y-1)*zbuffer->width+x] & 0x00FF0000) >> 16) + 
                        ((pixels[(y)*zbuffer->width+x] & 0x00FF0000) >> 16) +
                        ((pixels[(y+1)*zbuffer->width+x] & 0x00FF0000) >> 16));

                leftg = (((pixels[(y-1)*zbuffer->width+x] & 0x0000FF00) >> 8) + 
                        ((pixels[(y)*zbuffer->width+x] & 0x0000FF00) >> 8) +
                        ((pixels[(y+1)*zbuffer->width+x] & 0x0000FF00) >> 8));

                leftb = ((pixels[(y-1)*zbuffer->width+x] & 0x000000FF) + 
                        (pixels[(y)*zbuffer->width+x] & 0x000000FF) +
                        (pixels[(y+1)*zbuffer->width+x] & 0x000000FF));

                if (x >= 1) {
                    counter -= 3;
                }

                pixels[(y)*zbuffer->width+x] = red << 16 | green << 8 | blue;

            }
        }
    }

    for (int y = 0; y < buffer->height; ++y) {
        for (int x = 0; x < buffer->width; ++x) {
            uint32_t *zpixel = (uint32_t *)zbuffer->memory+y*buffer->width+x;
            uint32_t *bpixel = (uint32_t *)buffer->memory+y*buffer->width+x;
            uint8_t bred = (uint8_t)((*bpixel & 0x00FF0000) >> 16);
            uint8_t bgreen = (uint8_t)((*bpixel & 0x0000FF00) >> 8);
            uint8_t bblue = (uint8_t)(*bpixel & 0x000000FF);

            uint8_t zred = (uint8_t)((*zpixel & 0x00FF0000) >> 16);
            uint8_t zgreen = (uint8_t)((*zpixel & 0x0000FF00) >> 8);
            uint8_t zblue = (uint8_t)(*zpixel & 0x000000FF);

            bred + zred > 255 ? zred = 255 : zred += bred;
            bgreen + zgreen > 255 ? zgreen = 255 : zgreen += bgreen;
            bblue + zblue > 255 ? zblue = 255 : zblue += bblue;

            *bpixel = (zred << 16) | (zgreen << 8) | zblue;
        }
    }
}


static void renderer_draw(RendererState *render_state, 
                          GameCamera *camera, 
                          OffscreenBuffer *buffer, 
                          OffscreenBuffer *zbuffer) {

    renderer_transform_light_and_cull(render_state, camera, buffer->width, buffer->height);
    renderer_draw_triangles_filled(buffer, zbuffer, 
                                   render_state->vertex_list, 
                                   render_state->vertex_colors, 
                                   render_state->polygons_to_draw, 
                                   render_state->draw_count);
}

