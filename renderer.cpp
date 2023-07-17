#include <stdint.h>
#include <emmintrin.h>
#include <xmmintrin.h>
#include "renderer.h"
#if DEBUG_MODE
#define assert(expression) if(!(expression)) {*(int *)0 = 0;}
#else
#define assert(exppression)
#endif
#define WORLD_LIGHT 0

#define arraysize(array) (sizeof(array) / sizeof((array)[0]))

static V3 renderer_world_vertex_to_view(V3 world_pos, GameCamera *camera) {
    V3 result = v3_rotate_q4(world_pos-camera->pos,
                             (rotation_q4(-camera->theta_x, v3(1,0,0))*
                             rotation_q4(-camera->theta_y, v3(0,1,0))));

    return result;
}

static V3 renderer_vertex_invalid() {
    V3 result;
    *(uint32_t *)&result.x = 0xFFFFFFFF;
    *(uint32_t *)&result.y = 0xFFFFFFFF;
    *(uint32_t *)&result.z = 0xFFFFFFFF;

    return result;
}

static inline bool renderer_v3_is_invalid(V3 a, int buffer_width, int buffer_height){
    if ((a.x > buffer_width) || (a.y > buffer_height)) {return true;}
    if ((a.x < 0) || (a.y < 0)) {return true;}
    return false;
}

static Vertex4 vertex4_rotate(Vertex4 _vertex, Quaternion rotation) {
        Vertex4 result = {};

        Quaternion4 _q = {};
        _q.w = _mm_set1_ps(rotation.scalar);
        _q.vector.x = _mm_set1_ps(rotation.vector.x);
        _q.vector.y = _mm_set1_ps(rotation.vector.y);
        _q.vector.z = _mm_set1_ps(rotation.vector.z);

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

static inline float rgb_to_luminance(V3 rgb) {
    return v3_dot(rgb, v3(0.2126f, 0.7152f, 0.0722f));
}

static inline V3 triangle_get_attribute_from_buffer(uint32_t index, 
                                                    Vertex4 *buffer, 
                                                    uint32_t buffer_size) {
        assert(index/4 < buffer_size);
        V3 result = v3(((float *)&buffer[index/4].x)[index%4],
                      ((float *)&buffer[index/4].y)[index%4],
                      ((float *)&buffer[index/4].z)[index%4]);

        return result;
}


static inline void clamp_x(V3 *vertex, int buffer_width) {
    if (vertex->x < 0) {
        vertex->x = 0;
        return;
    }

    if (vertex->x >= buffer_width) {
        vertex->x = (float)(buffer_width - 1);
    }
}

static inline void clamp_y(V3 *vertex, int buffer_height) {
    if (vertex->y < 0) {
        vertex->y = 0;
        return;
    }

    if (vertex->y >= buffer_height) {
        vertex->y = (float)(buffer_height - 1);
    }
}

static inline uint32_t convert_v3_to_RGB(V3 vector) {
    //clamp to 1
    vector.x > 1 ? vector.x = 1: vector.x;
    vector.y > 1 ? vector.y = 1: vector.y;
    vector.z > 1 ? vector.z = 1: vector.z;

    vector.x < 0 ? vector.x = 0: vector.x;
    vector.y < 0 ? vector.y = 0: vector.y;
    vector.z < 0 ? vector.z = 0: vector.z;

    uint32_t result = (((uint8_t)roundf((0xFF*vector.x))) << 16) | 
                      (((uint8_t)roundf((0xFF*vector.y)))) << 8 | 
                      (uint8_t)roundf((0xFF*vector.z));
    return result;
}

static void renderer_draw_scanline(OffscreenBuffer *buffer,
                                   OffscreenBuffer *zbuffer,
                                   OffscreenBuffer *normal_buffer,
                                   V3 start,
                                   V3 end,
                                   V3 start_normal,
                                   V3 end_normal,
                                   V3 color,
                                   uint8_t light_volume_byte) {

    uint8_t bright = 0x0;

    if (color.x > 1 && color.y > 1 && color.z > 1) {
        bright = 0x01;
        color -= v3(1,1,1);
    }

    if (start.x > end.x) {
        v3_swap(&start, &end);
        v3_swap(&start_normal, &end_normal);
    }

    int x = floorf(start.x);
    int y = floorf(start.y);
    float z = start.z;

    int dx = floorf(end.x) - x ;
    float z_inc = (start.z / dx);
    V3 normal_inc = (end_normal - start_normal) / dx;
    uint32_t pixel_color = convert_v3_to_RGB(color);

    for (int i = 0; i <= abs(dx); ++i) {
        int offset = x*buffer->bytes_per_pixel + y*buffer->pitch;
        float *depth_value = (float *)((uint8_t *)zbuffer->memory + offset);
        uint32_t *pixel = (uint32_t *)((uint8_t *)buffer->memory + offset);
#if SHADING
        uint32_t *normal = (uint32_t *)((uint8_t *)normal_buffer->memory + offset);
#endif

        if ((z < *depth_value)) {
            pixel_color = ((bright) << 31) | pixel_color;

#if SHADING
            V3 normalised = v3_norm(start_normal);
            normalised.x = roundf((start_normal.x*0.5f+0.5f)*65535.0f);
            normalised.y = roundf((start_normal.y*0.5f+0.5f)*65535.0f);
            uint32_t normal_value = (((uint16_t)(normalised.x)) << 16) |
                                    (((uint16_t)(normalised.y)));
#endif

            *pixel = pixel_color;
            *depth_value = start.z;
#if SHADING
            *normal = normal_value;
#endif

        }
        x++;
        z += z_inc;
#if SHADING
        start_normal += normal_inc;
#endif
    }
}

static void renderer_draw_line(OffscreenBuffer *buffer,
                               OffscreenBuffer *zbuffer,
                               OffscreenBuffer *normal_buffer,
                               V3 start,
                               V3 end,
                               V3 start_normal,
                               V3 end_normal,
                               V3 color,
                               uint8_t light_volume_byte) {


#if 0
    clamp_x(&start, buffer->width);
    clamp_x(&end, buffer->width);
#endif

#if 0
    if (start.x > end.x) {
        v3_swap(&start, &end);
        v3_swap(&start_normal, &end_normal);
    }
#endif

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
        int offset = (int)((start.x))*buffer->bytes_per_pixel + 
            (int)((start.y))*buffer->pitch;
        float *depth_value = (float *)((uint8_t *)zbuffer->memory + offset);
        uint32_t *pixel = (uint32_t *)((uint8_t *)buffer->memory + offset);
#if SHADING
        uint32_t *normal = (uint32_t *)((uint8_t *)normal_buffer->memory + offset);
#endif

#if 0
        uint8_t light_volume_index = (light_volume_byte & (0b01111111));

        if (light_volume_index) {
            *pixel = (light_volume_index << 24) | *pixel;
            *pixel = 0xFFFFFFFF;
       }
#endif

        if ((start.z < *depth_value)) {
            uint32_t pixel_color = convert_v3_to_RGB(color);
            pixel_color = ((bright) << 31) | pixel_color;

            // normalise between 0 and 1
            //
#if SHADING
            V3 normalised = v3_norm(start_normal);
            normalised.x = roundf((start_normal.x*0.5f+0.5f)*65535.0f);
            normalised.y = roundf((start_normal.y*0.5f+0.5f)*65535.0f);
            uint32_t normal_value = (((uint16_t)(normalised.x)) << 16) | 
                                    (((uint16_t)(normalised.y)));
#endif



            *pixel = pixel_color;
#if SHADING
            *normal = normal_value;
#endif
            *depth_value = start.z;
        }

        start += start_inc;
#if SHADING
        start_normal += normal_inc;
#endif
    }
}

static bool early_depth_test(OffscreenBuffer *zbuffer,
                             V3 start,
                             V3 end) {

    float start_z = *(float *)(zbuffer->memory) + (int)roundf(start.x) + (int)roundf(start.y)*zbuffer->width;
    float end_z = *(float *)(zbuffer->memory) + (int)roundf(end.x) + (int)roundf(end.y)*zbuffer->width;

    if (start.z > start_z && end.z > end_z) {
        return false;
    }

    return true;
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
                                            V3 color,
                                            uint8_t light_volume_byte) {

    if (!early_depth_test) {
        return;
    }
#if 0
    clamp_y(&left, buffer->height);
    clamp_y(&right, buffer->height);
    clamp_y(&bottom, buffer->height);
#endif

    float dy = floorf(bottom.y)-floorf(right.y);

    V3 right_inc = v3((bottom.x-right.x)/dy, 1, -(bottom.z-right.z)/dy);
    V3 left_inc = v3((bottom.x-left.x)/dy, 1, -(bottom.z-left.z)/dy);
    V3 left_normal_inc = (bottom_normal - left_normal)/dy;
    V3 right_normal_inc = (bottom_normal - right_normal)/dy;


    for (int y = 0; y <= fabsf(dy); ++y) {
        renderer_draw_scanline(buffer, zbuffer, normal_buffer, left, right,
                               left_normal, right_normal, color, light_volume_byte);

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
                                               V3 color,
                                               uint8_t light_volume_byte) {

    if (!early_depth_test) {
        return;
    }
#if 0
    clamp_y(&left, buffer->height);
    clamp_y(&right, buffer->height);
    clamp_y(&top, buffer->height);
#endif
    float dy = floorf(top.y)-floorf(right.y);

    V3 left_inc = v3((top.x-left.x)/dy, 1, (top.z-left.z)/dy);
    V3 right_inc = v3((top.x-right.x)/dy, 1, (top.z-right.z)/dy);
    V3 right_normal_inc = (top_normal - right_normal)/dy;
    V3 left_normal_inc = (top_normal - left_normal)/dy;

    for (int y = 0; y <= fabsf(dy); ++y) {
        renderer_draw_scanline(buffer, zbuffer, normal_buffer, left, right,
                               left_normal, right_normal, color, light_volume_byte);

        left -= left_inc;
        right -= right_inc;
        left_normal -= left_normal_inc;
        right_normal -= right_normal_inc;
    }
}

static void draw_triangle(OffscreenBuffer *buffer, 
                          OffscreenBuffer *zbuffer, 
                          OffscreenBuffer *normal_buffer, 
                          V3 vert1,
                          V3 vert2,
                          V3 vert3,
                          V3 normal1,
                          V3 normal2,
                          V3 normal3,
                          V3 color,
                          uint8_t light_volume_byte) {

    if (renderer_v3_is_invalid(vert1, buffer->width, buffer->height) &&
        renderer_v3_is_invalid(vert2, buffer->width, buffer->height) &&
        renderer_v3_is_invalid(vert3, buffer->width, buffer->height)) {

        return;
    }

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
                    normal2, normal3, normal1, color, light_volume_byte);
        }
        else {
            renderer_draw_flat_top_triangle(buffer, zbuffer, normal_buffer, vert3, vert2, vert1, 
                    normal3, normal2, normal1, color, light_volume_byte);
        }

        return;
    }

    if (vert1.y == vert2.y) {
        if (vert1.x > vert2.x) {
            renderer_draw_flat_bottom_triangle(buffer, zbuffer, normal_buffer, vert2, vert1, vert3, 
                    normal2, normal1, normal3, color, light_volume_byte);
        }
        else {
            renderer_draw_flat_bottom_triangle(buffer, zbuffer, normal_buffer, vert1, vert2, vert3, 
                    normal1, normal2, normal3, color, light_volume_byte);
        }
        return;
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
                normal2, normal4, normal1, color, light_volume_byte);
        renderer_draw_flat_bottom_triangle(buffer, zbuffer, normal_buffer, vert2, vert4, vert3, 
                normal2, normal4, normal3, color, light_volume_byte);
    }

    else {
        renderer_draw_flat_top_triangle(buffer, zbuffer, normal_buffer, vert4, vert2, vert1, 
                normal4, normal2, normal1, color, light_volume_byte);

        renderer_draw_flat_bottom_triangle(buffer, zbuffer, normal_buffer, vert4, vert2, vert3, 
                normal4, normal2, normal3, color, light_volume_byte);
    }
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
    result.scale = scale;
    result.translation = translation;
    result.rotation = rotation;
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
    if (obj->mesh == nullptr) {
        return;
    }

    assert(render_state->vertex_count < arraysize(render_state->vertex_buffer));
    assert(obj->mesh->vert_count % 4 == 0);
    assert(obj->mesh->vertexn_count % 4 == 0);

	for (uint32_t vertex = 0; vertex < obj->mesh->vert_count; vertex += 4) {
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

static void renderer_render_obj_make_light_source(RendererState *render_state, RenderObj *obj, V3 position, V3 color, float attenuation) {
        render_state->light_sources[render_state->light_sources_count].position = position;
        render_state->light_sources[render_state->light_sources_count].color = color;
        render_state->light_sources[render_state->light_sources_count].attenuation = attenuation;
        render_state->light_sources[render_state->light_sources_count].obj = obj;
        ++render_state->light_sources_count;
}


static void renderer_draw(RendererState *render_state, 
                          GameCamera *camera, 
                          OffscreenBuffer *buffer, 
                          OffscreenBuffer *zbuffer,
                          OffscreenBuffer *normal_buffer,
                          OffscreenBuffer *postbuffer) {

    render_state->draw_count = 0;

    // clear buffers
    {
        uint32_t *color = (uint32_t *)buffer->memory;
        float *normal = (float *)normal_buffer->memory;
        float *depth = (float *)zbuffer->memory;
        for (int i = 0;i < buffer->height*buffer->width; ++i) {
            *color++ = 0;
            *depth++ = FLT_MAX;
            *normal++ = 0;
        }
    }

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
            *_out = vertex4_rotate(_v4, rotation);
        }

        for (uint32_t normal = 0; normal < render_state->normal_count; ++normal) {
            Vertex4 _n4 = render_state->normal_buffer[normal];
            Vertex4 *_out = &render_state->normal_out_buffer[normal];
#if WORLD_LIGHT
            *_out = _n4;
#else
            *_out = vertex4_rotate(_n4, rotation);
#endif

        }

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
            current.v1 = 0xFFFFFFFF;
            current.v2 = 0xFFFFFFFF;
            current.v3 = 0xFFFFFFFF;
        }
#endif
            render_state->polygons_to_draw[triangle] = current;
    }

    // change light sources to be 100 times the brightness for post processing
    for (uint32_t i = 0; i < render_state->light_sources_count; ++i) {
        LightSource light_source = render_state->light_sources[i]; 
        render_state->light_colors[i] = light_source.color;
#if WORLD_LIGHT
        render_state->light_positions[i] = light_source.position;
#else 
        render_state->light_positions[i] = renderer_world_vertex_to_view(light_source.position, camera);
#endif

        for (uint32_t j = light_source.obj->index_start; j < light_source.obj->index_end; ++j) {
            Triangle *polygon = &render_state->polygons_to_draw[j];
            polygon->color = v3(100,100,100)+polygon->color;
        }
    }

#if 0
    // sample the radius of the light sources at different x values to create a sphere
    // to then render to screen
    {
        uint32_t samples = 12;
        int circle_vertex_count = (samples*2+1);
        for (uint32_t light = 0; light < render_state->light_sources_count; ++light) {
            LightSource source = render_state->light_sources[light];
            source.position = renderer_world_vertex_to_view(source.position, camera);
            float circle_z = source.position.z;
            float min_intesity = 0.05f;
            float radius = 1/(min_intesity*source.attenuation);
            render_state->light_vertices[light*circle_vertex_count] = source.position;
            render_state->light_vertices[light*circle_vertex_count+1] = 
                v3(source.position.x-radius, source.position.y, circle_z);

            render_state->light_vertices[light*circle_vertex_count+samples+1] = 
                v3(source.position.x+radius, source.position.y, circle_z);

            for (uint32_t x = 1; x < samples; ++x) {
                float circle_x = f_lerp((float)samples, (float)(x), 
                        source.position.x-radius, 
                        source.position.x+radius);
                float circle_y =  
                    sqrtf((radius*radius) - (circle_x-source.position.x)*
                            (circle_x-source.position.x));
                render_state->light_vertices[light*circle_vertex_count+x+1] = v3(circle_x, source.position.y+circle_y, circle_z);
                render_state->light_vertices[(light+1)*(circle_vertex_count)-(x)] = 
                    v3(circle_x, source.position.y-circle_y, circle_z);
            }
        }

        // transform light volume vertices to screen
        for (uint32_t i = 0; i < render_state->light_sources_count*(samples*2+1); ++i) {
            V3 *position = &render_state->light_vertices[i];
            V3 result = {};
            float aspect_ratio = (float)buffer->width / (float)buffer->height;
            float length_x = tanf(camera->fov/2.0f)*camera->znear*aspect_ratio;
            float length_y = tanf(camera->fov/2.0f)*camera->znear;

            float normal_x = (position->x*(camera->znear/position->z))/length_x;
            float normal_y = (position->y*(camera->znear/position->z))/length_y;

            if (position->z > camera->zfar) {
                *position = renderer_vertex_invalid();
                continue;
            }

            if (position->z < camera->znear) {
                *position = renderer_vertex_invalid();
                continue;
        render_state->light_positions[i] = light_source.position;
            }

            position->x = (normal_x * (buffer->width/2.0f)) + buffer->width/2.0f;
            position->y = (-normal_y * (buffer->height/2.0f)) + buffer->height/2.0f;
        }
    }
#endif 

    // perspective projection and clipping
    {
        float aspect_ratio = (float)buffer->width / (float)buffer->height;
        for (uint32_t i = 0; i < render_state->vertex_count; ++i) {

            __m128 _vertex_x = render_state->vertex_out_buffer[i].x;
            __m128 _vertex_y = render_state->vertex_out_buffer[i].y;
            __m128 _vertex_z = render_state->vertex_out_buffer[i].z;

            __m128 _cam_znear = _mm_set1_ps(camera->znear);
            __m128 _cam_zfar = _mm_set1_ps(camera->zfar);

            __m128 _length_x = _mm_set1_ps(tanf(camera->fov/2.0f)*camera->znear*aspect_ratio);
            __m128 _length_y = _mm_set1_ps(tanf(camera->fov/2.0f)*camera->znear);

            __m128 _mask = _mm_cmplt_ps(_vertex_z, _cam_znear);
            _mask = _mm_or_ps(_mask, _mm_cmpgt_ps(_vertex_z, _cam_zfar));

            __m128 _normal_x = _mm_div_ps(_mm_mul_ps(_vertex_x, 
                        _mm_div_ps(_cam_znear, _vertex_z)),
                    _length_x);

            __m128 _normal_y = _mm_div_ps(_mm_mul_ps(_vertex_y, 
                        _mm_div_ps(_cam_znear, _vertex_z)),
                    _length_y);

#if 1
           _mask = _mm_or_ps(_mask, _mm_cmpgt_ps(_normal_y, _mm_set1_ps(1.0f)));
           _mask = _mm_or_ps(_mask, _mm_cmplt_ps(_normal_y, _mm_set1_ps(-1.0f)));
           _mask = _mm_or_ps(_mask, _mm_cmpgt_ps(_normal_x, _mm_set1_ps(1.0f)));
           _mask = _mm_or_ps(_mask, _mm_cmplt_ps(_normal_x, _mm_set1_ps(-1.0f)));
#endif


            __m128 _buffer_x = _mm_set1_ps(buffer->width/2.0f);
            __m128 _buffer_y = _mm_set1_ps(buffer->height/2.0f);

            __m128 _float_res_x = _mm_add_ps(_mm_mul_ps(_normal_x, _buffer_x), 
                    _buffer_x);
            __m128 _float_res_y = _mm_add_ps(_mm_mul_ps(_mm_mul_ps(_normal_y, 
                            _mm_set1_ps(-1.0f)), 
                        _buffer_y), 
                    _buffer_y);

            render_state->vertex_out_buffer[i].x = _mm_or_ps(_float_res_x, _mask);
            render_state->vertex_out_buffer[i].y = _mm_or_ps(_float_res_y, _mask);
            render_state->vertex_out_buffer[i].z = _mm_or_ps(_vertex_z, _mask);
        }
    }

    // Draw Triangles
    V3 white_point = v3(100,100,100);
    {
        Triangle *triangles = render_state->polygons_to_draw;
        Vertex4 *vertices = render_state->vertex_out_buffer;
        Vertex4 *normals = render_state->normal_out_buffer;

        for (uint32_t i = 0; i < render_state->polygon_count; ++i) {
            if (triangles[i].v1 == 0xFFFFFFFF) {
                continue;
            }

            V3 vert1 = triangle_get_attribute_from_buffer(triangles[i].v1, vertices, render_state->vertex_count);
            V3 vert2 = triangle_get_attribute_from_buffer(triangles[i].v2, vertices, render_state->vertex_count);
            V3 vert3 = triangle_get_attribute_from_buffer(triangles[i].v3, vertices, render_state->vertex_count);

            V3 normal1 = triangle_get_attribute_from_buffer(triangles[i].vn1, normals, render_state->normal_count);
            V3 normal2 = triangle_get_attribute_from_buffer(triangles[i].vn2, normals, render_state->normal_count);
            V3 normal3 = triangle_get_attribute_from_buffer(triangles[i].vn3, normals, render_state->normal_count);

            V3 color = triangles[i].color;

            // color correction & tone mapping
            if (color.x > 100 && color.y > 100 && color.z > 100) {
                color -= v3(100,100,100);
                color = 100.0f*color;
                float luminance = rgb_to_luminance(color);
                float white_l = rgb_to_luminance(white_point);
                float numerator = luminance * (1.0f + (luminance/(white_l * white_l)));
                float new_luminance = numerator / (1.0f + luminance);
                color = (new_luminance/luminance)*color;
                (color.x > 1) ? color.x = 1 : color.x;
                (color.y > 1) ? color.y = 1 : color.y;
                (color.z > 1) ? color.z = 1 : color.z;
                color += v3(1,1,1);
        }
            draw_triangle(buffer, zbuffer, normal_buffer, 
                          vert1, vert2, vert3, 
                          normal1, normal2, normal3, 
                          color, 0);
        }

    }
    
#if 0
    // draw light volumes
    {
        uint32_t samples = 12;
        uint32_t circle_index_count = (samples*2+1);
        for (uint32_t i = 0; i < render_state->light_sources_count; ++i) {
            V3 center = render_state->light_vertices[(i*circle_index_count)];
            V3 previous_vertex = render_state->light_vertices[(i*circle_index_count)+1];
#if 1
            for (uint32_t j = 2; j < circle_index_count; ++j) {
                V3 current_vertex = render_state->light_vertices[(i*circle_index_count)+j];
                draw_triangle(buffer, 
                              zbuffer, 
                              normal_buffer, 
                              center, previous_vertex, current_vertex, 
                              v3_zero(), v3_zero(), v3_zero(), v3(1,1,1), 
                              (uint8_t)i+1);
                previous_vertex = current_vertex;
            }
#endif

#if 1
            draw_triangle(buffer,
                         zbuffer,
                         normal_buffer,
                         center, 
                         render_state->light_vertices[(i*circle_index_count)+1],
                         render_state->light_vertices[((i+1)*circle_index_count)-1],
                         v3_zero(), v3_zero(), v3_zero(), v3(1,1,1),
                         (uint8_t)i+1);
#endif

        }
    }

#endif
    // Lighting Pass
#if SHADING
    {
        uint32_t *colors = (uint32_t *)buffer->memory;
        uint32_t *normals = (uint32_t *)normal_buffer->memory;
        float *depth = (float *)zbuffer->memory;
        V3 *light_positions = render_state->light_positions;
        V3 *light_colors = render_state->light_colors;

        for (int pixel = 0; pixel < buffer->width*buffer->height; pixel+=4) {

            bool should_skip = true;
            __m128i _depth_mask = {};
            for (int i = 0; i < 4; ++i) {
                if (!(depth[pixel+i] == FLT_MAX)) {
                    ((uint32_t *)&_depth_mask)[i] = 0xFFFFFFFF;
                    should_skip = false;
                }

                if (colors[pixel+i] & (1 << 31)) {
                    ((uint32_t *)&_depth_mask)[i] = 0;
                    should_skip = false;
                }

            }

            if (should_skip) {
                continue;
            }

            Vertex4 _extracted_norm = {};
            for (int i = 0; i < 4; ++i) {
                ((float *)&_extracted_norm.x)[i] = ((((normals[pixel+i] & 0xFFFF0000) >> 16)/65535.0f)-0.5f)*2;
                ((float *)&_extracted_norm.y)[i] = (((normals[pixel+i] & 0x0000FFFF)/65535.0f)-0.5f)*2;

            }

            __m128 _sum_sqrs = _mm_sub_ps(_mm_set1_ps(1.0f), _mm_add_ps(_mm_mul_ps(_extracted_norm.x, _extracted_norm.x),   
                        _mm_mul_ps(_extracted_norm.y, _extracted_norm.y)));

            _extracted_norm.z = _mm_mul_ps(_mm_set1_ps(-1.0f), (_mm_sqrt_ps(_sum_sqrs)));

            Vertex4 _extracted_pos = {};
            for (int i = 0; i < 4; ++i) {
                ((float *)&_extracted_pos.x)[i] = (float)((pixel+i) % 1280) + 0.5f; 
                ((float *)&_extracted_pos.y)[i] = (float)((pixel+i) / 1280) + 0.5f;
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
                _extracted_pos.y = _mm_div_ps(_mm_sub_ps(_extracted_pos.y, _buffer_height), (_buffer_height));
                _extracted_pos.y = _mm_mul_ps(_extracted_pos.y, _mm_set1_ps(-1.0f)); 

                _extracted_pos.x = _mm_div_ps(_mm_mul_ps(_extracted_pos.x, _length_x), 
                        _mm_div_ps(_znear, _extracted_pos.z));

                _extracted_pos.y = _mm_div_ps(_mm_mul_ps(_extracted_pos.y, _length_y), 
                        _mm_div_ps(_znear, _extracted_pos.z));

#if WORLD_LIGHT
                Quaternion inverse_view = q4_conj(rotation_q4(-camera->theta_x, v3(1,0,0))*
                                                  rotation_q4(-camera->theta_y, v3(0,1,0)));
                Vertex4 _camera = {};
                _camera.x = _mm_set1_ps(camera->pos.x);
                _camera.y = _mm_set1_ps(camera->pos.y);
                _camera.z = _mm_set1_ps(camera->pos.z);

                _extracted_pos = vertex4_rotate(_extracted_pos, inverse_view);
                _extracted_pos = vertex4_add(_extracted_pos, _camera);
#endif
            }



            Vertex4 _extracted_color {};
            for (int i = 0; i < 4; ++i) {
                ((float *)&_extracted_color.x)[i] = ((colors[pixel+i] & 0x00FF0000) >> 16)/255.0f;
                ((float *)&_extracted_color.y)[i] = ((colors[pixel+i] & 0x0000FF00) >> 8)/255.0f;
                ((float *)&_extracted_color.z)[i] = (colors[pixel+i] & 0x000000FF)/255.0f;
            }

            Vertex4 _color = {};
            _color.x = _mm_set1_ps(0.00f);
            _color.y = _mm_set1_ps(0.00f);
            _color.z = _mm_set1_ps(0.00f);

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
                __m128 _distance = _mm_sqrt_ps(vertex4_dot(_dl, _dl));

                __m128 _intensity = _mm_div_ps(vertex4_dot(vertex4_norm(_dl), vertex4_norm(_extracted_norm)),
                        _mm_mul_ps(_mm_set1_ps(render_state->light_sources[light].attenuation), 
                            _mm_sqrt_ps(vertex4_dot(_dl, _dl))));
                __m128 _above_zero_mask = _mm_cmpgt_ps(_intensity, _mm_setzero_ps());
                _intensity = _mm_and_ps(_intensity, _above_zero_mask);
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

            // unpack into the backbuffer
            {
                __m128 _one = _mm_set1_ps(1.0f);
                __m128 _zero = _mm_setzero_ps();

                //clamp between 0-1
                _color.x = _mm_min_ps(_mm_max_ps(_color.x, _zero), _one);
                _color.y = _mm_min_ps(_mm_max_ps(_color.y, _zero), _one);
                _color.z = _mm_min_ps(_mm_max_ps(_color.z, _zero), _one);

                __m128 _255 = _mm_set1_ps(255.0f);
                _color.x = _mm_mul_ps(_255, _color.x);
                _color.y = _mm_mul_ps(_255, _color.y);
                _color.z = _mm_mul_ps(_255, _color.z);

                __m128i _int_r = _mm_cvtps_epi32(_color.x);
                __m128i _int_g = _mm_cvtps_epi32(_color.y);
                __m128i _int_b = _mm_cvtps_epi32(_color.z);

                __m128i _out = _mm_or_si128(_mm_slli_epi32(_int_r, 16),
                        _mm_or_si128(_mm_slli_epi32(_int_g, 8), 
                            _int_b));

                __m128i _final_mask = _mm_and_si128(*(__m128i *)&colors[pixel], 
                                                    _mm_andnot_si128(_depth_mask, _mm_set1_epi32(0xFFFFFFFF)));
                *(__m128i *)&colors[pixel] = _mm_or_si128(_mm_and_si128(_out, _depth_mask), _final_mask);
            }
        }
    }
#endif

#if POST_FX
    // POST PROCESSING BLOOM EFFECT
    
    // zero the buffer
      
    {
        uint32_t *pixels = (uint32_t *)postbuffer->memory;
        for (int i = 0; i < postbuffer->width*postbuffer->height; ++i) {
            pixels[i] = 0; 
        }
    }

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
                        sumbright += ((in_pixels[(y*ratio+i)*buffer->width+x*ratio+j] & (1 << 31)) >> 24);
                    }
                }

                if (sumbright < (1 << 8)) {
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
