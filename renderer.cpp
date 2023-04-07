#include "particle_math.h"
#include <emmintrin.h>
#include <stdint.h>
#include <xmmintrin.h>

struct OffscreenBuffer {
    void *memory;
    int height;
    int width;
    int pitch;
    int bytes_per_pixel;
};

struct GameCamera {
    V3 pos;
    float width;
    float height;
    
    float znear;
    float zfar;
    float fov;
    float theta_x;
    float theta_y;
};

struct V2Screen {
    int x;
    int y;
};

struct V2Screen4 {
    __m128i x;
    __m128i y;
};

struct Vertex4 {
    __m128 x;
    __m128 y;
    __m128 z;
};

struct Vertex4Cube {
    Vertex4 vertices[2];
};

static void Vertex4_to_V2Screen(Vertex4 *in, 
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

static V3 transform_world_to_view(V3 world_pos, GameCamera *camera, 
                                  int buffer_width,
                                  int buffer_height) {
    V3 result;
    result = v3_rotate_q4(world_pos-camera->pos, 
                          (rotation_q4(-camera->theta_x, v3(1,0,0))*
                          rotation_q4(-camera->theta_y, v3(0,1,0))));

    return result;
}

static V2Screen transform_v3_to_screen(V3 world_pos, GameCamera *camera,
                                        int buffer_width, int buffer_height) {

    V2Screen result;

    V3 relative_pos = v3_rotate_q4(world_pos-camera->pos, 
            (rotation_q4(-camera->theta_x, v3(1,0,0))*
             rotation_q4(-camera->theta_y, v3(0,1,0))));

    float aspect_ratio = ((float)buffer_width/(float)buffer_height);
    float max_x = tanf(camera->fov/2.0f)*camera->znear * (relative_pos.z/camera->znear) * 
        aspect_ratio;

    float max_y = tanf(camera->fov/2.0f)*camera->znear * (relative_pos.z/camera->znear);

    if (relative_pos.z < camera->znear) {
        result.x = -1; 
        result.y = -1;
        return result;
    }

    if (relative_pos.z > camera->zfar) {
        result.x = -1; 
        result.y = -1;
        return result;
    }

    if (fabs(relative_pos.x) > max_x) {
        result.x = -1;
        result.y = -1;
        return result;
    }

    if (fabs(relative_pos.y) > max_y) {
        result.x = -1;
        result.y = -1;
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


static void transform_vertexes(V3 *in, GameCamera *camera, int count, 
                               int buffer_width, int buffer_height,
                               V2Screen *out) {

    for (int i = 0; i < count; ++i) {
        V2Screen *screen_pos = &out[i];

        V3 relative_pos = v3_rotate_q4(in[i]-camera->pos, 
                (rotation_q4(-camera->theta_x, v3(1,0,0))*
                 rotation_q4(-camera->theta_y, v3(0,1,0))));

    float aspect_ratio = ((float)buffer_width/(float)buffer_height);
    float max_x = tanf(camera->fov/2.0f)*camera->znear * (relative_pos.z/camera->znear) * 
        aspect_ratio;

    float max_y = tanf(camera->fov/2.0f)*camera->znear * (relative_pos.z/camera->znear);

    if (relative_pos.z < camera->znear) {
        out->x = -1; 
        out->y = -1;
    }

    if (relative_pos.z > camera->zfar) {
        out->x = -1; 
        out->y = -1;
    }

    if (fabs(relative_pos.x) > max_x) {
        out->x = -1;
        out->y = -1;
    }

    if (fabs(relative_pos.y) > max_y) {
        out->x = -1;
        out->y = -1;
    }


    float normal_x = (relative_pos.x*(camera->znear/relative_pos.z)) / 
                     ((tanf(camera->fov/2.0f)*camera->znear)*aspect_ratio);

    float normal_y = (relative_pos.y*(camera->znear/relative_pos.z)) / 
                     (tanf(camera->fov/2.0f)*camera->znear);
    
    screen_pos->x = (int)(normal_x * (buffer_width/2.0f) + buffer_width/2.0f);
    screen_pos->y = (int)(-normal_y * (buffer_height/2.0f) + buffer_height/2.0f);

    }
}

static inline V3 transform_mouse_to_world(V2Screen screen,
                                  GameCamera *camera,
                                  int buffer_width,
                                  int buffer_height) {

    V3 result = {};
    result.x = (((float)screen.x / (float)(buffer_width/2.0)) - 1) * 
               (float)camera->width/2.0f;

    result.y = (1 - ((float)screen.y / (float)(buffer_height/2.0))) * 
                (float)camera->height/2.0f;

    return result + camera->pos;
}

static void draw_background(OffscreenBuffer *buffer, GameCamera* camera) {
    uint8_t *row = (uint8_t *)buffer->memory;
    for (int y = 0;y < buffer->height; ++y) {
        uint32_t *pixel = (uint32_t *)row;

        for (int x =0;x < buffer->width; ++x) {
            *pixel++ = 0xFFFF00FF;
        }

        row += buffer->pitch;
    }
}

static void draw_line(OffscreenBuffer *buffer,
                 V2Screen start,
                 V2Screen end,
                 uint32_t color) {

    if (start.x == 0xFFFFFFFF) return;
    if (start.y == 0xFFFFFFFF) return;

    if (end.x == 0xFFFFFFFF) return;
    if (end.y == 0xFFFFFFFF) return;

    if (start.x < 0 || start.x > 1280) return;
    if (start.y < 0 || start.y > 720) return;

    if (end.x < 0 || end.x > 1280) return;
    if (end.y < 0 || end.y > 720) return;

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

static void transform_and_draw_line(OffscreenBuffer *buffer,
                 V3 v_start,
                 V3 v_end,
                 GameCamera *camera,
                 uint32_t color) {

    V2Screen start = transform_v3_to_screen(v_start, camera, buffer->width, buffer->height);
    V2Screen end = transform_v3_to_screen(v_end, camera, buffer->width, buffer->height);

    if (start.x < 0 || start.x > 1280) return;
    if (start.y < 0 || start.y > 720) return;

    if (end.x < 0 || end.x > 1280) return;
    if (end.y < 0 || end.y > 720) return;

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

static void render_square(OffscreenBuffer *buffer, 
                   V2Screen pixel_pos, 
                   int side_length, 
                   uint32_t color) {

    if (pixel_pos.x > buffer->width) {
        return;
    }
    if (pixel_pos.y > buffer->height) {
        return;
    }

    if (pixel_pos.x < 0) {
        return;
    }

    if (pixel_pos.y < 0) {
        return;
    }

    int start_x = pixel_pos.x - side_length/2;
    int start_y = pixel_pos.y - side_length/2;

    if (start_x < 0) {
        start_x = 0;
    }

    if (start_y < 0) {
        start_y = 0;
    }

    uint8_t *row = (uint8_t *)buffer->memory + 
        (start_x*buffer->bytes_per_pixel) + (start_y*buffer->pitch);

    for (int y = start_y; y < start_y+side_length; ++y) {
        if (y >= buffer->height) {
            break;
        }
        uint32_t *pixel = (uint32_t *)row;

        for (int x = start_x;x < start_x+side_length; ++x) {
                if (x >= buffer->width) {
                    break;
                }
                *pixel++ = color; 
            }

        row += buffer->pitch;
    }
}

static void draw_wire_cube(OffscreenBuffer *buffer, V3 pos, GameCamera *camera, uint32_t color) {
    transform_and_draw_line(buffer, v3(pos.x, pos.y, pos.z), v3(pos.x, pos.y+1, pos.z), camera, color);
    transform_and_draw_line(buffer, v3(pos.x, pos.y+1, pos.z), v3(pos.x, pos.y+1, pos.z+1), camera, color);
    transform_and_draw_line(buffer, v3(pos.x, pos.y+1, pos.z+1), v3(pos.x, pos.y, pos.z+1), camera, color);
    transform_and_draw_line(buffer, v3(pos.x, pos.y+1, pos.z+1), v3(pos.x, pos.y, pos.z+1), camera, color);
    transform_and_draw_line(buffer, v3(pos.x, pos.y, pos.z+1), v3(pos.x+1, pos.y, pos.z+1), camera, color);
    transform_and_draw_line(buffer, v3(pos.x+1, pos.y, pos.z+1), v3(pos.x+1, pos.y, pos.z), camera, color);
    transform_and_draw_line(buffer, v3(pos.x+1, pos.y, pos.z), v3(pos.x+1, pos.y+1, pos.z), camera, color);
    transform_and_draw_line(buffer, v3(pos.x+1, pos.y+1, pos.z), v3(pos.x+1, pos.y+1, pos.z+1), camera, color);
    transform_and_draw_line(buffer, v3(pos.x+1, pos.y+1, pos.z+1), v3(pos.x, pos.y+1, pos.z+1), camera, color);
    transform_and_draw_line(buffer, v3(pos.x, pos.y, pos.z), v3(pos.x+1, pos.y, pos.z), camera, color);
    transform_and_draw_line(buffer, v3(pos.x, pos.y, pos.z), v3(pos.x, pos.y, pos.z+1), camera, color);
    transform_and_draw_line(buffer, v3(pos.x+1, pos.y+1, pos.z+1), v3(pos.x+1, pos.y, pos.z+1), camera, color);
    transform_and_draw_line(buffer, v3(pos.x, pos.y+1, pos.z), v3(pos.x+1, pos.y+1, pos.z), camera, color);
}
