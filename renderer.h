#include "particle_math.h"
#include <emmintrin.h>
#include <stdint.h>

#define BLACK 0;
#define WHITE 0xFFFFFFFF
#define RED 0xFFFF0000
#define GREEN 0xFF00FF00
#define BLUE 0xFF0000FF
#define YELLOW 0xFFFFFF00
#define PINK 0xFFFF00FF
#define CYAN 0xFF00FFFF

struct OffscreenBuffer {
    void *memory;
    int height;
    int width;
    int bytes_per_pixel;
    int pitch;
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

struct V3Screen {
    int x;
    int y;
    float z;
};

struct V3Screen4 {
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

// indexes into a vertex array;
struct Triangle {
    int v1;
    int v2;
    int v3;
    uint32_t color;
};

struct CubeMesh {
    Triangle triangles[12];
};

static void renderer_vertex4_to_v2screen(Vertex4 *in,
                                GameCamera *camera,
                                int screen_width,
                                int screen_height,
                                int count,
                                V3Screen4 *out);

static V3 renderer_world_vertex_to_view(V3 world_pos, 
                                                  GameCamera *camera, 
                                                  int buffer_width, 
                                                  int buffer_height);

static bool renderer_v3_should_clip(V3 pos, GameCamera *camera, float aspect_ratio);

static V3Screen renderer_world_vertex_to_screen(V3 world_pos, 
                                                          GameCamera *camera, 
                                                          int buffer_width, 
                                                          int buffer_height);

static void renderer_world_vertices_to_screen(V3 *in, 
                                                        GameCamera *camera, 
                                                        int count, 
                                                        int buffer_width, 
                                                        int buffer_height, 
                                                        V3Screen *out);

static bool renderer_check_v2screen_invalid(V3Screen screen, 
                                            int screen_width, 
                                            int screen_height);

static void renderer_draw_background(OffscreenBuffer *buffer);
static void renderer_draw_line(OffscreenBuffer *buffer, 
                               V3Screen start, 
                               V3Screen end, 
                               uint32_t color);

static void renderer_transform_and_draw_line(OffscreenBuffer *buffer, 
                                             V3 v_start, 
                                             V3 v_end, 
                                             GameCamera *camera, 
                                             uint32_t color);

static void renderer_draw_triangle_wireframe(OffscreenBuffer *buffer, 
                                             V3Screen v1,
                                             V3Screen v2,
                                             V3Screen v3,
                                             uint32_t color);

static void renderer_v2screen4_draw_triangle_wireframe(OffscreenBuffer *buffer, 
                                                       V3Screen4 *vertices, 
                                                       Triangle *triangles, 
                                                       uint32_t color, 
                                                       int count);

static void renderer_draw_flat_top_triangle(OffscreenBuffer *buffer, 
                                            V3Screen v1, 
                                            V3Screen v2, 
                                            V3Screen v3, 
                                            uint32_t color);

static void renderer_draw_flat_bottom_triangle(OffscreenBuffer *buffer, 
                                               V3Screen v1,
                                               V3Screen v2,
                                               V3Screen v3,
                                               uint32_t color);

static void renderer_draw_triangles_filled(OffscreenBuffer *buffer,
                                           V3Screen *vertices,
                                           Triangle *triangles,
                                           uint32_t *colors,
                                           int count);

static void renderer_v2screen4_draw_triangle_filled(OffscreenBuffer *buffer,
                                                    V3Screen4 *vertices,
                                                    Triangle *triangles,
                                                    uint32_t *colors,
                                                    int colors_size,
                                                    int tricount);
