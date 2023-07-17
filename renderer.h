#include "particle_math.h"
#include <emmintrin.h>
#include <stdint.h>

#define BLACK 0
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

// indexes into a vertex array;
struct Triangle {
    int v1;
    int v2;
    int v3;

    int vn1;
    int vn2;
    int vn3;

    V3 color;
};

struct Mesh {
    V3 vertices[65536*2];
    uint32_t vert_count;
    
    V3 vertexn[65536*2];
    uint32_t vertexn_count;

    Triangle polygons[65536];
    uint32_t poly_count;
    float min_y;
    float max_y;
};

struct RenderObj {
    uint32_t vstart;
    uint32_t vend;
    uint32_t nstart;
    uint32_t nend;
    uint32_t index_start;
    uint32_t index_end;
    Mesh *mesh;
    V3 scale;
    V3 translation;
    Quaternion rotation;
    V3 color;
};

struct LightSource {
    RenderObj *obj;
    V3 position;
    V3 color;
    float attenuation;
};

struct RendererState {
    uint32_t vertex_count;
    Vertex4 vertex_buffer[65536*25];
    Vertex4 vertex_out_buffer[65536*25];

    uint32_t normal_count;
    Vertex4 normal_buffer[65536*25];
    Vertex4 normal_out_buffer[65536*25];

    uint32_t polygon_count;
    Triangle polygons[65536*50];

    uint32_t draw_count;
    Triangle polygons_to_draw[65536*50];

    uint32_t light_sources_count;
    LightSource light_sources[128];
    V3 light_positions[128];
    V3 light_colors[128];
    V3 light_vertices[128*(12*2+1)];
};
