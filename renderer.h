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

    int vn1;
    int vn2;
    int vn3;
};

struct Mesh {
    V3 vertices[65536*2];
    int vert_count;
    
    V3 vertexn[65536*2];
    int vertexn_count;

    Triangle polygons[65536];
    int poly_count;
    float min_y;
    float max_y;
};

struct RenderObj {
    int index;
    Mesh *mesh;
    V3 color;
};

struct LightSource {
    RenderObj *obj;
    V3 position;
    V3 color;
    float falloff;
};

struct RendererState {
    int vertex_count;
    V3 vertex_list[65536*100];
    V3 vertex_colors[65536*100];

    int vertexn_count;
    V3 vertexn_list[65536*100];

    int polygon_count;
    Triangle polygons[65536*50];

    int draw_count;
    Triangle polygons_to_draw[65536*30];

    int light_sources_count;
    LightSource light_sources[255];
};
