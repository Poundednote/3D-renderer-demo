#include <stdint.h>
#include <windows.h>
#include <windowsx.h>
#include <winuser.h>
#include <stdio.h>

#include "particle_fun.h"
#include "particle_fun.cpp"


struct VideoOffscreenBuffer {
    BITMAPINFO info = {};
    void *bitmap_memory; 
    int height;
    int width;  
    int pitch;
    int bytes_per_pixel;
};

static bool RUNNING = false;
static VideoOffscreenBuffer vid_buf;

void win32_stretch_buffer_to_window(HWND window, HDC device_context,
                                    VideoOffscreenBuffer buffer) {


    RECT rect;
    GetClientRect(window, &rect);
    int width = rect.right - rect.left;
    int height = rect.bottom - rect.top;

    StretchDIBits(device_context, 
                  0, 0, width, height, 
                  0, 0, buffer.width, buffer.height,
                  buffer.bitmap_memory, 
                  &buffer.info, DIB_RGB_COLORS, SRCCOPY);
}


inline uint32_t safe_truncate_uint64(uint64_t value) {
    assert(value <= 0xFFFFFFFF);
    uint32_t result = (uint32_t)value;
    return result;
}

static LRESULT CALLBACK win32_window_proc(HWND window, UINT message, 
                                 WPARAM w_param, LPARAM l_param) {

    LRESULT result = 0;
    switch (message) {
        case WM_CLOSE:
            RUNNING = false;
            OutputDebugStringA("WM_CLOSE\n");
            break;

        case WM_PAINT: {
            OutputDebugStringA("WM_PAINT\n");
            PAINTSTRUCT paint;
            HDC dc_paint = BeginPaint(window, &paint);
            win32_stretch_buffer_to_window(window,  dc_paint, vid_buf);

            EndPaint(window, &paint);
        }

        default:
            result = DefWindowProcA(window, message, w_param, l_param);
            break;
    }
    
    return result;
}

int CALLBACK WinMain(HINSTANCE instance, HINSTANCE prev_instance, 
            LPSTR cmd_line, int show_cmd) {
    
    LARGE_INTEGER perf_counter_freq;
    QueryPerformanceFrequency(&perf_counter_freq);
    double time_elapsed_for_count = 1.0f/(double)perf_counter_freq.QuadPart;

    WNDCLASS window_class = {};
    window_class.style = CS_OWNDC|CS_HREDRAW|CS_VREDRAW; 
    window_class.lpfnWndProc = win32_window_proc;
    window_class.hInstance = instance;
    window_class.lpszClassName = "PhysicsEngine";

    if (!(RegisterClass(&window_class))) {
        OutputDebugStringA("CANT REGISTER WINDOW\n");
        return 0;
    }
    
    HWND window = CreateWindowEx(NULL, window_class.lpszClassName,
                                 "Physics Engine", 
                                 WS_OVERLAPPEDWINDOW|WS_VISIBLE,
                                 CW_USEDEFAULT,
                                 CW_USEDEFAULT,
                                 1280,
                                 720,
                                 NULL, NULL, instance, NULL);
    if (!window) {
        OutputDebugStringA("ERROR CREATING WINDOW\n");
        return 0;
    }

    HDC dc_window = GetDC(window);
    
    vid_buf.height = 720;
    vid_buf.width = 1280;
    vid_buf.bytes_per_pixel = 4;
    vid_buf.pitch = vid_buf.width*vid_buf.bytes_per_pixel;
    vid_buf.info.bmiHeader.biSize = sizeof(vid_buf.info.bmiHeader);
    vid_buf.info.bmiHeader.biWidth = vid_buf.width;
    vid_buf.info.bmiHeader.biHeight = -vid_buf.height;
    vid_buf.info.bmiHeader.biPlanes = 1;
    vid_buf.info.bmiHeader.biBitCount = 32;
    vid_buf.info.bmiHeader.biCompression = BI_RGB;
    int bitmap_memory_size = vid_buf.width * vid_buf.height * 
                             vid_buf.bytes_per_pixel;

    vid_buf.bitmap_memory = VirtualAlloc(0, bitmap_memory_size, 
                                         MEM_RESERVE|MEM_COMMIT, 
                                         PAGE_READWRITE);

    GameMemory game_memory = {};
    game_memory.permanent_stoarage_size = megabytes(64);
    game_memory.permanent_storage = VirtualAlloc(
                0, game_memory.permanent_stoarage_size, 
                MEM_RESERVE|MEM_COMMIT, PAGE_READWRITE);

    OffscreenBuffer game_buffer = {};
    game_buffer.memory = vid_buf.bitmap_memory;
    game_buffer.width = vid_buf.width;
    game_buffer.height = vid_buf.height;
    game_buffer.pitch = vid_buf.pitch;
    game_buffer.bytes_per_pixel = vid_buf.bytes_per_pixel;

    OffscreenBuffer zbuffer = {};
    zbuffer.width = game_buffer.width;
    zbuffer.height = game_buffer.height;
    zbuffer.pitch = game_buffer.pitch;
    zbuffer.bytes_per_pixel = game_buffer.bytes_per_pixel;
    
    uint32_t zbuffer_size = zbuffer.width*zbuffer.height*zbuffer.bytes_per_pixel;
    zbuffer.memory = VirtualAlloc(0, zbuffer_size, MEM_RESERVE|MEM_COMMIT, PAGE_READWRITE);
    
    LARGE_INTEGER previous_counter; 
    QueryPerformanceCounter(&previous_counter);
    GameInput input = {}; 

    RUNNING = true;
    while (RUNNING) {

        LARGE_INTEGER start_counter = previous_counter;

        MSG msg;

        HDESK desktop = OpenInputDesktop(0, 0, DELETE|DESKTOP_SWITCHDESKTOP);
        SetThreadDesktop(desktop);
        
        POINT cursor;
        GetCursorPos(&cursor);
        ScreenToClient(window, &cursor);
        input.mouse_pos.x = cursor.x;
        input.mouse_pos.y = cursor.y;

        input.mouse_lbutton_click = false;
        input.action = false;

        static float time = 0;
        static uint32_t count;

        while (PeekMessageA(&msg, NULL, 0, 0, PM_REMOVE)) {

            switch (msg.message) {
            
                case WM_QUIT: 
                    RUNNING = false;
                    break;

                case WM_KEYDOWN: {
                    WPARAM key = msg.wParam;  
                    bool was_down = (msg.lParam >> 30);
                    switch (key) {
                        case VK_RETURN: {
                            OutputDebugStringA("ENTTER\n");
                            input.action = true;
                            if (was_down) {
                                OutputDebugStringA("WASDOWN\n");
                            }
                            else {
                                OutputDebugStringA("first press\n");
                            }
                            break;
                        }

                        case VK_LEFT: {
                            input.camleft = true;
                            break;
                        }

                        case VK_RIGHT: {
                            input.camright = true;
                            break;
                        }

                        case VK_UP: {
                            input.camup = true;
                            break;
                        }

                        case VK_DOWN: {
                            input.camdown = true;
                            break;
                        }

                        case 'Q': {
                            input.camzoomout = true;
                            break;
                        }
                        case 'E': {
                            input.camzoomin = true;
                            break;
                        }

                        case 'A': {
                            input.left = true;
                            break;
                        }

                        case 'D': {
                            input.right = true;
                            break;
                        }

                        case 'W': {
                            input.up = true;
                            break;
                        }

                        case 'S': {
                            input.down = true;
                            break;
                        }
                    }
                    break;

                } break;

                case WM_KEYUP: {
                    WPARAM key = msg.wParam;  
                    bool was_down = (msg.lParam >> 30);
                    switch (key) {

                        case VK_LEFT: {
                            input.camleft = false;
                            OutputDebugStringA("leftup");
                            break;
                        }

                        case VK_RIGHT: {
                            input.camright = false;
                            OutputDebugStringA("rightUp");
                            break;
                        }

                        case VK_UP: {
                            input.camup = false;
                            break;
                        }

                        case VK_DOWN: {
                            input.camdown = false;
                            break;
                        }

                        case 'Q': {
                            input.camzoomout = false;
                            break;
                        }
                        case 'E': {
                            input.camzoomin = false;
                            break;
                        }

                        case 'A': {
                            input.left = false;
                            break;
                        }

                        case 'D': {
                            input.right = false;
                            break;
                        }

                        case 'W': {
                            input.up = false;
                            break;
                        }

                        case 'S': {
                            input.down = false;
                            break;
                        }

                    }
                } break;

                case WM_LBUTTONDOWN: {
                    input.mouse_lbutton_down = true;
                    input.start_click_pos.x = GET_X_LPARAM(msg.lParam);
                    input.start_click_pos.y = GET_Y_LPARAM(msg.lParam);
                    OutputDebugStringA("LBUTTONDOWN\n");
                    break;

                }

                case WM_LBUTTONUP: {
                    if (input.mouse_lclickdrag) {
                        input.mouse_lclickdrag = false;
                    }
                    else {
                        input.mouse_lbutton_down = false;
                        input.mouse_lbutton_click = true;
                    }
                    OutputDebugStringA("LBUTTONUP\n");
                    break;
                }

                case WM_MOUSEMOVE: {
                    // if leftbutton held down
                    if (msg.wParam & MK_LBUTTON) {
                        input.mouse_lclickdrag = true;
                    }

                   break;
                }

                default: 
                    TranslateMessage(&msg);
                    DispatchMessageA(&msg);
            }
        }

        game_update_and_render(&game_memory, 
                               &game_buffer, 
                               &zbuffer,
                               &input);

        if (!vid_buf.bitmap_memory) {
            OutputDebugStringA("lost bitmap_memory\n");
        }
        win32_stretch_buffer_to_window(window, dc_window, vid_buf);

        LARGE_INTEGER end_counter; 
        QueryPerformanceCounter(&end_counter);

        double time_elapsed_for_frame = ((double)end_counter.QuadPart - 
                                         (double)start_counter.QuadPart) *
                                        time_elapsed_for_count;
#if 1
        while (time_elapsed_for_frame < 1.0f/GAME_UPDATE_HZ) {
            time_elapsed_for_frame = (end_counter.QuadPart - 
                                      start_counter.QuadPart) *
                                     time_elapsed_for_count;

            QueryPerformanceCounter(&end_counter);
        }

#endif
        previous_counter = end_counter;

#if DEBUG_MODE
        char frame_debug_buffer[256];
        ++count;
        snprintf(frame_debug_buffer, 
                 sizeof(frame_debug_buffer), 
                 "time elapsed for frame: %f, game loops: %u\n",
                 time_elapsed_for_frame, count);

        OutputDebugStringA(frame_debug_buffer);
#endif
    }
    return 0;
}
