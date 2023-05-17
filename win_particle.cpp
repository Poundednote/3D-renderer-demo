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

static void PlatformFreeFile(void *filebuffer) {
    if (filebuffer) {
        VirtualFree(filebuffer, 0, MEM_RELEASE);
        filebuffer = 0;
    }
}

static ReadFileResult PlatformReadFile(char *filepath) {
    ReadFileResult result = {};
    HANDLE fhandle = CreateFile(filepath, GENERIC_READ, FILE_SHARE_READ, 
                                0, OPEN_EXISTING, 0, 0);
    
    LARGE_INTEGER filesize64;
    if (!GetFileSizeEx(fhandle, &filesize64)) {
        // 
    }
    DWORD filesize32 = safe_truncate_uint64(filesize64.QuadPart);
    result.file = VirtualAlloc(0, filesize32, MEM_RESERVE|MEM_COMMIT, PAGE_READWRITE);
    if (!result.file) {
        result.size = 0;
    }

    DWORD bytesread;
    if (ReadFile(fhandle, result.file, filesize32, &bytesread, 0) && 
            (filesize32 == bytesread)) {
        result.size = filesize32;
    }

    else {
        PlatformFreeFile(result.file); 
        result.size = 0;
    }
    CloseHandle(fhandle);
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
            break;
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

    game_memory.transient_storage_size = gigabytes(1);
    game_memory.transient_storage = VirtualAlloc(
            0, game_memory.transient_storage_size, 
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

    OffscreenBuffer postfx_buffer = {};
    postfx_buffer.width = game_buffer.width;
    postfx_buffer.height = game_buffer.height;
    postfx_buffer.pitch = game_buffer.pitch;
    postfx_buffer.bytes_per_pixel = game_buffer.bytes_per_pixel;
    
    uint32_t postfx_buffer_size = postfx_buffer.width*postfx_buffer.height*postfx_buffer.bytes_per_pixel;
    postfx_buffer.memory = VirtualAlloc(0, postfx_buffer_size, MEM_RESERVE|MEM_COMMIT, PAGE_READWRITE);
    
    GameInput input = {}; 
    bool window_in_focus = true;
    RECT old_clip;
    GetClipCursor(&old_clip);
    LARGE_INTEGER previous_counter; 
    QueryPerformanceCounter(&previous_counter);
    RUNNING = true;
    while (RUNNING) {
        if (GetActiveWindow() != window) {
            window_in_focus = false;
        }

        else {
            window_in_focus = true;
        }

        LARGE_INTEGER start_counter = previous_counter;

        MSG msg;

        input.mouse_pos.x = 0;
        input.mouse_pos.y = 0;
        input.mouse_lbutton_click = false;
        input.action = false;

        static float time = 0;
        static uint32_t count;
        RECT client_rect;
        GetClientRect(window, &client_rect);
        int mid_x = (client_rect.left + client_rect.right) / 2;
        int mid_y = (client_rect.bottom + client_rect.top) / 2;
        POINT center = {};
        center.x = mid_x;
        center.y = mid_y;
        ClientToScreen(window, &center);
        RECT center_clip;
        center_clip.left = center.x-2;
        center_clip.right = center.x+2;
        center_clip.top = center.y-2;
        center_clip.bottom = center.y+2;
        if (window_in_focus) {
            SetCursorPos(center.x, center.y);
            ClipCursor(&center_clip);
        }

        else {
            ClipCursor(&old_clip);
        }
        
        ScreenToClient(window, &center);

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
                    if (window_in_focus) {
                        input.mouse_pos.x = GET_X_LPARAM(msg.lParam)-center.x;
                        input.mouse_pos.y = GET_Y_LPARAM(msg.lParam)-center.y;
                    }
                    OutputDebugString("MOUSEMOVE\n");
                }

                default: 
                    TranslateMessage(&msg);
                    DispatchMessageA(&msg);
            }
        }

        game_update_and_render(&game_memory, 
                               &game_buffer, 
                               &zbuffer,
                               &postfx_buffer,
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
