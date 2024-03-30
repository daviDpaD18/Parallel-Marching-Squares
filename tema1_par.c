// Author: APD team, except where source was noted

#include "helpers.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <pthread.h>


#define CONTOUR_CONFIG_COUNT    16
#define FILENAME_MAX_SIZE       50
#define STEP                    8
#define SIGMA                   200
#define RESCALE_X               2048
#define RESCALE_Y               2048



#define min(a, b) ((a) < (b) ? (a) : (b))
#define CLAMP(v, min, max) if(v < min) { v = min; } else if(v > max) { v = max; }


struct struct_thread
{   

    long id;
    int P;
    ppm_image *image;
    ppm_image **contour_map;
    ppm_image *scaled_image;
    int step_x;
    int step_y;
    unsigned char **grid;
    pthread_barrier_t *barrier;
};

// Creates a map between the binary configuration (e.g. 0110_2) and the corresponding pixels
// that need to be set on the output image. An array is used for this map since the keys are
// binary numbers in 0-15. Contour images are located in the './contours' directory.
ppm_image **init_contour_map() {
    ppm_image **map = (ppm_image **)malloc(CONTOUR_CONFIG_COUNT * sizeof(ppm_image *));
    if (!map) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        char filename[FILENAME_MAX_SIZE];
        sprintf(filename, "./contours/%d.ppm", i);
        map[i] = read_ppm(filename);
    }

    return map;
}

// Updates a particular section of an image with the corresponding contour pixels.
// Used to create the complete contour image.
void update_image(ppm_image *image, ppm_image *contour, int x, int y) {
    for (int i = 0; i < contour->x; i++) {
        for (int j = 0; j < contour->y; j++) {
            int contour_pixel_index = contour->x * i + j;
            int image_pixel_index = (x + i) * image->y + y + j;

            image->data[image_pixel_index].red = contour->data[contour_pixel_index].red;
            image->data[image_pixel_index].green = contour->data[contour_pixel_index].green;
            image->data[image_pixel_index].blue = contour->data[contour_pixel_index].blue;
        }
    }
}

// Corresponds to step 1 of the marching squares algorithm, which focuses on sampling the image.
// Builds a p x q grid of points with values which can be either 0 or 1, depending on how the
// pixel values compare to the `sigma` reference value. The points are taken at equal distances
// in the original image, based on the `step_x` and `step_y` arguments.
unsigned char **sample_grid(ppm_image *image, int step_x, int step_y, unsigned char sigma) {
    int p = image->x / step_x;
    int q = image->y / step_y;

    unsigned char **grid = (unsigned char **)malloc((p + 1) * sizeof(unsigned char*));
    if (!grid) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i <= p; i++) {
        grid[i] = (unsigned char *)malloc((q + 1) * sizeof(unsigned char));
        if (!grid[i]) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
    }

    for (int i = 0; i < p; i++) {
        for (int j = 0; j < q; j++) {
            ppm_pixel curr_pixel = image->data[i * step_x * image->y + j * step_y];

            unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

            if (curr_color > sigma) {
                grid[i][j] = 0;
            } else {
                grid[i][j] = 1;
            }
        }
    }
    grid[p][q] = 0;

    // last sample points have no neighbors below / to the right, so we use pixels on the
    // last row / column of the input image for them
    for (int i = 0; i < p; i++) {
        ppm_pixel curr_pixel = image->data[i * step_x * image->y + image->x - 1];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > sigma) {
            grid[i][q] = 0;
        } else {
            grid[i][q] = 1;
        }
    }
    for (int j = 0; j < q; j++) {
        ppm_pixel curr_pixel = image->data[(image->x - 1) * image->y + j * step_y];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > sigma) {
            grid[p][j] = 0;
        } else {
            grid[p][j] = 1;
        }
    }

    return grid;
}

// Corresponds to step 2 of the marching squares algorithm, which focuses on identifying the
// type of contour which corresponds to each subgrid. It determines the binary value of each
// sample fragment of the original image and replaces the pixels in the original image with
// the pixels of the corresponding contour image accordingly.
void march(long id, int P ,ppm_image *image, unsigned char **grid, ppm_image **contour_map, int step_x, int step_y) {
    int p = image->x / step_x;
    int q = image->y / step_y;

    int start2 = id * (double)(p) / P;
    int end2 = min((id + 1) * (double)(p)/ P , p);


    for (int i = start2; i < end2; i++) {
        for (int j = 0; j < q; j++) {
            unsigned char k = 8 * grid[i][j] + 4 * grid[i][j + 1] + 2 * grid[i + 1][j + 1] + 1 * grid[i + 1][j];
            update_image(image, contour_map[k], i * step_x, j * step_y);
        }
    }
}

// Calls `free` method on the utilized resources.
void free_resources(ppm_image *image, ppm_image **contour_map, unsigned char **grid, int step_x) {
    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        free(contour_map[i]->data);
        free(contour_map[i]);
    }
    free(contour_map);

    for (int i = 0; i <= image->x / step_x; i++) {
        free(grid[i]);
    }
    free(grid);

    free(image->data);
    free(image);
}

void rescale_image(long id, int P ,ppm_image *image, ppm_image *new_image) {
    uint8_t sample[3];

    
    int start = id * (double)(2048) / P;
    int end = min((id + 1) * (double)(2048) / P, 2048);
    printf("%d %d\n", start, end);
    

    // use bicubic interpolation for scaling
    for (int i = start; i < end  ; i++) {
        for (int j = 0; j < 2048; j++) {
            float u = (float)i / (float)(2048 - 1);
            float v = (float)j / (float)(2048 - 1);
            sample_bicubic(image, u, v, sample);

            new_image->data[i * new_image->y + j].red = sample[0];
            new_image->data[i * new_image->y + j].green = sample[1];
            new_image->data[i * new_image->y + j].blue = sample[2];
        }
    }

    

 
}

void *thread_function(void *arg)    {
     uint8_t sample[3];

    struct struct_thread *struct_thread = (struct struct_thread *)arg;
    int step_x = struct_thread->step_x;
    int step_y = struct_thread->step_y;
    int i;
    int id = struct_thread->id;
    int start1, end1;
    int a;

    
    if(struct_thread->image->x <= RESCALE_X && struct_thread->image->y <= RESCALE_Y) {
        
        // imaginea nu trebuie rescalata  deci se face sampling direct 
        int p = struct_thread->image->x / step_x;
            int q = struct_thread->image->y / step_y;

            int start2 = id * (double)(p) / struct_thread->P;
            int end2 = min((id + 1) * (double)(p)/ struct_thread->P , p);

            for(i = start2; i < end2; i++) {
                for (int j = 0; j < q; j++) {
                ppm_pixel curr_pixel = struct_thread->image->data[i * step_x * struct_thread->image->y + j * step_y];

                unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

                if (curr_color > SIGMA) {
                    struct_thread->grid[i][j] = 0;
                } else {
                    struct_thread->grid[i][j] = 1;
                }
            }
            }
            struct_thread->grid[p][q] = 0;

            // last sample points have no neighbors below / to the right, so we use pixels on the
            // last row / column of the input image for them

        
            for(i = start2; i < end2; i++) {
                ppm_pixel curr_pixel = struct_thread->image->data[i * step_x * struct_thread->image->y + struct_thread->image->x - 1];

                unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

                if (curr_color > SIGMA) {
                    struct_thread->grid[i][q] = 0;
                } else {
                    struct_thread->grid[i][q] = 1;
            }
            }

            int start2_noneigh = id * (double)(q) / struct_thread->P;
            int end2_noneigh = min((id + 1) * (double)(q)/ struct_thread->P , q);
            for(i = start2_noneigh; i < end2_noneigh; i++) {
                ppm_pixel curr_pixel = struct_thread->image->data[(struct_thread->image->x - 1) * struct_thread->image->y + i * step_y];

                unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

                if (curr_color > SIGMA) {
                    struct_thread->grid[p][i] = 0;
                } else {
                    struct_thread->grid[p][i] = 1;
            }
            }
        
        //march the squares
        march(id, struct_thread->P, struct_thread->image, struct_thread->grid, struct_thread->contour_map, step_x, step_y);
        
    }
        else {
        
        rescale_image(id, struct_thread->P, struct_thread->image, struct_thread->scaled_image);
        pthread_barrier_wait(struct_thread->barrier);
        
        
        
        
        int p = struct_thread->scaled_image->x / step_x;
        int q = struct_thread->scaled_image->y / step_y;

        int start2 = id * (double)(p) / struct_thread->P;
        int end2 = min((id + 1) * (double)(p)/ struct_thread->P , p);

        for(i = start2; i < end2; i++) {
            for (int j = 0; j < q; j++) {
                ppm_pixel curr_pixel = struct_thread->scaled_image->data[i * step_x * struct_thread->scaled_image->y + j * step_y];

                unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

                if (curr_color > SIGMA) {
                    struct_thread->grid[i][j] = 0;
                } else {
                    struct_thread->grid[i][j] = 1;
                }
            }
        }
        struct_thread->grid[p][q] = 0;

        // last sample points have no neighbors below / to the right, so we use pixels on the
        // last row / column of the input image for them

        
        for(i = start2; i < end2; i++) {
            ppm_pixel curr_pixel = struct_thread->scaled_image->data[i * step_x * struct_thread->scaled_image->y + struct_thread->scaled_image->x - 1];

            unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

            if (curr_color > SIGMA) {
                struct_thread->grid[i][q] = 0;
            } else {
                struct_thread->grid[i][q] = 1;
            }
        }

        int start2_noneigh = id * (double)(q) / struct_thread->P;
        int end2_noneigh = min((id + 1) * (double)(q)/ struct_thread->P , q);
        for(i = start2_noneigh; i < end2_noneigh; i++) {
            ppm_pixel curr_pixel = struct_thread->scaled_image->data[(struct_thread->scaled_image->x - 1) * struct_thread->scaled_image->y + i * step_y];

            unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

            if (curr_color > SIGMA) {
               struct_thread->grid[p][i] = 0;
            } else {
                struct_thread->grid[p][i] = 1;
            }
        }
        

        //march the squares
        march(id, struct_thread->P, struct_thread->scaled_image, struct_thread->grid, struct_thread->contour_map, step_x, step_y);
       
            
        }
    
    
    
    pthread_exit(NULL);

}   

int main(int argc, char *argv[]) {
    if (argc < 4) {
        fprintf(stderr, "Usage: ./tema1 <in_file> <out_file> <P>\n");
        return 1;
    }
    
    pthread_barrier_t barrier;

    int nr_threads = atoi(argv[3]);
    pthread_t tid[nr_threads];
    long i;
    int a;
    int r;

    a = pthread_barrier_init(&barrier, NULL, nr_threads);
    if(a != 0) {
        printf("Error creating barrier\n");
        return 1;
    }

    
    int step_x = STEP;
    int step_y = STEP;

    // 0. Initialize contour map

    ppm_image *image = read_ppm(argv[1]);
    ppm_image **contour_map = init_contour_map();
    ppm_image *new_image = (ppm_image*)malloc(sizeof(ppm_image));
    if (!new_image) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }
    new_image->x = RESCALE_X;
    new_image->y = RESCALE_Y;

    new_image->data = (ppm_pixel*)malloc(new_image->x * new_image->y * sizeof(ppm_pixel));
    if (!new_image) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    int p = min(2048, image->x) / step_x;
    int q = min(2048, image->y) / step_y;
    unsigned char **grid = (unsigned char **)malloc((p + 1) * sizeof(unsigned char*));
    if (!grid) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i <= p; i++) {
        grid[i] = (unsigned char *)malloc((q + 1) * sizeof(unsigned char));
        if (!grid[i]) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
    }
     
   struct struct_thread struct_thread[nr_threads] ;

    for(i = 0; i < nr_threads; i++) {
        struct_thread[i].id = i;
        struct_thread[i].P = nr_threads;
        struct_thread[i].image = image;
        struct_thread[i].contour_map = contour_map;
        struct_thread[i].scaled_image = new_image;
        struct_thread[i].scaled_image->x = RESCALE_X;
        struct_thread[i].scaled_image->y = RESCALE_Y;
        struct_thread[i].step_x = step_x;
        struct_thread[i].step_y = step_y;
        struct_thread[i].grid = grid;
        struct_thread[i].barrier = &barrier;

       r = pthread_create(&(tid[i]), NULL, thread_function, &(struct_thread[i]));
        if(r != 0) {
            printf("Error creating thread\n");
            return 1;
        }
    }


    for(i = 0; i < nr_threads; i++) {
        r = pthread_join(tid[i], NULL);
        if(r != 0) {
            printf("Error joining thread\n");
            return 1;
        }
    }

    // 4. Write output
    
    if(image->x <= RESCALE_X && image->y <= RESCALE_Y) {
        write_ppm(image, argv[2]);
        
    }
    else {
        
     
        write_ppm(new_image, argv[2]);
        
    }
        

    //free_resources(image, contour_map, struct_thread->grid, struct_thread->step_x);
    //free(struct_thread->scaled_image->data);
    //free(struct_thread->scaled_image);
    pthread_exit(NULL);
    pthread_barrier_destroy(&barrier);

    return 0;





}
