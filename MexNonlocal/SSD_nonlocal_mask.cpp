/* Distance de similarité sous hypothèse de bruit Cauchy */

#include "matrix.h"
#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*
 * similarity.cpp
 *
 * Auteur : Fabien Pierre
 * Date : 20/06/2017
 *
 *
 * Compilation Linux with Matlab linking:
 * mex  -largeArrayDims cauchy_similarity.cpp -output cauchy_similarity
 *
 * input : image gray scale single precision NxM
 * gamma (parameter of Cauchy estimator)
 * K (number of retained patches)
 * patch radius
 * search window radius
 *
 * output : single precision matrix, size NxMxK
 * for i,j, output(i,j,:) is the vector of the K nearest patch center
 *
 */


class FLOATMAP { public:
    int w, h;
    float *data;
    FLOATMAP(int w_, int h_) :w(w_), h(h_) { data = new float[w*h]; }
    FLOATMAP(int w_, int h_, float *data_) :w(w_), h(h_), data(data_) {}
    ~FLOATMAP() {  }
    float *operator[](int y) { return &data[y*w]; }
};


/*Global variables*/
int patch_w;
int W_radius;
int patch_radius;

int K;

/*Define the Cauchy similarity. */
float ssd_similarity(float & a, float & b){
    float ret = (a-b)*(a-b);
    return ret;
}


/* Measure distance between 2 patches
 * with center point (ax, ay) and (bx, by)
 ****Based on Cauchy similarity
 */
float dist(FLOATMAP *a, int ax, int ay, int bx, int by) {
    float ans = 0.0f;
    // Scaling factors
    float scalex = (patch_w^4)*(2.5*255.0/a->w)*(2.5*255.0/a->w);
    float scaley = (patch_w^4)*(2.5*255.0/a->h)*(2.5*255.0/a->h);
    
    /* compare patches*/
    for (int dy = -patch_radius; dy <= patch_radius; dy++) {
        
        int ady= (ay+dy < a->h)? ay+dy: 2*a->h-(ay+dy)-1;
        int bdy= (by+dy < a->h)? by+dy: 2*a->h-(by+dy)-1;
        ady= (ady < 0)? -ady-1 : ady;
        bdy= (bdy < 0)? -bdy-1 : bdy;
        
        for (int dx = -patch_radius; dx <= patch_radius; dx++) {
            int adx = (ax+dx < a->w) ? ax+dx :  2*a->w-(ax+dx)-1;
            int bdx = (bx+dx < a->w) ? bx+dx :  2*a->w-(bx+dx)-1;
            adx = (adx < 0)? -adx-1 : adx;
            bdx = (bdx < 0)? -bdx-1 : bdx;
            
            float ac = (*a)[ady][adx];
            float bc = (*a)[bdy][bdx];
            
            ans += ssd_similarity(ac, bc);
        }
    }
    // Modified distance
    ans += ((ax-bx)*(ax-bx))*scalex + ((ay-by)*(ay-by))*scaley;
    return ans;
}

float * distances;

int compare (const void * b, const void * a)
{
    float to_test =  distances[*(int*)a] - distances[*(int*)b] ;
    
    if(to_test>0)
        return -1;
    else if(to_test<0)
        return 1;
    else
        return 0;
}

/* search of the K neearest neighbors */
void KNN(FLOATMAP *a, FLOATMAP *mask,  float * k_best_centers) {
    
    int width  = a->w ;
    int height = a->h ;
    
    /*initialize variables*/
    int pixel_arround_in_window = (2*W_radius+1)*(2*W_radius+1);
    distances = (float*) malloc(pixel_arround_in_window*sizeof(float));
    int * indexes = (int*) malloc(pixel_arround_in_window*sizeof(int));
    int * idbx = (int*) malloc(pixel_arround_in_window*sizeof(int));
    int * idby = (int*) malloc(pixel_arround_in_window*sizeof(int));
    
    
    for(int i=0; i<pixel_arround_in_window; i++){
        distances[i]=0.0f;
        indexes[i]=i;
    }
    
    int aew = a->w;
    int aeh = a->h;
    
    for (int ay = 0; ay < aeh; ay ++) {
        for (int ax = 0; ax < aew; ax ++) {
            
            if ((*mask)[ay][ax]>0.5){
                /*Set similarities to zero*/
                for(int i=0; i< pixel_arround_in_window; i++){
                    distances[i]=0.0f;
                    indexes[i]=i;
                }
                
                /*define an accumulator*/
                int last_pix_in_windows = 0;
                
                /*define limit of the search window*/
                int ystart_windows = ay-W_radius < 0 ? 0 : ay-W_radius;
                int yend_windows = ay+W_radius >= aeh ? aeh-1 : ay+W_radius;
                int xstart_windows = ax-W_radius < 0 ? 0 : ax-W_radius;
                int xend_windows = ax+W_radius >= aew ? aew-1 : ax+W_radius;
                
                /*computation of distances*/
                for  (int by = ystart_windows; by <= yend_windows; by ++) {
                    for (int bx = xstart_windows; bx <= xend_windows; bx ++) {
                        distances[last_pix_in_windows]=
                                dist(a,  ax,  ay,  bx,  by);
                        idbx[last_pix_in_windows]= bx;
                        idby[last_pix_in_windows]= by;
                        last_pix_in_windows++;
                    }
                }
                
                
                /*sort idexes*/
                qsort (indexes, last_pix_in_windows, sizeof(int), compare);
                
                /*copy the K best patch centers.*/
                for(int id_k = 0; id_k < K; id_k ++){
                    k_best_centers[ax+ ay*width + id_k*height*width] =
                            distances[indexes[id_k+1]];
                    k_best_centers[height*width*K+ ax+ay*width + id_k*height*width] =
                            idbx[indexes[id_k+1]];
                    k_best_centers[2*height*width*K+ ax+ay*width + id_k*height*width] =
                            idby[indexes[id_k+1]];
                    
                }
            }
        }
    }
    
    free(distances);
    free(indexes);
    free(idbx);
    free(idby);
    
    return;
}




void mexFunction (int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    /* Vérification du nombre d'arguments d'entrée :
     * image input , gamma, K, patch width, search window radius
     **/
    if (nrhs != 5)
    {
        mexErrMsgIdAndTxt("KNNS:nrhs","5 inputs required");
    }
    
    if (nlhs != 1)
    {
        mexErrMsgIdAndTxt("KNNS:nlhs","1 output required");
    }
    
    // mexPrintf("\n--------------------------\n");
    // mexPrintf("Call KNN similarity Cauchy\n");
    // mexPrintf("--------------------------\n");
    
// Gestion entree-sortie
    
    const size_t * dims = (const size_t*)(mxGetDimensions(prhs[0]));
    int dim = mxGetNumberOfDimensions(prhs[0]);
    
    int wA= (int)dims[0];
    int hA = (int)dims[1];
    
    
    patch_w = (int)(mxGetScalar(prhs[2]));
    patch_radius = (patch_w-1)/2;
    
    W_radius = (int)(mxGetScalar(prhs[3]));
    K = (int)(mxGetScalar(prhs[1]));
    
    /* penser à vérifier en entrée que */
    if (((W_radius+1)*(W_radius+1))-1 < K)
    {
        mexWarnMsgIdAndTxt("KNNS:checkParam",
                "K is too large in comparison of window radius.\nFrom now K=%d\n",
                ((W_radius+1)*(W_radius+1))-1);
        K=((W_radius+1)*(W_radius+1))-1;
    }
    
    const int n_dim_result= 4;
    size_t* dims_result = (size_t *)(malloc(sizeof(size_t)*n_dim_result));
    dims_result[0]= (size_t)wA;
    dims_result[1]= (size_t)hA;
    dims_result[2]= (size_t)K;
    dims_result[3]= (size_t)3.0;
    
    plhs[0]= mxCreateNumericArray(n_dim_result, dims_result ,
            mxSINGLE_CLASS, mxREAL     ) ;
    
    FLOATMAP *A = new FLOATMAP (wA,hA,(float*)(mxGetPr(prhs[0])));
    FLOATMAP *Mask = new FLOATMAP (wA,hA,(float*)(mxGetPr(prhs[4])));
    
    // mexPrintf("size image:\t %d x %d \ngamma=\t\t %lf \n", wA,hA,gamma_param);
    // mexPrintf("K= \t\t\t%d \npatch width=\t\t%d\n",  K,  2*patch_radius+1);
    // mexPrintf("search window width=\t%d\n", 2*W_radius+1  );
    // mexPrintf("--------------------------\n\n");
    
    KNN(A, Mask,  (float*)(mxGetPr(plhs[0])));
    
    delete A;
    delete Mask;
    free(dims_result);
    
    return;
}