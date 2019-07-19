#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <algorithm>
#include "mex.h"

using namespace std;

/* Input Arguments */
#define	X_size_fin			prhs[0]
#define	Y_size_fin          prhs[1]
#define	UV_diameter         prhs[2]
#define	Im_in_remap         prhs[3]
#define Alpha_min           prhs[4]
#define Alpha_max           prhs[5]
#define Depth_res           prhs[6]
#define Orient              prhs[7]

/* Output Arguments */
#define	Depth               plhs[0]
#define Depth_var           plhs[1]
#define Depth_var_sum       plhs[2]
#define Refocus             plhs[3]

double alpha_min;
double alpha_max;
int depth_res;
double alpha_step;
int uv_diameter;
int uv_radius;
int uv_size;
int width;
int height;
int pixelNum;
int remap_height;
int remap_width;
int remap_pixelNum;
double ratio = 1;          // weight of refocus cue to correspondence cue
int bsz = 5;               // spatial window size of bilateral filter for depth cost
double sigma_range = 0.01; // range sigma of bilateral filter
float* refocus;            // the refocused image

int checkX(int x) {
    if (x < 0) return 0;
    else if (x >= width) return width-1;
    else return x;
}

int checkY(int y) {
    if (y < 0) return 0;
    else if (y >= height) return height-1;
    else return y;
}

// compute mean and variance
double getVar(double array[], int size, double* mean)
{
    *mean = 0;
    int num = 0;
    for (int i = 0; i < size; i++) {
        if (array[i] >= 0) {
            *mean += array[i];
            num++;
        }
    }
    *mean /= num;
    double var = 0;
    for (int i = 0; i < size; i++) {
        if (array[i] >= 0) {
            var += (array[i]-*mean) * (array[i]-*mean);
        }
    }
    if (num == 1) var = INT_MAX;
    else var /= (num-1);
    return var;
}

// interpolate pixels to gather the angular patch for a particular spatial pixel
void remapping(double im_in_remap[], int x, int y, int alpha_num, double remap[])
{
    double alpha = alpha_min + alpha_step * alpha_num;
    for (int i = -uv_radius; i <= uv_radius; i++) {
        for (int j = -uv_radius; j <= uv_radius; j++) {
            float y_ind   = (float)i*alpha + y;
            float x_ind   = (float)j*alpha + x;
            
            int x_floor = floor(x_ind);
            int y_floor = floor(y_ind);
            
            int x_1     = checkX(x_floor);
            int y_1     = checkY(y_floor);
            int x_2     = checkX(x_floor+1);
            int y_2     = checkY(y_floor+1);
            
            float x_1_w   = 1-(x_ind-x_floor);
            float x_2_w   = 1-x_1_w;
            float y_1_w   = 1-(y_ind-y_floor);
            float y_2_w   = 1-y_1_w;
            
            int x_1_index = j+uv_radius + (x_1)*uv_diameter;
            int y_1_index = i+uv_radius + (y_1)*uv_diameter;
            int x_2_index = j+uv_radius + (x_2)*uv_diameter;
            int y_2_index = i+uv_radius + (y_2)*uv_diameter;
            
            for (int c = 0; c < 3; c++) {
                remap[(i+uv_radius)*uv_diameter+(j+uv_radius)+c*uv_size] =
                y_1_w * x_1_w * im_in_remap[y_1_index+x_1_index*remap_height+c*remap_pixelNum] +
                y_2_w * x_1_w * im_in_remap[y_2_index+x_1_index*remap_height+c*remap_pixelNum] +
                y_1_w * x_2_w * im_in_remap[y_1_index+x_2_index*remap_height+c*remap_pixelNum] +
                y_2_w * x_2_w * im_in_remap[y_2_index+x_2_index*remap_height+c*remap_pixelNum];
            }
        }
    }
}

// compute the depth cost for a pixel
void computeDepthCost(int depth[], float depth_cost[], double im_in_remap[], double im_pinhole[], double orient[], int x, int y)
{
    double* remap    = new double[uv_size*3];
    double* remap_p1 = new double[uv_size*3];
    double* remap_p2 = new double[uv_size*3];
    double mean_p1[3], mean_p2[3];
    
    for (int alpha_num = 0; alpha_num < depth_res; alpha_num++) {
        remapping(im_in_remap, x, y, alpha_num, remap);
        
        if (orient[y+x*height] < -100*M_PI/180+0.1) { // no edge, apply traditional photo-consistency
            double var1 = getVar(remap, uv_size, mean_p1) + getVar(&remap[uv_size], uv_size, &mean_p1[1]) + getVar(&remap[uv_size*2], uv_size, &mean_p1[2]);
            for (int c = 0; c < 3; c++) {
                refocus[y + x*height + c*pixelNum + alpha_num*3*pixelNum] = mean_p1[c];
                refocus[y + x*height + c*pixelNum + alpha_num*3*pixelNum + depth_res*pixelNum*3] = mean_p1[c];
            }
            depth_cost[y + x*height + alpha_num*pixelNum] = sqrt(var1);
            depth_cost[y + x*height + alpha_num*pixelNum + depth_res*pixelNum] = sqrt(var1);
            continue;
        }
        
        // edge exists, compute the cost for both halves of the angular patch
        // generate the two half patches
        double theta = orient[y+x*height];
        for (int c = 0; c < 3; c++) {
            for (int i = 0; i < uv_diameter; i++) {
                for (int j = 0; j < uv_diameter; j++) {
                    remap_p1[i*uv_diameter + j + c*uv_size] = (uv_radius-i-sin(theta+M_PI/2) > (tan(theta) * (j-cos(theta+M_PI/2)-uv_radius)))?
                                                                remap[i*uv_diameter + j + c*uv_size] : -1;
                    remap_p2[i*uv_diameter + j + c*uv_size] = (uv_radius-i+sin(theta+M_PI/2) < (tan(theta) * (j+cos(theta+M_PI/2)-uv_radius)))?
                                                                remap[i*uv_diameter + j + c*uv_size] : -1;
                }
            }
        }
        // compute both the correspondence cue: patch variance, and refocus cue: patch mean
        double var1 = getVar(remap_p1, uv_size, mean_p1) + getVar(&remap_p1[uv_size], uv_size, &mean_p1[1]) + getVar(&remap_p1[uv_size*2], uv_size, &mean_p1[2]);
        double var2 = getVar(remap_p2, uv_size, mean_p2) + getVar(&remap_p2[uv_size], uv_size, &mean_p2[1]) + getVar(&remap_p2[uv_size*2], uv_size, &mean_p2[2]);
        
        // color consistency constraint
        double px1[3], px2[3], dif1 = 0, dif2 = 0;
        for (int c = 0; c < 3; c++) {
            px1[c] = im_pinhole[checkY((int)floor(y-2*sin(theta+M_PI/2)+0.5)) + height*checkX((int)floor(x+2*cos(theta+M_PI/2)+0.5)) + c*pixelNum];
            px2[c] = im_pinhole[checkY((int)floor(y+2*sin(theta+M_PI/2)+0.5)) + height*checkX((int)floor(x-2*cos(theta+M_PI/2)+0.5)) + c*pixelNum];
        }
        for (int c = 0; c < 3; c++) {
            dif1 += pow(mean_p1[c]-px1[c], 2) + pow(mean_p2[c]-px2[c], 2);
            dif2 += pow(mean_p1[c]-px2[c], 2) + pow(mean_p2[c]-px1[c], 2);
        }
        
        if (dif1 < dif2 + 0.05) { // pass color consistency check, apply the cost computed above
            for (int c = 0; c < 3; c++) {
                refocus[y + x*height + c*pixelNum + alpha_num*3*pixelNum] = var1<var2? mean_p1[c]:mean_p2[c];
                refocus[y + x*height + c*pixelNum + alpha_num*3*pixelNum + depth_res*pixelNum*3] = var1<var2? mean_p2[c]:mean_p1[c];
            }
            depth_cost[y + x*height + alpha_num*pixelNum] = min(sqrt(var1), sqrt(var2));
            depth_cost[y + x*height + alpha_num*pixelNum + depth_res*pixelNum] = max(sqrt(var1), sqrt(var2));
        }
        else { // fail color consistency check, apply largest cost
            for (int c = 0; c < 3; c++) {
                refocus[y + x*height + c*pixelNum + alpha_num*3*pixelNum] = -1;
                refocus[y + x*height + c*pixelNum + alpha_num*3*pixelNum + depth_res*pixelNum*3] = -1;
            }
            depth_cost[y + x*height + alpha_num*pixelNum] = 1;
            depth_cost[y + x*height + alpha_num*pixelNum + depth_res*pixelNum] = 1;
        }
    }
    delete[] remap;
    delete[] remap_p1;
    delete[] remap_p2;
}

// apply bilateral filter on the depth cost (may be optimized for better efficiency)
void costFilter(int depth[], float depth_cost[], float depth_cost_sum[], double im_pinhole[])
{
    for (int i = 0; i < depth_res; i++) {
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                double sumVal_c = 0, sumVal_f = 0;
                double weight_sum_c = 0, weight_sum_f = 0;
                for (int u = -bsz/2; u <= bsz/2; u++) {
                    for (int v = -bsz/2; v <= bsz/2; v++) {
                        int xu = checkX(x+u);
                        int yv = checkY(y+v);
                        
                        // compute bilateral weight
                        double color_dif = 0;
                        for (int c = 0; c < 3; c++)
                            color_dif += pow(im_pinhole[y+x*height+c*pixelNum]-im_pinhole[yv+xu*height+c*pixelNum], 2);
                        double weight = exp(-color_dif / (2*pow(sigma_range, 2)));
                        
                        // correspondence cost
                        sumVal_c += weight * depth_cost[xu*height+yv+i*pixelNum];
                        weight_sum_c += weight;
                        
                        // refocus cost
                        double gra = 0;
                        for (int c = 0; c < 3; c++) gra += fabs(refocus[yv + xu*height + c*pixelNum + i*3*pixelNum] - im_pinhole[yv + xu*height + c*pixelNum]);
                        sumVal_f += weight * gra;
                        weight_sum_f += weight;
                    }
                }
                sumVal_c = weight_sum_c != 0? sumVal_c / weight_sum_c : 1;
                sumVal_f = weight_sum_f != 0? sumVal_f / weight_sum_f : 1;
                depth_cost_sum[y+x*height+i*pixelNum] = sumVal_c + ratio*sumVal_f;
            }
        }
    }
}

// get the depth with the lowest depth cost
int getMinIdx(int depth[], float depth_cost_sum[], double im_pinhole[], int x, int y)
{
    int minIdx;
    double minVal = INT_MAX;
    for (int i = 0; i < depth_res; i++) {
        double sumVal = depth_cost_sum[y+x*height+i*pixelNum];
        if (sumVal < minVal) {
            minVal = sumVal;
            minIdx = i;
        }
    }
    return minIdx;
}

// the main depth estimation function
void depthEstimate(double* im_in_remap, double* im_pinhole, int* depth, float* depth_cost, float* depth_cost_sum, double* orient)
{
    // compute depth cost for each depth
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            computeDepthCost(depth, depth_cost, im_in_remap, im_pinhole, orient, x, y);
        }
    }
    
    // apply a bilateral filter on the depth cost
    costFilter(depth, depth_cost, depth_cost_sum, im_pinhole);
    
    // get the depth with the lowest depth cost
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            depth[x*height+y] = getMinIdx(depth, depth_cost_sum, im_pinhole, x, y);
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray* prhs[])
{
	/* Check for proper number of arguments */
	if (nrhs != 8)
		mexErrMsgTxt("Eight input arguments required.");
	else if (nlhs > 4)
		mexErrMsgTxt("Too many output arguments.");
    
    
	/* Assign pointers to the various parameters */
	double* x_size_fin_pt   = mxGetPr(X_size_fin);
	double* y_size_fin_pt   = mxGetPr(Y_size_fin);
	double* uv_diameter_pt  = mxGetPr(UV_diameter);
    double* im_in_remap_pt  = mxGetPr(Im_in_remap);
    double* alpha_min_pt    = mxGetPr(Alpha_min);
    double* alpha_max_pt    = mxGetPr(Alpha_max);
    double* depth_res_pt    = mxGetPr(Depth_res);
    double* orient_pt       = mxGetPr(Orient);
    
    
    // Assign parameters
    alpha_min = *alpha_min_pt;
    alpha_max = *alpha_max_pt;
    depth_res = *depth_res_pt;
    alpha_step = (alpha_max-alpha_min) / (depth_res-1);
    
    uv_diameter = *uv_diameter_pt;
    uv_radius = uv_diameter/2;
    uv_size = uv_diameter * uv_diameter;
    
    width = *x_size_fin_pt;
    height = *y_size_fin_pt;
    pixelNum = width*height;
    
    remap_height = height * uv_diameter;
    remap_width  = width  * uv_diameter;
    remap_pixelNum = remap_height * remap_width;
    
    /* Create the output matrix */
    Depth = mxCreateNumericMatrix(height, width, mxINT32_CLASS, mxREAL);
    int* depth_pt = (int*)mxGetData(Depth);
    
    const mwSize dims[] = { height, width, 2 };
    Depth_var = mxCreateNumericArray(3, dims, mxSINGLE_CLASS, mxREAL);
    float* depth_var_pt = (float*)mxGetPr(Depth_var);
    
    const mwSize dims2[] = { height, width, depth_res };
    Depth_var_sum = mxCreateNumericArray(3, dims2, mxSINGLE_CLASS, mxREAL);
    float* depth_var_sum_pt = (float*)mxGetPr(Depth_var_sum);
    
    const mwSize dimsIm[] = { height, width, 3, 2 };
    Refocus = mxCreateNumericArray(4, dimsIm, mxSINGLE_CLASS, mxREAL);
    float* refocus_pt = (float*)mxGetPr(Refocus);

    
    // Assign matrices
    refocus = new float[pixelNum*3*depth_res*2];
    float* depth_cost = new float[pixelNum*depth_res*2];
    double* im_pinhole = new double[pixelNum*3];
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            for (int c = 0; c < 3; c++) {
                im_pinhole[y + x*height + c*pixelNum] =
                    im_in_remap_pt[y*uv_diameter+uv_radius + (x*uv_diameter+uv_radius)*remap_height + c*remap_pixelNum];
            }
        }
    }

    /* Do the actual computations in a subroutine */
    depthEstimate(im_in_remap_pt, im_pinhole, depth_pt, depth_cost, depth_var_sum_pt, orient_pt);
    
    
    // Copy result to output matrix
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            depth_var_pt[y+x*height] = depth_cost[y + x*height + depth_pt[y+x*height]*pixelNum];
            depth_var_pt[y+x*height+pixelNum] = depth_cost[y + x*height + depth_pt[y+x*height]*pixelNum + depth_res*pixelNum];
            
            for (int c = 0; c < 3; c++) {
                refocus_pt[y + x*height + c*pixelNum] = refocus[y + x*height + c*pixelNum + 3*depth_pt[y+x*height]*pixelNum];
                refocus_pt[y + x*height + c*pixelNum + 3*pixelNum] = refocus[y + x*height + c*pixelNum + 3*depth_pt[y+x*height]*pixelNum + 3*depth_res*pixelNum];
            }
        }
    }
    
    delete[] refocus;
    delete[] depth_cost;
    delete[] im_pinhole;
}