/********************************************************
 * Kernels to be optimized for the CS:APP Performance Lab
 ********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "defs.h"

/* 
 * Please fill in the following team struct 
 */
team_t team = {
    "shape_of_june",              /* Team name */

    "Junhyung Kim",     /* First member full name */
    "kyle14916@naver.com",  /* First member email address */

    "",                   /* Second member full name (leave blank if none) */
    "" /* Second member email addr (leave blank if none) */
};

/***************
 * ROTATE KERNEL
 ***************/

/******************************************************
 * Your different versions of the rotate kernel go here
 ******************************************************/

/* 
 * naive_rotate - The naive baseline version of rotate 
 */
char naive_rotate_descr[] = "naive_rotate: Naive baseline implementation";
void naive_rotate(int dim, pixel *src, pixel *dst) 
{
    int i, j;

    for (i = 0; i < dim; i++)
	for (j = 0; j < dim; j++)
	    dst[RIDX(dim-1-j, i, dim)] = src[RIDX(i, j, dim)];
}

/* 
 * rotate - Your current working version of rotate
 * IMPORTANT: This is the version you will be graded on
 */
char rotate_descr[] = "rotate: Current working version";
void rotate(int dim, pixel *src, pixel *dst) 
{
    int i, j;

	for (j = 0; j < dim; j++)
    for (i = 0; i < dim; i++)
	    dst[(dim-1-j) * dim + i] = src[i * dim + j];

    return;
}

char rotate0_descr[] = "rotate: version0";
void rotate0(int dim, pixel *src, pixel *dst) 
{
    int i, j;

	for (j = 0; j < dim; j++)
    for (i = 0; i < dim; i++)
	    dst[RIDX(dim-1-j, i, dim)] = src[RIDX(i, j, dim)];

    return;
}

char rotate1_descr[] = "rotate: version1";
void rotate1(int dim, pixel *src, pixel *dst) 
{
    int i, j;

	for (j = 0; j < dim; j++)
    for (i = 0; i < dim; i++)
	    dst[(dim-1-j) * dim + i] = src[i * dim + j];

    return;
}

char rotate2_descr[] = "rotate: version2";
void rotate2(int dim, pixel *src, pixel *dst) 
{
    int i, j;
    int size = dim*dim;
    pixel tmpdst[size];

    for (i = 0; i < dim; i++){
	for (j = 0; j < dim; j++)
	    tmpdst[(dim-1-j) * dim + i] = src[i * dim + j];
    }

	for(i=0;i<size;i++){
        dst[i] = tmpdst[i];
    }

    return;
}

char rotate3_descr[] = "rotate: version3";
void rotate3(int dim, pixel *src, pixel *dst) 
{
    int i, j;

    int tdst = (dim-1)*dim, tsrc = 0;
	for(j=0;j<dim;j++){
        for(i=0;i<dim;i++){
            dst[tdst++] = src[tsrc];
            tsrc += dim;       
        }
        tdst -= 2*dim;
        tsrc -= dim*dim-1;
    }
    return;
}

char rotate4_descr[] = "rotate: version4";
void rotate4(int dim, pixel *src, pixel *dst) 
{
    int i, j;
    int t1,t2;
    pixel tmpdst[dim*dim];

    for(i=0;i<dim;i++){
        for(j=i;j<dim;j++){
            t1 = i*dim+j;
            t2 = j*dim+i;
            tmpdst[t1] = src[t2];
            tmpdst[t2] = src[t1];
        }
    }

    for(i=0;i<dim/2;i++){
        for(j=0;j<dim;j++){
            t1 = i*dim+j;
            t2 = (dim-1-i)*dim+j;
            dst[t1] = tmpdst[t2];
            dst[t2] = tmpdst[t1];
        }
    }

    return;
}

/*********************************************************************
 * register_rotate_functions - Register all of your different versions
 *     of the rotate kernel with the driver by calling the
 *     add_rotate_function() for each test function. When you run the
 *     driver program, it will test and report the performance of each
 *     registered test function.  
 *********************************************************************/

void register_rotate_functions() 
{
    add_rotate_function(&naive_rotate, naive_rotate_descr);   
    add_rotate_function(&rotate, rotate_descr);     

    /* ... Register additional test functions here */
    // add_rotate_function(&rotate0, rotate0_descr);
    // add_rotate_function(&rotate1, rotate1_descr);   
    // add_rotate_function(&rotate2, rotate2_descr);   
    // add_rotate_function(&rotate3, rotate3_descr); 
    // add_rotate_function(&rotate4, rotate4_descr); 
}


/***************
 * SMOOTH KERNEL
 **************/

/***************************************************************
 * Various typedefs and helper functions for the smooth function
 * You may modify these any way you like.
 **************************************************************/

/* A struct used to compute averaged pixel value */
typedef struct {
    int red;
    int green;
    int blue;
    int num;
} pixel_sum;

/* Compute min and max of two integers, respectively */
static int min(int a, int b) { return (a < b ? a : b); }
static int max(int a, int b) { return (a > b ? a : b); }

/* 
 * initialize_pixel_sum - Initializes all fields of sum to 0 
 */
static void initialize_pixel_sum(pixel_sum *sum) 
{
    sum->red = sum->green = sum->blue = 0;
    sum->num = 0;
    return;
}

/* 
 * accumulate_sum - Accumulates field values of p in corresponding 
 * fields of sum 
 */
static void accumulate_sum(pixel_sum *sum, pixel p) 
{
    sum->red += (int) p.red;
    sum->green += (int) p.green;
    sum->blue += (int) p.blue;
    sum->num++;
    return;
}

/* 
 * assign_sum_to_pixel - Computes averaged pixel value in current_pixel 
 */
static void assign_sum_to_pixel(pixel *current_pixel, pixel_sum sum) 
{
    current_pixel->red = (unsigned short) (sum.red/sum.num);
    current_pixel->green = (unsigned short) (sum.green/sum.num);
    current_pixel->blue = (unsigned short) (sum.blue/sum.num);
    return;
}

/* 
 * avg - Returns averaged pixel value at (i,j) 
 */
static pixel avg(int dim, int i, int j, pixel *src) 
{
    int ii, jj;
    pixel_sum sum;
    pixel current_pixel;

    initialize_pixel_sum(&sum);
    for(ii = max(i-1, 0); ii <= min(i+1, dim-1); ii++) 
	for(jj = max(j-1, 0); jj <= min(j+1, dim-1); jj++)
	    accumulate_sum(&sum, src[RIDX(ii, jj, dim)]);

    assign_sum_to_pixel(&current_pixel, sum);
    return current_pixel;
}

/******************************************************
 * Your different versions of the smooth kernel go here
 ******************************************************/

/*
 * naive_smooth - The naive baseline version of smooth 
 */
char naive_smooth_descr[] = "naive_smooth: Naive baseline implementation";
void naive_smooth(int dim, pixel *src, pixel *dst) 
{
    int i, j;

    for (i = 0; i < dim; i++)
	for (j = 0; j < dim; j++)
	    dst[RIDX(i, j, dim)] = avg(dim, i, j, src);
}

/*
 * smooth - Your current working version of smooth. 
 * IMPORTANT: This is the version you will be graded on
 */


char smooth_descr[] = "smooth: Current working version";
void smooth(int dim, pixel *src, pixel *dst) 
{
    int i,j,cur;
    // corners
    cur = 0;
    dst[cur].red = (src[cur].red+src[cur+1].red+src[cur+dim].red+src[cur+dim+1].red)/4;
    dst[cur].green = (src[cur].green+src[cur+1].green+src[cur+dim].green+src[cur+dim+1].green)/4;
    dst[cur].blue = (src[cur].blue+src[cur+1].blue+src[cur+dim].blue+src[cur+dim+1].blue)/4;

    cur = dim-1;
    dst[cur].red = (src[cur].red+src[cur-1].red+src[cur+dim].red+src[cur+dim-1].red)/4;
    dst[cur].green = (src[cur].green+src[cur-1].green+src[cur+dim].green+src[cur+dim-1].green)/4;
    dst[cur].blue = (src[cur].blue+src[cur-1].blue+src[cur+dim].blue+src[cur+dim-1].blue)/4;

    cur = dim*(dim-1);
    dst[cur].red = (src[cur].red+src[cur+1].red+src[cur-dim].red+src[cur-dim+1].red)/4;
    dst[cur].green = (src[cur].green+src[cur+1].green+src[cur-dim].green+src[cur-dim+1].green)/4;
    dst[cur].blue = (src[cur].blue+src[cur+1].blue+src[cur-dim].blue+src[cur-dim+1].blue)/4;

    cur = dim*dim-1;
    dst[cur].red = (src[cur].red+src[cur-1].red+src[cur-dim].red+src[cur-dim-1].red)/4;
    dst[cur].green = (src[cur].green+src[cur-1].green+src[cur-dim].green+src[cur-dim-1].green)/4;
    dst[cur].blue = (src[cur].blue+src[cur-1].blue+src[cur-dim].blue+src[cur-dim-1].blue)/4;

    // edges
    for(cur=1;cur<dim-1;cur++){
        dst[cur].red = (src[cur-1].red+src[cur].red+src[cur+1].red+src[cur+dim-1].red+src[cur+dim].red+src[cur+dim+1].red)/6;
        dst[cur].green = (src[cur-1].green+src[cur].green+src[cur+1].green+src[cur+dim-1].green+src[cur+dim].green+src[cur+dim+1].green)/6;
        dst[cur].blue = (src[cur-1].blue+src[cur].blue+src[cur+1].blue+src[cur+dim-1].blue+src[cur+dim].blue+src[cur+dim+1].blue)/6;
    }

    for(cur=dim*(dim-1)+1;cur<dim*dim-1;cur++){
        dst[cur].red = (src[cur-1].red+src[cur].red+src[cur+1].red+src[cur-dim-1].red+src[cur-dim].red+src[cur-dim+1].red)/6;
        dst[cur].green = (src[cur-1].green+src[cur].green+src[cur+1].green+src[cur-dim-1].green+src[cur-dim].green+src[cur-dim+1].green)/6;
        dst[cur].blue = (src[cur-1].blue+src[cur].blue+src[cur+1].blue+src[cur-dim-1].blue+src[cur-dim].blue+src[cur-dim+1].blue)/6;
    }

    for(cur=dim;cur<dim*(dim-1);cur+=dim){
        dst[cur].red = (src[cur-dim].red+src[cur-dim+1].red+src[cur].red+src[cur+1].red+src[cur+dim].red+src[cur+dim+1].red)/6;
        dst[cur].green = (src[cur-dim].green+src[cur-dim+1].green+src[cur].green+src[cur+1].green+src[cur+dim].green+src[cur+dim+1].green)/6;
        dst[cur].blue = (src[cur-dim].blue+src[cur-dim+1].blue+src[cur].blue+src[cur+1].blue+src[cur+dim].blue+src[cur+dim+1].blue)/6;
    }

    for(cur=2*dim-1;cur<dim*(dim-1);cur+=dim){
        dst[cur].red = (src[cur-dim-1].red+src[cur-dim].red+src[cur-1].red+src[cur].red+src[cur+dim-1].red+src[cur+dim].red)/6;
        dst[cur].green = (src[cur-dim-1].green+src[cur-dim].green+src[cur-1].green+src[cur].green+src[cur+dim-1].green+src[cur+dim].green)/6;
        dst[cur].blue = (src[cur-dim-1].blue+src[cur-dim].blue+src[cur-1].blue+src[cur].blue+src[cur+dim-1].blue+src[cur+dim].blue)/6;
    }

    // inners
    for(i=1;i<dim-1;i++){
        for(j=1;j<dim-1;j++){
            cur = i*dim+j;
            dst[cur].red = (src[cur-dim-1].red+src[cur-dim].red+src[cur-dim+1].red+
                            src[cur-1].red+src[cur].red+src[cur+1].red+
                            src[cur+dim-1].red+src[cur+dim].red+src[cur+dim+1].red)/9;
            dst[cur].green = (src[cur-dim-1].green+src[cur-dim].green+src[cur-dim+1].green+
                            src[cur-1].green+src[cur].green+src[cur+1].green+
                            src[cur+dim-1].green+src[cur+dim].green+src[cur+dim+1].green)/9;
            dst[cur].blue = (src[cur-dim-1].blue+src[cur-dim].blue+src[cur-dim+1].blue+
                            src[cur-1].blue+src[cur].blue+src[cur+1].blue+
                            src[cur+dim-1].blue+src[cur+dim].blue+src[cur+dim+1].blue)/9;                                                        
        }
    }


    return;
}

char smooth0_descr[] = "smooth: version0";
void smooth0(int dim, pixel *src, pixel *dst) 
{
    int i,j;
    int ii,jj;
    pixel_sum tmpdst[dim*dim];

    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++) initialize_pixel_sum(&tmpdst[i*dim+j]);
    }
    
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            for(ii = max(i-1, 0); ii <= min(i+1, dim-1); ii++){
                for(jj = max(j-1, 0); jj <= min(j+1, dim-1); jj++){
                    accumulate_sum(&tmpdst[ii*dim+jj],src[i*dim+j]);
                }
            }
        }
    }

    int size = dim*dim;
    for(i=0;i<size;i++) assign_sum_to_pixel(&dst[i],tmpdst[i]);
    return;
}

char smooth1_descr[] = "smooth: version1";
void smooth1(int dim, pixel *src, pixel *dst) 
{
    int i,j,ii,jj,ni,nj;
    int ind[3] = {-1,0,1};
    
    pixel_sum tmp;

    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            initialize_pixel_sum(&tmp);
            for(ii=0;ii<3;ii++){
                ni = i+ind[ii];
                if(ni<0 || ni>=dim) continue;
                for(jj=0;jj<3;jj++){
                    nj = j+ind[jj];
                    if(nj<0 || nj>=dim) continue;
                    accumulate_sum(&tmp,src[ni*dim+nj]);
                }
            }
            assign_sum_to_pixel(&dst[i*dim+j],tmp);
        }
    }

    return;
}

char smooth2_descr[] = "smooth: version2";
void smooth2(int dim, pixel *src, pixel *dst) 
{
    int i,j,ii,jj;
    int n = 4;

    for(i=0;i<dim;i+=n){
        for(j=0;j<dim;j+=n){
            for(ii=i;ii<i+n;ii++){
                for(jj=j;jj<j+n;jj++){
                    dst[ii*dim+jj] = avg(dim,ii,jj,src);
                }
            }
        }
    }

    return;
}

char smoot3_descr[] = "smooth: Current working version";
void smooth3(int dim, pixel *src, pixel *dst) 
{
    int i,j,k;
    long long dp[dim+1][dim+1][3];
    int bi,si,bj,sj;
    for(i=0;i<dim+1;i++){
        for(k=0;k<3;k++){
            dp[0][i][k] = 0;
            dp[i][0][k] = 0;
        }
    }

    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            dp[i+1][j+1][0] = (dp[i+1][j][0] + dp[i][j+1][0] - dp[i][j][0] +
            (long long) src[i*dim+j].blue);
            dp[i+1][j+1][1] = (dp[i+1][j][1] + dp[i][j+1][1] - dp[i][j][1] +
            (long long) src[i*dim+j].green);
            dp[i+1][j+1][2] = (dp[i+1][j][2] + dp[i][j+1][2] - dp[i][j][2] +
            (long long) src[i*dim+j].red);
        }
    }

    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            si = max(0,i-1);
            sj = max(0,j-1);
            bi = min(dim,i+2);
            bj = min(dim,j+2);
            dst[i*dim+j].blue = (unsigned short)((dp[bi][bj][0] - dp[bi][sj][0] - 
            dp[si][bj][0] + dp[si][sj][0])/((bi-si)*(bj-sj)));
            dst[i*dim+j].green = (unsigned short)((dp[bi][bj][1] - dp[bi][sj][1] - 
            dp[si][bj][1] + dp[si][sj][1])/((bi-si)*(bj-sj)));
            dst[i*dim+j].red = (unsigned short)((dp[bi][bj][2] - dp[bi][sj][2] - 
            dp[si][bj][2] + dp[si][sj][2])/((bi-si)*(bj-sj)));
        }
    }
    return;
}


/********************************************************************* 
 * register_smooth_functions - Register all of your different versions
 *     of the smooth kernel with the driver by calling the
 *     add_smooth_function() for each test function.  When you run the
 *     driver program, it will test and report the performance of each
 *     registered test function.  
 *********************************************************************/

void register_smooth_functions() {
    add_smooth_function(&naive_smooth, naive_smooth_descr);
    add_smooth_function(&smooth, smooth_descr);

    /* ... Register additional test functions here */
    // add_smooth_function(&smooth0, smooth0_descr);
    // add_smooth_function(&smooth1, smooth1_descr);
    // add_smooth_function(&smooth2, smooth2_descr);
    // add_smooth_function(&smooth3, smooth3_descr);
}
