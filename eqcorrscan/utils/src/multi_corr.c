/*
 * =====================================================================================
 *
 *       Filename:  multi_corr.c
 *
 *        Purpose:  Routines for computing cross-correlations
 *
 *        Created:  03/07/17 02:25:07
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (Calum Chamberlain),
 *   Organization:  EQcorrscan
 *      Copyright:  EQcorrscan developers.
 *        License:  GNU Lesser General Public License, Version 3
 *                  (https://www.gnu.org/copyleft/lesser.html)
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Prototypes
int normxcorr_fftw_1d(float *signala, int a_len, float *signalb, int b_len, float *ncc, int N);

int xcorr_fftw_1d(float *signala, int a_len, float *signalb, int b_len, float *ncc, int N);

int run_std_mean(int a_len, float *signalb, int b_len, float *run_std, float *run_mean);

int xcorr (float *signala, int a_len, float *signalb, int b_len, float *ccc);

int multi_corr (float *templates, int template_len, int n_templates, float *image, int image_len, float *ccc);

int multi_normalise(float *ccc, int ccc_len, float *image, float *norm_sum, int template_len, int n);

// Functions
int run_std_mean(int a_len, float *signalb, int b_len, float *run_std, float *run_mean){
    int i;
	double sum = 0.0, mean, stdev, old_mean, var=0.0, new_samp, old_samp;

	for (i=0; i < a_len; ++i){
		sum += (double) signalb[i];
	}
	mean = sum / a_len;

	// Compute starting standard deviation
	for (i=0; i < a_len; ++i){
		var += pow(signalb[i] - mean, 2) / (a_len);
	}
	stdev = sqrt(var);

	run_std[0] = (float) stdev;
	run_mean[0] = (float) mean;
	for(i = 1; i < b_len; ++i){
	    new_samp = signalb[i + a_len - 1];
	    old_samp = signalb[i - 1];
		old_mean = mean;
		mean = mean + (new_samp - old_samp) / a_len;
		var += (new_samp - old_samp) * (new_samp - mean + old_samp - old_mean) / (a_len);
		stdev = sqrt(var);
        run_mean[i] = (float) mean;
        run_std[i] = (float) stdev;
	}
	return 0;
}

int normxcorr_fftw_1d(float *signala, int a_len, float *signalb, int b_len,
				  float *ncc, int N){
  /*
  Purpose: compute frequency domain normalised cross-correlation of real data using fftw
  Author: Calum J. Chamberlain
  Date: 12/06/2017
  Args:
	signala:  Template signal
	a_len:    Length of signala
	signalb:  Image signal (to scan through)
	b_len:    Length of signalb
	ncc:      Output for cross-correlation - should be pointer to memory -
			  must be b_len - a_len + 1 long
	N:        Size for fft
  */
    printf("A hollow shell of a function\n");
	return 0;
}


int xcorr_fftw_1d(float *signala, int a_len, float *signalb, int b_len,
				  float *ncc, int N){
  /*
  Purpose: compute frequency domain cross-correlation of real data using fftw
  Author: Calum J. Chamberlain

  Note: This is NOT normalised

  Date: 12/06/2017
  Args:
	signala:  Template signal
	a_len:    Length of signala
	signalb:  Image signal (to scan through)
	b_len:    Length of signalb
	ncc:      Output for cross-correlation - should be pointer to memory -
			  must be b_len - a_len + 1 long
	N:        Size for fft
  */
    printf("A hollow shell of a function\n");

	return 0;
}


int xcorr(float *signala, int a_len, float *signalb, int b_len, float *ccc){
	int p, k;
	int steps = b_len - a_len + 1;
	float numerator, denom;
	float auto_a = 0.0, auto_b = 0.0;

	for(p = 0; p < a_len; ++p){
		auto_a += signala[p] * signala[p];
	}
	for(k = 0; k < steps; ++k){
		numerator = 0.0;
		auto_b = 0.0;
		for(p = 0; p < a_len; ++p){
			numerator += signala[p] * signalb[p + k];
		}
		for(p = 0; p < a_len; ++p){
			auto_b += signalb[p + k] * signalb[p + k];
		}
		denom = sqrtf(auto_a * auto_b);
		ccc[k] = numerator / denom;
	}
	return 0;
}


int multi_corr(float *templates, int template_len, int n_templates, float *image, int image_len, float *ccc){
	int i;
	#pragma omp parallel for
	for (i = 0; i < n_templates; ++i){
		xcorr(&templates[template_len * i], template_len, image, image_len, &ccc[(image_len - template_len) * i]);
	}
	return 0;
}


int multi_normalise(float *ccc, int ccc_len, float *image, float *norm_sum, int template_len, int n)
{
	int i, k;
	float mean, old_mean, std, sum=0.0, var=0.0;
	float acceptedDiff = 0.00000001;

	// Compute starting mean, will update this
	for (i=0; i < template_len; ++i)
	{
		sum += image[i];
	}
	mean = sum / template_len;

	// Compute starting standard deviation, will update this
	for (i=0; i < template_len; ++i)
	{
		var += powf(image[i] - mean, 2);
	}
	std = sqrtf(var / (template_len));
	if (std < acceptedDiff)
	{
		for (k=0; k<n; ++k)
		{
			ccc[k] = 0.0;
		}
	}
	else
	{
		for (k=0; k<n; ++k)
		{
			ccc[k] = (ccc[k] - norm_sum[k] * mean ) / std;
		}
	}
	if (ccc[0] != ccc[0])
	{
		ccc[0] = 0.0;
	}
	// Loop through updating as we go
	for(i=1; i<ccc_len; ++i)
	{
		old_mean = mean;
		mean = mean + (image[i + template_len - 1] - image[i - 1]) / template_len;
		// Don't know of a usefully accurate and efficient method :(
		var += (image[i + template_len - 1] - image[i - 1]) * (image[i + template_len -1] - mean + image[i - 1] - old_mean) / template_len;
		std = sqrtf(var / template_len);
		if (std < acceptedDiff)
		{
			for (k=0; k<n; ++k)
			{
				ccc[(i * n) + k] = 0.0;
			}
		}
		else
		{
			// Normalize by the std and mean
			for (k=0; k<n; ++k)
			{
				ccc[(i * n) + k] = (ccc[(i * n) + k] - norm_sum[k] * mean ) / std;
			}
		}
		// Convert nans
		if (ccc[i] != ccc[i])
		{
			ccc[i] = 0.0;
		}
	}
	return 0;
}

