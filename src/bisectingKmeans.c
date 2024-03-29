/*
    Kalign - a multiple sequence alignment program

    Copyright 2006, 2019 Timo Lassmann

    This file is part of kalign.

    Kalign is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

*/

#include <xmmintrin.h>

#include <float.h>
#include "msa.h"
#include "bpm.h"
#include "bisectingKmeans.h"

#include "euclidean_dist.h"
#include "sequence_distance.h"

#include "pick_anchor.h"

struct node{
        struct node* left;
        struct node* right;
        int* samples;
        int num_samples;
        int id;
};


struct kmeans_result{
        int* sl;
        int* sr;
        int nl;
        int nr;
        float score;
};

static struct kmeans_result* alloc_kmeans_result(int num_samples);
static void free_kmeans_results(struct kmeans_result* k);

struct node* upgma(float **dm,int* samples, int numseq);
struct node* alloc_node(void);

int merge_clusters(struct node*n, struct msa* msa, int threshold);
int test_for_merge(struct node* n, struct msa* msa, int threshold);

int label_internal(struct node*n, int label);
void printTree(struct node* curr,int depth);
struct node* bisecting_kmeans(struct msa* msa, struct node* n, float** dm,int* samples,int numseq, int num_anchors,int num_samples,struct rng_state* rng);


int build_tree_kmeans(struct msa* msa)
{
        //struct drand48_data randBuffer;
        struct node* root = NULL;
        float** dm = NULL;
        int* tree = NULL;
        int* samples = NULL;
        int* anchors = NULL;
        int num_anchors;
        int numseq;

        int i;

        struct rng_state* rng;


        rng = init_rng(42);
        ASSERT(msa != NULL, "No alignment.");
        //ASSERT(param != NULL, "No input parameters.");



        // nbeed to alloc tree.... (?)


        numseq = msa->numseq;

        DECLARE_TIMER(timer);
        /* pick anchors . */
        LOG_MSG("Calculating pairwise distances");
        START_TIMER(timer);
        RUNP(anchors = pick_anchor(msa, &num_anchors));

        RUNP(dm = d_estimation(msa, anchors, num_anchors,0));//les,int pair)

        STOP_TIMER(timer);

        LOG_MSG("Done in %f sec.", GET_TIMING(timer));
        LOG_MSG("Using %d anchors.", num_anchors);
        MFREE(anchors);

        MMALLOC(samples, sizeof(int)* numseq);
        for(i = 0; i < numseq;i++){
                samples[i] = i;
        }



        //RUNP(root = alloc_node());

        START_TIMER(timer);

        LOG_MSG("Building guide tree.");

        //need RNG!
        RUNP(root = bisecting_kmeans(msa,root, dm, samples, numseq, num_anchors, numseq, rng));
        STOP_TIMER(timer);

        LOG_MSG("Done in %f sec.", GET_TIMING(timer));

        LOG_MSG("First Merge.");
        START_TIMER(timer);
        //merge_clusters(root, msa, 4);
        STOP_TIMER(timer);

        LOG_MSG("Done in %f sec.", GET_TIMING(timer));

        LOG_MSG("Second Merge.");
        START_TIMER(timer);
        //merge_clusters(root, msa, 4);
        STOP_TIMER(timer);

        LOG_MSG("Done in %f sec.", GET_TIMING(timer));

        exit(0);
        MFREE(root);
        for(i =0 ; i < msa->numseq;i++){
                _mm_free(dm[i]);
        }
        MFREE(dm);

        return OK;
ERROR:
        return FAIL;
}

int merge_clusters(struct node*n, struct msa* msa, int threshold)
{
        RUN(test_for_merge(n,msa,threshold));

        if(n->left){
                return merge_clusters(n->left,msa,threshold);
        }
        RUN(test_for_merge(n,msa,threshold));

        if(n->right){
                return merge_clusters(n->right,msa,threshold);
        }
        RUN(test_for_merge(n,msa,threshold));
        return OK;
ERROR:
        return FAIL;
}

int test_for_merge(struct node* n, struct msa* msa, int threshold)
{
        int* samples_a;
        int* samples_b;
        uint8_t* a;
        uint8_t* b;
        int len_a,len_b;

        int num_a,num_b;
        uint8_t d;
        int i,j;


        if(n->left && n->right){
                if(n->left->num_samples && n->right->num_samples){
                        samples_a = n->left->samples;
                        num_a = n->left->num_samples;
                        samples_b = n->right->samples;
                        num_b = n->right->num_samples;
                        LOG_MSG("Merging: %d and %d", num_a,num_b);
                        for(i = 0; i < num_a;i++){
                                a = msa->sequences[samples_a[i]]->s;
                                len_a = msa->sequences[samples_a[i]]->len;
                                for(j = 0;j < num_b;j++){
                                        b = msa->sequences[samples_b[j]]->s;
                                        len_b = msa->sequences[samples_b[j]]->len;
                                        d = MACRO_MAX(
                                                bpm_256(a,b,len_a,len_b),
                                                bpm_256(b,a,len_b,len_a)
                                                );
                                        if(d <= threshold){
                                                LOG_MSG("We are merging because of:");
                                                LOG_MSG("%s",msa->sequences[samples_a[i]]->seq);
                                                LOG_MSG("%s",msa->sequences[samples_b[j]]->seq);
                                                LOG_MSG("%d\tedits",d);


                                                n->num_samples = num_a+num_b;
                                                MMALLOC(n->samples, sizeof(int) * n->num_samples);
                                                for(i = 0; i < num_a;i++){
                                                        n->samples[i] = samples_a[i];
                                                }
                                                for(j = 0;j < num_b;j++){
                                                        n->samples[j+num_a] = samples_b[j];
                                                }

                                                MFREE(n->left->samples);
                                                MFREE(n->right->samples);
                                                MFREE(n->left);
                                                MFREE(n->right);
                                                n->left = NULL;
                                                n->right = NULL;
                                                i = num_a + 1;
                                                j = num_b + 1;
                                                break;
                                        }
                                }
                        }

                }
        }

        return OK;
ERROR:
        return FAIL;
}




struct node* bisecting_kmeans(struct msa* msa, struct node* n, float** dm,int* samples,int numseq, int num_anchors,int num_samples,struct rng_state* rng)
{
        struct kmeans_result* res_tmp = NULL;
        struct kmeans_result* best = NULL;
        struct kmeans_result* res_ptr = NULL;

        struct node* tmp = NULL;

        int tries = 50;
        int t_iter;
        int r;
        int* sl = NULL;
        int* sr = NULL;
        int num_l,num_r;
        float* w = NULL;
        float* wl = NULL;
        float* wr = NULL;
        float* cl = NULL;

        float* cr = NULL;
        float dl = 0.0f;
        float dr = 0.0f;
        float score;
        int i,j,s;
        int num_var;

        int stop = 0;

        if(num_samples < 1000){
                tmp = alloc_node();
                tmp->samples = samples;
                tmp->num_samples = num_samples;
                return tmp;
        }

        num_var = num_anchors / 8;
        if( num_anchors%8){
                num_var++;
        }
        num_var = num_var << 3;


        wr = _mm_malloc(sizeof(float) * num_var,32);
        wl = _mm_malloc(sizeof(float) * num_var,32);
        cr = _mm_malloc(sizeof(float) * num_var,32);
        cl = _mm_malloc(sizeof(float) * num_var,32);


        RUNP(best = alloc_kmeans_result(num_samples));
        RUNP(res_tmp = alloc_kmeans_result(num_samples));

        best->score = FLT_MAX;

        for(t_iter = 0;t_iter < tries;t_iter++){
                res_tmp->score = FLT_MAX;
                sl = res_tmp->sl;
                sr = res_tmp->sr;

                w = _mm_malloc(sizeof(float) * num_var,32);
                for(i = 0; i < num_var;i++){
                        w[i] = 0.0f;
                        wr[i] = 0.0f;
                        wl[i] = 0.0f;
                        cr[i] = 0.0f;
                        cl[i] = 0.0f;
                }
                for(i = 0; i < num_samples;i++){
                        s = samples[i];
                        for(j = 0; j < num_anchors;j++){
                                w[j] += dm[s][j];
                        }
                }

                for(j = 0; j < num_anchors;j++){
                        w[j] /= num_samples;
                }
                r = tl_random_int(rng  , num_samples);


                s = samples[r];
                //LOG_MSG("Selected %d\n",s);
                for(j = 0; j < num_anchors;j++){
                        cl[j] = dm[s][j];
                }

                for(j = 0; j < num_anchors;j++){
                        cr[j] = w[j] - (cl[j] - w[j]);
                        //      fprintf(stdout,"%f %f  %f\n", cl[j],cr[j],w[j]);
                }

                _mm_free(w);

                /* check if cr == cl - we have identical sequences  */
                s = 0;
                for(j = 0; j < num_anchors;j++){
                        if(fabsf(cl[j]-cr[j]) >  1.0E-6){
                                s = 1;
                                break;
                        }

                }

                if(!s){

                        tmp = alloc_node();
                        tmp->samples = samples;
                        tmp->num_samples = num_samples;
                        return tmp;
                }else{
                        w = NULL;
                        while(1){
                                stop++;
                                if(stop == 10000){
                                        ERROR_MSG("Failed.");
                                }
                                num_l = 0;
                                num_r = 0;

                                for(i = 0; i < num_anchors;i++){

                                        wr[i] = 0.0f;
                                        wl[i] = 0.0f;
                                }
                                score = 0.0f;
                                for(i = 0; i < num_samples;i++){
                                        s = samples[i];
#ifdef HAVE_AVX2
                                        edist_256(dm[s], cl, num_anchors, &dl);
                                        edist_256(dm[s], cr, num_anchors, &dr);
#else
                                        edist_serial(dm[s], cl, num_anchors, &dl);
                                        edist_serial(dm[s], cr, num_anchors, &dr);
#endif
                                        score += MACRO_MIN(dl,dr);

                                        if(dr < dl){
                                                w = wr;
                                                sr[num_r] = s;
                                                num_r++;
                                        }else{
                                                w = wl;
                                                sl[num_l] = s;
                                                num_l++;
                                        }

                                        for(j = 0; j < num_anchors;j++){
                                                w[j] += dm[s][j];
                                        }

                                }

                                for(j = 0; j < num_anchors;j++){
                                        wl[j] /= num_l;
                                        wr[j] /= num_r;
                                }

                                s = 0;

                                for(j = 0; j < num_anchors;j++){
                                        if(wl[j] != cl[j]){
                                                s = 1;
                                                break;
                                        }
                                        if(wr[j] != cr[j]){
                                                s = 1;
                                                break;

                                        }
                                }
                                if(s){

                                        w = cl;
                                        cl = wl;
                                        wl = w;

                                        w = cr;
                                        cr = wr;
                                        wr = w;
                                }else{
                                        break;
                                }
                        }
                }
                res_tmp->nl =  num_l;
                res_tmp->nr =  num_r;
                res_tmp->score = score;

                if(res_tmp->score < best->score){
                        res_ptr = res_tmp;
                        res_tmp = best;
                        best = res_ptr;
                }
        }
        free_kmeans_results(res_tmp);

        sl = best->sl;
        sr = best->sr;

        num_l = best->nl;

        num_r = best->nr;

        MFREE(best);

        _mm_free(wr);
        _mm_free(wl);
        _mm_free(cr);
        _mm_free(cl);

        n = alloc_node();
        n->samples = samples;
        n->num_samples = num_samples;
        //LOG_MSG("%d left\n%d right\n", num_l,num_r);
        RUNP(n->left = bisecting_kmeans(msa,n->left, dm, sl, numseq, num_anchors, num_l,rng));

        RUNP(n->right = bisecting_kmeans(msa,n->right, dm, sr, numseq, num_anchors, num_r,rng));

        return n;
ERROR:
        return NULL;
}

struct node* upgma(float **dm,int* samples, int numseq)
{
        struct node** tree = NULL;
        struct node* tmp = NULL;

        int i,j;
        int *as = NULL;
        float max;
        int node_a = 0;
        int node_b = 0;
        int cnode = numseq;
        int numprofiles;


        numprofiles = (numseq << 1) - 1;

        MMALLOC(as,sizeof(int)*numseq);
        for (i = numseq; i--;){
                as[i] = i+1;
        }

        MMALLOC(tree,sizeof(struct node*)*numseq);
        for (i = 0;i < numseq;i++){
                tree[i] = NULL;
                tree[i] = alloc_node();
                tree[i]->id = samples[i];
        }

        while (cnode != numprofiles){
                max = FLT_MAX;
                for (i = 0;i < numseq-1; i++){
                        if (as[i]){
                                for ( j = i + 1;j < numseq;j++){
                                        if (as[j]){
                                                if (dm[i][j] < max){
                                                        max = dm[i][j];
                                                        node_a = i;
                                                        node_b = j;
                                                }
                                        }
                                }
                        }
                }
                tmp = NULL;
                tmp = alloc_node();
                tmp->left = tree[node_a];
                tmp->right = tree[node_b];


                tree[node_a] = tmp;
                tree[node_b] = NULL;

                /*deactivate  sequences to be joined*/
                as[node_a] = cnode+1;
                as[node_b] = 0;
                cnode++;

                /*calculate new distances*/
                for (j = numseq;j--;){
                        if (j != node_b){
                                dm[node_a][j] = (dm[node_a][j] + dm[node_b][j])*0.5f;

                        }
                }
                dm[node_a][node_a] = 0.0f;
                for (j = numseq;j--;){
                        dm[j][node_a] = dm[node_a][j];
                }
        }
        tmp = tree[node_a];
        MFREE(tree);
        MFREE(as);
        return tmp;
ERROR:
        return NULL;
}

struct node* alloc_node(void)
{
        struct node* n = NULL;
        MMALLOC(n, sizeof(struct node));
        n->left = NULL;
        n->right = NULL;
        n->samples = NULL;
        n->num_samples = 0;
        n->id = -1;
        return n;
ERROR:
        return NULL;
}


struct kmeans_result* alloc_kmeans_result(int num_samples)
{
        struct kmeans_result* k = NULL;
        ASSERT(num_samples != 0, "No samples???");

        MMALLOC(k, sizeof(struct kmeans_result));

        k->nl = 0;
        k->nr = 0;
        k->sl = NULL;
        k->sr = NULL;
        MMALLOC(k->sl, sizeof(int) * num_samples);
        MMALLOC(k->sr, sizeof(int) * num_samples);
        k->score = FLT_MAX;
        return k;
ERROR:
        free_kmeans_results(k);
        return NULL;
}

void free_kmeans_results(struct kmeans_result* k)
{
        if(k){
                if(k->sl){
                        MFREE(k->sl);
                }
                if(k->sr){
                        MFREE(k->sr);
                }
                MFREE(k);
        }
}
