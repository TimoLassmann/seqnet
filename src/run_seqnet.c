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

#include "global.h"
#include "msa.h"
#include "parameters.h"
#include "bpm.h"
#include <getopt.h>
#include "alphabet.h"

#include "matrix_io.h"

#define OPT_T_UNIQUE 1
#define OPT_T_TOTAL 2

#define OPT_SHOWW 5

int run_seqnet(struct parameters* param);

int print_seqnet_header(void);
int print_seqnet_help(int argc, char * argv[]);
int print_seqnet_warranty(void);
int print_AVX_warning(void);
static int calc_diff(struct msa* msa, uint8_t* seq_a,int len_a,  int i);

static int compare_seq_based_on_count(const void *a, const void *b);

int print_seqnet_help(int argc, char * argv[])
{
        const char usage[] = " -i <seq file> -o <out prefix> ";
        fprintf(stdout,"\nUsage: %s %s\n\n",basename(argv[0]) ,usage);
        fprintf(stdout,"Options:\n\n");

        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--threshold","Number of edits." ,"[1]"  );

        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--mintotal","Minimum number of sequences to form a cluster." ,"[0]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--minuniq","Minimum number of unique sequences to make up a cluster." ,"[NA]"  );

        fprintf(stdout,"\n");

        return OK;
}

int print_seqnet_warranty(void)
{
        fprintf(stdout,"Here is the Disclaimer of Warranty section of the GNU General Public License (GPL):\n");
        fprintf(stdout,"\n");
        fprintf(stdout,"15. Disclaimer of Warranty.\n");
        fprintf(stdout,"THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY\n");
        fprintf(stdout,"APPLICABLE LAW.  EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT\n");
        fprintf(stdout,"HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM \"AS IS\" WITHOUT WARRANTY\n");
        fprintf(stdout,"OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO,\n");
        fprintf(stdout,"THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR\n");
        fprintf(stdout,"PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM\n");
        fprintf(stdout,"IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF\n");
        fprintf(stdout,"ALL NECESSARY SERVICING, REPAIR OR CORRECTION.\n");
        fprintf(stdout,"\n");
        fprintf(stdout,"A complete copy of the GPL can be found in the \"COPYING\" file.\n");
        return OK;
}

int print_seqnet_header(void)
{
        fprintf(stdout,"\n");
        fprintf(stdout,"Seqnet (%s)\n", PACKAGE_VERSION);
        fprintf(stdout,"\n");
        fprintf(stdout,"Copyright (C) 2019 Timo Lassmann\n");
        fprintf(stdout,"\n");
        fprintf(stdout,"This program comes with ABSOLUTELY NO WARRANTY; for details type:\n");
        fprintf(stdout,"`seqnet -showw'.\n");
        fprintf(stdout,"This is free software, and you are welcome to redistribute it\n");
        fprintf(stdout,"under certain conditions; consult the COPYING file for details.\n");
        fprintf(stdout,"\n");


        return OK;
}

int print_AVX_warning(void)
{
        fprintf(stdout,"\n");
        fprintf(stdout,"WARNING: AVX2 instruction set not found!\n");
        fprintf(stdout,"         Seqnet will not run optimally.\n");
        fprintf(stdout,"\n");

        return OK;
}


int main(int argc, char *argv[])
{
        int c;
        int showw = 0;
        struct parameters* param = NULL;

        RUNP(param = init_param());
        param->threshold = 1;
        while (1){
                static struct option long_options[] ={
                        {"showw", 0,0,OPT_SHOWW },
                        {"input",  required_argument, 0, 'i'},
                        {"infile",  required_argument, 0, 'i'},
                        {"in",  required_argument, 0, 'i'},
                        {"threshold",  required_argument, 0, 't'},
                        {"mintotal",  required_argument, 0, OPT_T_TOTAL},
                        {"minuniq",  required_argument, 0, OPT_T_UNIQUE},
                        {"output",  required_argument, 0, 'o'},
                        {"outfile",  required_argument, 0, 'o'},
                        {"out",  required_argument, 0, 'o'},
                        {"help",   no_argument,0,'h'},
                        {"quiet",  0, 0, 'q'},
                        {0, 0, 0, 0}
                };

                int option_index = 0;

                c = getopt_long_only (argc, argv,"i:o:t:hq",long_options, &option_index);

                /* Detect the end of the options. */
                if (c == -1){
                        break;
                }
                switch(c) {
                case OPT_SHOWW:
                        showw = 1;
                        break;
                case OPT_T_TOTAL:
                        param->t_total = atof(optarg);
                        break;
                case OPT_T_UNIQUE :
                        param->t_unique = atof(optarg);
                        break;

                case 'h':
                        param->help_flag = 1;
                        break;
                case 'i':
                        param->num_infiles =1;
                        MMALLOC(param->infile, sizeof(char*));
                        param->infile[0] = optarg;

                        break;
                case 't':
                        param->threshold = atoi(optarg);
                        break;
                case 'o':
                        param->outfile = optarg;
                        break;
                case '?':
                        free_parameters(param);
                        exit(1);
                        break;
                default:

                        abort ();
                }
        }

        print_seqnet_header();
#ifndef HAVE_AVX2
        RUN(print_AVX_warning());
#endif

        if(showw){
                print_seqnet_warranty();
                free_parameters(param);
                return EXIT_SUCCESS;

        }
        if(param->help_flag){
                RUN(print_seqnet_help(argc, argv));
                free_parameters(param);
                return EXIT_SUCCESS;
        }

        if (optind < argc){
                c = param->num_infiles;
                param->num_infiles += argc-optind;
                MREALLOC(param->infile, sizeof(char*) * param->num_infiles);
                while (optind < argc){
                        param->infile[c] =  argv[optind++];
                        c++;
                }
        }


        if (param->num_infiles == 0){
                RUN(print_seqnet_help(argc, argv));
                LOG_MSG("No infiles");
                return EXIT_SUCCESS;
        }
        if (param->num_infiles >= 2){
                LOG_MSG("Too many input files:");
                for(c = 0; c < param->num_infiles;c++){
                        LOG_MSG("  %s",param->infile[c]);
                }
                LOG_MSG("Version %s only accepts one input file!\n", PACKAGE_VERSION);
                return EXIT_SUCCESS;
        }


        log_command_line(argc, argv);

        RUN(run_seqnet(param));
        free_parameters(param);
        return EXIT_SUCCESS;
ERROR:
        free_parameters(param);
        return EXIT_FAILURE;
}


int run_seqnet(struct parameters* param)
{
        struct msa* msa = NULL;
        FILE* f_ptr = NULL;

        int i,j;


        uint8_t d;
        uint8_t* seq_a;
        uint8_t* seq_b;
        int len_a,len_b;
        int num_threads = 8;
        char* tmp = NULL;
        char* buffer = NULL;
        char* t1;
        char* t2;
        int total_count;
        int count;
        int max_name_len;
        DECLARE_TIMER(t1);
        /* Step 1: read all input sequences & figure out output  */
        START_TIMER(t1);
        for(i = 0; i < param->num_infiles;i++){
                RUNP(msa = read_input(param->infile[i],msa));
        }
        STOP_TIMER(t1);
        LOG_MSG("Detected: %d sequences in %f sec.", msa->numseq,GET_TIMING(t1));
        max_name_len = 0;
        for(i = 0; i < msa->numseq;i++){
                if(max_name_len < msa->sequences[i]->name_len){
                        max_name_len = msa->sequences[i]->name_len;
                }
        }
        LOG_MSG("Longest name: %d",max_name_len);
        MMALLOC(buffer, sizeof(char) * (max_name_len));

        total_count = 0;
        //tmp = mal
        /* Get sequence counts; from my sample : ; delimited naming scheme   */
        for(i = 0; i < msa->numseq;i++){
                strncpy(buffer, msa->sequences[i]->name, max_name_len);
                tmp= buffer;
                //tmp = msa->sequences[i]->name;
                count = 0;
                t1  = strchr(tmp, ';');
                while (t1 != NULL){
                        /* String to scan is in string..token */
                        *t1++ = '\0';
                        //printf("a = %s\n", tmp);
                        t2 = strtok(tmp, ":");
                        j = 0;
                        while (t2 != NULL){
                                if(j & 1){
                                        count += atoi(t2);
                                }
                                t2 = strtok(NULL, ":");
                                j++;
                        }
                        tmp = t1;
                        t1 = strchr(tmp, ';');
                }
                t2 = strtok(tmp, ":");
                j = 0;
                while (t2 != NULL){
                        if(j & 1){
                                count += atoi(t2);
                        }
                        t2 = strtok(NULL, ":");
                        j++;
                }
                msa->sequences[i]->count = count;
                total_count += count;
        }


        qsort(msa->sequences, msa->numseq, sizeof(struct msa_seq* ),compare_seq_based_on_count);


        //>640_RS:1;mTCR_215-RS-1:1;mTCR_412-RS-5:1
        /* Kind of essential */
#ifdef HAVE_AVX2
        set_broadcast_mask();
#endif


        int num_clu = 1;
        int* seq_in_clu = NULL;
        int num_seq_in_clu = 0;
        int left = msa->numseq;
        int counts_in_clu;


        MMALLOC(seq_in_clu, sizeof(int) * msa->numseq);

        while(1){
                /* select seed  */
                j = -1;
                for(i = 0; i < msa->numseq;i++){
                        if(msa->sequences[i]->cluster == 0){
                                j = i;
                                break;
                        }
                }
                if(j == -1){
                        LOG_MSG("Quitting");
                        break;
                }
                seq_a = msa->sequences[j]->s;
                len_a = msa->sequences[j]->len;

                num_seq_in_clu =0;
                counts_in_clu = 0;
                for(i = j; i < msa->numseq;i++){
                        //START_TIMER(t1);
                        if(msa->sequences[i]->cluster == 0 ){
                                seq_b = msa->sequences[i]->s;
                                len_b = msa->sequences[i]->len;

                                d = MACRO_MAX(
                                        bpm_256(seq_a,seq_b,len_a,len_b),
                                        bpm_256(seq_b,seq_a,len_b,len_a)
                                        );

                                if(d <= param->threshold){
                                        seq_in_clu[num_seq_in_clu] = i;
                                        num_seq_in_clu++;
                                        counts_in_clu += msa->sequences[i]->count;
                                }

                                //calc_diff(msa, seq_a, len_a, i);
                                //STOP_TIMER(t1);
                        }

                }
                /* shall I print out the sequences?  */
                if(num_seq_in_clu >= param->t_unique && counts_in_clu >= param->t_total){
                        fprintf(stdout,"CLUSTER%d: %d unique %d total number of sequences\n",num_clu, num_seq_in_clu, counts_in_clu);

                        snprintf(buffer, max_name_len,"%s_cluster%d_t%d_u%d.fa",param->outfile, num_clu,counts_in_clu, num_seq_in_clu);
                        f_ptr = fopen(buffer,"w");
                        for(i = 0; i < num_seq_in_clu;i++){
                                j = seq_in_clu[i];
                                fprintf(f_ptr,">%s\n%s\n", msa->sequences[j]->name,msa->sequences[j]->seq);
                                //                left--;
                                //j = seq_in_clu[i];
                                //msa->sequences[j]->cluster = num_clu;
                                //fprintf(stdout,"%d\t%s\n",msa->sequences[j]->count,msa->sequences[j]->seq);
                                //msa->sequences[j]->count = 0;
                        }
                        //fprintf(stdout,"%d remaining\n",left);

                        fclose(f_ptr);
                        num_clu++;
                }
                for(i = 0; i < num_seq_in_clu;i++){
                        left--;
                        j = seq_in_clu[i];
                        msa->sequences[j]->cluster = num_clu;
                        //fprintf(stdout,"%d\t%s\n",msa->sequences[j]->count,msa->sequences[j]->seq);
                        //msa->sequences[j]->count = 0;
                }

                //if(num_clu == 10){
                //break;
                //}
        }

        MFREE(seq_in_clu);
        MFREE(buffer);

        free_msa(msa);
        /* If we just want to reformat end here */
        return OK;
ERROR:
        return FAIL;
}



int calc_diff(struct msa* msa, uint8_t* seq_a,int len_a,  int i)
{
        uint8_t* seq_b;
        int len_b;

        int j;
        uint8_t d;
        for(j = 0;j < msa->numseq;j++){
                seq_b = msa->sequences[j]->s;
                len_b = msa->sequences[j]->len;
                d =    bpm_256(seq_a,seq_b,len_a,len_b);
                //lk_dist[d]++;
        }
        return OK;
}



int compare_seq_based_on_count(const void *a, const void *b)
{
        struct msa_seq* const *one = a;
        struct msa_seq* const *two = b;

        if((*one)->count  < (*two)->count){
                return 1;
        }
        if((*one)->count  ==  (*two)->count){
                return 0;
        }
        return -1;

/* integer comparison: returns negative if b > a
   and positive if a > b */
}
