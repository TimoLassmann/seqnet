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

#define OPT_SHOWW 5

int run_seqnet(struct parameters* param);

int print_seqnet_header(void);
int print_seqnet_help(int argc, char * argv[]);
int print_seqnet_warranty(void);
int print_AVX_warning(void);
static int calc_diff(struct msa* msa, uint8_t* seq_a,int len_a,  int i);

int print_seqnet_help(int argc, char * argv[])
{
        const char usage[] = " -i <seq file> -o <out aln> ";
        fprintf(stdout,"\nUsage: %s %s\n\n",basename(argv[0]) ,usage);
        fprintf(stdout,"Options:\n\n");

        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--format","Output format." ,"[Fasta]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--reformat","Reformat existing alignment." ,"[NA]"  );
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
        fprintf(stdout,"Copyright (C) 2006,2019 Timo Lassmann\n");
        fprintf(stdout,"\n");
        fprintf(stdout,"This program comes with ABSOLUTELY NO WARRANTY; for details type:\n");
        fprintf(stdout,"`kalign -showw'.\n");
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

        while (1){
                static struct option long_options[] ={
                        {"showw", 0,0,OPT_SHOWW },
                        {"input",  required_argument, 0, 'i'},
                        {"infile",  required_argument, 0, 'i'},
                        {"in",  required_argument, 0, 'i'},
                        {"output",  required_argument, 0, 'o'},
                        {"outfile",  required_argument, 0, 'o'},
                        {"out",  required_argument, 0, 'o'},
                        {"help",   no_argument,0,'h'},
                        {"quiet",  0, 0, 'q'},
                        {0, 0, 0, 0}
                };

                int option_index = 0;

                c = getopt_long_only (argc, argv,"i:o:hq",long_options, &option_index);

                /* Detect the end of the options. */
                if (c == -1){
                        break;
                }
                switch(c) {
                case OPT_SHOWW:
                        showw = 1;
                        break;
                case 'h':
                        param->help_flag = 1;
                        break;
                case 'i':
                        param->num_infiles =1;
                        MMALLOC(param->infile, sizeof(char*));
                        param->infile[0] = optarg;

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

        int i,j;

        int k_dist[101];
        uint8_t d;
        uint8_t* seq_a;
        uint8_t* seq_b;
        int len_a,len_b;
        int num_threads = 8;
        DECLARE_TIMER(t1);
        /* Step 1: read all input sequences & figure out output  */
        START_TIMER(t1);
        for(i = 0; i < param->num_infiles;i++){
                RUNP(msa = read_input(param->infile[i],msa));
        }
        STOP_TIMER(t1);
        LOG_MSG("Detected: %d sequences in %f sec.", msa->numseq,GET_TIMING(t1));

        RUN(build_tree_kmeans(msa));

        for(i = 0; i < 101;i++){

                k_dist[i] = 0;
        }


#ifdef HAVE_OPENMP
#pragma omp parallel shared(num_threads,msa) private(i)

        {
#pragma omp for schedule(dynamic) nowait
#endif
        for(i = 0; i < msa->numseq;i++){
                //START_TIMER(t1);

                seq_a = msa->sequences[i]->s;
                len_a = msa->sequences[i]->len;
                calc_diff(msa, seq_a, len_a, i);
                //STOP_TIMER(t1);
                LOG_MSG("%d in %f sec.", i,1.0);//GET_TIMING(t1));

        }

        //thr_pool_wait(pool);
#ifdef HAVE_OPENMP
        }
#endif



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
                d = bpm_256(seq_a,seq_b,len_a,len_b);
                //lk_dist[d]++;
        }
        return OK;
}
